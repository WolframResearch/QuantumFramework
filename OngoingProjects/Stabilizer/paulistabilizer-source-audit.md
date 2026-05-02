# `PauliStabilizer.m` — Deep Audit

> Anchor: `Wolfram/QuantumFramework 1.6.5`, SHA `cbbe9368`. Source: `QuantumFramework/Kernel/PauliStabilizer.m` (494 LOC). Last modified 2024-05-29.
>
> Purpose: surface what works today vs. what needs development, **with bias toward symbolic capability**, before mapping requirements against the 33-paper SDK literature in `Planning for research/Stabilizer/`.

## 0. Anchors and conventions

- All file:line citations refer to `QuantumFramework/Kernel/PauliStabilizer.m` unless noted.
- "AG" = Aaronson–Gottesman tableau formalism (Phys. Rev. A 70 052328, [arxiv:quant-ph/0406196](https://arxiv.org/abs/quant-ph/0406196)).
- "Tableau" in this file means the **rank-3 array** `Dimensions == {2, n, 2n}` — first axis splits X/Z bits, second axis is qubit index (n qubits), third axis is row index (2n rows: n destabilizers + n stabilizers, in that order).
- "Signs" is a length-`2n` integer list `∈ {-1, +1}` that tracks the ±1 phase per row.
- `phase ∈ {0, 1}` is `(1 - sign) / 2` — i.e., 0 ↔ +1, 1 ↔ -1. Both forms are exposed (line 207).

## 1. Surface API: data structure and invariants

`PauliStabilizer[<|"Signs" -> {±1...}, "Tableau" -> array|>]` (line 13).

Validity is enforced by `PauliStabilizerQ` and `PauliTableauQ`:
- `PauliTableauQ[t]` (line 12): `ArrayQ[t, 3, MatchQ[0 | 1 | -1]] && MatchQ[Dimensions[t], {2, n_, m_}]`. **Note**: tableau entries are restricted to **`{0, 1, -1}` integers** — no symbolic tableau values.
- `PauliStabilizerQ` (line 13): additionally requires `Length[signs] == Dimensions[tableau][[3]]` and `signs ∈ {-1, +1}` literals.

A valid object always has **2n rows** (Stinespring-doubled: destabilizers first, stabilizers second). Single-row constructors auto-pad with the destabilizer half (lines 102-107).

**Properties contract** (line 17):
```
{"Qubits", "Generators", "Matrix", "Phase", "TableauForm", "State", "Operator", "Circuit"}
```
This list is **incomplete vs. what's actually defined**: at least 25 more property names dispatch through the property handlers (see §3 below).

## 2. Constructors

| # | Signature | Line | Path |
|---|---|---|---|
| 2.1 | `PauliStabilizer[qs_QuantumState]` | 40-48 | Pauli-expectation tomography → `PauliStabilizerTableau` (line 22) |
| 2.2 | `PauliStabilizer[<\|"Phase" -> ...\|>]` | 51 | Convert `Phase` to `Signs = 1 - 2 phase` |
| 2.3 | `PauliStabilizer[<\|"Matrix" -> mat\|>]` | 52 | Reshape flat 2n×2n binary matrix back to rank-3 tableau |
| 2.4 | `PauliStabilizer[<\|"Tableau" -> t\|>]` | 53 | Auto-fill `Signs` with `+1`s |
| 2.5 | `PauliStabilizer[qo_QuantumOperator]` | 88-92 | `QuantumOperatorTableau` (line 79) — conjugate {X,Z} of each qubit by `qo` and read tableau via `PauliRow` (line 58) |
| 2.6 | `PauliStabilizer[qco_QuantumCircuitOperator]` | 94-97 | `Fold` over `qco["NormalOperators"]`, applying each Clifford gate to a freshly constructed n-qubit `\|0…0⟩` stabilizer |
| 2.7 | `PauliStabilizer[basis]` (rank-3 ±1/0 array) | 99-100 | Resolve sign from each column's union |
| 2.8 | `PauliStabilizer[tableau]` and `[sign, tableau]` | 101-107 | Pad with destabilizer half |
| 2.9 | `PauliStabilizer[{"XZZXI", ...}]` | 109-119 | **String constructor** — Pauli-string list with optional ± prefix. Pattern: `Repeated["-" \| "+", {0, 1}] ~~ ("I" \| "X" \| "Y" \| "Z") ...`. **Qubits only.** |
| 2.10 | `PauliStabilizer[stabStrings, destabStrings]` | 121-129 | Both halves explicit |
| 2.11 | `PauliStabilizer[q_Integer]` | 131 | q-qubit `\|0…0⟩` register |
| 2.12 | `PauliStabilizer["5QubitCode" \| ...]` | 133-140 | 9 named codes (see §2.13) |
| 2.13 | `PauliStabilizer["Random", n]` | 184 | `RandomClifford[n]` — Mallows distribution sampler (line 165) |
| 2.14 | `PauliStabilizer[shortcut_String]` | 187 | Routes through `QuantumCircuitOperator[shortcut]` then 2.6 |

### 2.13 — Named codes catalog (line 9)
```
$PauliStabilizerNames = {"5QubitCode", "5QubitCode1", "SteaneCode", "7QubitCode", "7QubitCode1",
                         "SteaneCode", "SteaneCode1", "9QubitCode", "9QubitCode1", "Random"}
```
**🐛 Duplicate**: `"SteaneCode"` appears twice in the list (positions 3 and 6). Cosmetic — runtime dispatches via patterns at lines 136-137 which use `"7QubitCode" | "SteaneCode"` etc., so duplicate keys in the catalog are harmless but list looks malformed.

The `*1` variants invert the sign of the last stabilizer (logical `|1⟩_L` instead of `|0⟩_L`).

### 2.1 — State → Stabilizer is *expensive*
`PauliStabilizerTableau` (line 22) computes `⟨v|σ|v⟩` for **all 4ⁿ Pauli strings** of length n (`Tuples[Range[0, 3], n]` at line 25), runs them through `ResourceFunction["RowSpace"]`, and partitions the result. For n ≥ 8 this is very slow and memory-heavy; for n ≥ 12 it's effectively unusable.

This is the **only path** for `QuantumState → PauliStabilizer`. There is no compiled or sparse fallback. Worth flagging because users will reach for it.

## 3. Properties / accessors

The property handler `ps[prop_String]` (line 192) does direct association lookup if `prop` is already a stored key (`"Signs"`, `"Tableau"`); otherwise pattern-defined property rules dispatch.

### Direct (line 192-210)
- `"Qubits"`/`"Qudits"` (line 194)
- `"GeneratorCount"` (line 195) — **this returns `n`, not `2n`** (it's `Dimensions[Tableau][[3]] / 2`)
- `"X"`/`"Z"` (lines 208-209) — full 2n-row X-bits / Z-bits
- `"Phase"` (line 207) — `(1 - Signs) / 2`

### Stabilizer / destabilizer split
| Property | Line | Notes |
|---|---|---|
| `"Stabilizer"` / `"StabilizerTableau"` | 198 | Last n rows |
| `"StabilizerX"` / `"StabilizerZ"` | 199-200 | First / second axis of `"Stabilizer"` |
| `"StabilizerSigns"` | 197 | Last n entries of `Signs` |
| `"Destabilizer"` / `"DestabilizerTableau"` | 203 | First n rows |
| `"DestabilizerX"` / `"DestabilizerZ"` | 204-205 | |
| `"DestabilizerSigns"` | 202 | First n entries of `Signs` |

### Flat 2n×2n binary representations
- `"Matrix"` (line 211) — `2n × 2n` binary matrix in `[X | Z]` block layout. This is what the `Matrix` constructor consumes.
- `"p"` (line 212) — `Diagonal[M . Ω . Mᵀ]` where Ω is the symplectic form. Used in `Dagger`.
- `"TableauPhase"` (line 214) — `Matrix` augmented with a phase column.

### Display-form properties (§9)
`"PauliForm"`, `"Generators"`, `"Stabilizers"`, `"Destabilizers"`, `"PauliStrings"`, `"PauliSymbols"`, `"StabilizerTableauForm"`, `"TableauForm"` (lines 469-475).

## 4. Gate update rules (Clifford generators)

All gate updates are direct tableau mutations — no matrix multiplication, no expand. Each returns a *new* `PauliStabilizer`.

| Gate | Property | Line | Update rule |
|---|---|---|---|
| **H** | `ps["H", j]` | 216 | Swap `Tableau[[1,j]]` ↔ `Tableau[[2,j]]`; flip sign of any row where `x[j] == z[j] == 1` |
| **S** | `ps["S", j]` | 222 | `Z[j] ← Z[j] ⊕ X[j]`; flip sign on rows where `x[j] == z[j] == 1` |
| **S†** | `ps[SuperDagger["S"], j]` | 228 | Same Z-update; sign flip on `x[j] == 1, z[j] == 0` |
| **CNOT/CX** | `ps["CNOT", j, k]` or `["CX", j, k]` | 234 | `X[k] ← X[k] ⊕ X[j]`; `Z[j] ← Z[j] ⊕ Z[k]`; sign rule per Aaronson-Gottesman |
| **X** | `ps["X", j]` | 243 | Phase ← Phase ⊕ Z[j] (tableau unchanged) |
| **Y** | `ps["Y", j]` | 244 | Phase ← Phase ⊕ X[j] ⊕ Z[j] |
| **Z** | `ps["Z", j]` | 245 | Phase ← Phase ⊕ X[j] |
| **CZ** | `ps["CZ", j, k]` | 246 | Decomposed: `H_k . CNOT_jk . H_k` |
| **SWAP** | `ps["SWAP", j, k]` | 247 | Decomposed: `PermuteQudits[Cycles[{{j,k}}]]` |
| **V** = √X | `ps["V", j]` | 318 | Defined as `PauliStabilizer[{"-Y"}, {"X"}]` applied at qubit j |
| **V†** | `ps[SuperDagger["V"], j]` | 319 | Defined as `PauliStabilizer[{"Y"}, {"X"}]` applied at qubit j |
| **Permute** | `ps["Permute", perm]` | 249 | Permutes 2n rows |
| **PermuteQudits** | `ps["PermuteQudits", perm]` | 256 | Permutes columns (n qubits) |
| **Dagger / Inverse** | `ps["Dagger"]` or `["Inverse"]` | 263 | `Inverse[Matrix, Modulus -> 2]` for tableau; phase recomputed |
| **PadLeft / PadRight** | `ps["PadRight", n]` | 276-277 | Tensor with identity stabilizer to grow to n qubits |

**Gate name aliases**: `"CNOT" | "CX"` (line 234). Note: `"CY"` is **NOT** a Clifford generator here. There is no built-in CY rule.

**Order argument convention** (line 327): `ps[op -> order]` flattens to `ps[op, Sequence @@ Flatten[{order}]]`. So `ps["H" -> 3]` works; `ps["CNOT" -> {1, 3}]` works.

## 5. Non-Clifford gates: the symbolic boundary

The file contains three "non-Clifford" entry points. These do *not* extend stabilizer simulation — they leave the formalism, but the kernel does not flag this clearly to the user.

### 5.1 `ps["P"[phase], j]` (line 321)
```mathematica
ps_PauliStabilizer["P"[phase_], j_Integer] := With[{c = Exp[I phase / 2]},
    (1 + c) / 2 ps + (1 - c) / 2 ps["Z", j]
]
```
This returns a **symbolic linear combination** of two `PauliStabilizer` objects with complex coefficients. The result is **not a `PauliStabilizer`** — it's `<expr> + <expr>` at the top level.

What that means in practice:
- The result is **not closed under further gate application** — `(...)["H", j]` does not distribute correctly.
- `["State"]`, `["Circuit"]`, etc. on the outer expression error.
- This is the **only place** in the file where symbolic coefficients appear, and it doesn't propagate.

So `"P"[θ]` is a token gesture, not real symbolic-phase support. Useful only as a one-shot expansion for hand calculations.

### 5.2 `ps["T", j]` (line 324) and `ps[SuperDagger["T"], j]` (line 325)
Both alias `"P"[Pi/2]` and `"P"[-Pi/2]` respectively — same superposition return type, same lack of further composability.

### 5.3 `PauliStabilizerApply` (line 426) for non-Clifford gates
```mathematica
PauliStabilizerApply[qco_QuantumCircuitOperator, qs] := Fold[
    #1[Replace[#2, {"C", gate : "NOT" | "X" | "Z" -> t_, c_, _} :> "C" <> gate -> Join[c, t]]] &,
    ...
]
```
This is the dispatcher for `qco[qs, Method -> "Stabilizer"]`. It only reformats `"C"` controls. **Any operator in `qco` that doesn't match a defined update rule on `PauliStabilizer` will return *unevaluated*** — `ps[someUnknownGate, ...]` simply doesn't reduce. There's no error, no fallback. Silent failure.

This is brittle. A proper SDK needs to either (a) raise on non-Clifford, or (b) fall back to a dense path.

## 6. Measurement (Aaronson–Gottesman)

### 6.1 The `g` and `rowsum` primitives
- `g[x1, z1, x2, z2]` (line 280) is the standard AG symplectic phase function.
- `rowsum[{r, t}, h, i]` (line 282) replaces row h with row i ⊕ row h, updating sign according to `g`.
- Exposed as `ps["RowSum", h, i]` (line 288).

### 6.2 `ps["Measure", a]` / `ps["M", a]` — single qubit Pauli-Z (line 291)
Returns an **`Association`**:
- **Deterministic case** (no stabilizer anticommutes with Z_a): single key, e.g. `<|0 -> ps_post|>`. Outcome computed via fold of `rowsum` to a "scratch row", outcome = `(1 - last_sign) / 2`.
- **Non-deterministic case** (some stabilizer anticommutes): two keys, `<|0 -> ps0, 1 -> ps1|>`. The anticommuting stabilizer becomes the outcome's new generator.

### 6.3 Multi-qubit measurement (line 314)
`ps["M", {q1, q2, ...}]` — recursively measures each qubit. Result is a **flat Association keyed by tuples** of outcomes:
```
<| {0, 0} -> ps_00, {0, 1} -> ps_01, {1, 0} -> ps_10, {1, 1} -> ps_11 |>
```
With deterministic qubits, the corresponding tuple positions collapse to a single value.

### 6.4 Convenience overloads (lines 329-332)
- `ps["M", q1, q2, ...]` — variadic
- `ps[q1, q2, ...]` — positional (same as `["M", ...]`)
- `ps[]` — measure all qubits

### 6.5 What's missing
- **Pauli-X / Pauli-Y / arbitrary Pauli-string measurement** — only Z-basis on a single qubit. To measure `XZZXI` (e.g., a 5-qubit-code stabilizer) the user must rotate first. This is a major limitation for QEC stabilizer-syndrome workflows.
- **No expectation value `⟨P⟩` for arbitrary Pauli `P`** — the AG framework supports this in O(n²); not exposed.
- **Sampling** — `ps["M", q]` returns conditional outcomes, not a sampled outcome. The user has to do `RandomChoice[Keys[result] -> Values[result]]` themselves. There is no built-in batched sampler — Stim's headline feature ([arxiv:2103.02202](https://arxiv.org/abs/2103.02202)) is not implemented.

## 7. Composition and tensor product

### 7.1 `left[right]` — composition (line 389)
`left_PauliStabilizer[right_PauliStabilizer]` → applies `left` after `right` (right-to-left). Implementation: pad both to common width, multiply tableaux mod 2, accumulate phase via `phaseLookup` (a precomputed sparse rank-4 array, line 380).

This is the symplectic-mod-2 multiplication algorithm. **Performance**: `O(n³)` for n qubits (matrix multiply dominated). No specialization for sparse tableaux.

### 7.2 `QuantumTensorProduct[a, b]` (line 408)
Block-diagonal merge of destabilizers and stabilizers. Linear in `n_a + n_b`.

### 7.3 What's missing
- **Inner product `⟨ψ|φ⟩` between stabilizer states** — paper [GarMarCro12, arxiv:1210.6646](https://arxiv.org/abs/1210.6646) gives an O(n³) algorithm. Currently the only way to compute it is `ps_a["State"]["Conjugate"] @ ps_b["State"]`, which materializes 2ⁿ-dim vectors.
- **Stabilizer-frame composition** ([GarMar15, arxiv:1712.03554](https://arxiv.org/abs/1712.03554)) — not implemented.

## 8. Conversions: PauliStabilizer ↔ {QuantumState, QuantumOperator, QuantumCircuitOperator}

### 8.1 `ps["State"]` (line 335)
Computes the `+1`-eigenspace of `(I - σ_a)` for each stabilizer `σ_a`, takes a normalized total of the null space.

**Cost**: builds a `2ⁿ × 2ⁿ` identity matrix per stabilizer, materializes Pauli tensors via `kroneckerProduct @@@ Map[PauliMatrix, ...]`. Memory: `O(4ⁿ)`. **Practical limit**: n ≤ ~10. Above that, the conversion explodes.

**Typo**: the property name `"QuantumSttate"` (with double t) is included as an alias on line 335 — likely an unintentional typo.

### 8.2 `ps["Circuit"]` (line 344)
Greedy AG canonicalization: for each qubit, set destabilizer-X column to a single 1, then zero out the destabilizer-Z column, then zero out the stabilizer-X column. Produces a Clifford circuit whose dagger equals `ps`. Standard, well-tested algorithm. **Performance**: `O(n³)` gates in the worst case; circuit length is **not minimized** — papers [Reid24, arxiv:2404.19408](https://arxiv.org/abs/2404.19408), [Winderl23, arxiv:2309.08972](https://arxiv.org/abs/2309.08972) propose better synthesis.

### 8.3 `ps["Operator"]` (line 375)
Builds the circuit, then converts to a `QuantumOperator`. Materializes the matrix — same `O(4ⁿ)` memory blow-up as `["State"]` but for unitaries.

### 8.4 UpValues (lines 432-436)
- `qo_QuantumOperator[ps_]` → `PauliStabilizerApply[QuantumCircuitOperator[qo], ps]`
- `qmo_QuantumMeasurementOperator[ps_]` → likewise
- `QuantumState[ps]` → `ps["State"]`
- `QuantumCircuitOperator[ps]` → `ps["Circuit"]`
- `QuantumOperator[ps]` → `ps["Circuit"]["QuantumOperator"]`

These are the integration points with the rest of QF.

## 9. Display forms

- `PauliForm[ps, n]` (line 441) — string list like `{"XZZXI", "IXZZX", ...}` with optional ± prefix.
- `TableauForm[ps, n]` (line 453) — bordered grid display.
- TraditionalForm summary (line 480) — defers to `ToBoxes[ps["State"], TraditionalForm]`. **This means TraditionalForm of a stabilizer with n > ~10 will hang or OOM** (it materializes the state vector).
- Default summary box (line 484) — guarded: only shows `PauliForm`/`TableauForm` if `GeneratorCount < 32`; otherwise shows just qubit/generator counts.

## 10. Circuit-pipeline integration: `PauliStabilizerApply`

`qco[qs, Method -> "Stabilizer"]` routes here ([function-index.md:157](function-index.md:157)).

### 10.1 The dispatcher (line 426)
```
PauliStabilizerApply[qco_, qs : _QuantumState | _PauliStabilizer : Automatic] := Fold[
    #1[Replace[#2, {"C", gate : "NOT" | "X" | "Z" -> t_, c_, _} :> "C" <> gate -> Join[c, t]]] &,
    Replace[qs, {Automatic :> PauliStabilizer[qco["Arity"]], s_QuantumState :> PauliStabilizer[s]}],
    QuantumShortcut[qco]
]
```
Key facts:
- If `qs` is `Automatic`, starts from `|0…0⟩`.
- If `qs` is a `QuantumState`, **converts via the 4ⁿ tomography path** (§2.1) — silent expensive operation.
- Each circuit element is run through `QuantumShortcut`, which returns a label like `"H" -> 1` or `"CNOT" -> {1, 2}`. Only labels that match a defined `ps[label, ...]` rule reduce.
- **No precheck for non-Clifford gates** — passes them through; result either stays unevaluated (silent failure) or routes to `"P"[θ]`/`"T"` and dumps to a sum (§5).

## 11. Random Clifford sampling

`RandomClifford[n]` (line 165) implements the Bravyi–Maslov / Koenig–Smolin Mallows-distribution sampler ([KoeSmo14, arxiv:1406.2170](https://arxiv.org/abs/1406.2170)). Needs `O(n²)` random bits, `O(n³)` matrix work for binary inverse. Phase is randomized via `RandomInteger[1, 2n]`.

Exposed as `PauliStabilizer["Random", n]` (line 184). Default n=5.

`SampleMallows[n]` (line 147) is the underlying permutation-with-Mallows-bias sampler.

## 12. ✅ What works today

### 12.1 Core qubit Clifford simulation
- All standard Clifford generators: H, S, S†, X, Y, Z, CNOT, CZ, SWAP, V, V† (lines 216-319).
- Symplectic-mod-2 composition with proper phase tracking (line 389).
- Tensor product (line 408).
- Inverse / dagger (line 263).
- Single-qubit Z-basis measurement, deterministic and non-deterministic, returning a conditional Association (line 291).
- Multi-qubit Z-basis measurement returning nested Association (line 314).
- AG canonical decomposition into a circuit (line 344).

### 12.2 Constructors
- From explicit `Signs` + `Tableau` Association.
- From a Pauli-string list with optional ± per row (line 109).
- From explicit destabilizer-and-stabilizer string halves (line 121).
- From `q_Integer` for `|0…0⟩` (line 131).
- From a `QuantumState` (expensive — see §2.1 caveat).
- From a `QuantumOperator` if its input/output orders match.
- From a `QuantumCircuitOperator`.
- 9 named codes: `"5QubitCode"`, `"5QubitCode1"`, `"7QubitCode"` (= `"SteaneCode"`), `"7QubitCode1"`, `"9QubitCode"`, `"9QubitCode1"` (lines 133-140).
- `"Random"` Clifford via Mallows distribution.

### 12.3 Conversions
- `QuantumState[ps]`, `QuantumOperator[ps]`, `QuantumCircuitOperator[ps]` UpValues all wired (lines 432-436).
- Used as `Method -> "Stabilizer"` in `qco[qs, ...]` (`PauliStabilizerApply`, line 426).

### 12.4 Display
- `PauliForm` (string list).
- `TableauForm` (grid).
- Full SummaryBox with smart truncation at 32+ generators.
- `MakeBoxes[..., TraditionalForm]` returns the materialized state (line 480).

## 13. ❌ What needs development

### 13.1 Symbolic capability — the priority area

This is what the user flagged as the focus. The current state:

| Symbolic axis | Status | Where this hurts |
|---|---|---|
| **Symbolic Pauli coefficients in tableau** | ❌ Hard-coded `{0, 1, -1}` ints (`PauliTableauQ`, line 12) | Cannot carry symbolic weights or parameter expressions through the tableau |
| **Symbolic signs** | ❌ Hard-coded `{-1, +1}` ints (line 13) | No `±s` for parameter `s`; no `Exp[I φ]` global phase tracking |
| **Symbolic phase tracking** | ❌ `Phase` is `{0, 1}` only (line 207) | Cannot express `e^{iφ_k}` per row; needed for tracking T-gate accumulation, parameterized circuits |
| **`P[θ]` / `T` / `T†`** | ⚠️ Returns symbolic *sum of stabilizers* (line 321), but the sum is **not closed** under further gates | Cannot stay in the formalism after a single T |
| **Parametric Clifford** | ❌ No `Rx[θ]`, `Ry[θ]`, `Rz[θ]` even at θ ∈ {π/2, π, 3π/2} as named gates | User has to compose H/S/X manually |
| **Symbolic measurement outcome** | ⚠️ Outcomes are `{0, 1}` only — no `Bernoulli[p]` form | Cannot do "expected stabilizer" calculations |
| **Phase-symbolization (à la SymPhase)** | ❌ Not implemented | [arxiv:2311.03906](https://arxiv.org/abs/2311.03906) shows this is the right way to get symbolic outcomes for Pauli noise channels |

**The pivotal design choice**: the current file commits to ints in `PauliTableauQ` (line 12) and `signs` (line 13). To add symbolic support, that contract has to be loosened *and* every gate update rule (lines 216-319), the rowsum function (line 282), the composition path (line 389), and the random sampler (line 165) need to handle symbolic entries. This is a non-trivial refactor.

### 13.2 Qubit-only restriction

- `PauliStabilizer[stabStrings]` (line 109) only accepts `I/X/Y/Z`.
- `["State"]` (line 335) hard-codes `PauliMatrix[0..3]`.
- Named codes (lines 133-140) are all qubit codes.

But `PauliRow` (line 58) takes a `d` parameter and the inner functions would generalize. The qudit extension is described in [HosDehMoo04, arxiv:quant-ph/0408190](https://arxiv.org/abs/quant-ph/0408190) (already on disk) and [Beaudrap11, arxiv:1102.3354](https://arxiv.org/abs/1102.3354) (already on disk).

### 13.3 Measurement gaps

- **Pauli-string measurement** (e.g., `ps["M", "XZZXI"]`) — not supported. Critical for syndrome extraction. Would be a thin wrapper over the existing single-qubit Z-measure plus rotations.
- **Expectation value** `⟨P⟩` — not exposed.
- **Sampling** — no batched/repeated-shots simulator. Stim's [arxiv:2103.02202](https://arxiv.org/abs/2103.02202) headline feature is reference-sample + Pauli-frame propagation; would need a fundamentally different code path.

### 13.4 Performance / scale
- `["State"]` and `["Operator"]` materialize 2ⁿ×2ⁿ matrices. Practical limit ~n=10.
- `PauliStabilizer[qs_QuantumState]` runs 4ⁿ Pauli expectation values. Practical limit ~n=8.
- AG canonicalization (`["Circuit"]`, line 344) is O(n³) gates. No optimization passes; doesn't use [Winderl23](https://arxiv.org/abs/2309.08972) or [Reid24](https://arxiv.org/abs/2404.19408).
- Composition (line 389) is dense matrix multiply mod 2. No bit-packed SIMD path à la Stim.
- No caching of property values (e.g., `["Matrix"]` recomputed every call).

### 13.5 Robustness gaps
- `PauliStabilizerApply` (line 426) silently passes through unknown gates.
- `ps["Operator"]` errors if `qo` mixes input/output orders (constructor 2.5 has the precondition `Sort[qo["OutputOrder"]] == Sort[qo["InputOrder"]]`).
- `Tests/PauliStabilizer.wlt` does **not exist** — only `Tests/QuantumDistance.wlt`.

### 13.6 Documentation gaps
- **No `Usage.m` entry** — `PauliStabilizer` has no `?PauliStabilizer` message.
- **No documentation page** under `Documentation/.../Symbols/`.
- **The 3 demo notebooks (`Stabilizers.nb`, `StabilizerError.nb`, `PauliGroupTheory.nb`) reference `PauliStabilizer` zero times** — they demonstrate stabilizer concepts via raw `QuantumOperator["XZZXI"]` and group theory tools instead. Verified by grep over the converted markdown.
- The Properties contract (line 17) lists 8 properties; the actual API exposes 25+.

### 13.7 Cosmetic / typos
- `"QuantumSttate"` typo (line 335) — should be `"QuantumState"`.
- `$PauliStabilizerNames` (line 9) lists `"SteaneCode"` twice.
- `_PauliStabilizer["Properties"]` (line 17) is incomplete — out of sync with the actual property dispatch.

### 13.8 Features absent vs. the SDK literature

Mapping our 28 fetched papers against the kernel:

| Capability | Paper | Status here |
|---|---|---|
| Tableau simulation | [AarGot04](https://arxiv.org/abs/quant-ph/0406196) | ✅ Implemented |
| Mallows random Clifford | [KoeSmo14](https://arxiv.org/abs/1406.2170) | ✅ Implemented |
| Bit-packed SIMD + Pauli-frame sampling | [Gid21 Stim](https://arxiv.org/abs/2103.02202) | ❌ |
| Graph-state representation | [AndBri05](https://arxiv.org/abs/quant-ph/0504117) | ❌ |
| Inner product algorithm | [GarMarCro12](https://arxiv.org/abs/1210.6646) | ❌ |
| Stabilizer frames | [GarMar15](https://arxiv.org/abs/1712.03554) | ❌ |
| Symbolic phase | [FangYing23 SymPhase](https://arxiv.org/abs/2311.03906) | ❌ |
| GF(2) algebraic specifications | [DehMoo03](https://arxiv.org/abs/quant-ph/0304125) | Partial — `Matrix` form exists, no high-level API |
| Better synthesis from tableau | [Reid24](https://arxiv.org/abs/2404.19408), [Winderl23](https://arxiv.org/abs/2309.08972) | ❌ |
| Pauli tracking / frame propagation | [RuhDev25](https://arxiv.org/abs/2405.03970), [Paler14](https://arxiv.org/abs/1401.5872) | ❌ |
| Qudit stabilizers | [HosDehMoo04](https://arxiv.org/abs/quant-ph/0408190), [Beaudrap11](https://arxiv.org/abs/1102.3354) | ❌ |
| Projective Clifford / qudit programming | [PayWin24a](https://arxiv.org/abs/2407.16801), [WinPay24b](https://arxiv.org/abs/2407.16861) | ❌ |
| Symbolic Pauli arithmetic engine | [Mueller26 PauliEngine](https://arxiv.org/abs/2601.02233) | ❌ |
| Fast classical specifications interconversion | [deSilSalYin23](https://arxiv.org/abs/2311.10357) | ❌ |
| Stabilizer-to-Boolean-linear-algebra reduction | [Yashin25](https://arxiv.org/abs/2504.14101) | ❌ |
| Discrete Wigner derivation of AG | [KocHuaLov17](https://arxiv.org/abs/1703.04630) | ❌ |
| Graphical Pauli-measurement rules | [EllEasCav08](https://arxiv.org/abs/0806.2651), [PatGuh26](https://arxiv.org/abs/2312.02377) | ❌ |
| Pedagogical / introductory | [Got97](https://arxiv.org/abs/quant-ph/9705052), [Got00](https://arxiv.org/abs/quant-ph/0004072), [Got09](https://arxiv.org/abs/0904.2557), [MonPar23](https://arxiv.org/abs/2309.11793), [Biswas24](https://arxiv.org/abs/2405.13590) | (Reference material, not features) |

## 14. Open questions for the SDK design pass

When we read the 28 papers in detail, these are the questions to pin down:

1. **Symbolic representation choice** — extend `PauliTableauQ` to accept symbolic entries, *or* introduce a parallel symbolic class (e.g., `SymbolicPauliStabilizer`) that lowers to numeric on demand?
2. **Phase representation** — keep `{0, 1}` ints for compatibility, or move to `Exp[I θ]` symbolic with simplification rules?
3. **Non-Clifford strategy** — magic-state injection, sum-over-Cliffords decomposition (Aaronson-Gottesman §VII), or pure Pauli-frame + reference sample (Stim style)?
4. **Qudit roadmap** — promote `PauliRow`'s `d` parameter through the API, or keep a separate `QuditStabilizer` head?
5. **Performance target** — match Stim on 1000+ qubit pure-Clifford circuits, or stay in the Wolfram-symbolic regime where small-n + symbolic parameters are the use case?
6. **Test contract** — what's the minimum suite to lock down the existing implementation before refactoring? `Tests/PauliStabilizer.wlt` doesn't exist.

A read of [PatGuh26 / 2312.02377](https://arxiv.org/abs/2312.02377) is probably the highest-yield single paper for designing the API, since it covers the full set of stabilizer manipulations (not just simulation) with explicit graphical rules.
