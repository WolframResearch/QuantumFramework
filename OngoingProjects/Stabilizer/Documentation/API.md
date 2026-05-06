# Stabilizer subsystem — API reference

> Function-by-function reference for the 6 top-level public symbols + the method-grade operations dispatched on `PauliStabilizer` (4 of which were Phase-3 top-level symbols, demoted to methods in Phase 6 — see [§API consolidation in ROADMAP.md](ROADMAP.md#a0--api-consolidation-2026-05-06)). Every code example is verified by [`verify-API.wls`](verify-API.wls) — `wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-API.wls`.

## Public API at a glance

| Symbol | Phase | Purpose |
|---|---|---|
| [`PauliStabilizer`](#paulistabilizer) | 1 | Stabilizer-state head: tableau-encoded n-qubit state |
| [`StabilizerFrame`](#stabilizerframe) | 4 | Superposition of stabilizer states (for non-Clifford) |
| [`StabilizerInnerProduct`](#stabilizerinnerproduct) | 4 | `<ψ\|φ>` for two stabilizer/frame states |
| [`StabilizerExpectation`](#stabilizerexpectation) | 4 | `<ψ\|P\|ψ>` for an arbitrary Pauli string |
| [`GraphState`](#graphstate) | 5 | Graph-state representation (AndBri05) |
| [`LocalComplement`](#localcomplement) | 5 | Local complementation on a graph or graph state |

### Method-grade operations on `PauliStabilizer` (Phase 6, formerly top-level)

| Method | Old top-level form | Purpose |
|---|---|---|
| [`PauliStabilizer["Random", n]`](#paulistabilizername_string--named-stabilizer-code) | `RandomClifford[n]` | Uniformly random n-qubit Clifford state (Mallows sampler) |
| [`ps["SymbolicMeasure", q]`](#methods-symbolic-measurement) | `StabilizerMeasure[ps, q]` | Symbolic Z-basis measurement (allocates fresh outcome symbol) |
| [`ps["SubstituteOutcomes", rules]`](#methods-symbolic-measurement) | `SubstituteOutcomes[ps, rules]` | Plug concrete values into measurement-outcome symbols |
| [`ps["SampleOutcomes", n]`](#methods-symbolic-measurement) | `SampleOutcomes[ps, n]` | Random samples by independent symbol substitution |

**Companion:** [`synthesis-implementation.md`](synthesis-implementation.md) walks through synthesis §1–§11 by capability. [`ROADMAP.md`](ROADMAP.md) tracks the partial / deferred / known-bug items.

**Re-verify:** `wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-API.wls`

---

# PauliStabilizer

The atomic stabilizer-state object. Encodes an n-qubit state as a `(2n × 2n)` symplectic matrix plus a length-`2n` sign vector, packaged as `PauliStabilizer[<|"Tableau" -> ..., "Signs" -> ...|>]`.

## Constructors

### `PauliStabilizer[stabStrings_List]` — from a list of Pauli strings

Build a stabilizer state from a list of Pauli strings (one per stabilizer). Optional `+` / `-` prefix sets the sign of each stabilizer.

```wolfram
PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"]
```
```
{"XX", "ZZ"}
```

```wolfram
PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"]
```
```
{-1, 1}
```

### `PauliStabilizer[stabStrings, destabStrings]` — explicit stab + destab halves

```wolfram
With[{ps = PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]},
    {ps["Stabilizers"], ps["Destabilizers"]}
]
```
```
{{"XX", "ZZ"}, {"IX", "ZI"}}
```

### `PauliStabilizer[n_Integer]` — n-qubit `|0...0⟩` register

```wolfram
With[{ps = PauliStabilizer[3]},
    {ps["Qubits"], ps["Stabilizers"], ps["Destabilizers"]}
]
```
```
{3, {"ZII", "IZI", "IIZ"}, {"XII", "IXI", "IIX"}}
```

### `PauliStabilizer[name_String]` — named stabilizer code

Available names: `"5QubitCode"`, `"5QubitCode1"`, `"SteaneCode"` (= `"7QubitCode"`), `"SteaneCode1"` (= `"7QubitCode1"`), `"9QubitCode"`, `"9QubitCode1"`, `"Random"`. The `*1` suffix denotes the logical `|1_L⟩` codeword.

```wolfram
With[{ps = PauliStabilizer["5QubitCode"]},
    {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"]}
]
```
```
{5, 5, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"}}
```

```wolfram
With[{ps = PauliStabilizer["SteaneCode"]},
    {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"][[7]]}
]
```
```
{7, 7, "XXXXXXX"}
```

### `PauliStabilizer[qco_QuantumCircuitOperator]` — apply a Clifford circuit to `|0...0⟩`

Folds the circuit's normal operators over the n-qubit register. Each operator must be a Clifford gate; non-Clifford gates trigger `PauliStabilizer::nonclifford`.

```wolfram
PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"]
```
```
{"XX", "ZZ"}
```

### `PauliStabilizer[qo_QuantumOperator]` — Clifford gate to tableau

Conjugate the qubit-aligned Pauli generators by `qo`; extract the resulting tableau.

```wolfram
PauliStabilizer[QuantumOperator["H", 1]]["Stabilizers"]
```
```
{"X"}
```

### `PauliStabilizer["Random", n_Integer]` — uniformly-random n-qubit Clifford state

Bravyi–Maslov / Koenig–Smolin Mallows sampler. See also [`RandomClifford`](#randomclifford).

```wolfram
SeedRandom[20260430];
With[{ps = PauliStabilizer["Random", 3]},
    {ps["Qubits"], ps["Stabilizers"]}
]
```
```
{3, {"-ZXZ", "-ZZX", "XXX"}}
```

## Properties

Access via `ps[propName]`. Full list via `ps["Properties"]`.

### Shape / structural

```wolfram
With[{ps = PauliStabilizer["5QubitCode"]},
    <|
        "Qubits"           -> ps["Qubits"],
        "GeneratorCount"   -> ps["GeneratorCount"],
        "Tableau-shape"    -> Dimensions[ps["Tableau"]],
        "Matrix-shape"     -> Dimensions[ps["Matrix"]],
        "Signs-length"     -> Length[ps["Signs"]]
    |>
]
```
```
<|"Qubits" -> 5, "GeneratorCount" -> 5,
  "Tableau-shape" -> {2, 5, 10}, "Matrix-shape" -> {10, 10},
  "Signs-length" -> 10|>
```

The tableau is shape `{2, n, 2n}` (X/Z block · qubit · row); the matrix is the flattened `[X|Z]` form of shape `{2n, 2n}`.

### Stabilizer / destabilizer split

```wolfram
With[{ps = PauliStabilizer["5QubitCode"]},
    <|
        "Stabilizers"        -> ps["Stabilizers"],
        "Destabilizers"      -> ps["Destabilizers"],
        "StabilizerSigns"    -> ps["StabilizerSigns"],
        "DestabilizerSigns"  -> ps["DestabilizerSigns"]
    |>
]
```
```
<|"Stabilizers" -> {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ", "XXXXX"},
  "Destabilizers" -> {"ZXXZI", "IZXXZ", "ZIZXX", "XZIZX", "ZZZZZ"},
  "StabilizerSigns" -> {1, 1, 1, 1, 1},
  "DestabilizerSigns" -> {1, 1, 1, 1, 1}|>
```

### Bit views

```wolfram
With[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    <|"X" -> ps["X"], "Z" -> ps["Z"], "Phase" -> ps["Phase"]|>
]
```
```
<|"X" -> {{0, 0, 1, 0}, {0, 1, 1, 0}},
  "Z" -> {{1, 0, 0, 1}, {0, 0, 0, 1}},
  "Phase" -> {0, 0, 0, 0}|>
```

`ps["X"]` and `ps["Z"]` are shape `{n_qubits, 2*GeneratorCount}` — rows are qubits, columns are tableau rows (destabilizers first, then stabilizers). `ps["Phase"] = (1 - ps["Signs"]) / 2`.

### Pauli string list

```wolfram
With[{ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    ps["PauliStrings"]
]
```
```
{"XX", "ZZ", "ZI", "IX"}
```

Concatenated `[Stabilizers, Destabilizers]`.

### Full property list

| Property | Returns |
|---|---|
| `"Qubits"`, `"Qudits"` | n |
| `"GeneratorCount"` | n (= rank of stabilizer group) |
| `"Tableau"` | rank-3 binary array, shape `{2, n, 2n}` |
| `"Matrix"` | 2n×2n binary matrix in `[X\|Z]` block layout |
| `"Signs"` | length-2n list of ±1 (or symbolic in Phase 3+) |
| `"Phase"` | `(1 - Signs) / 2` |
| `"X"`, `"Z"` | X-bits and Z-bits (qubits × rows) |
| `"Stabilizer"` / `"StabilizerTableau"` | last n columns of Tableau |
| `"StabilizerX"`, `"StabilizerZ"` | X / Z bits of the n stabilizer rows |
| `"StabilizerSigns"` | last n signs |
| `"Destabilizer"` etc. | analogous for the first n rows |
| `"PauliForm"`, `"Generators"`, `"Stabilizers"` | string list of stabilizers (with sign prefix) |
| `"Destabilizers"` | string list |
| `"PauliStrings"` | concatenated `Stabilizers ++ Destabilizers` |
| `"PauliSymbols"` | same but as `±symbol` entries |
| `"TableauForm"`, `"StabilizerTableauForm"` | grid display |
| `"State"` / `"QuantumState"` | materialized `QuantumState` (cost `2ⁿ`) |
| `"Circuit"` / `"QuantumCircuitOperator"` | AG-greedy synthesized Clifford circuit |
| `"Operator"` / `"QuantumOperator"` | matrix-form `QuantumOperator` (cost `4ⁿ`) |
| `"GlobalPhase"` | optional complex unit set when constructed from `QuantumState` / `QuantumOperator`; multiplies `["State"]` and `["QuantumOperator"]` outputs (Phase 5c, default 1) |
| `"Properties"` | this list |

## Methods (Clifford gates)

Each method returns a new `PauliStabilizer` (Clifford gates) or a `StabilizerFrame` (non-Clifford).

| Method | Signature | Effect |
|---|---|---|
| `"H"` | `ps["H", q]` | Hadamard on qubit `q` |
| `"S"` | `ps["S", q]` | π/4 phase gate |
| `SuperDagger["S"]` | `ps[SuperDagger["S"], q]` | S-dagger |
| `"X"`, `"Y"`, `"Z"` | `ps["X", q]` | Pauli on qubit q (only flips phase) |
| `"CNOT"` / `"CX"` | `ps["CNOT", c, t]` | controlled-NOT |
| `"CZ"` | `ps["CZ", c, t]` | controlled-Z |
| `"SWAP"` | `ps["SWAP", a, b]` | swap |
| `"V"` | `ps["V", q]` | √X |
| `SuperDagger["V"]` | `ps[SuperDagger["V"], q]` | √X-dagger |
| `"Permute"`, `"PermuteQudits"` | `ps["PermuteQudits", Cycles[{...}]]` | row / column permutation |
| `"Dagger"` / `"Inverse"` | `ps["Dagger"]` | tableau inverse (⚠ Phase 1 known bug at `Dagger ∘ Dagger` — see [ROADMAP §A.3](ROADMAP.md)) |
| `"PadLeft"`, `"PadRight"` | `ps["PadRight", n]` | tensor identity to grow to n qubits |
| `"P"[θ]` | `ps["P"[θ], q]` | non-Clifford rotation; returns `StabilizerFrame` |
| `"T"` | `ps["T", q]` | alias for `"P"[Pi/2]` |
| `SuperDagger["T"]` | `ps[SuperDagger["T"], q]` | alias for `"P"[-Pi/2]` |

### Examples

**H takes Z to X (Heisenberg conjugation):**
```wolfram
PauliStabilizer[1]["H", 1]["Stabilizers"]
```
```
{"X"}
```

**S takes X to Y:**
```wolfram
PauliStabilizer[1]["H", 1]["S", 1]["Stabilizers"]
```
```
{"Y"}
```

**CNOT(1,2) on H|0⟩⊗|0⟩ = Bell state:**
```wolfram
Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
```
```
{"XX", "ZZ"}
```

**CZ(1,2) on H|0⟩⊗|0⟩:**
```wolfram
Sort @ PauliStabilizer[2]["H", 1]["CZ", 1, 2]["Stabilizers"]
```
```
{"IZ", "XZ"}
```

**SWAP(1,2) on H|0⟩⊗|0⟩:**
```wolfram
Sort @ PauliStabilizer[2]["H", 1]["SWAP", 1, 2]["Stabilizers"]
```
```
{"IX", "ZI"}
```

**X gate flips signs of stabilizers with Z-content at q (no tableau change):**
```wolfram
With[{ps5 = PauliStabilizer["5QubitCode"]},
    ps5["X", 1]["StabilizerSigns"]
]
```
```
{1, 1, 1, -1, 1}
```

`X_1` anticommutes with the 4th stabilizer `ZXIXZ` (because it has Z at qubit 1).

**Z gate similarly:**
```wolfram
With[{ps5 = PauliStabilizer["5QubitCode"]},
    ps5["Z", 3]["StabilizerSigns"]
]
```
```
{1, 1, -1, 1, -1}
```

`Z_3` anticommutes with stabilizers 3 (`XIXZZ`, has X at q3) and 5 (`XXXXX`, has X at q3).

**S†, V, V†:**
```wolfram
With[{ps = PauliStabilizer[1]["H", 1]},   (* |+>, stab X *)
    <|
        "Sdag-on-X" -> ps[SuperDagger["S"], 1]["Stabilizers"],
        "V-on-X"    -> ps["V", 1]["Stabilizers"],
        "Vdag-on-X" -> ps[SuperDagger["V"], 1]["Stabilizers"]
    |>
]
```
```
<|"Sdag-on-X" -> {"-Y"}, "V-on-X" -> {"X"}, "Vdag-on-X" -> {"X"}|>
```

S†: `S†XS = -Y`. V (= √X): `V X V† = X` (X commutes with √X).

**Convenience syntax `op -> order` flattens to `op, args`:**
```wolfram
With[{ps = PauliStabilizer[2]["H" -> 1]["CNOT" -> {1, 2}]},
    Sort @ ps["Stabilizers"]
]
```
```
{"XX", "ZZ"}
```

**T returns a `StabilizerFrame` (non-Clifford boundary):**
```wolfram
Head @ PauliStabilizer[1]["T", 1]
```
```
StabilizerFrame
```

The original `Plus` return from Phase 1 was migrated to `StabilizerFrame` in Phase 4. See [`StabilizerFrame`](#stabilizerframe).

## Methods (Measurement)

### Z-basis on a single qubit: `ps["M", q]` or `ps[q]`

Returns an `Association` keyed by outcome bit (0 or 1) and valued by the post-measurement `PauliStabilizer`.

**Deterministic:**
```wolfram
PauliStabilizer[1]["M", 1]
```
```
<|0 -> PauliStabilizer[<|"Signs" -> {1, 1}, "Tableau" -> {{{1, 0}}, {{0, 1}}}|>]|>
```

**Random:**
```wolfram
Sort @ Keys @ PauliStabilizer[1]["H", 1]["M", 1]
```
```
{0, 1}
```

### Pauli-string measurement: `ps["M", "XZZXI"]`

Measures an arbitrary Pauli observable on the stabilizer state. Same Association return type. For the 5-qubit code state, all 5 stabilizers are deterministic with outcome bit 0:

```wolfram
With[{ps5 = PauliStabilizer["5QubitCode"]},
    Keys @ ps5["M", "XZZXI"]
]
```
```
{0}
```

### Multi-qubit: `ps["M", {q1, q2, ...}]`

Returns Association keyed by outcome tuples:

```wolfram
With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    Sort @ Keys @ psBell["M", {1, 2}]
]
```
```
{{0, 0}, {1, 1}}
```

Bell ZZ correlation: outcomes are perfectly correlated.

## Composition + tensor product

**Compose: `ps1[ps2]`** — applies ps1 after ps2 (right-to-left).

```wolfram
Module[{psH = PauliStabilizer[1]["H", 1]["S", 1]["H", 1]},
    psH["Stabilizers"]
]
```
```
{"-Y"}
```

**Tensor product: `QuantumTensorProduct[ps1, ps2]`**

```wolfram
With[{a = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]],
      b = PauliStabilizer[1]},
    QuantumTensorProduct[a, b]["Stabilizers"]
]
```
```
{"XXI", "ZZI", "IIZ"}
```

## Conversions

### `ps["State"]` — materialize to `QuantumState` (cost `2ⁿ`)

```wolfram
With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    Normal @ psBell["State"]["StateVector"]
]
```
```
{1/Sqrt[2], 0, 0, 1/Sqrt[2]}
```

### `ps["Circuit"]` — synthesize Clifford circuit from tableau (AG greedy)

The circuit's *Hermitian conjugate* applied to `|0...0⟩` reproduces `ps`. Note that `PauliStabilizer[ps["Circuit"]]` does **not** generally produce string-equal stabilizers (different generating set of the same group), but the materialized state vector matches up to global phase:

```wolfram
Module[{psBell, circ, fromCircuit, vec1, vec2},
    psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
    circ = psBell["Circuit"];
    fromCircuit = PauliStabilizer[circ];
    vec1 = N @ Normal @ psBell["State"]["StateVector"];
    vec2 = N @ Normal @ fromCircuit["State"]["StateVector"];
    Chop @ Abs[Conjugate[vec1] . vec2] - 1
]
```
```
-2.220446049250313*^-16
```

Numerical zero — the round-trip is correct up to global phase. For optimized synthesis methods (Reid24, Winderl23), see [ROADMAP §A.5](ROADMAP.md).

### `ps["Operator"]` — full unitary matrix (cost `4ⁿ`)

Returns a `QuantumOperator` whose action on `|0...0⟩` reproduces the state. Practical limit `n ≤ ~10`.

Returns a `QuantumOperator` whose action on `|0...0⟩` reproduces the state. Practical limit `n ≤ ~10`. The round-trip contract is documented in the next section.

## Round-trip contract

Phase 5c gives `PauliStabilizer` a phase-aware contract for the first hop (constructor → accessor) and an explicit up-to-phase contract for subsequent gate updates.

| Path | Equality |
|---|---|
| `PauliStabilizer[qo_QuantumOperator]["QuantumOperator"]` vs `qo` | **exact** (`===`) |
| `PauliStabilizer[qs_QuantumState]["State"]` vs `qs` | **exact** (`===`) |
| `PauliStabilizer[qs]["gate", q]["State"]` vs `gate @ qs` | up to global phase |
| `PauliStabilizer[gate @ qs]["State"]` vs `gate @ qs` | **exact** — the escape hatch for gate-update phase loss |
| `PauliStabilizer[ps["Circuit"]]["State"]` vs `ps["State"]` | up to global phase (no source state, only tableau) |

### First hop is exact (Phase 5c)

The Y gate is the canary: AG-decomposition recovers `Z·X = i·Y`, but the constructor captures the missing `−i` factor under `"GlobalPhase"`, so `["QuantumOperator"]` returns `Y` exactly:

```wolfram
With[{ps = PauliStabilizer[QuantumOperator["Y"]]},
    <|"GlobalPhase" -> ps["GlobalPhase"], "matches" -> ps["QuantumOperator"]["Matrix"] === QuantumOperator["Y"]["Matrix"]|>
]
```
```
<|"GlobalPhase" -> -I, "matches" -> True|>
```

State construction is the same:

```wolfram
With[{ps = PauliStabilizer[QuantumState[{0, 1}]]},
    <|"GlobalPhase" -> ps["GlobalPhase"], "matches" -> ps["State"]["StateVector"] === QuantumState[{0, 1}]["StateVector"]|>
]
```
```
<|"GlobalPhase" -> 1, "matches" -> True|>
```

### Gate-update path is up to phase (ROADMAP §A.9)

Clifford gate updates do **not** propagate `"GlobalPhase"`. Recovering the new phase exactly would require materializing the state at every gate update, which is `O(2ⁿ)` and defeats the AG `O(n²)` advantage. Documented as an inherent trade-off:

```wolfram
With[{
    actual   = PauliStabilizer[QuantumState[{0, 1}]]["Y", 1]["State"]["StateVector"] // Normal,
    expected = (QuantumOperator["Y"][QuantumState[{0, 1}]])["StateVector"] // Normal
},
    <|"actual" -> actual, "expected" -> expected, "abs-overlap" -> Quiet @ Simplify[Abs[Conjugate[actual] . expected]]|>
]
```
```
<|"actual" -> {1, 0}, "expected" -> {-I, 0}, "abs-overlap" -> 1|>
```

The actual is `\|0⟩` and the expected is `−i\|0⟩`; they differ by a phase but `\|⟨a\|b⟩\| = 1`.

### Escape hatch for exact equality after a gate

Re-construct the `PauliStabilizer` from the post-gate state:

```wolfram
With[{
    qs = QuantumState[{0, 1}],
    gate = QuantumOperator["Y"]
},
    PauliStabilizer[gate[qs]]["State"]["StateVector"] // Normal // Simplify
]
```
```
{-I, 0}
```

This re-tomographs the Y-rotated state, captures the new `GlobalPhase`, and `["State"]` returns `−i\|0⟩` exactly.

For the full test surface see [TIER 1.4a/1.4b/1.4c/1.4d in `Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) and the cross-module probes in [`Tests/Roundtrips.wlt`](../../../Tests/Roundtrips.wlt). Root-cause + design rationale: [`post-mortem-phase-5c.md`](post-mortem-phase-5c.md).

## Integration with QuantumFramework

```wolfram
QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"]
```
```
{"XX", "ZZ"}
```

```wolfram
QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"]
```
```
{"XXX", "ZZI", "IZZ"}
```

UpValues are wired so `QuantumOperator[gate, q] @ ps` works:

```wolfram
With[{ps = PauliStabilizer[1]},
    (QuantumOperator["H", 1] @ ps)["Stabilizers"]
]
```
```
{"X"}
```

## Messages

- `PauliStabilizer::nonclifford` — emitted when a circuit contains a gate that doesn't match any tableau-update rule (and isn't `P[θ]`/`T`/`T†`). Returns the input state unchanged so the user can recover.
- `PauliStabilizer::tdeprecated` — historical message announcing the Phase 4 migration of `P[θ]/T/T†` from `Plus` to `StabilizerFrame`.

## Methods (Symbolic measurement)

Phase 3 SymPhase methods (FangYing23 §3, arxiv:2311.03906). All three were Phase-3 top-level public symbols (`StabilizerMeasure`, `SubstituteOutcomes`, `SampleOutcomes`); demoted in Phase 6 to methods on `PauliStabilizer` per the API consolidation.

### `ps["SymbolicMeasure", q_Integer]` and `ps["SymbolicMeasure", qudits_List]`

Symbolic Z-basis measurement. Allocates a fresh `\[FormalS][k]` symbol per non-deterministic measurement; returns one `PauliStabilizer` instead of an `Association`.

**Deterministic case** (no anticommuting stabilizer): returns the post-state directly with no fresh symbol allocated.

```wolfram
PauliStabilizer[1]["SymbolicMeasure", 1]["Stabilizers"]
```
```
{"Z"}
```

**Non-deterministic case** (one anticommuting stabilizer): allocates a fresh `\[FormalS][k]` symbol and stamps it into the appropriate phase position.

```wolfram
Module[{psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
    <|
        "Head"          -> Head[psSym],
        "Phase"         -> psSym["Phase"],
        "FreshSymbols"  -> DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity]
    |>
]
```
```
<|"Head" -> PauliStabilizer, "Phase" -> {0, \[FormalS][1]},
  "FreshSymbols" -> {\[FormalS][1]}|>
```

**Phase 3 known limitation:** when a deterministic measurement follows a prior symbolic one (e.g. Bell ZZ correlation), the deterministic outcome polynomial is computed correctly by AG but is **not** stamped into the post-state's signs. Tracked in [ROADMAP §A.4](ROADMAP.md). Fix lives in Phase 4 `StabilizerFrame` integration.

### `ps["SubstituteOutcomes", rules]`

Replace measurement-outcome symbols with concrete 0/1 values, reducing signs back to {-1, 1} via `Mod 2`. Substituting a symbol with 0 reproduces the outcome-0 branch; with 1 reproduces outcome-1:

```wolfram
Module[{psSym, sym, ps0, ps1Reg},
    psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
    sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
    ps0 = psSym["SubstituteOutcomes", sym -> 0];
    ps1Reg = psSym["SubstituteOutcomes", sym -> 1];
    <|
        "outcome 0 stabilizers" -> ps0["Stabilizers"],
        "outcome 1 stabilizers" -> ps1Reg["Stabilizers"]
    |>
]
```
```
<|"outcome 0 stabilizers" -> {"Z"}, "outcome 1 stabilizers" -> {"-Z"}|>
```

### `ps["SampleOutcomes"]` and `ps["SampleOutcomes", n_Integer]`

Draw `n` random samples by independently substituting each outcome symbol with a uniformly-random 0 or 1. Single-arg form returns one sample; with `n`, returns a list.

```wolfram
SeedRandom[42];
With[{psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]},
    #["Stabilizers"] & /@ psSym["SampleOutcomes", 10]
]
```
```
{{"Z"}, {"-Z"}, {"-Z"}, {"Z"}, {"-Z"}, {"-Z"}, {"Z"}, {"-Z"}, {"Z"}, {"Z"}}
```

## Methods (Random Clifford)

Uniformly random sampler from the n-qubit Clifford group via the Bravyi–Maslov / Koenig–Smolin Mallows-distribution algorithm. The standalone `RandomClifford[n]` was demoted in Phase 6; reach the same sampler via the named-pattern constructor:

```wolfram
SeedRandom[42];
With[{r = PauliStabilizer["Random", 3]},
    <|"Qubits" -> r["Qubits"], "Stabilizers" -> r["Stabilizers"]|>
]
```
```
<|"Qubits" -> 3, "Stabilizers" -> {"-XIY", "XYY", "-XII"}|>
```

**Cardinality** (KoeSmo14 Eq 2): `|C_n| = 2^(n² + 2n) · Π_{j=1}^n (4^j − 1)`. For `n = 1`, 24 elements. For `n = 2`, 11520. For `n = 3`, ~9.29 × 10⁷.

**Uniformity smoke test:** 200 samples on n=1 hit most of the small number of distinct stabilizer states.

```wolfram
SeedRandom[42];
Length @ DeleteDuplicates @ Table[
    With[{r = PauliStabilizer["Random", 1]}, {r["Stabilizers"], r["Signs"]}],
    {200}
]
```
```
12
```

See KoeSmo14 §3.2 (arxiv:1406.2170) and [ROADMAP §B.5](ROADMAP.md) (`IndexClifford` inverse map, deferred).

10 samples, mix of outcome-0 (`{"Z"}`) and outcome-1 (`{"-Z"}`) — approximately 50/50 split.

---

# StabilizerFrame

Phase 4 superposition-of-stabilizer-states head: `Σ_i c_i |s_i⟩` with (possibly symbolic) coefficients `c_i` and stabilizer states `|s_i⟩`. Produced by non-Clifford gates (`P[θ]`, `T`, `T†`) on a `PauliStabilizer`. Closes under further Clifford gates (which distribute over components) and under further non-Clifford gates (which double the component count).

## Constructors

### From a list of {coefficient, PauliStabilizer} pairs

```wolfram
With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
    With[{fH = f["H", 1]},
        <|"Length" -> fH["Length"],
          "Stabilizers-each" -> #["Stabilizers"] & /@ fH["Stabilizers"]|>
    ]
]
```
```
<|"Length" -> 2, "Stabilizers-each" -> {"Z"}|>
```

### From a single PauliStabilizer (coefficient = 1)

```wolfram
With[{f = StabilizerFrame[PauliStabilizer[1]]},
    <|"Length" -> f["Length"], "Components" -> f["Components"]|>
]
```
```
<|"Length" -> 1,
  "Components" -> {{1, PauliStabilizer[<|"Signs" -> {1, 1},
                                          "Tableau" -> {{{1, 0}}, {{0, 1}}}|>]}}|>
```

### From a non-Clifford gate

```wolfram
With[{psT = PauliStabilizer[1]["T", 1]},
    <|
        "Head"         -> Head[psT],
        "Length"       -> psT["Length"],
        "Coefficients" -> psT["Coefficients"]
    |>
]
```
```
<|"Head" -> StabilizerFrame, "Length" -> 2,
  "Coefficients" -> {(1 + E^((I/4)*Pi))/2, (1 - E^((I/4)*Pi))/2}|>
```

## Properties

| Property | Returns |
|---|---|
| `"Components"` | List of `{coeff, PauliStabilizer}` pairs |
| `"Coefficients"` | List of coefficients only |
| `"Stabilizers"` | List of underlying `PauliStabilizer`s |
| `"Length"` | Number of components |
| `"Qubits"`, `"Qudits"` | Qubit count of components (assumed identical) |
| `"GeneratorCount"` | Generator count of components |
| `"StateVector"` | Materialized state vector (cost `2ⁿ`) |
| `"State"` | Materialized `QuantumState` |

## Methods

### Clifford gates distribute over components

```wolfram
With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
    With[{fH = f["H", 1]},
        <|"Length" -> fH["Length"],
          "Stabilizers-each" -> #["Stabilizers"] & /@ fH["Stabilizers"]|>
    ]
]
```
```
<|"Length" -> 2, "Stabilizers-each" -> {"Z"}|>
```

(Both components had `H` applied; output frame still has 2 components.)

### Non-Clifford gates double component count

```wolfram
With[{psT2 = PauliStabilizer[1]["H", 1]["T", 1]["T", 1]},
    {Head[psT2], psT2["Length"]}
]
```
```
{StabilizerFrame, 4}
```

`T ∘ T` on `|+⟩`: each `T` doubles the frame, ending at 4 components.

### Materialization

`T|0⟩ = |0⟩` (eigenstate). The frame materializes correctly:

```wolfram
Module[{psT = PauliStabilizer[1]["T", 1], vec},
    vec = Normal @ psT["StateVector"];
    Chop @ N @ FullSimplify[vec - {1, 0}]
]
```
```
{0, 0}
```

After `FullSimplify`, the diff is exactly zero — confirming `T|0⟩ = |0⟩`.

## Arithmetic upvalues

- `Plus` of two frames: concatenates their components (`StabilizerFrame /: Plus[a, b]`).
- `Times` by a scalar: scales all coefficients (`StabilizerFrame /: Times[c, f]`).
- `Equal`: structural equality on components.

## See also

- [`StabilizerInnerProduct`](#stabilizerinnerproduct) — works for both `PauliStabilizer` and `StabilizerFrame`.
- GarMar15 §3 (arxiv:1712.03554) — Quipu stabilizer frames.

---

# StabilizerInnerProduct

Compute `⟨ψ|φ⟩` for two stabilizer states (or stabilizer frames). Phase 4 v1 uses direct vector materialization (cost `2ⁿ`); the closed-form `O(n³)` GarMarCro12 algorithm is on [ROADMAP §A.1](ROADMAP.md).

## Signatures

```wolfram
StabilizerInnerProduct[psi_PauliStabilizer, phi_PauliStabilizer]
StabilizerInnerProduct[fA_StabilizerFrame, fB_StabilizerFrame]
StabilizerInnerProduct[psi_PauliStabilizer, fB_StabilizerFrame]
StabilizerInnerProduct[fA_StabilizerFrame, phi_PauliStabilizer]
```

## Examples

```wolfram
With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    StabilizerInnerProduct[psBell, psBell]
]
```
```
1
```

```wolfram
Module[{psPhiPlus, psPhiMinus},
    psPhiPlus = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
    psPhiMinus = psPhiPlus["Z", 1];
    Chop @ N @ StabilizerInnerProduct[psPhiPlus, psPhiMinus]
]
```
```
0
```

(`|Φ+⟩` and `|Φ−⟩` are orthogonal Bell states.)

**Mixed PauliStabilizer/StabilizerFrame:**

```wolfram
With[{psT = PauliStabilizer[1]["T", 1], ps0 = PauliStabilizer[1]},
    FullSimplify @ StabilizerInnerProduct[psT, ps0]
]
```
```
1
```

`⟨T|0⟩|0⟩ = ⟨0|T|0⟩ = 1` since `T|0⟩ = |0⟩`.

---

# StabilizerExpectation

Compute `⟨ψ|P|ψ⟩` for an arbitrary Pauli string P. Returns ±1 for stabilizer-group elements, 0 for anticommuting Paulis, and the exact expectation via direct vector for Paulis in `N(S) \ S`.

## Signature

```wolfram
StabilizerExpectation[ps_PauliStabilizer, pauliString_String]
```

The `pauliString` is `[+/-][I|X|Y|Z]+` (n characters, with optional sign prefix).

## Examples

**Stabilizer-group elements give +1:**

```wolfram
With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    <|
        "<XX>" -> StabilizerExpectation[psBell, "XX"],
        "<ZZ>" -> StabilizerExpectation[psBell, "ZZ"],
        "<YY>" -> StabilizerExpectation[psBell, "YY"]
    |>
]
```
```
<|"<XX>" -> 1, "<ZZ>" -> 1, "<YY>" -> -1|>
```

`⟨XX⟩ = ⟨ZZ⟩ = +1` because both are in the stabilizer group `⟨XX, ZZ⟩`. **`⟨YY⟩ = −1`** because `YY = (iXZ)⊗(iXZ) = i² · XX · ZZ = −1 · XX · ZZ` — the i-factor matters! This case currently uses the direct-vector fallback; the AG-closed-form i-factor tracking is on [ROADMAP §A.2](ROADMAP.md).

**Anticommuting Paulis give 0:**

```wolfram
With[{psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]},
    StabilizerExpectation[psBell, "XI"]
]
```
```
0
```

**5Q stabilizer measurements all give +1:**

```wolfram
With[{ps5 = PauliStabilizer["5QubitCode"]},
    StabilizerExpectation[ps5, #] & /@ {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}
]
```
```
{1, 1, 1, 1}
```

## Messages

- `StabilizerExpectation::dim` — emitted when the Pauli string length doesn't match the qubit count of the input.

---

# GraphState

Phase 5 graph-state representation. `GraphState[<|"Graph" -> g_Graph, "VOPs" -> {0,...,0}|>]`. The stabilizer at vertex `i` is `K_i = X_i ⊗ Π_{j ∈ N(i)} Z_j` (AndBri05 Eq 1).

## Constructors

### From a Graph

VOPs default to identity (all zeros).

```wolfram
With[{gs = GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]},
    <|
        "Vertices"    -> gs["VertexCount"],
        "Edges"       -> gs["EdgeCount"],
        "Stabilizers" -> gs["Stabilizers"]
    |>
]
```
```
<|"Vertices" -> 3, "Edges" -> 2, "Stabilizers" -> {"XZI", "ZXZ", "IZX"}|>
```

The linear cluster `1—2—3` has stabilizers `K_1 = XZI`, `K_2 = ZXZ`, `K_3 = IZX`.

### Cluster-5

```wolfram
GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"]
```
```
{"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"}
```

### From a graph-form PauliStabilizer

`GraphState[ps]` works when `ps` already has graph-form stabilizers (each generator is `X_i` at one position with `Z` at neighbors).

## Properties

| Property | Returns |
|---|---|
| `"Graph"` | the underlying `Graph` |
| `"VOPs"` | vertex-operator indices (Phase 5 v1: all 0 = identity) |
| `"Vertices"` | `VertexList[graph]` |
| `"Edges"` | `EdgeList[graph]` |
| `"VertexCount"` | n |
| `"EdgeCount"` | edge count |
| `"Qubits"` / `"Qudits"` | n |
| `"AdjacencyMatrix"` | dense adjacency matrix |
| `"Stabilizers"` | string list of stabilizer generators (each `K_i`) |
| `"PauliStabilizer"` | converts to a `PauliStabilizer` head |

## Conversion to PauliStabilizer

```wolfram
With[{gs = GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]]},
    gs["PauliStabilizer"]["Stabilizers"]
]
```
```
{"XZI", "ZXZ", "IZX"}
```

## See also

- [`LocalComplement`](#localcomplement) — graph operation that transforms graph states by a local Clifford.
- [ROADMAP §A.6](ROADMAP.md) — VOP tracking under LC (deferred).
- AndBri05 §2 (arxiv:quant-ph/0504117).

---

# LocalComplement

Phase 5 local-complementation operation on a `Graph` or `GraphState`. Toggles all edges among the neighbors of a vertex. Phase 5 v1 does not yet update VOPs ([ROADMAP §A.6](ROADMAP.md)).

## Signatures

```wolfram
LocalComplement[g_Graph, v]
LocalComplement[gs_GraphState, v]
```

## Examples

**LC at the center of a star turns the leaves into a clique:**

```wolfram
With[{g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}]},
    Sort @ EdgeList @ LocalComplement[g, 1]
]
```
```
{UndirectedEdge[1, 2], UndirectedEdge[1, 3], UndirectedEdge[1, 4],
 UndirectedEdge[2, 3], UndirectedEdge[2, 4], UndirectedEdge[3, 4]}
```

The star `1—{2,3,4}` becomes `K_4`-minus-the-{2,3,4}-clique-attached-to-1: original edges `{1-2, 1-3, 1-4}` plus three new edges `{2-3, 2-4, 3-4}`.

**LC is involutive:**

```wolfram
With[{g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]},
    Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
]
```
```
True
```

`LC ∘ LC = id` for any vertex (AndBri05 Definition 1).

## See also

- [`GraphState`](#graphstate) — the graph-state head.
- AndBri05 Theorem 1 (LC preserves entanglement spectrum).
- PatGuh26 §3.2 (arxiv:2312.02377) — graphical rules.
- [ROADMAP §A.6](ROADMAP.md) — VOP tracking.
- [ROADMAP §B.3](ROADMAP.md) — 24-element LocalClifford table.

---

# Re-verification

To re-run all 49 code blocks and confirm the embedded outputs are still correct:

```bash
wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-API.wls
```

Drift between this document's `output` blocks and the verifier's printed output is a regression signal.

---

# Quick reference card

```wolfram
(* Constructor *)
ps = PauliStabilizer["5QubitCode"];           (* named code *)
ps = PauliStabilizer[{"XX", "ZZ"}];           (* string list *)
ps = PauliStabilizer[3];                       (* |0...0> register *)
ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];

(* Random *)
ps = RandomClifford[3];

(* Properties *)
ps["Qubits"]; ps["Stabilizers"]; ps["Tableau"]; ps["Matrix"]; ps["Phase"];

(* Clifford gates *)
ps["H", q]; ps["S", q]; ps["X", q]; ps["CNOT", c, t]; ps["CZ", c, t];

(* Non-Clifford -> StabilizerFrame *)
psT = ps["T", q];   (* StabilizerFrame head *)
psT["StateVector"]; (* materialize *)

(* Measurement *)
ps["M", q];           (* Z-basis: <|0 -> ps0, 1 -> ps1|> *)
ps["M", "XZZXI"];     (* arbitrary Pauli string *)
ps["M", {1, 2, 3}];   (* multi-qubit *)

(* Symbolic measurement (Phase 3) *)
psSym = StabilizerMeasure[ps, q];
ps0 = SubstituteOutcomes[psSym, \[FormalS][1] -> 0];
samples = SampleOutcomes[psSym, 100];

(* Inner product / expectation (Phase 4) *)
StabilizerInnerProduct[psA, psB];
StabilizerExpectation[ps, "XYZIX"];

(* Graph state (Phase 5) *)
gs = GraphState[Graph[Range[5], ...]];
gs["Stabilizers"]; gs["PauliStabilizer"];
LocalComplement[gs, vertex];

(* QuantumFramework integration *)
QuantumCircuitOperator[gates][Method -> "Stabilizer"];
```
