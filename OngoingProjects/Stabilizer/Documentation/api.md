# Stabilizer subsystem — API reference

> Function-by-function reference for the **5 top-level public symbols** plus the method-grade operations dispatched on `PauliStabilizer`, `StabilizerFrame`, and `CliffordChannel`. Reflects the kernel state after Phases 1–8.3 + the WL-idiom audit (2026-05-06). For the full phase log and deferred items, see [roadmap.md](roadmap.md). For the Phase 5c design lesson on global-phase contracts, see [postmortem.md](postmortem.md).

## Public surface at a glance

| Symbol | Phase | Purpose |
|---|---|---|
| [`PauliStabilizer`](#paulistabilizer) | 1 | Tableau-encoded n-qubit stabilizer state |
| [`StabilizerFrame`](#stabilizerframe) | 4 | Superposition `Σᵢ cᵢ \|sᵢ⟩` of stabilizer states (non-Clifford boundary) |
| [`GraphState`](#graphstate) | 5 | Graph-state representation (AndBri05) |
| [`LocalComplement`](#localcomplement) | 5 | Local complementation on a graph or graph state |
| [`CliffordChannel`](#cliffordchannel) | 8 | Choi-tableau Clifford channel `[U_A \| U_B \| c]` (Yashin25 §2.3) |

### Method-grade operations (Phase 6 / 6.5 demoted from top-level)

| Method | Old top-level form (pre-Phase 6) | Purpose |
|---|---|---|
| [`PauliStabilizer["Random", n]`](#named-pattern-random-n) | `RandomClifford[n]` | Uniformly random n-qubit Clifford state (Mallows sampler) |
| [`ps["SymbolicMeasure", q]`](#symbolic-measurement-fangying23) | `StabilizerMeasure[ps, q]` | Symbolic Z-basis measurement |
| [`ps["SubstituteOutcomes", rules]`](#symbolic-measurement-fangying23) | `SubstituteOutcomes[ps, rules]` | Plug concrete values into outcome symbols |
| [`ps["SampleOutcomes", n]`](#symbolic-measurement-fangying23) | `SampleOutcomes[ps, n]` | Random outcome samples |
| [`ps["InnerProduct", other]`](#inner-product-and-expectation) | `StabilizerInnerProduct[ps, other]` | `⟨ps\|other⟩` (PS or Frame on either side) |
| [`ps["Expectation", pauli]`](#inner-product-and-expectation) | `StabilizerExpectation[ps, pauli]` | `⟨ps\|P\|ps⟩` for a Pauli string |

### Verification

The 684-test suite under [`Tests/Stabilizer/`](../../../Tests/Stabilizer/) is the authoritative oracle:

```bash
wolframscript -f Tests/RunTests.wls Stabilizer
```

The breadth-first coverage matrix is [`Tests/Stabilizer/AuditMatrix.wlt`](../../../Tests/Stabilizer/AuditMatrix.wlt) (155 tests, 15 TIERs). Run it first when reviewing a PR that touches the public surface.

---

# PauliStabilizer

The atomic stabilizer-state object. Encodes an n-qubit state as a `(2n × 2n)` symplectic matrix plus a length-`2n` sign vector, packaged as `PauliStabilizer[<|"Tableau" -> rank3binarray, "Signs" -> {±1, ...}|>]`. The optional `"GlobalPhase"` key (Phase 5c) carries the complex unit that the AG decomposition drops, so that constructor → accessor round-trips exactly on the first hop.

## Constructors

### Empty form: `PauliStabilizer[]`

Delegates to `PauliStabilizer[1]` (single qubit `|0⟩`).

### Integer form: `PauliStabilizer[n]` — n-qubit `|0…0⟩` register

```wolfram
PauliStabilizer[3]["Stabilizers"]
(* {"ZII", "IZI", "IIZ"} *)
```

`n = 0` auto-pads to one qubit (`Max[n, 1]` in [Constructors.m:241](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m)).

### List of Pauli strings: `PauliStabilizer[{stab1, stab2, ...}]`

Each string is `[+|-]?[IXYZ]+`. Sign prefix optional.

```wolfram
PauliStabilizer[{"-XX", "+ZZ"}]["StabilizerSigns"]
(* {-1, 1} *)
```

Destabilizers are auto-padded via the Reverse rule ([Constructors.m:203](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m)). The result satisfies stabilizer-stabilizer commutation but **not** the full AG canonical pairing — circuit-built fixtures do.

### Stabilizers + destabilizers: `PauliStabilizer[stabs, destabs]`

Both lists are Pauli strings.

```wolfram
With[{ps = PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]},
    {ps["Stabilizers"], ps["Destabilizers"]}
]
(* {{"XX", "ZZ"}, {"IX", "ZI"}} *)
```

### Named codes: `PauliStabilizer[name_String]`

Available names (see [`$PauliStabilizerNames`](../../../QuantumFramework/Kernel/Stabilizer/NamedCodes.m)):

| Name | Aliases | What |
|---|---|---|
| `"5QubitCode"` | — | `\|0_L⟩` of `[[5,1,3]]` (Got97 §3.5) |
| `"5QubitCode1"` | — | `\|1_L⟩` (last stabilizer flipped) |
| `"7QubitCode"` | `"SteaneCode"` | `\|0_L⟩` of `[[7,1,3]]` (Got00 §4) |
| `"7QubitCode1"` | `"SteaneCode1"` | `\|1_L⟩` |
| `"9QubitCode"` | — | `\|0_L⟩` of `[[9,1,3]]` (Shor) |
| `"9QubitCode1"` | — | `\|1_L⟩` |
| `"Random"` | — | See [Random named pattern](#named-pattern-random-n) below |

### Named pattern: `PauliStabilizer["Random", n]`

Bravyi–Maslov / Koenig–Smolin Mallows sampler ([RandomClifford.m](../../../QuantumFramework/Kernel/Stabilizer/RandomClifford.m); KoeSmo14 §3.2). Cardinality `\|C_n\| = 2^(n²+2n) Π(4ʲ−1)` (24 for n=1, 11520 for n=2, ~9.29×10⁷ for n=3).

```wolfram
SeedRandom[20260507];
PauliStabilizer["Random", 3]["Stabilizers"]
```

### From a `QuantumCircuitOperator`: `PauliStabilizer[qco]`

Folds `qco`'s normal operators over the `|0…0⟩` register. Each gate must be Clifford; non-Clifford gates emit `PauliStabilizer::nonclifford` and the state is returned unchanged.

```wolfram
PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Stabilizers"]
(* {"XX", "ZZ"} *)
```

### From a `QuantumOperator`: `PauliStabilizer[qo, n_Integer : 1]`

Conjugates the qubit-aligned Pauli generators by `qo`; extracts the resulting tableau. The optional second arg pads the result to `Max[n, qo's input order]` qubits. **Captures `"GlobalPhase"`** so that `["QuantumOperator"]` returns the input exactly.

```wolfram
PauliStabilizer[QuantumOperator["Y"]]["GlobalPhase"]
(* -I *)

PauliStabilizer[QuantumOperator["H", 1], 3]["Qubits"]
(* 3  -- pad to 3 qubits *)
```

Requires `Sort[qo["OutputOrder"]] == Sort[qo["InputOrder"]]` (square operator).

### From a `QuantumState`: `PauliStabilizer[qs]`

Runs 4ⁿ Pauli-expectation tomography ([Constructors.m:88](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m)) and picks n linearly-independent generators by greedy `𝔽₂` rank-growth. Practical for `n ≤ ~8`. **Captures `"GlobalPhase"`**.

```wolfram
PauliStabilizer[QuantumState[{0, 1}]]["Stabilizers"]
(* {"-Z"} *)
```

### Association forms

| Form | Effect |
|---|---|
| `PauliStabilizer[<\|"Tableau" -> t, "Signs" -> s\|>]` | Direct construction; signs and tableau set |
| `PauliStabilizer[<\|"Tableau" -> t\|>]` | Default signs `{1, …, 1}` |
| `PauliStabilizer[<\|"Phase" -> ph, "Tableau" -> t\|>]` | Phase form (`Signs = 1 - 2 ph`) |
| `PauliStabilizer[<\|"Matrix" -> m, "Phase" -> ph\|>]` | Symplectic-block layout `[X\|Z]` reshaped into the rank-3 tableau |
| `PauliStabilizer[<\|…, "GlobalPhase" -> z\|>]` | Multiplies `["State"]` and `["QuantumOperator"]` outputs by complex unit `z` |

### Tableau-only positional forms

For interoperating with raw bit arrays (e.g., when porting from another stabilizer library):

| Form | Effect |
|---|---|
| `PauliStabilizer[t_?PauliTableauQ]` | Single rank-3 `{0,1}` tableau, shape `{2, n, m}`; default sign +1, paired with auto-padded destabilizer half via the Reverse rule. |
| `PauliStabilizer[sign : -1\|1, t_?PauliTableauQ]` | One sign for the first row + the tableau. |
| `PauliStabilizer[signs : {(-1\|1) ...}, t_?PauliTableauQ]` | Per-row signs (PadRight'd to the tableau's row count). |
| `PauliStabilizer[basis_]` | Rank-3 array of `{0, 1, -1}` of shape `{2, n, m}`. The sign is read from each row by `Union @@@ #` (must collapse to `{0, 1} \| {1}` → `+1`, `{-1, 0} \| {-1}` → `-1`). Internally calls the signs+tableau form. |

### Shortcut form: `PauliStabilizer["H" -> 1]` (or list / string / `name -> qudit`)

Any `_String | (_String -> _ ? orderQ | _Integer) | _List` shortcut is forwarded to `PauliStabilizer[QuantumCircuitOperator[shortcut]]`. Equivalent to writing the circuit, then taking its tableau.

## Properties

Access via `ps[propName]`. Source-of-truth: [Properties.m](../../../QuantumFramework/Kernel/Stabilizer/Properties.m). The full key list is `_PauliStabilizer["Properties"]` (40 entries, set in [PauliStabilizer.m:55](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m)). If you see a stale short list at runtime, run `Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]]; PacletDirectoryLoad["…/QuantumFramework"]` first to clear an old paclet install.

### Shape

| Property | Returns |
|---|---|
| `"Qubits"` / `"Qudits"` | n |
| `"GeneratorCount"` | n (rank of the stabilizer group) |
| `"Tableau"` | rank-3 binary array, shape `{2, n, 2n}` (X-block, qubit, row) |
| `"Matrix"` | `2n × 2n` binary matrix in `[X \| Z]` block layout |

### Signs / phase

| Property | Returns |
|---|---|
| `"Signs"` | length-2n list of ±1 (or symbolic Phase-3 polynomials) |
| `"Phase"` | `(1 - Signs) / 2` |
| `"StabilizerSigns"` | last n signs |
| `"DestabilizerSigns"` | first n signs |

### Bit views

| Property | Returns |
|---|---|
| `"X"` | shape `{n, 2n}`: rows are qubits, columns are tableau rows (destab first) |
| `"Z"` | analogous Z-bits |
| `"StabilizerX"`, `"StabilizerZ"` | last n columns of `X` / `Z` |
| `"DestabilizerX"`, `"DestabilizerZ"` | first n columns |
| `"p"` | symplectic-product diagonal (length 2n) |
| `"TableauPhase"` | concatenation `[Matrix \| Phase]` |

### Stabilizer / destabilizer slices

| Property | Returns |
|---|---|
| `"Stabilizer"` / `"StabilizerTableau"` | last n columns of Tableau, shape `{2, n, n}` |
| `"Destabilizer"` / `"DestabilizerTableau"` | first n columns of Tableau |

### Pauli string lists

| Property | Returns |
|---|---|
| `"Stabilizers"` / `"PauliForm"` / `"Generators"` | string list of stabilizers (with `-` prefix where appropriate) |
| `"Destabilizers"` | string list of destabilizers |
| `"PauliStrings"` | concatenation `Stabilizers ++ Destabilizers` |
| `"PauliSymbols"` | same but with `±symbol` form (for `MakeBoxes`) |

Optional 2nd arg `n` truncates the list to the first `n` entries: `ps["Stabilizers", 3]`.

### Display

| Property | Returns |
|---|---|
| `"TableauForm"` | `Row[…Grid…]` for the full tableau |
| `"StabilizerTableauForm"` | grid for stabilizer half only |

### Materialization (cost ≥ 2ⁿ)

| Property | Cost | Returns |
|---|---|---|
| `"State"` / `"QuantumState"` | `O(2ⁿ)` | `QuantumState` (multiplied by `"GlobalPhase"`) |
| `"Operator"` / `"QuantumOperator"` | `O(4ⁿ)` | `QuantumOperator` |
| `"Circuit"` / `"QuantumCircuit"` / `"QuantumCircuitOperator"` | `O(n³)` gates | `QuantumCircuitOperator` whose dagger applied to `\|0…0⟩` reproduces ps |

### Phase 5c global-phase contract

| Path | Equality | Notes |
|---|---|---|
| `PauliStabilizer[qo]["QuantumOperator"]` vs `qo` | **exact** (`===`) | `"GlobalPhase"` captured |
| `PauliStabilizer[qs]["State"]` vs `qs` | **exact** (`===`) | `"GlobalPhase"` captured |
| `PauliStabilizer[qs]["gate", q]["State"]` vs `gate[qs]` | **up to phase** | Gate updates do not propagate `"GlobalPhase"` ([ROADMAP §A.9](roadmap.md)) |
| `PauliStabilizer[gate[qs]]["State"]` vs `gate[qs]` | **exact** | Re-tomograph after the gate to recover phase |
| `PauliStabilizer[ps["Circuit"]]["State"]` vs `ps["State"]` | **up to phase** | No source state, only tableau |

The Y row is the canary: `Y = i X Z` → AG-decomposed circuit `Z·X = i·Y`, so the constructor records `"GlobalPhase" -> -I` and `["QuantumOperator"]` returns Y exactly.

## Methods (Clifford gates)

Each gate returns a new `PauliStabilizer`. Source: [GateUpdates.m](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m).

| Method | Signature | Effect |
|---|---|---|
| `"H"` | `ps["H", q]` | Hadamard on qubit q |
| `"S"` | `ps["S", q]` | π/4 phase gate |
| `SuperDagger["S"]` | `ps[SuperDagger["S"], q]` | S-dagger |
| `"X"`, `"Y"`, `"Z"` | `ps["X", q]` etc. | Pauli (only flips signs) |
| `"V"` | `ps["V", q]` | √X |
| `SuperDagger["V"]` | `ps[SuperDagger["V"], q]` | √X-dagger |
| `"CNOT"` / `"CX"` | `ps["CNOT", c, t]` | Controlled-NOT |
| `"CZ"` | `ps["CZ", c, t]` | Controlled-Z (= H_t · CNOT · H_t) |
| `"SWAP"` | `ps["SWAP", a, b]` | Swap (= column permutation) |
| `"Permute"` | `ps["Permute", Cycles[…]]` | Permute tableau rows |
| `"PermuteQudits"` | `ps["PermuteQudits", Cycles[…]]` | Permute tableau columns |
| `"PadRight"` | `ps["PadRight", n]` | Tensor identity to grow to n qubits |
| `"PadLeft"` | `ps["PadLeft", n]` | Tensor identity from the left |
| `"Dagger"` / `"Inverse"` | `ps["Dagger"]` | Tableau inverse over `𝔽₂` |
| `"P"[θ]` | `ps["P"[θ], q]` | Non-Clifford rotation; returns `StabilizerFrame` |
| `"T"` | `ps["T", q]` | Alias for `"P"[π/2]` |
| `SuperDagger["T"]` | `ps[SuperDagger["T"], q]` | Alias for `"P"[-π/2]` |

### Convenience syntax: `op -> order`

```wolfram
ps[gate -> q]      (* same as ps[gate, q] *)
ps[gate -> {c, t}] (* same as ps[gate, c, t] for two-qubit gates *)
```

## Methods (Measurement)

### Implicit measurement shortcuts

These are syntactic sugar that all rewrite to `ps["M", …]` ([PauliStabilizer.m:82-85](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m)):

| Form | Equivalent to |
|---|---|
| `ps[q1, q2, …]` (positional integers) | `ps["M", {q1, q2, …}]` |
| `ps[{q1, q2, …}]` | `ps["M", {q1, q2, …}]` |
| `ps[]` (no args) | `ps["M", Range[ps["Qudits"]]]` (measure every qudit) |

```wolfram
PauliStabilizer[3][1, 2, 3]   (* same as ps["M", {1, 2, 3}] *)
PauliStabilizer[3][]          (* measure all 3 qudits *)
```

### Z-basis: `ps["M", q]` or `ps[q]`

Returns an `Association` keyed by outcome bit (0 or 1) and valued by the post-measurement `PauliStabilizer`.

```wolfram
PauliStabilizer[1]["M", 1]
(* <|0 -> PauliStabilizer[<|"Signs" -> {1, 1}, "Tableau" -> {{{1, 0}}, {{0, 1}}}|>]|> *)
```

Single key = deterministic; two keys = non-deterministic.

### Pauli string: `ps["M", "XZZXI"]`

Measure an arbitrary Pauli observable. Same Association return type. Source: [PauliMeasure.m](../../../QuantumFramework/Kernel/Stabilizer/PauliMeasure.m).

```wolfram
PauliStabilizer["5QubitCode"]["M", "XZZXI"]
(* <|0 -> PauliStabilizer[…]|>  -- deterministic, eigenvalue +1 *)
```

### Multi-qubit: `ps["M", {q1, q2, …}]`

Returns Association keyed by outcome tuples.

```wolfram
PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}}["M", {1, 2}]
(* <|{0, 0} -> ..., {1, 1} -> ...|>  -- Bell ZZ correlation *)
```

### AG primitive: `ps["RowSum", h, i]`

Low-level Aaronson-Gottesman row-multiplication primitive ([Measurement.m:40](../../../QuantumFramework/Kernel/Stabilizer/Measurement.m); AarGot04 §3 Lemma 3). Replaces tableau row `h` with row `i` XOR row `h`, tracking the per-qubit phase via the `agPhase` g-function. Used internally by `["M", a]` measurement updates; exposed for testing or for users who need to manipulate the tableau directly.

```wolfram
PauliStabilizer[2]["H", 1]["RowSum", 4, 3]   (* combine stabilizer rows 3 and 4 *)
```

### Symbolic measurement (FangYing23)

Demoted from top-level `StabilizerMeasure` / `SubstituteOutcomes` / `SampleOutcomes` in Phase 6. Source: [SymbolicMeasure.m](../../../QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m).

#### `ps["SymbolicMeasure", q]` / `ps["SymbolicMeasure", {q1, q2, …}]`

Allocates a fresh `\[FormalS][k]` symbol per non-deterministic measurement; returns a single `PauliStabilizer` instead of an Association. Deterministic measurements get a concrete sign with no fresh symbol.

```wolfram
PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1]["Phase"]
(* {0, \[FormalS][1]}  -- second sign-bit stamped with the symbol *)
```

#### `ps["SubstituteOutcomes", rules]`

Replace outcome symbols with concrete 0/1 values; reduce signs back to `{-1, 1}` via `Mod 2`.

```wolfram
psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
psSym["SubstituteOutcomes", \[FormalS][1] -> 0]["Stabilizers"]
(* {"Z"} *)
```

#### `ps["SampleOutcomes"]` / `ps["SampleOutcomes", n]`

Draw n random samples by independently substituting each outcome symbol with a uniform 0/1.

```wolfram
SeedRandom[42];
psSym = PauliStabilizer[1]["H", 1]["SymbolicMeasure", 1];
#["Stabilizers"] & /@ psSym["SampleOutcomes", 5]
(* {{"Z"}, {"-Z"}, {"-Z"}, {"Z"}, {"-Z"}} *)
```

**Phase 3 known limitation:** when a deterministic measurement follows a prior symbolic one (e.g. Bell ZZ correlation), the deterministic outcome polynomial is computed correctly but is not stamped into the post-state's signs. Tracked in [ROADMAP §A.4](roadmap.md).

## Inner product and expectation

Phase 6.5 demoted from top-level `StabilizerInnerProduct` / `StabilizerExpectation`. Source: [InnerProduct.m](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m).

### `ps["InnerProduct", other]`

`other` is a `PauliStabilizer` or a `StabilizerFrame`. Phase 4 v1 uses direct vector materialization (`O(2ⁿ)`); the closed-form `O(n³)` GarMarCro12 algorithm is on [ROADMAP §A.1](roadmap.md).

```wolfram
psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
psBell["InnerProduct", psBell]
(* 1 *)

psBell["InnerProduct", psBell["Z", 1]]
(* 0  -- |Φ+>, |Φ-> are orthogonal *)
```

### `ps["Expectation", pauliString]`

`⟨ps\|P\|ps⟩` for a Pauli string `P`. Returns ±1 for stabilizer-group elements, 0 for anticommuting Paulis, exact value via direct vector for `P ∈ N(S) \ S`.

```wolfram
psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
{psBell["Expectation", "XX"], psBell["Expectation", "ZZ"], psBell["Expectation", "YY"]}
(* {1, 1, -1}  -- YY = -XX·ZZ because Y = iXZ *)
```

The `YY = -1` case currently uses the direct-vector fallback because the AG closed-form `i`-factor tracking ([ROADMAP §A.2](roadmap.md)) is not yet implemented.

Message: `PauliStabilizer::expectationdim` if string length doesn't match qubit count.

## Composition and tensor product

### `ps1[ps2]` — apply ps1 after ps2

Sequential composition. Source: [Compose.m](../../../QuantumFramework/Kernel/Stabilizer/Compose.m).

### `QuantumTensorProduct[ps1, ps2]`

Block-diagonal merge of destabilizers and stabilizers. Linear in `n_a + n_b`.

```wolfram
With[{a = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}},
      b = PauliStabilizer[1]},
    QuantumTensorProduct[a, b]["Stabilizers"]
]
(* {"XXI", "ZZI", "IIZ"} *)
```

## Connections to other QuantumFramework heads

PauliStabilizer participates in 9 cross-head dispatch points:

| Pattern | Type | Defined in | Effect |
|---|---|---|---|
| `QuantumState[ps]` | UpValue | [Conversions.m:145](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) | Same as `ps["State"]` |
| `QuantumOperator[ps]` | DownValue | [Conversions.m:147](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) | Same as `ps["Circuit"]["QuantumOperator"]` |
| `QuantumCircuitOperator[ps]` | DownValue | [Conversions.m:146](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) | Same as `ps["Circuit"]` |
| `qo_QuantumOperator[ps]` | UpValue | [Conversions.m:144](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) | `PauliStabilizerApply[QuantumCircuitOperator[qo], ps]` |
| `qmo_QuantumMeasurementOperator[ps]` | UpValue | [HybridInterop.m:193](../../../QuantumFramework/Kernel/Stabilizer/HybridInterop.m) | Pauli fast path or `::nonpaulibasis` fallback (see [Hybrid interop](#hybrid-interop-phase-7)) |
| `qc_QuantumChannel[ps]` | UpValue | [HybridInterop.m:260](../../../QuantumFramework/Kernel/Stabilizer/HybridInterop.m) | Tableau-mixture for named Pauli channels, dense fallback otherwise |
| `cc_CliffordChannel[ps]` | DownValue | [CliffordChannel.m:451](../../../QuantumFramework/Kernel/Stabilizer/CliffordChannel.m) | Identity / state-prep / dim-matched composition (see [CliffordChannel](#cliffordchannel)) |
| `QuantumTensorProduct[ps_a, ps_b]` | DownValue | [Compose.m:57](../../../QuantumFramework/Kernel/Stabilizer/Compose.m) | Block-diagonal tensor merge |
| `QuantumCircuitOperator[…][Method -> "Stabilizer"]` | Method dispatch | (host paclet) | Routes through `PauliStabilizerApply` |

```wolfram
QuantumState[ps]                                   (* via UpValue *)
QuantumCircuitOperator[…][Method -> "Stabilizer"]  (* via Method *)
QuantumOperator["H", 1] @ ps                       (* via UpValue *)
```

## Internals (PackageScope)

These are not part of the public surface but are referenced from tests, documentation, or the design rationale. Access via `Wolfram``QuantumFramework``PackageScope``…`.

| Symbol | Role |
|---|---|
| `PauliStabilizerQ[expr]` | Predicate. Matches `PauliStabilizer[<\|"Signs" -> List, "Tableau" -> rank-3 array of ±1/0\|>]` with consistent dimensions. Phase 3 loosened to accept symbolic signs. |
| `ConcretePauliStabilizerQ[expr]` | Strict predicate: signs are concrete `(-1\|1)` integers (no symbolic phase). Used to gate the HybridInterop fast paths and `cc[ps]` evolution. |
| `PauliTableauQ[t]` | Predicate. Matches a rank-3 binary array of shape `{2, n, m}` with `0`/`1`/`-1` entries. |
| `PauliStabilizerApply[qco, qs : Automatic \| _QuantumState \| _PauliStabilizer : Automatic]` | The folder over a `QuantumCircuitOperator`'s normal operators that applies each gate via `state[gate]` and bails to `Plus` / `::nonclifford` on non-Clifford gates. Entry point for `Method -> "Stabilizer"`. |
| `FromFullTableau[t]` | Construct a `PauliStabilizer` from a `(2n) × (2n+1)` AG-format matrix `[Matrix \| Phase]`. |
| `PauliRow[mat, n, d]` | Decode a Pauli matrix back to `{xbits, zbits, phase}` (used by `QuantumOperatorTableau`). |
| `QuantumOperatorTableau[qo]` | Build the AG tableau from a `QuantumOperator` by conjugating the qubit-aligned `X` and `Z` generators. |
| `agPhase[x1, z1, x2, z2]` | AG g-function (AarGot04 §3 Lemma 3). The integer i-power of `P(src) · P(dst)` per qubit. |
| `rowsum[{r, t}, h, i]` | Row-multiplication primitive used by `["M", a]` and `["RowSum", h, i]`. |

## Messages

| Message | Fired when |
|---|---|
| `PauliStabilizer::nonclifford` | A circuit contains a gate with no tableau-update rule and that isn't `P[θ] / T / T†`. State returned unchanged. |
| `PauliStabilizer::tdeprecated` | Historical: announces the Phase 4 migration of `P[θ] / T / T†` from `Plus` to `StabilizerFrame`. |
| `PauliStabilizer::expectationdim` | `ps["Expectation", P]` with `Length[P] ≠ Qubits`. |
| `PauliStabilizer::nonpaulibasis` | A non-Pauli QMO or non-Pauli channel acts on a stabilizer; falls back to legacy circuit path. |

---

# StabilizerFrame

Phase 4 superposition of stabilizer states: `Σᵢ cᵢ \|sᵢ⟩` with possibly symbolic coefficients. Produced by non-Clifford gates (`P[θ]`, `T`, `T†`) on a `PauliStabilizer`. Closes under further Clifford gates (which distribute over components) and under further non-Clifford gates (which double the component count). Source: [StabilizerFrame.m](../../../QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m).

Predicate (PackageScope): `Wolfram``QuantumFramework``PackageScope``StabilizerFrameQ`.

## Constructors

### From a list of `{coefficient, PauliStabilizer}` pairs

```wolfram
StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]
```

### From a single PauliStabilizer (coefficient 1)

```wolfram
StabilizerFrame[PauliStabilizer[1]]
```

### From a non-Clifford gate

```wolfram
PauliStabilizer[1]["T", 1]
(* StabilizerFrame[…] of length 2 *)
```

### From an Association

```wolfram
StabilizerFrame[<|"Components" -> {{c1, ps1}, {c2, ps2}, …}|>]
```

## Properties

| Property | Returns |
|---|---|
| `"Components"` | List of `{coefficient, PauliStabilizer}` pairs |
| `"Coefficients"` | List of coefficients only |
| `"Stabilizers"` | List of underlying `PauliStabilizer`s |
| `"Length"` | Number of components |
| `"Qubits"` / `"Qudits"` | Qubit count of components (assumed identical) |
| `"GeneratorCount"` | Generator count of components |
| `"StateVector"` | Materialized state vector (`O(2ⁿ)`) |
| `"State"` | Materialized `QuantumState` |

## Methods

### Clifford gates distribute over components

```wolfram
With[{f = StabilizerFrame[{{1, PauliStabilizer[1]}, {1, PauliStabilizer[1]["H", 1]}}]},
    f["H", 1]["Length"]
]
(* 2  -- Same component count, each had H applied *)
```

### Non-Clifford gates double the count

```wolfram
PauliStabilizer[1]["H", 1]["T", 1]["T", 1]["Length"]
(* 4 *)
```

### Inner product (Phase 6.5)

```wolfram
frame["InnerProduct", other]
(* same as ps["InnerProduct", other] *)
```

## Arithmetic UpValues

| Form | Effect |
|---|---|
| `frame1 + frame2` | Concatenate component lists |
| `c * frame` (scalar c) | Scale all coefficients |
| `frame1 == frame2` | Structural equality on components |

## Materialization

`T\|0⟩ = \|0⟩` (eigenstate). The frame materializes correctly:

```wolfram
Module[{psT = PauliStabilizer[1]["T", 1], vec},
    vec = Normal @ psT["StateVector"];
    Chop @ N @ FullSimplify[vec - {1, 0}]
]
(* {0, 0} *)
```

---

# GraphState

Phase 5 graph-state representation (AndBri05 §2). Each vertex `i` carries the stabilizer `K_i = X_i ⊗ Π_{j ∈ N(i)} Z_j`. Source: [GraphState.m](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m).

Predicate (PackageScope): `Wolfram``QuantumFramework``PackageScope``GraphStateQ`.

## Constructors

### From a `Graph`

VOPs default to identity (all zeros).

```wolfram
gs = GraphState[Graph[Range[3], {1 \[UndirectedEdge] 2, 2 \[UndirectedEdge] 3}]];
gs["Stabilizers"]
(* {"XZI", "ZXZ", "IZX"} *)
```

### From a graph-form `PauliStabilizer`

```wolfram
GraphState[PauliStabilizer[{"XZI", "ZXZ", "IZX"}]]
(* GraphState[<|"Graph" -> linear cluster, "VOPs" -> {0, 0, 0}|>] *)
```

Phase 5 v1 only handles graph-form input (each generator is `X_i` at one position with `Z` at neighbors).

### Association form

```wolfram
GraphState[<|"Graph" -> g, "VOPs" -> {0, …, 0}|>]
```

VOP indices > 0 are reserved for the 24-element `LocalCliffordGroup` (deferred — [ROADMAP §B.3](roadmap.md)).

## Properties

| Property | Returns |
|---|---|
| `"Graph"` | Underlying `Graph` |
| `"VOPs"` | Vertex-operator indices (Phase 5 v1: all 0 = identity) |
| `"Vertices"` | `VertexList[graph]` |
| `"Edges"` | `EdgeList[graph]` |
| `"VertexCount"` | n |
| `"EdgeCount"` | edge count |
| `"Qubits"` / `"Qudits"` | n |
| `"AdjacencyMatrix"` | Dense adjacency matrix |
| `"Stabilizers"` | String list of stabilizer generators |
| `"PauliStabilizer"` | Convert to `PauliStabilizer` head |

## Examples

### Linear cluster on 5 vertices

```wolfram
GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]]["Stabilizers"]
(* {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"} *)
```

### Round-trip with cluster-state circuit

```wolfram
g = Graph[Range[4], Table[i \[UndirectedEdge] (i + 1), {i, 3}]];
ps = GraphState[g]["PauliStabilizer"];
psFromCirc = PauliStabilizer[QuantumCircuitOperator[
    Join[Table["H" -> i, {i, 4}], Table["CZ" -> {i, i + 1}, {i, 3}]]
]];
ps["Stabilizers"] === psFromCirc["Stabilizers"]
(* True *)
```

---

# LocalComplement

Phase 5 local-complementation operation (AndBri05 Definition 1). Toggles all edges among the neighbors of a vertex. Source: [GraphState.m:135](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m).

## Signatures

```wolfram
LocalComplement[g_Graph, v]
LocalComplement[gs_GraphState, v]
```

## Examples

### Star → K₄ minus the star

```wolfram
g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}];
Sort @ EdgeList @ LocalComplement[g, 1]
(* All 6 edges of K_4 *)
```

### Involutivity

```wolfram
g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]];
Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
(* True *)
```

Phase 5 v1 returns the new graph but does **not** update VOPs (deferred to [ROADMAP §A.6](roadmap.md)).

---

# CliffordChannel

Phase 8 Choi-tableau representation of an arbitrary Clifford channel `Φ_{A→B}` per Yashin25 §2.3 (arxiv:2504.14101). Encoded as `[U_A \| U_B \| c]` where:

- `U_A` is a `k × 2|A|` bit matrix on the input system (or `{}` if `|A| = 0`)
- `U_B` is a `k × 2|B|` bit matrix on the output system
- `c` is a length-`k` bit vector (signs)
- `k` ≤ rank of the Choi state

Each row is a Pauli superoperator `π(u_A \| u_B \| c)[ρ_A] = (-1)^c · 2^|A| · Tr[ρ_A P(u_A)] · P(u_B)`, and the channel is `Φ[ρ] = (1 / 2^{|A|+|B|}) Σᵢ π(rowᵢ)[ρ]`.

Source: [CliffordChannel.m](../../../QuantumFramework/Kernel/Stabilizer/CliffordChannel.m).

Predicate (PackageScope): `Wolfram``QuantumFramework``PackageScope``CliffordChannelQ`.

## Constructors

### Idempotent: `CliffordChannel[cc]`

If the argument is already a `CliffordChannel`, returns it unchanged. Used to make wrapper code idempotent.

### Identity channel: `CliffordChannel["Identity", n]`

```wolfram
ccI = CliffordChannel["Identity", 1];
{ccI["InputQubits"], ccI["OutputQubits"], ccI["Rank"]}
(* {1, 1, 2} *)
```

### From a `PauliStabilizer`: `CliffordChannel[ps]`

State-preparation channel (`|A| = 0`, rank = `|B|`). The Choi tableau is the state's stabilizer half.

```wolfram
With[{cc = CliffordChannel[PauliStabilizer[2]]},
    {cc["InputQubits"], cc["OutputQubits"], cc["Rank"]}
]
(* {0, 2, 2} *)
```

### From a `QuantumChannel`: `CliffordChannel[qc]`

Phase 8.2 v1 detects deterministic single-Pauli channels by Label (e.g., `"X"`, `"-XX"`); stochastic Pauli channels emit `CliffordChannel::stochastic`.

### Association form

```wolfram
CliffordChannel[<|
    "UA"           -> uA,    (* k x 2|A| bit matrix or {} *)
    "UB"           -> uB,    (* k x 2|B| bit matrix *)
    "c"            -> c,     (* length-k bit vector *)
    "InputQubits"  -> nA,
    "OutputQubits" -> nB,
    "Source"       -> tag    (* informational *)
|>]
```

## Properties

| Property | Returns |
|---|---|
| `"UA"`, `"UB"`, `"c"` | The three components of the Choi tableau |
| `"InputQubits"` | nA |
| `"OutputQubits"` | nB |
| `"Rank"` | k = `Length[c]` |
| `"Tableau"` | Concatenation `[UA \| UB \| c]` |
| `"Source"` | Provenance tag (`"Identity"`, `"PauliStabilizer"`, `"Composition"`, …) |

Predicate (PackageScope): `Wolfram``QuantumFramework``PackageScope``CliffordChannelQ`.

## Methods

### Composition: `cc1[cc2]`

Apply `cc2` first, then `cc1`. The composition algorithm (Phase 8.2/8.3, Yashin25 §3.2/§3.3):

1. Stack `cc2.UB` (k1 rows) above `cc1.UA` (k2 rows) into a `(k1+k2) × 2nB` matrix.
2. Compute the **left null space** `λ s.t. λ · stack = 0` via `NullSpace[Transpose[stack], Modulus -> 2]`.
3. For each kernel vector `λ = (λ_A, λ_C)`, the composition row is:
   - `u_A_new = λ_A · cc2.UA mod 2`
   - `u_C_new = λ_C · cc1.UB mod 2`
   - `c_new = (λ · cConcat) ⊕ ((rowSumPhase + contractionPhase) / 2)`
4. Drop trivial rows; deduplicate.

The **rowSumPhase** tracks the cumulative AG g-function (`agPhase`) phase across F₂-summing rows on the A, both B, and C sides; the **contractionPhase** adds `(-1)^{Σ_q x_q z_q}` for the combined u_B Pauli (Y_B contributes -1 since `Y^T = -Y`). Without these, the simple c-bit XOR gives the wrong sign for Y-bearing combined u_B (e.g., `cc_S² = cc_Z` requires the contraction sign for the X→-X row).

```wolfram
ccS = CliffordChannel[<|"UA" -> {{1, 0}, {0, 1}}, "UB" -> {{1, 1}, {0, 1}},
                        "c" -> {0, 0}, "InputQubits" -> 1, "OutputQubits" -> 1,
                        "Source" -> "S"|>];
ccS3 = ccS[ccS[ccS]];
{ccS3["UA"], ccS3["UB"], ccS3["c"]}
(* {{{0, 1}, {1, 0}}, {{0, 1}, {1, 1}}, {0, 1}}  -- Z->Z (c=0), X->-Y (c=1) *)
```

### State evolution: `cc[ps]`

Apply the channel to a `PauliStabilizer` state. Three recognized cases:

- **Identity channel** (`Source -> "Identity"`, dimensions match): return ps unchanged.
- **State-prep channel** (`nA == 0`): return the state encoded by cc.
- **Dim-matched channel**: build `CliffordChannel[ps]`, compose, convert back to `PauliStabilizer`.

```wolfram
ccS2 = ccS[ccS];   (* = cc_Z *)
ccS2[PauliStabilizer[1]["X", 1]]["Stabilizers"]
(* {"-Z"}  -- Z|1> = -|1>, with Phase 8.3 phase tracking flowing through *)
```

Falls back to dense materialization with `CliffordChannel::stateevol` for unrecognized cases.

## Helpers (PackageScope)

| Helper | Purpose |
|---|---|
| `stabilizerRowSumAGPhase[U, lambda, n]` | Cumulative AG i-power mod 4 for F₂-summing Pauli rows |
| `stabilizerContractionPhase[uB, n]` | `(-1)^{Σ_q x_q z_q}` for the |Φ⁺⟩_BB' contraction sign |
| `cliffordChannelToPauliStabilizer[cc]` | State-channel → PauliStabilizer (uses string-list constructor) |

## Messages

| Message | Fired when |
|---|---|
| `CliffordChannel::dimmismatch` | Composition dimension mismatch (`cc2.OutputQubits ≠ cc1.InputQubits`) |
| `CliffordChannel::stochastic` | `CliffordChannel[qc]` for a stochastic channel (rank > 1 Kraus). Phase 8.2 only handles deterministic single-Pauli channels. |
| `CliffordChannel::stateevol` | `cc[ps]` for an unrecognized case |

---

# Hybrid interop (Phase 7)

Cross-head dispatch so that `QuantumMeasurementOperator` / `QuantumChannel` consume a `PauliStabilizer` or `StabilizerFrame` natively, without forcing the caller to materialize via `ps["State"]` (which costs `O(2ⁿ)`). Source: [HybridInterop.m](../../../QuantumFramework/Kernel/Stabilizer/HybridInterop.m).

## `qmo[ps_PauliStabilizer]`

UpValue on `PauliStabilizer`. Dispatches by the QMO's operator label.

| Label form | Path | Cost |
|---|---|---|
| Pauli string `[+-]?[IXYZ]+` | `ps["M", pauli]` AG fast path | `O(n²)` |
| `Times[-1, Superscript[X\|Y\|Z\|I, CircleTimes[m]]]` (Phase 7.3) | `"-XXXX…"` then AG fast path | `O(n²)` |
| `Superscript[X\|Y\|Z\|I, CircleTimes[m]]` (Phase 7.3) | `"XXXX…"` then AG fast path | `O(n²)` |
| `Times[-1, str]` for Pauli string `str` (Phase 7.3) | `"-" <> str` then AG fast path | `O(n²)` |
| QMO matrix matches a Pauli string for `n ≤ 4` (Phase 7.4) | `stabilizerPauliFromMatrix`, then AG fast path | `O(4ⁿ)` for detection, `O(n²)` after |
| Other (e.g. computational basis) | `PauliStabilizer::nonpaulibasis` info, fall back to `PauliStabilizerApply[QuantumCircuitOperator[qmo], ps]` | dense |

```wolfram
psBell = PauliStabilizer @ QuantumCircuitOperator @ {"H" -> 1, "CNOT" -> {1, 2}};
QuantumMeasurementOperator["XX", {1, 2}][psBell]
(* <|0 -> ...|>  -- routed through ps["M", "XX"], no fallback *)

QuantumMeasurementOperator[QuantumOperator[KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]], {1, 2}][psBell]
(* Phase 7.4 matrix-iteration: detects "XZ", routes through fast path *)
```

The Phase 7.4 cap is `Wolfram``QuantumFramework``PackageScope``$stabilizerPauliMatrixSearchMaxQubits` (default 4); override via `Block`.

## `qmo[sf_StabilizerFrame]`

UpValue on `StabilizerFrame`. Currently always materializes the frame (`sf["State"]`) and applies the QMO on the materialized state. Frame-native fast path is [Phase 7.5 deferred](roadmap.md).

## `qc[ps_PauliStabilizer]`

UpValue on `PauliStabilizer`. Detects the four named single-qubit Pauli channels (`BitFlip`, `PhaseFlip`, `BitPhaseFlip`, `Depolarizing`) and returns a probabilistic-mixture list `{{prob, ps_after_pauli}, …}`:

| Channel | Branches |
|---|---|
| `BitFlip[p]` | `{(1-p, ps), (p, ps[X, q])}` |
| `PhaseFlip[p]` | `{(1-p, ps), (p, ps[Z, q])}` |
| `BitPhaseFlip[p]` | `{(1-p, ps), (p, ps[Y, q])}` |
| `Depolarizing[p]` | `{(1-3p/4, ps), (p/4, ps[X, q]), (p/4, ps[Y, q]), (p/4, ps[Z, q])}` |

Cost: `O(n)` per branch.

Non-Clifford channels (`AmplitudeDamping`, `PhaseDamping`, `GeneralizedAmplitudeDamping`, `ResetError`, multi-qubit named channels) emit `PauliStabilizer::nonpaulibasis` and fall back to dense materialization (`qc[ps["State"]]`).

```wolfram
QuantumChannel["BitFlip"[1/3], {1}][PauliStabilizer[1]]
(* {{2/3, PauliStabilizer[…]}, {1/3, PauliStabilizer[…]}} *)
```

## `qc[sf_StabilizerFrame]`

Same fallback as `qmo[sf]` — currently materializes; deferred to a frame-native path.

## Helpers (PackageScope)

| Symbol | Role |
|---|---|
| `stabilizerPauliLabelFromQMO[qmo]` | Map a QMO's operator label to a Pauli string (with optional `-` prefix) or `Missing[…]`. Implements the four label-form detectors + the matrix-iteration fallback. |
| `stabilizerPauliFromMatrix[mat, n]` | Phase 7.4 matrix-iteration detector. Iterates `Tuples[{"I","X","Y","Z"}, n]` against `mat` and returns the first Pauli string match (with sign), or `Missing["TooManyQubits"]` if `n > $stabilizerPauliMatrixSearchMaxQubits`. |
| `stabilizerPauliMatrixFromString[str]` | Internal helper: build the `2ⁿ × 2ⁿ` matrix for a (signed) Pauli string, handling the n=1 single-qubit case. |
| `stabilizerCliffordChannelMixture[qc]` | Map a `QuantumChannel` Label (e.g. `"BitFlip"[p]`) to a list of `{prob, paulistring, qubit}` triples, or `Missing["MultiQubitChannel"\|"NotClifford"]`. |
| `$stabilizerPauliMatrixSearchMaxQubits` | Cap on Phase 7.4 matrix iteration. Default 4 (= 512 candidates × 4ⁿ size). User-tunable via `Block`. |

---

# Cross-package fixtures

## Stim

Live fixtures in [`Tests/Stabilizer/fixtures/stim_fixtures.json`](../../../Tests/Stabilizer/fixtures/stim_fixtures.json), generated by [`generate_stim_fixtures.py`](../../../Tests/Stabilizer/fixtures/generate_stim_fixtures.py). Tested in [`Tests/Stabilizer/CrossPackage_Stim.wlt`](../../../Tests/Stabilizer/CrossPackage_Stim.wlt) (44 tests).

Stim uses `_` for I and explicit `+` / `-` signs (e.g., `"+_XX"`); the WLT files normalize via:

```wolfram
stimNormalize[s_String] := StringReplace[
    StringReplace[s, "_" -> "I"],
    StartOfString ~~ "+" -> ""
]
```

## QuantumClifford.jl

Hand-coded canonical results (no live import). Source repo at [OngoingProjects/Stabilizer/External Packages/QuantumClifford.jl](../../External%20Packages/QuantumClifford.jl/). Tested in [`Tests/Stabilizer/CrossPackage_QuantumClifford.wlt`](../../../Tests/Stabilizer/CrossPackage_QuantumClifford.wlt) (15 tests).

Key cross-validated identity: AG g-function per qubit (`agPhase(0, 1, 1, 0) = 1`), so `Z*X = i Y` and three-qubit `ZZZ * XXX = i^3 YYY = -i YYY`, matching QC.jl `prodphase`.

---

# Quick reference card

```wolfram
(* Construction *)
ps = PauliStabilizer[3];                                       (* |000> *)
ps = PauliStabilizer["5QubitCode"];                            (* named code *)
ps = PauliStabilizer[{"XX", "ZZ"}];                            (* string list *)
ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];

(* Random *)
ps = PauliStabilizer["Random", 3];

(* Properties *)
ps["Qubits"]; ps["Stabilizers"]; ps["Tableau"]; ps["Matrix"]; ps["Phase"];
ps["State"]; ps["QuantumOperator"]; ps["Circuit"];

(* Clifford gates *)
ps["H", q]; ps["S", q]; ps[SuperDagger["S"], q];
ps["X", q]; ps["Y", q]; ps["Z", q]; ps["V", q]; ps[SuperDagger["V"], q];
ps["CNOT", c, t]; ps["CZ", c, t]; ps["SWAP", a, b];
ps["Permute", Cycles[…]]; ps["PadRight", n]; ps["Dagger"];

(* Non-Clifford (returns StabilizerFrame) *)
ps["T", q]; ps["P"[\[Theta]], q];

(* Measurement *)
ps["M", q];                  (* Z basis: <|0 -> ps0, 1 -> ps1|> *)
ps["M", "XZZXI"];            (* arbitrary Pauli string *)
ps["M", {1, 2, 3}];          (* multi-qubit *)

(* Symbolic measurement (Phase 3) *)
psSym = ps["SymbolicMeasure", q];
ps0 = psSym["SubstituteOutcomes", \[FormalS][1] -> 0];
samples = psSym["SampleOutcomes", 100];

(* Inner product / expectation *)
ps["InnerProduct", other];
ps["Expectation", "XYZIX"];

(* Graph state *)
gs = GraphState[Graph[Range[5], …]];
gs["Stabilizers"]; gs["PauliStabilizer"];
LocalComplement[gs, vertex];

(* Clifford channel (Phase 8) *)
ccI = CliffordChannel["Identity", n];
ccS = CliffordChannel[<|"UA" -> …, "UB" -> …, "c" -> …,
                        "InputQubits" -> 1, "OutputQubits" -> 1|>];
ccI[ccS];                    (* compose *)
ccS[ps];                     (* state evolution *)

(* Hybrid interop (Phase 7) *)
QuantumMeasurementOperator["ZZ", {1, 2}][ps];                  (* AG fast path *)
QuantumChannel["BitFlip"[p], {1}][ps];                         (* tableau-mixture *)

(* QuantumFramework integration *)
QuantumCircuitOperator[gates][Method -> "Stabilizer"];

(* PackageScope tunables *)
Wolfram`QuantumFramework`PackageScope`$stabilizerPauliMatrixSearchMaxQubits  (* default 4 *)
```

---

# References

| Tag | Paper | Used in |
|---|---|---|
| AarGot04 | Aaronson, Gottesman, "Improved simulation of stabilizer circuits" (arxiv:quant-ph/0406196) | Tableau, gate updates, measurement (`agPhase`, `rowsum`) |
| AndBri05 | Anders, Briegel, "Fast simulation of stabilizer circuits using a graph-state representation" (arxiv:quant-ph/0504117) | GraphState, LocalComplement |
| FangYing23 | Fang, Ying, "SymPhase: symbolic phase representation for stabilizer simulation" (arxiv:2311.03906) | SymbolicMeasure / SubstituteOutcomes / SampleOutcomes |
| GarMar15 | García, Markov, "Simulation of quantum circuits via stabilizer frames" (arxiv:1712.03554) | StabilizerFrame |
| GarMarCro12 | García, Markov, Cross, "Efficient inner-product algorithm for stabilizer states" (arxiv:1210.6646) | InnerProduct (closed-form deferred to ROADMAP §A.1) |
| KoeSmo14 | Koenig, Smolin, "How to efficiently select an arbitrary Clifford group element" (arxiv:1406.2170) | RandomClifford (Mallows sampler) |
| Yashin25 | Yashin, "Choi-tableau formalism for Clifford channels" (arxiv:2504.14101) | CliffordChannel (composition + phase tracking) |
| Got97 | Gottesman, "Stabilizer codes and quantum error correction" (PhD thesis) | 5-qubit code |
| Got00 | Gottesman, "An introduction to quantum error correction and fault-tolerant computation" | Steane / Shor codes |
