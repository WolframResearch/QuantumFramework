# Synthesis → Implementation: working WL code for every section of `package-design-synthesis.md`

> **Companion document.** This walks through `package-design-synthesis.md` §1–§11 and shows the Wolfram Language code that implements each capability in the QuantumFramework Stabilizer subsystem (Phases 1–5). Every code block is verified by `wolframscript`; the output shown is captured directly from a verifier run.

## Provenance

| Item | Value |
|---|---|
| Synthesis source | [`OngoingProjects/Stabilizer/package-design-synthesis.md`](../package-design-synthesis.md) (28 papers) |
| Kernel | [`QuantumFramework/Kernel/Stabilizer/`](../../../QuantumFramework/Kernel/Stabilizer/) (15 files, ~1500 LOC) |
| Verifier | [`verify-synthesis-implementation.wls`](verify-synthesis-implementation.wls) (re-runs every block) |
| Branch | `stabilizer-phases-1-4` |
| Test suite | [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) — 8 tiers, 185 tests, 100% passing |
| Run command | `wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-synthesis-implementation.wls` |
| Generated | 2026-04-30 (moved to OngoingProjects/Stabilizer/ on 2026-05-02) |

## Coverage status

| § | Capability | Phase | Symbol(s) | Status |
|---|---|---|---|---|
| 1.1 | Pauli string (symplectic) | 1 | `PauliStabilizer`, internal `PauliRow` | ✅ |
| 1.2 | Tableau (with destabilizers) | 1 | `PauliStabilizer` | ✅ |
| 1.3 | Graph state | 5 | `GraphState`, `LocalComplement` | ✅ |
| 1.4 | Quadratic-form triple | — | — | ⏸ deferred (DehMoo03 Theorem 3) |
| 1.5 | Stabilizer frame | 4 | `StabilizerFrame` | ✅ |
| 1.6 | Clifford channel (Choi) | — | — | ⏸ deferred (Yashin25) |
| 2.1 | Clifford gate updates | 1 | `ps[gate, q]` | ✅ |
| 2.2 | 24-element local Clifford group | — | — | ⏸ deferred (AndBri05 §2 footnote) |
| 2.3 | Local complementation | 5 | `LocalComplement` | ✅ (no VOP tracking yet) |
| 2.4 | Z-basis & Pauli-string measurement | 1, 4 | `ps["M", q]`, `ps["M", "XZZXI"]` | ✅ |
| 2.5 | Inner products & expectation | 4 | `StabilizerInnerProduct`, `StabilizerExpectation` | ✅ (direct vector; closed-form TODO) |
| 2.6 | Distance / nearest neighbors | — | — | ⏸ deferred (GarMarCro12 §5) |
| 2.7 | Counting / enumeration | 1 | (formulas at test-fixture level) | ✅ |
| 2.8 | Random Clifford | 2 | `RandomClifford` | ✅ |
| 2.9 | Canonical forms | 1 | `ps["Circuit"]` (AG only) | partial |
| 2.10 | Synthesis from a tableau | 1 | `ps["Circuit"]` | partial (no Reid24/Winderl23) |
| 2.11 | Pauli tracking | — | — | ⏸ deferred (Paler14/RuhDev25) |
| 2.12 | QEC code extraction | 1 | named codes + syndromes | ✅ |
| 3.1 | Symbolic Pauli arithmetic | 3 | loosened predicates | ✅ |
| 3.2 | Symbolic measurement | 3 | `StabilizerMeasure`, `SubstituteOutcomes`, `SampleOutcomes` | ✅ |
| 3.3 | Symbolic Clifford parameters | — | — | ⏸ deferred (Mueller26) |
| 3.4 | Symbolic dynamical Lie algebras | — | — | ⏸ deferred |
| 3.5 | Symbolic interconversion | 1, 4 | partial | partial |
| 3.6 | Symbolic qudits | — | — | ⏸ v2 |
| 4-6 | User-facing menus | — | (tables) | mostly ✅ |
| 7-11 | Discussion / architecture / priorities | — | — | (narrative) |

**Legend:** ✅ working, partial, ⏸ deferred (with paper anchor).

---

## §1 — Core data structures

### §1.1 — Pauli string

> **Synthesis** (`package-design-synthesis.md:46-58`): "The atomic object. Three coexisting representations, all interconvertible: symbolic tensor; symplectic bit-vector pair `(x, z) ∈ 𝔽₂ⁿ × 𝔽₂ⁿ` with `(0,0)=I, (1,0)=X, (0,1)=Z, (1,1)=Y`; phase-tracked symplectic `(ε, δ, x, z)`."

**QF kernel** (Phase 1, [Stabilizer/PauliStabilizer.m](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m), [Constructors.m](../../../QuantumFramework/Kernel/Stabilizer/Constructors.m), [Conversions.m:`PauliRow`](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m)). The `PauliStabilizer[<|"Tableau" -> ..., "Signs" -> ...|>]` head stores symplectic bits in the rank-3 `Tableau` array (shape `{2, n_qubits, 2*GeneratorCount}` — first axis splits X/Z) and phase via `Signs ∈ {-1, +1}` or `Phase = (1 - Signs)/2 ∈ {0, 1}`.

**Code:**
```wolfram
ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
Association[
    "Stabilizers"  -> ps["Stabilizers"],
    "X-bits-shape" -> Dimensions[ps["X"]],
    "X-bits"       -> ps["X"],
    "Z-bits"       -> ps["Z"],
    "Phase"        -> ps["Phase"]
]
```

**Verified output** (block `1.1-symplectic-encoding`):
```
<|"Stabilizers" -> {"XX", "ZZ"},
  "X-bits-shape" -> {2, 4},
  "X-bits" -> {{0, 0, 1, 0}, {0, 1, 1, 0}},
  "Z-bits" -> {{1, 0, 0, 1}, {0, 0, 0, 1}},
  "Phase" -> {0, 0, 0, 0}|>
```

The `X-bits` and `Z-bits` arrays are shape `{n_qubits, 2*GeneratorCount}` — rows are qubits, columns are tableau rows (destabilizers first, then stabilizers). For Bell state, last 2 columns are stabilizers `XX` and `ZZ`: column 3 = X-bits `{1, 1}` + Z-bits `{0, 0}` = `XX`; column 4 = X-bits `{0, 0}` + Z-bits `{1, 1}` = `ZZ`.

**Pauli multiplication via the symplectic group.** The reason the symplectic encoding is useful: Pauli composition `P · Q` reduces to `BitXor` on the `(x, z)` bit pairs, while the prefactor (the `±1` or `±i`) accumulates separately in the phase. Concretely, `X · Z = −i Y` at the matrix level corresponds to `BitXor[(1, 0), (0, 1)] = (1, 1) = bits(Y)` at the symplectic level — the `−i` lives in the phase, not in the bits.

**Code:**
```wolfram
Module[{xMat, zMat, yMat, prodMat},
    xMat = Normal @ QuantumOperator["X"]["Matrix"];
    zMat = Normal @ QuantumOperator["Z"]["Matrix"];
    yMat = Normal @ QuantumOperator["Y"]["Matrix"];
    prodMat = xMat . zMat;
    <|
        "X . Z (matrix)"          -> prodMat,
        "matches -i * Y"          -> (prodMat === -I * yMat),
        "BitXor[bits X, bits Z]"  -> BitXor[{1, 0}, {0, 1}],
        "bits Y"                  -> {1, 1},
        "symplectic match"        -> (BitXor[{1, 0}, {0, 1}] === {1, 1})
    |>
]
```

**Verified output** (block `1.1-multiplication-via-symplectic`):
```
<|"X . Z (matrix)" -> {{0, -1}, {1, 0}},
  "matches -i * Y" -> True,
  "BitXor[bits X, bits Z]" -> {1, 1},
  "bits Y" -> {1, 1},
  "symplectic match" -> True|>
```

The bit-XOR rule generalizes to multi-qubit Paulis: tensor-product bits XOR component-wise, and the phase cocycle (Yashin25 §2.3) tracks `i` factors per pair. This is why every Clifford gate update in `Stabilizer/GateUpdates.m` is some flavor of `BitXor` on the `Tableau` rows.

**Cross-reference:** AarGot04 §2; Mueller26 §2 (binary symplectic representation); Yashin25 §2.3 (phase cocycle). Internal helper `PauliRow` at [Stabilizer/Conversions.m](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m).


### §1.2 — Tableau (extended / improved with destabilizers)

> **Synthesis** (`:60-74`): "Aaronson–Gottesman tableau (AarGot04 §3): a `2n × (2n+1)` binary matrix … Top half = destabilizers (added to halve measurement cost from `O(n³)` to `O(n²)`). Bottom = stabilizer generators."

**QF kernel** (Phase 1). The tableau is shape `{2, n, 2n}`. Stabilizer rows are the last `n` of the third axis; destabilizer rows are the first `n`. Property accessors split the two halves: `ps["Stabilizer"]`, `ps["Destabilizer"]`, `ps["Matrix"]`.

**Code:**
```wolfram
ps = PauliStabilizer["5QubitCode"];
<|
    "Qubits"             -> ps["Qubits"],
    "GeneratorCount"     -> ps["GeneratorCount"],
    "TableauDimensions"  -> Dimensions[ps["Tableau"]],
    "MatrixDimensions"   -> Dimensions[ps["Matrix"]],
    "DestabilizerExample" -> ps["Destabilizers"][[1]],
    "StabilizerExample"   -> ps["Stabilizers"][[1]]
|>
```

**Verified output** (block `1.2-tableau-shape`):
```
<|"Qubits" -> 5, "GeneratorCount" -> 5, "TableauDimensions" -> {2, 5, 10},
  "MatrixDimensions" -> {10, 10},
  "DestabilizerExample" -> "ZXXZI",
  "StabilizerExample" -> "XZZXI"|>
```

**AG invariant verification** (synthesis §1.2 line 74): `R J R^T = J` over `𝔽₂` where `J` is the symplectic form (block off-diagonal identity).

**Code:**
```wolfram
ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
n = 2;
omega = ArrayFlatten[{
    {ConstantArray[0, {2, 2}], IdentityMatrix[2]},
    {IdentityMatrix[2], ConstantArray[0, {2, 2}]}
}];
With[{m = ps["Matrix"]}, Mod[m . omega . Transpose[m] - omega, 2]]
```

**Verified output** (block `1.2-AG-symplectic-invariant`):
```
{{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
```

The all-zeros matrix confirms `M Ω Mᵀ ≡ Ω (mod 2)`. This is `Tier 2` of `Tests/PauliStabilizer.wlt` and is checked for Bell, GHZ-3, GHZ-5, Cluster-5, 5Q, Steane, Shor, and `RandomClifford[4]`.

**Cross-reference:** AarGot04 Prop 2; Yashin25 §2.3.


### §1.3 — Graph state representation (Anders & Briegel)

> **Synthesis** (`:76-84`): "Every stabilizer state is local-Clifford-equivalent to a graph state (AndBri05 §2). Stored as: an adjacency list of the graph G; a list of n vertex operators (VOPs), each one of the 24 single-qubit Clifford operators."

**QF kernel** (Phase 5, [Stabilizer/GraphState.m](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m)). `GraphState[<|"Graph" -> g_Graph, "VOPs" -> {0, ..., 0}|>]`. Phase 5 v1 supports identity VOPs only (richer 24-element table deferred — see synthesis §2.2).

The stabilizer at vertex `i` is `K_i = X_i ⊗ Π_{j ∈ N(i)} Z_j` (AndBri05 Eq 1).

**Code:**
```wolfram
gs = GraphState[Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]]];
<|
    "Vertices"    -> gs["VertexCount"],
    "Edges"       -> gs["EdgeCount"],
    "Stabilizers" -> gs["Stabilizers"]
|>
```

**Verified output** (block `1.3-graph-state-stabilizers`):
```
<|"Vertices" -> 5, "Edges" -> 4,
  "Stabilizers" -> {"XZIII", "ZXZII", "IZXZI", "IIZXZ", "IIIZX"}|>
```

Linear cluster on 5 vertices: K_i has X at qubit `i` and Z at neighbors `i-1, i+1` (with boundary conditions).

**Round-trip with cluster-state circuit:**
```wolfram
n = 4;
gs = GraphState[Graph[Range[n], Table[i \[UndirectedEdge] (i + 1), {i, n - 1}]]];
ps = gs["PauliStabilizer"];
psFromCirc = PauliStabilizer[QuantumCircuitOperator[
    Join[Table["H" -> i, {i, n}], Table["CZ" -> {i, i + 1}, {i, n - 1}]]
]];
ps["Stabilizers"] === psFromCirc["Stabilizers"]
```

**Verified output** (block `1.3-graph-state-equals-cluster-circuit`):
```
True
```

The graph-state stabilizer extraction matches the circuit-based construction (H to all qubits, CZ between adjacent qubits).

**Cross-reference:** AndBri05 §2 (graph + 24-VOP encoding); van den Nest, Dehaene, De Moor 2004 (LC-equivalence theorem).


### §1.4 — Quadratic-form triple (Dehaene–De Moor)

> **Synthesis** (`:86-90`): "`|s⟩ ∝ Σ_{z ∈ V} (-1)^{Q(z)} i^{ℓ(z)} |z + z₀⟩` with V ⊂ 𝔽₂ⁿ vector subspace … `Q` quadratic, `ℓ` linear … the form most natural for symbolic computation."

**Status:** ⏸ **deferred to Phase 6+ / v2.** The DehMoo03 quadratic-form representation is mathematically the most compact and best for closed-form symbolic computation, but extracting `(V, Q, ℓ)` from a tableau requires Gaussian elimination over `𝔽₂` plus a separate quadratic-form recovery (deSilSalYin23 §3 Theorem 2.5).

For Phase 1–5, the kernel uses the tableau form directly (synthesis §1.2). The quadratic-form is reachable through `ps["State"]["StateVector"]` for `n ≤ 8`, but the structured `(V, Q, ℓ)` form is not yet a first-class output.

**Cross-reference:** DehMoo03 §4 Theorem 3; HosDehMoo04 §V Theorem 1; deSilSalYin23 §3.


### §1.5 — Stabilizer frame (Quipu)

> **Synthesis** (`:92-98`): "A list of stabilizer-state generators sharing a global phase (GarMar15 §3). Used to represent superpositions of stabilizer states — i.e. arbitrary states with bounded stabilizer rank. Critical for handling non-Clifford gates (Toffoli, T) symbolically."

**QF kernel** (Phase 4, [Stabilizer/StabilizerFrame.m](../../../QuantumFramework/Kernel/Stabilizer/StabilizerFrame.m)). `StabilizerFrame[<|"Components" -> {{c_i, ps_i}, ...}|>]` represents `Σ_i c_i |s_i⟩`. Closes under Clifford gates (which distribute over components) and under non-Clifford `P[θ] / T / T†` (which double the component count).

**Code (T gate produces a 2-component frame):**
```wolfram
psT = PauliStabilizer[1]["T", 1];
<|
    "Head"          -> Head[psT],
    "Length"        -> psT["Length"],
    "Coefficients"  -> psT["Coefficients"]
|>
```

**Verified output** (block `1.5-stabilizer-frame-from-T-gate`):
```
<|"Head" -> StabilizerFrame, "Length" -> 2,
  "Coefficients" -> {(1 + E^((I/4)*Pi))/2, (1 - E^((I/4)*Pi))/2}|>
```

**Closure under further Clifford:**
```wolfram
psPlus = PauliStabilizer[1]["H", 1];
psTH = psPlus["T", 1]["H", 1];
{Head[psTH], psTH["Length"]}
```

**Verified output** (block `1.5-frame-closes-under-Clifford`):
```
{StabilizerFrame, 2}
```

After T followed by H, still a `StabilizerFrame` with 2 components — Clifford gates distribute over the components without doubling the frame size.

**Materialization (T|0⟩ = |0⟩, eigenstate):**
```wolfram
psT = PauliStabilizer[1]["T", 1];
vec = Normal @ psT["StateVector"];
Chop @ N @ FullSimplify[vec - {1, 0}]
```

**Verified output** (block `1.5-frame-materialization-T-on-zero`):
```
{0, 0}
```

The materialized state vector `(1+e^(iπ/4))/2 * |0⟩ + (1-e^(iπ/4))/2 * |0⟩ = |0⟩` after simplification — confirming `T|0⟩ = |0⟩`.

**Cross-reference:** GarMar15 §3 (Quipu stabilizer frames).


### §1.6 — Clifford channel (Yashin)

> **Synthesis** (`:100-104`): "Choi-state stabilizer tableau `[U_A | U_B | c]` describing arbitrary stabilizer operations — including measurements, dephasing, qubit discarding, mixed-state preparations (Yashin25 §2.3). Composition becomes 'find a basis of the intersection of two vector subspaces'."

**Status:** ⏸ **deferred to Phase 6+ / v2.** The Yashin25 Choi-tableau formalism unifies pure states, mixed states, measurements, and post-selection into a single `[U_A | U_B | c]` matrix. Composition via vector-space intersection (Gaussian elimination) is mathematically clean but requires a new `CliffordChannel` head with measurement-record propagation.

For Phase 1–5, the kernel handles Clifford channels through the existing `QuantumChannel` infrastructure (Stinespring dilation form), routed via `PauliStabilizerApply` for the Clifford subset.

**Cross-reference:** Yashin25 §2.3, §3.2.


---

## §2 — What the package must compute

### §2.1 — Tableau-update rules for Clifford gates

> **Synthesis** (`:110-118`): "For gate `U ∈ {H, S, CNOT}` applied to qubit `a` (or `a → b` for CNOT), update each tableau row …"

**QF kernel** (Phase 1, [Stabilizer/GateUpdates.m](../../../QuantumFramework/Kernel/Stabilizer/GateUpdates.m)). All Clifford generators implemented: H, S, S†, X, Y, Z, CNOT, CX (alias), CZ, SWAP, V, V†.

**Code (Heisenberg conjugation: H takes Z to X):**
```wolfram
PauliStabilizer[1]["H", 1]["Stabilizers"]
```

**Verified output** (block `2.1-Heisenberg-H-on-Z`):
```
{"X"}
```

Starting from `|0⟩` (stabilizer Z), applying H gives stabilizer `X` — confirming `H Z H† = X`.

**Code (CNOT takes X⊗I to X⊗X):**
```wolfram
Sort @ PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
```

**Verified output** (block `2.1-Heisenberg-CNOT-XI-to-XX`):
```
{"XX", "ZZ"}
```

H|0⟩⊗|0⟩ has stabilizers `{XI, IZ}`; after CNOT(1,2), stabilizers are `{XX, ZZ}` — the Bell state.

**Cross-reference:** AarGot04 §3; Biswas24 §3 (pedagogical derivation); PatGuh26 §3.3 (Karnaugh-map derivation).


### §2.2 — The 24-element local Clifford group

> **Synthesis** (`:126-134`): "every single-qubit Clifford acts on Pauli operators as a `3! × ± = 24`-fold action on `{X, Y, Z}` (AndBri05 §2). Encode as a `24 × 24` multiplication table (lookup) plus a 24-entry decomposition table."

**Status:** ⏸ **deferred.** Phase 5 v1 supports only the identity VOP (index 0) in `GraphState`. The `LocalClifford[]` 24-element catalog and the 24×24 multiplication table are required for the richer graph-state algorithms (e.g., propagating Clifford updates without unrolling to circuits).

**Cross-reference:** AndBri05 §2 footnote, And05 (the supplementary "24-VOP" table).


### §2.3 — Local complementation

> **Synthesis** (`:136-140`): "Given a vertex `a` in graph `G`, `LocalComplement[G, a]` complements all edges among `a`'s neighbors (AndBri05 Def 1)."

**QF kernel** (Phase 5, [Stabilizer/GraphState.m](../../../QuantumFramework/Kernel/Stabilizer/GraphState.m)). `LocalComplement[g, v]` for `g_Graph` toggles edges among `AdjacencyList[g, v]`. Also accepts `GraphState` input (passes VOPs through unchanged — VOP tracking deferred).

**Code (LC at the center of a star turns it into the complete graph on the leaves):**
```wolfram
g = Graph[{1, 2, 3, 4}, {1 \[UndirectedEdge] 2, 1 \[UndirectedEdge] 3, 1 \[UndirectedEdge] 4}];
Sort @ EdgeList @ LocalComplement[g, 1]
```

**Verified output** (block `2.3-local-complement-star-to-wheel`):
```
{UndirectedEdge[1, 2], UndirectedEdge[1, 3], UndirectedEdge[1, 4],
 UndirectedEdge[2, 3], UndirectedEdge[2, 4], UndirectedEdge[3, 4]}
```

The edges among `{2, 3, 4}` (none originally) get toggled to all 3 edges, giving K₄ on the leaves plus the star.

**Code (LC is involutive):**
```wolfram
g = Graph[Range[5], Table[i \[UndirectedEdge] (i + 1), {i, 4}]];
Sort @ EdgeList @ LocalComplement[LocalComplement[g, 3], 3] === Sort @ EdgeList[g]
```

**Verified output** (block `2.3-local-complement-involutive`):
```
True
```

`LC ∘ LC = id` for any vertex.

**Cross-reference:** AndBri05 Def 1, Theorem 1 (LC preserves entanglement spectrum); PatGuh26 §3.2.


### §2.4 — Measurement (random and deterministic outcomes)

> **Synthesis** (`:142-160`): "Z-basis measurement of qubit `a` on a tableau state … Random (some `p ∈ {n+1, …, 2n}` has `x_{pa} = 1`) … Deterministic (no such `p`) … `X` and `Y` measurements: precondition with `H_a` or `SH_a`, measure `Z_a`, postcondition."

**QF kernel** (Phase 1 + Phase 4). Single Z-basis at [Stabilizer/Measurement.m](../../../QuantumFramework/Kernel/Stabilizer/Measurement.m); arbitrary Pauli string at [Stabilizer/PauliMeasure.m](../../../QuantumFramework/Kernel/Stabilizer/PauliMeasure.m). Result is an `Association` keyed by outcome bit.

**Z-basis measurement (Bell, qubit 1) — non-deterministic:**
```wolfram
psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
Sort @ Keys @ psBell["M", 1]
```

**Verified output** (block `2.4-Z-measurement-Bell`):
```
{0, 1}
```

Two keys → non-deterministic 50/50 outcome.

**Pauli-string measurement (5-qubit code stabilizers — deterministic):**
```wolfram
ps5Q = PauliStabilizer["5QubitCode"];
{
    Sort @ Keys @ ps5Q["M", "XZZXI"],
    Sort @ Keys @ ps5Q["M", "IXZZX"],
    Sort @ Keys @ ps5Q["M", "XIXZZ"],
    Sort @ Keys @ ps5Q["M", "ZXIXZ"]
}
```

**Verified output** (block `2.4-Pauli-string-measurement-5Q`):
```
{{0}, {0}, {0}, {0}}
```

All 4 stabilizer measurements on `|0_L⟩` are deterministic with outcome bit 0 (eigenvalue +1) — the defining property of the encoded code state.

**Cross-reference:** AarGot04 §3 (measurement-update rules with `rowsum` primitive); PatGuh26 §3.4 (rotation-based Pauli measurement).


### §2.5 — Inner products and expectation values

> **Synthesis** (`:162-166`): "`StabilizerInnerProduct[ψ, φ]`: zero if the stabilizer groups have a Pauli with opposite signs; otherwise `2^(-s/2)` where `s` is the minimal symmetric difference of generators (GarMarCro12 §3, `O(n³)` algorithm)."

**QF kernel** (Phase 4, [Stabilizer/InnerProduct.m](../../../QuantumFramework/Kernel/Stabilizer/InnerProduct.m)). `StabilizerInnerProduct[ψ, φ]` computes `⟨ψ|φ⟩` and `StabilizerExpectation[ps, "XZZXI"]` returns `⟨ψ|P|ψ⟩`. Phase 4 v1 uses **direct vector materialization** (cost `2ⁿ`); the `O(n³)` GarMarCro12 closed-form is a TODO for Phase 6+.

**Code (Bell self inner product):**
```wolfram
psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
StabilizerInnerProduct[psBell, psBell]
```

**Verified output** (block `2.5-inner-product-self`):
```
1
```

**Code (orthogonal Bell states):**
```wolfram
psPhiPlus = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
psPhiMinus = psPhiPlus["Z", 1];   (* |\[CapitalPhi]+\[RightAngleBracket] -> |\[CapitalPhi]-\[RightAngleBracket] *)
Chop @ N @ StabilizerInnerProduct[psPhiPlus, psPhiMinus]
```

**Verified output** (block `2.5-inner-product-orthogonal`):
```
0
```

**Code (Bell expectation values):**
```wolfram
psBell = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]];
<|
    "<XX>" -> StabilizerExpectation[psBell, "XX"],
    "<ZZ>" -> StabilizerExpectation[psBell, "ZZ"],
    "<YY>" -> StabilizerExpectation[psBell, "YY"],
    "<XI>" -> StabilizerExpectation[psBell, "XI"]
|>
```

**Verified output** (block `2.5-stabilizer-expectation-on-Bell`):
```
<|"<XX>" -> 1, "<ZZ>" -> 1, "<YY>" -> -1, "<XI>" -> 0|>
```

`⟨XX⟩ = ⟨ZZ⟩ = +1` because both are stabilizers. `⟨YY⟩ = -1` because `YY = (iXZ)⊗(iXZ) = -1·XX·ZZ` (the i-factor matters! — recovered correctly via the direct-vector fallback). `⟨XI⟩ = 0` because `XI` anticommutes with `ZZ`.

**Cross-reference:** GarMarCro12 §3; deSilSalYin23 §3 (faster algorithms via the quadratic-form path).


### §2.6 — Distance / nearest neighbors

> **Synthesis** (`:168-170`): "For an n-qubit stabilizer state, there are exactly `4(2ⁿ - 1)` nearest-neighbor stabilizer states with `|⟨ψ|φ⟩| = 2^(-1/2)` (GarMarCro12 §5)."

**Status:** ⏸ **deferred.** Useful as a sanity check + building block for stabilizer rank algorithms; not on the v1 critical path.

**Cross-reference:** GarMarCro12 §5.


### §2.7 — Counting and enumeration

> **Synthesis** (`:172-178`): "Number of n-qubit stabilizer states: `N(n) = 2ⁿ Π(2^(n-k)+1)` (AarGot04 Prop 1). `|C_n| = 2^(n²+2n) Π(4ʲ-1)` (KoeSmo14 Eq 2)."

**QF kernel:** these are formula-level checks; embedded in the test fixtures (Tier 2.4 of `Tests/PauliStabilizer.wlt`).

**Code:**
```wolfram
(* N(n) = 2^n * Product[2^(n-k) + 1, {k, 0, n-1}] *)
Function[n, 2^n Product[2^(n - k) + 1, {k, 0, n - 1}]] /@ Range[1, 4]
```

**Verified output** (block `2.7-stabilizer-state-count`):
```
{6, 60, 1080, 36720}
```

`N(1) = 6, N(2) = 60, N(3) = 1080, N(4) = 36720`.

**Code:**
```wolfram
(* |C_n| = 2^(n^2 + 2n) * Product[4^j - 1, {j, 1, n}] *)
Function[n, 2^(n^2 + 2 n) Product[4^j - 1, {j, 1, n}]] /@ Range[1, 3]
```

**Verified output** (block `2.7-clifford-group-order`):
```
{24, 11520, 92897280}
```

`|C_1| = 24, |C_2| = 11520, |C_3| ≈ 9.29 × 10⁷`.

**Cross-reference:** AarGot04 Prop 1; KoeSmo14 §1.1, Eq 2.


### §2.8 — Random Clifford / random stabilizer state

> **Synthesis** (`:180-188`): "`RandomClifford[n]` and `RandomStabilizerState[n]`: KoeSmo14 §3.2 gives an `O(n³)` algorithm via symplectic transvections that maps `Range[CliffordGroupOrder[n]] → Sp(2n, 𝔽₂)` bijectively."

**QF kernel** (Phase 2, [Stabilizer/RandomClifford.m](../../../QuantumFramework/Kernel/Stabilizer/RandomClifford.m)). `RandomClifford[n]` implements the Bravyi-Maslov / Koenig-Smolin Mallows-distribution sampler. Promoted to public symbol in Phase 2.

**Code:**
```wolfram
SeedRandom[20260430];
With[{ps = RandomClifford[3]},
    <|"Qubits" -> ps["Qubits"], "Stabilizers" -> ps["Stabilizers"]|>
]
```

**Verified output** (block `2.8-random-clifford`):
```
<|"Qubits" -> 3, "Stabilizers" -> {"-ZXZ", "-ZZX", "XXX"}|>
```

**Code (uniformity smoke test — 200 samples on n=1):**
```wolfram
SeedRandom[12345];
Length @ DeleteDuplicates @ Table[
    With[{r = RandomClifford[1]}, {r["Stabilizers"], r["Signs"]}],
    {200}
]
```

**Verified output** (block `2.8-random-clifford-uniformity`):
```
12
```

12 distinct {stab strings, sign combinations} appeared in 200 samples. Note: this is fewer than 24 because `RandomClifford` only returns a *state* (one of 6 single-qubit stabilizer states modulo signs); the full Clifford gate group has 24 elements but maps to 6 distinct *stabilizer states* (modulo phase).

**Cross-reference:** KoeSmo14 §3.2; Bravyi-Maslov 2020 (the same algorithm).


### §2.9, §2.10 — Canonical forms and synthesis from a tableau

> **Synthesis** (`:190-205`): "Aaronson-Gottesman canonical form (`H-C-P-C-P-C-H-P-C-P-C`) … Garcia-Markov-Cross canonical form (`H-C-CZ-P-H`) … Standard form `H_s` for QEC encoders … `CliffordSynthesize[tableau, ConnectivityGraph -> g, GateSet -> {...}]`."

**Status:** ⚠️ **partial.** [Stabilizer/Conversions.m](../../../QuantumFramework/Kernel/Stabilizer/Conversions.m) contains the AG greedy canonicalization (`ps["Circuit"]`) which produces a Clifford circuit whose dagger equals `ps`. It is **not** minimized (audit §8.2: "circuit length is **not minimized** — Reid24, Winderl23 propose better").

**Deferred:**
- `Method -> "GarciaMarkov"` (GarMarCro12 H-C-CZ-P-H form) — useful for inner-product computation.
- `Method -> "Winderl"` (Winderl23 Algorithm 1) — connectivity-aware Steiner-tree-based pivot.
- `Method -> "Reid"` (Reid24) — empirical 2Q-gate minimization.

**Cross-reference:** AarGot04 §5 Theorem 8; GarMarCro12 §3; Reid24; Winderl23 §IV Algorithm 1; MonPar23 §IV (encoder synthesis from `H_s`).


### §2.11 — Pauli tracking

> **Synthesis** (`:207-217`): "Given a `Clifford ∘ Measurement ∘ Clifford ∘ ...` circuit, propagate Pauli corrections without applying them on hardware (Paler14, RuhDev25). Per qubit, the Pauli 'frame' is one of `{I, X, Z, XZ}`."

**Status:** ⏸ **deferred.** Required for the QF MBQC pipeline. Phase 5 v1 does not include `PauliTrack[circuit, frame]` or the `PauliFrame[]` data type.

**Cross-reference:** Paler14 §V (algorithm), RuhDev25 §3 (modern library + MBQC scheduling).


### §2.12 — QEC code extraction

> **Synthesis** (`:219-229`): "Given a stabilizer subgroup, return: `[[n, k, d]]` parameters; logical `X̄, Z̄` operators; encoder circuit; syndrome decoder; pretty-printed stabilizer table."

**QF kernel.** Named codes provided (Phase 1, [Stabilizer/NamedCodes.m](../../../QuantumFramework/Kernel/Stabilizer/NamedCodes.m)): `"5QubitCode"`, `"5QubitCode1"`, `"SteaneCode"` (= `"7QubitCode"`), `"7QubitCode1"`, `"SteaneCode1"`, `"9QubitCode"`, `"9QubitCode1"`, `"Random"`. Syndrome extraction via `StabilizerExpectation` (commute/anticommute test) or `ps["M", "XZZXI"]` (Pauli-string measurement). Code distance via direct enumeration (Tier 3 of `Tests/PauliStabilizer.wlt`).

**Code (5Q syndrome uniqueness — defining property of `[[5,1,3]]`):**
```wolfram
n = 5;
sympIP[v1_, v2_] := Mod[v1[[;; n]] . v2[[n + 1 ;; 2 n]] + v2[[;; n]] . v1[[n + 1 ;; 2 n]], 2];
gens = (Module[{xs, zs},
    {xs, zs} = Transpose @ Replace[Characters[#],
        {"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}}, {1}];
    Join[xs, zs]
]) & /@ Take[PauliStabilizer["5QubitCode"]["Stabilizers"], 4];
errors = Flatten[Table[
    With[{e = ConstantArray[0, 2 n]},
        Switch[op,
            "X", ReplacePart[e, i -> 1],
            "Y", ReplacePart[e, {i -> 1, n + i -> 1}],
            "Z", ReplacePart[e, n + i -> 1]
        ]
    ],
    {i, n}, {op, {"X", "Y", "Z"}}], 1];
syndromes = Table[Table[sympIP[err, g], {g, gens}], {err, errors}];
<|"NumErrors" -> Length[errors], "DistinctSyndromes" -> Length @ Union[syndromes]|>
```

**Verified output** (block `2.12-5Q-syndromes`):
```
<|"NumErrors" -> 15, "DistinctSyndromes" -> 15|>
```

All 15 single-qubit Pauli errors (X_i, Y_i, Z_i for i ∈ 1..5) produce distinct 4-bit syndromes — the defining property of an `[[5, 1, 3]]` code (corrects 1 arbitrary error).

**Code distance d=3** is verified directly in `Tests/PauliStabilizer.wlt` `Tier 3.B` via enumeration of `N(S) \ S` (1024 Pauli vectors filtered to commuting normalizer minus the 16-element stabilizer subgroup, min weight = 3).

**Cross-reference:** Got97 §3.5 (5-qubit cyclic code); Got00 §4 (Steane CSS construction); MonPar23 §IV (encoder synthesis); the Quantum Singleton bound `n - k ≥ 2(d-1)` is saturated for `[[5,1,3]]` (5 - 1 = 2(3-1) = 4).


---

## §3 — Symbolic operations

### §3.1 — Symbolic Pauli arithmetic

> **Synthesis** (`:236-244`): "`PauliMultiply[X_1 Z_2, Y_1 Z_2]` should return `i Z_1` symbolically … Commutators: `[P, Q] = 0` iff the symplectic inner product `Σ_i (x_i z'_i + x'_i z_i) mod 2 = 0` (Got98 §2)."

**QF kernel** (Phase 3). The predicate `PauliStabilizerQ` was loosened to accept symbolic signs (Phase 3); `BitXor` propagates symbolic phases through Clifford gate updates without modification (audit §G Risk 2 verified). Concrete-only paths use `ConcretePauliStabilizerQ` ([Stabilizer/PauliStabilizer.m:33-41](../../../QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m)).

The full Mueller26 `PauliEngine` symbolic-coefficients-on-both-sides commutator engine is **deferred**.

**Cross-reference:** Mueller26 §3 (PauliEngine algorithms); audit §3.1.


### §3.2 — Symbolic phases & symbolic measurement outcomes

> **Synthesis** (`:246-258`): "FangYing23 §3 (SymPhase): represent the sign vector `r⃗` in the tableau as bit-vectors over `𝔽₂^(n_s + 1)`, where `n_s` = number of fresh symbols introduced so far."

**QF kernel** (Phase 3, [Stabilizer/SymbolicMeasure.m](../../../QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m)). Three new public symbols:

- `StabilizerMeasure[ps, q]` — single Z-basis measurement; allocates fresh `\[FormalS][k]` symbol for non-deterministic case.
- `SubstituteOutcomes[ps, rules]` — replace measurement-outcome symbols with concrete 0/1.
- `SampleOutcomes[ps, n]` — n random samples by independently substituting each symbol.

**Code (StabilizerMeasure allocates fresh symbol):**
```wolfram
psPlus = PauliStabilizer[1]["H", 1];
psSym = StabilizerMeasure[psPlus, 1];
syms = DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
<|
    "Head"          -> Head[psSym],
    "FreshSymbols"  -> syms,
    "Phase"         -> psSym["Phase"]
|>
```

**Verified output** (block `3.2-symbolic-measurement-allocates-fresh-symbol`):
```
<|"Head" -> PauliStabilizer, "FreshSymbols" -> {\[FormalS][1]},
  "Phase" -> {0, \[FormalS][1]}|>
```

The `\[FormalS][1]` symbol stamps the second phase entry (= the stabilizer Z's sign-bit). `SubstituteOutcomes[psSym, \[FormalS][1] -> 0]` returns the outcome-0 branch; `-> 1` returns the outcome-1 branch.

**Code (substitute-outcomes round-trip vs regular `["M"]`):**
```wolfram
psPlus = PauliStabilizer[1]["H", 1];
psSym = StabilizerMeasure[psPlus, 1];
sym = First @ DeleteDuplicates @ Cases[psSym["Phase"], _\[FormalS], Infinity];
psSub0 = SubstituteOutcomes[psSym, sym -> 0];
psSub1 = SubstituteOutcomes[psSym, sym -> 1];
ps0 = psPlus["M", 1][0];
ps1 = psPlus["M", 1][1];
<|
    "outcome 0 stabilizers match" -> (psSub0["Stabilizers"] === ps0["Stabilizers"]),
    "outcome 1 stabilizers match" -> (psSub1["Stabilizers"] === ps1["Stabilizers"])
|>
```

**Verified output** (block `3.2-substitute-outcomes-roundtrip`):
```
<|"outcome 0 stabilizers match" -> True,
  "outcome 1 stabilizers match" -> True|>
```

**Phase 3 known limitation** (locked down in `Tier 6 KNOWN LIMITATIONS` of `Tests/PauliStabilizer.wlt`): when a deterministic measurement follows a prior symbolic one (e.g. Bell ZZ correlation), the deterministic outcome polynomial is computed correctly by the AG algorithm but is not stamped into the post-state's signs. The post-state IS physically correct, but `SampleOutcomes` cannot directly recover the deterministic outcome from stabilizer signs alone. Phase 4's `StabilizerFrame` adds the outcome-record machinery needed to fix this; the `Phase3-LIMITATION-DeterministicOutcomeNotStamped` test will flip from passing to failing once that is implemented.

**Cross-reference:** FangYing23 §3 SymPhase (arxiv:2311.03906).


### §3.3 — Symbolic Clifford parameters

> **Synthesis** (`:260-270`): "For variational quantum algorithms and ansatz design, the Clifford gates depend on parameters `θ` that are symbolic until fixed. Mueller26 PauliEngine demonstrates this for Pauli rotations `e^(-iθP/2)`."

**Status:** ⏸ **deferred.** Parametric Clifford rotations + parameter-shift-rule gradients (Mueller26 §3 Generator Gradients) are required for VQE/QAOA differentiability but not yet implemented.

**Cross-reference:** Mueller26 §3; Schuld 2019 (parameter-shift rule).


### §3.4 — Symbolic dynamical Lie algebras

> **Synthesis** (`:272-282`): "Mueller26 §3 builds the DLA `g = ⟨iG⟩_Lie` for a parameterized circuit by computing nested commutators of Pauli strings until closure."

**Status:** ⏸ **deferred.**

**Cross-reference:** Mueller26 §3.


### §3.5 — Symbolic interconversion between representations

> **Synthesis** (`:284-292`): "`StabilizerStateFromAmplitudes[vec]` … `CliffordTableauFromMatrix[U]` → stabilizer tableau, in `O(n·2ⁿ)` time without reading every matrix entry (deSilSalYin23 §4.3)."

**Status:** ⚠️ **partial.** The `O(2ⁿ)`/`O(4ⁿ)` direct-conversion paths exist (constructor `PauliStabilizer[qs_QuantumState]` runs the 4ⁿ Pauli-expectation tomography). The deSilSalYin23 closed-form algorithms (`O(n·2ⁿ)` and faster) are deferred.

**Cross-reference:** deSilSalYin23 §3 (10 fast interconversion algorithms).


### §3.6 — Symbolic qudits

> **Synthesis** (`:294-303`): "HosDehMoo04 and Beaudrap11 generalize everything to `d`-dimensional qudits using `ℤ_d` (or `ℤ_(2d)` when `d` even). WinPay24b condensed encodings unifies the parity treatment."

**Status:** ⏸ **deferred to v2.** This is the synthesis's biggest "comparative advantage" claim ("the only one in the literature handling all `d` in a single uniform code path") but requires a substantial qudit-aware refactor of the Pauli arithmetic. Not on the v1 critical path.

**Cross-reference:** Beaudrap11 §2-3 (Weyl-operator linearization); HosDehMoo04 §III-IV; WinPay24b (condensed encoding for even `d`).


---

## §4–§6 — User-facing menus

### §4 — Stabilizer-state menu (`:307-326`)

| Function | Status | Notes |
|---|---|---|
| `StabilizerStateQ[ψ]` | partial | use `PauliStabilizerQ`/`MatchQ[..., _PauliStabilizer]` |
| `StabilizerCheckMatrix[ψ]` | ✅ | `ps["Matrix"]` |
| `StabilizerGenerators[ψ]` | ✅ | `ps["Stabilizers"]` |
| `Destabilizers[ψ]` | ✅ | `ps["Destabilizers"]` |
| `StabilizerTableau[ψ]` | ✅ | `ps["Tableau"]` |
| `StabilizerToGraph[ψ]` | ✅ | `GraphState[ps]` (graph-form input only) |
| `GraphToStabilizer[g, vops]` | ✅ | `GraphState[g]["PauliStabilizer"]` (identity VOPs) |
| `StabilizerToQuadraticForm[ψ]` | ⏸ | DehMoo03 §4 |
| `StabilizerInnerProduct[ψ, φ]` | ✅ | direct vector |
| `StabilizerDistance[ψ, φ]` | ⏸ | quantum Singleton |
| `StabilizerEntanglement[ψ, partition]` | partial | via Schmidt rank of materialized state |
| `StabilizerNearestNeighbors[ψ]` | ⏸ | GarMarCro12 §5 |
| `LocalCliffordEquivalent[ψ, φ]` | ⏸ | LC-equivalence test |
| `LocalComplement[g, vertex]` | ✅ | graph only; VOP tracking deferred |

### §5 — Clifford-operations menu (`:330-348`)

| Function | Status | Notes |
|---|---|---|
| `CliffordOperatorQ[U]` | ⏸ | deSilSalYin23 §4 |
| `CliffordTableau[U]` | partial | via `PauliStabilizer[qo_QuantumOperator]` |
| `CliffordMatrix[T]` | ✅ | `ps["Operator"]` (cost 4ⁿ) |
| `CliffordCompose[T1, T2]` | ✅ | `ps1[ps2]` |
| `CliffordInverse[T]` | partial | `ps["Dagger"]` (latent infinite-recursion bug — see §3.G of plan) |
| `CliffordSymplectic[T]` | ✅ | `ps["Matrix"]` |
| `CliffordCircuit[T, gateSet, connectivity]` | partial | `ps["Circuit"]` (AG only, no connectivity) |
| `CliffordCanonicalForm[T, "AaronsonGottesman"]` | ✅ | `ps["Circuit"]` |
| `RandomClifford[n]` | ✅ | Phase 2 |
| `IndexClifford[T]` | ⏸ | KoeSmo14 §3.3 inverse map |
| `ParametricCliffordRotation[P, θ]` | ⏸ | Mueller26 |
| `CliffordChannel[circuit]` | ⏸ | Yashin25 |
| `ChannelCompose[Φ, Ψ]` | ⏸ | Yashin25 vector-space intersection |
| `PauliTrack[circuit, frame]` | ⏸ | Paler14, RuhDev25 |

### §6 — Specifically *symbolic* operations (`:352-365`)

| Function | Status |
|---|---|
| `SymbolicMeasurementOutcome[tableau, qubit]` | ✅ `StabilizerMeasure` |
| `SubstituteOutcomes[symbolic_state, rules]` | ✅ `SubstituteOutcomes` |
| `SymbolicPauliCommutator[P, Q]` | ⏸ Mueller26 |
| `DynamicalLieAlgebra[generators]` | ⏸ Mueller26 |
| `LieClosure[generators]` | ⏸ |
| `StructureConstants[basis]` | ⏸ |
| `StabilizerEntropy[ψ, n]` | partial via Schmidt |
| `MagicMonotone[ρ]` | ⏸ |
| `ParametricCircuit[gates, params]` | ⏸ |
| `GradientPauli[expectation, parameter]` | ⏸ Mueller26, parameter-shift |
| `StabilizerRankDecomposition[ρ]` | partial via `StabilizerFrame` |


---

## §7 — Specific things that surprised me

The synthesis (`:368-386`) lists 8 surprises. Phase 1–5 implementations and deferrals:

- **Yashin25's contribution is the cleanest.** Single Boolean matrix `[U_A | U_B | c]` for state, channel, measurement, post-selection. ⏸ **Not yet adopted**; the kernel still has separate code paths for `PauliStabilizer` and `QuantumChannel`. Phase 6+ recommendation: introduce `CliffordChannel` head with Yashin25 contraction semantics.

- **AndBri05 implementation cost.** "1400 lines, no discrepancies after 4×10⁶ ops." Phase 1–5 kernel is ~1500 LOC with 185 verification tests; the same correctness oracle (cross-check `["State"]` for small n) is exercised by Tier 3.

- **Patil & Guha cluster-state rule book.** Graph rewrites for X/Y/Z measurements on cluster states. ⚠️ Pauli-string measurement IS implemented (Phase 4), but the explicit graph-rewrite forms (PatGuh26 §3, §4) are not yet a separate `ClusterStateRules.wl` module.

- **Beaudrap11's Weyl-operator trick.** ⏸ Required for the qudit unification (§3.6, deferred to v2).

- **deSilSalYin23 §4.3 tableau-from-unitary in O(n·2ⁿ).** ⏸ Not yet implemented; Phase 1 uses the `4ⁿ` tomography path (`PauliStabilizerTableau`).

- **MonPar23 §IV encoder synthesis from `H_q`.** ⏸ Useful for arbitrary stabilizer codes; Phase 5 does not include `EncoderCircuit[code]`.

- **GarMar15's stabilizer frames.** ✅ Phase 4 `StabilizerFrame`.

- **KoeSmo14 §3.3 inverse map.** ⏸ Deferred (see §5 menu).

- **No paper covers symbolic commutators for parametric coefs on both sides.** This is a Wolfram research opportunity. ⏸ Not addressed.


## §8 — Suggested architecture for the package

The synthesis (`:391-447`) proposes an 11-module sub-directory layout. Phase 1 chose a **flatter 10-file** structure (revised to **15 files** by Phase 5):

```
QuantumFramework/Kernel/Stabilizer/
├── PauliStabilizer.m       Predicate, dispatcher, Properties contract
├── Constructors.m          ~14 constructor patterns
├── NamedCodes.m            $PauliStabilizerNames + named codes
├── GateUpdates.m           Clifford gates + non-Clifford boundary
├── Measurement.m           Z-basis measurement + AG rowsum primitive
├── Properties.m            Property dispatch + display-form properties
├── Conversions.m           PauliRow, State, Operator, Circuit synthesis
├── Compose.m               Symplectic multiplication + tensor product
├── Synthesis.m             (merged into Conversions.m for v1)
├── RandomClifford.m        Mallows distribution sampler
├── Formatting.m            PauliForm, TableauForm, MakeBoxes
├── SymbolicMeasure.m       Phase 3: StabilizerMeasure / Substitute / Sample
├── StabilizerFrame.m       Phase 4: superpositions of stabilizer states
├── InnerProduct.m          Phase 4: StabilizerInnerProduct + Expectation
├── PauliMeasure.m          Phase 4: ps["M", "XZZXI"] arbitrary Pauli
└── GraphState.m            Phase 5: GraphState + LocalComplement
```

The synthesis's `Pauli/`, `Tableau/`, `Frame/`, `Synthesis/`, etc. nested-subdirectories were deferred — that's a v2 promotion candidate when the subsystem grows past ~3000 LOC.


## §9 — What *not* to do

The synthesis's 7 don'ts (`:451-458`) are honored:

- **Don't try to outrun Stim/Qiskit on raw qubit counts.** ✅ Phase 1–5 explicitly does not bit-pack or SIMD-optimize. The TODO in `Stabilizer/Compose.m` for an opt-in bit-packed path is documented but deferred.
- **Don't write loops when broadcasting / `BitXor` over `SparseArray`s suffices.** ✅ All Phase 1 gate updates use `MapIndexed` + `BitXor`; the `Sum`-as-loop in `rowsum` was replaced with `Total @ Table` (Phase 1 vectorization).
- **Don't store global phases of stabilizer states.** ⚠️ **Partially relaxed in Phase 5c.** The kernel still tracks `Phase ∈ {0, 1}` per row for the tableau itself. *In addition*, when a `PauliStabilizer` is constructed from a `QuantumState` or `QuantumOperator`, the constructor now records the overall complex phase that the AG decomposition drops under a `"GlobalPhase"` association key, so that `["State"]` / `["QuantumOperator"]` round-trip the input exactly. Tableau-derived `PauliStabilizer`s (no source state/operator) leave `GlobalPhase` unset (default 1). Gate updates do not yet propagate `GlobalPhase` — see [ROADMAP §A.9](ROADMAP.md).
- **Don't reimplement the qubit-only case for qudits.** ⏸ Phase 1–5 is qubit-only by design; v2 unification deferred.
- **Don't expose Karnaugh-map-derived gate rules to the end user.** ✅ Internal documentation only.
- **Don't try to symbolically diagonalize 2ⁿ × 2ⁿ matrices.** ✅ Phase 4's direct-vector inner-product fallback is gated by user choice (TODO Phase 5+ for closed-form replacement).
- **Don't separate `pure stabilizer state` and `Clifford channel` into different code paths.** ⚠️ **VIOLATED** in v1. Phase 6+ should adopt the Yashin25 unified `CliffordChannel` representation.


## §10 — Direct citations driving each design choice

The synthesis's citation map (`:464-476`) is preserved in the audit document. Phase-by-phase:

- **Phase 1** — AarGot04 (tableau), KoeSmo14 (random Clifford via Mallows).
- **Phase 2** — minor-only; PacletInfo + Usage hygiene.
- **Phase 3** — FangYing23 (SymPhase symbolic measurement).
- **Phase 4** — GarMar15 (stabilizer frames), GarMarCro12 (closed-form inner product, partial — direct-vector fallback in v1).
- **Phase 5** — AndBri05 (graph state, local complementation).

**Awaiting:** Yashin25 (channel tableau, Phase 6+), PatGuh26 (cluster-state rule book, Phase 6+), Mueller26 (PauliEngine symbolic gradients, v2), HosDehMoo04 / Beaudrap11 / WinPay24b (qudits, v2), Reid24 / Winderl23 (hardware-aware synthesis, Phase 6+), Paler14 / RuhDev25 (Pauli tracking, Phase 6+), MonPar23 (encoder synthesis, Phase 6+), deSilSalYin23 (interconversion algorithms, Phase 6+), KocHuaLov17 (Wigner-tableau equivalence, qudits-related v2), PayWin24a (lambda calculus type system, conceptual reference).


## §11 — Final priorities for v1

The synthesis's 10 v1 priorities (`:482-494`) status:

1. **Pauli string + symplectic + symbolic phase** (Mueller26 + FangYing23) — ✅ Phase 1 + 3.
2. **Tableau with destabilizers + AG update rules** (AarGot04) — ✅ Phase 1.
3. **Graph-state + 24-VOP** (AndBri05) — ⚠️ Graph state ✅ (Phase 5); 24-VOP table ⏸.
4. **Quadratic form triple + 6 interconversion algorithms** (DehMoo03 / deSilSalYin23) — ⏸.
5. **Clifford channels via Choi tableau** (Yashin25) — ⏸ deferred.
6. **Random Clifford** (KoeSmo14) — ✅ Phase 2.
7. **Inner product** (GarMarCro12 / deSilSalYin23) — ⚠️ direct-vector ✅ (Phase 4); closed-form ⏸.
8. **Stabilizer frames** (GarMar15) — ✅ Phase 4.
9. **Hardware-aware synthesis** (Winderl23) — ⏸.
10. **Pauli tracking + MBQC scheduling** (RuhDev25) — ⏸.

**v1 score:** 5 ✅ + 2 ⚠️ partial + 3 ⏸ deferred. The Phase 6+ followup work is well-scoped and citation-anchored.


---

## Integration smoke tests

The `Method -> "Stabilizer"` path in `QuantumCircuitOperator` (the load-bearing user-facing API) routes through `PauliStabilizerApply` → the new Stabilizer subsystem.

**Code:**
```wolfram
QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}][Method -> "Stabilizer"]["Stabilizers"]
```

**Verified output** (block `integration-method-stabilizer`):
```
{"XX", "ZZ"}
```

**Code (named circuits):**
```wolfram
QuantumCircuitOperator["GHZ"[3]][Method -> "Stabilizer"]["Stabilizers"]
```

**Verified output** (block `integration-named-circuit-GHZ`):
```
{"XXX", "ZZI", "IZZ"}
```

Both match the expected stabilizer generators — Bell `{XX, ZZ}` and `GHZ_3 {XXX, ZZI, IZZ}` (Got97 §2.2).

---

## Re-verification

To re-run all 26 code blocks and confirm the embedded outputs are still correct:

```bash
wolframscript -f OngoingProjects/Stabilizer/Documentation/verify-synthesis-implementation.wls
```

The verifier loads the local paclet, runs every `block[id, expr]`, and prints `=== id === \n <output>` for each. Diff this output against the `Verified output` blocks in this document to detect drift.

---

## Companion files

All files now live under [`OngoingProjects/Stabilizer/`](..) (moved 2026-05-02 from the previously-gitignored `audit/Stabilizer/` plus `Documentation/Stabilizer/`):

| File | Role |
|---|---|
| [`Documentation/API.md`](API.md) | Per-function API reference (10 public symbols, 49 verified code blocks) |
| [`Documentation/synthesis-implementation.md`](synthesis-implementation.md) | This document — capability tour through synthesis §1–§11 |
| [`Documentation/ROADMAP.md`](ROADMAP.md) | 21 partial / deferred / known-bug items with concrete next steps |
| [`Documentation/verify-synthesis-implementation.wls`](verify-synthesis-implementation.wls) | Runnable verifier for every code block in this document |
| [`Documentation/verify-API.wls`](verify-API.wls) | Runnable verifier for `API.md` |
| [`package-design-synthesis.md`](../package-design-synthesis.md) | Source synthesis distilled from 28 papers |
| [`paulistabilizer-source-audit.md`](../paulistabilizer-source-audit.md) | Line-by-line audit of the pre-Phase-1 monolith |
| [`external-packages-audit.md`](../external-packages-audit.md) | QuantumClifford.jl + Stim reference |
| [`paper-bibliography.md`](../paper-bibliography.md) | 33-paper catalog |
| [`paper-fetch-report.md`](../paper-fetch-report.md) | arXiv fetch provenance |
| [`tex/`](../tex/) | Extracted TeX sources for the 28 papers |
| [`External Packages/`](../External%20Packages/) | QuantumClifford.jl + Stim source code |

**Other:**
- [`Tests/PauliStabilizer.wlt`](../../../Tests/PauliStabilizer.wlt) — 8-tier, 185-test suite (100% passing).
- [`QuantumFramework/Kernel/Stabilizer/`](../../../QuantumFramework/Kernel/Stabilizer/) — kernel implementation (15 files).
- Plan: `/Users/mohammadb/.claude/plans/audit-this-users-mohammadb-documents-git-robust-russell.md` (local-only).
