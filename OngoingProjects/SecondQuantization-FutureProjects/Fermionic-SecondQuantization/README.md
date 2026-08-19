# Fermionic Second Quantization — Development Plan

**Owner:** Bruno Tenorio · **Reviewer:** TBD
**Target:** extend ``Wolfram`QuantumFramework`SecondQuantization` `` with a fermionic layer that mirrors the existing bosonic layer, symbol for symbol.
**Status:** design locked in this document; no kernel code written yet.
**Supersedes:** the one-line stub *"Fermionic Second Quantization (extra)"* in [`../PlanningIdeas.nb`](../PlanningIdeas.nb).

---

## 1. Goal

The `SecondQuantization` context currently implements one algebra: bosons. It does so on two levels at once — a **symbolic** level (`BosonicNormalOrder`, `BosonicVEV`, `FieldVariables`, `ToBosonicOperator`) built on WL's noncommutative algebra, and a **numeric** level (`AnnihilationOperator`, `FockState`, `CovarianceMatrix`) built on truncated `QuantumOperator` matrices, with `ToBosonicOperator` as the bridge between them.

The goal is to add fermions on both levels, with the same bridge, so that a user can:

- derive a fermionic operator identity symbolically (CAR normal ordering, Wick contraction),
- realize it as a concrete `QuantumOperator` on `n` qubits via Jordan–Wigner,
- evolve or measure it with the rest of QuantumFramework (`QuantumEvolve`, `QuantumMeasurementOperator`, `QuantumCircuitOperator`),
- and do all of it in one kernel, including the **mixed** boson–fermion case that no other tool handles symbolically.

The mixed case is the lead capability. Everything else is table stakes that OpenFermion, Qiskit Nature, and PennyLane already provide in Python.

---

## 2. Feasibility — verified in-kernel, not assumed

The load-bearing question was whether WL 15's `NonCommutativePolynomialReduce` can carry canonical **anti**commutation relations plus nilpotency, the way `grobnerNormalOrder` in [BosonicAlgebra.m](../../../QuantumFramework/Kernel/SecondQuantization/BosonicAlgebra.m) carries bosonic ones. It can. All results below were reproduced with `wolframscript`, kernel 15.0.0, paclet 2.1.1.

**Single mode.** Relations `{c, c†} = 1`, `c² = 0`, `c†² = 0`:

```wl
c = \[FormalC]; d = SuperDagger[\[FormalC]];
rels = {Anticommutator[c, d, NonCommutativeMultiply] - 1, c ** c, d ** d};
alg = NonCommutativeAlgebra[<|"Generators" -> {{d}, {c}}|>];
NonCommutativePolynomialReduce[#, rels, alg][[2]] & /@ {c ** d ** c, c ** d ** c ** d, d ** c ** d}
```

| Input | Reduced | Hand check |
|---|---|---|
| `c c† c` | `c` | `(1 - c†c)c = c` ✓ |
| `c c† c c†` | `1 - c†c` | `(1 - c†c)² = 1 - c†c` ✓ |
| `c† c c†` | `c†` | ✓ |

**Two modes**, full CAR set (six anticommutators plus four nilpotency relations), generators graded as `{{d1}, {d2}, {c1}, {c2}}`:

| Input | Reduced | Note |
|---|---|---|
| `c1 c2†` | `-c2† c1` | sign across distinct modes ✓ |
| `c1† c2 c2† c1` | `n1 + c1†c2†c1c2` | matches hand reduction ✓ |
| `[n1, n2]` | `0` | non-trivial consistency check ✓ |

**Numeric side.** Jordan–Wigner assembled from primitives QuantumFramework *already* ships, with no new kernel code:

```wl
sm[k_] := QuantumOperator[{{0, 1}, {0, 0}}, {k}]
jw[k_] := Fold[#1 @ #2 &, sm[k], QuantumOperator["Z", {#}] & /@ Range[k - 1]]
```

At `n = 4` this satisfies every CAR relation exactly: `{c1, c1†} = I`, `{c1, c2†} = 0`, `{c2, c3} = 0`, `c1² = 0`.

**A pitfall confirmed at the same time.** Reaching for the existing bosonic operator at truncation 2 does *not* give a fermion — `AnnihilationOperator[2, {1}]` and `AnnihilationOperator[2, {2}]` yield `{a1, a2†}` with maximum absolute entry `2`, not `0`. Truncating a boson to two levels reproduces the nilpotency but not the anticommutation, because the JW string is missing. This belongs in the tutorial as an explicit warning.

**Four conclusions for the plan.**

1. The symbolic layer needs **no new algebra engine** — it is a relations table plus a generator-ordering convention, mirroring `BosonicRelations` / `grobnerNormalOrder`.
2. The numeric layer needs **no new tensor machinery** — `Fold[#1 @ #2 &, ...]` over disjoint orders is the whole construction.
3. `Anticommutator` is a real ``System` `` builtin, spelled with a **lowercase c**, and takes a third argument for the product. See [Anticommutator-Postmortem.md](../../Courses/QuantumInDiscreteSpace/Anticommutator-Postmortem.md) — this exact spelling has already cost one round of confabulation; do not re-derive it from memory.
4. When an algebra is built with `"Generators"`, `NonCommutativePolynomialReduce` **ignores** an explicit variable list in the third argument and emits `::gvars`. The bosonic `grobnerNormalOrder` currently passes both in its no-scalars branch. The fermionic mirror should pass generators *only* through the algebra, and the bosonic branch should be cleaned up to match.

---

## 3. What exists today (the mirror table)

Every fermionic symbol below has a bosonic counterpart already in the paclet. This is the core design constraint: **a user who knows the bosonic layer should be able to guess the fermionic one.**

| Bosonic (shipped) | File | Fermionic (proposed) | Semantic difference |
|---|---|---|---|
| `AnnihilationOperator[size, order]` | `Operators.m` | `FermionicAnnihilationOperator[k, n]` | no truncation; `n` = mode count, dim `2^n`, JW string attached |
| `QuadratureOperators[size]` | `Operators.m` | `MajoranaOperators[n]` | the `γ` pairs are the fermionic quadratures |
| `DisplacementOperator`, `SqueezeOperator` | `Operators.m` | *(deferred, §8)* | Grassmann-valued; no truncated-matrix analogue |
| `BeamSplitterOperator` | `Operators.m` | `FermionicBeamSplitter` | Givens rotation / mode mixing |
| `PhaseShiftOperator[θ]` | `Operators.m` | reuse — `e^{iθn}` is statistics-agnostic | none |
| `FockState[{n1, n2, …}]` | `States.m` | `OccupationState[{0, 1, …}]` | occupations ∈ {0,1}; JW sign convention fixed by mode order |
| `CoherentState`, `CatState` | `States.m` | `SlaterDeterminantState[U]` | the fermionic Gaussian pure state |
| `ThermalState[nbar]` | `States.m` | `FermionicThermalState[nbar, n]` | Fermi–Dirac, not Bose–Einstein |
| `BosonicRelations[vars]` | `BosonicAlgebra.m` | `FermionicRelations[vars]` | anticommutators **plus nilpotency** — strictly more relations |
| `BosonicNormalOrder`, `BosonicAntinormalOrder` | `BosonicAlgebra.m` | `FermionicNormalOrder`, `FermionicAntinormalOrder` | signs; `ladderSwap` needs a fermionic analogue |
| `BosonicVEV`, `BosonicMatrixElement` | `BosonicAlgebra.m` | `FermionicVEV`, `FermionicMatrixElement` | Wick's theorem with signed pairings |
| `FieldVariables[var, labels]` | `Utils.m` | `FermionicFieldVariables[…]` | must be *distinguishable* by statistics — see §5 |
| `ToBosonicOperator[expr]` | `Utils.m` | `ToFermionicOperator[expr]` | realizes via JW rather than per-mode truncation |
| `CovarianceMatrix[state]` | `Utils.m` | `FermionicCovarianceMatrix[state]` | Majorana covariance: antisymmetric, not symplectic |
| `G2Coherence[state]` | `Utils.m` | reuse, or `FermionicG2` | `g²(0) = 0` exactly — antibunching |
| `WignerRepresentation`, `HusimiQRepresentation` | `PhaseSpaceRepresentations.m` | *(deferred, §8)* | requires Grassmann phase space |
| `$FockSize`, `SetFockSpaceSize` | `Utils.m` | **no global** — see §5 | mode count changes the Hilbert dimension; a silent default is a bug factory |

Three pieces have **no** bosonic counterpart and are genuinely new:

- **`JordanWignerTransform` / `BravyiKitaevTransform` / `ParityTransform`** — the bridge from fermionic operators to `QuantumOperator`, `PauliStabilizer`, and `QuantumCircuitOperator`. This is what connects the module to the rest of QuantumFramework.
- **A statistics-aware variable registry** — so one expression can hold both `\[FormalA]` (bosonic) and `\[FormalC]` (fermionic) and be reduced correctly.
- **`FermionicHamiltonian[h, V]`** — build `Σ hᵢⱼ cᵢ†cⱼ + Σ Vᵢⱼₖₗ cᵢ†cⱼ†cₖcₗ` from coefficient tensors, the entry point every electronic-structure workflow needs.

---

## 4. Competitive landscape

⚠️ **Every claim in this table is from recall and must be re-verified against current sources before it appears in any paper or announcement.** The bosonic README's §2 facts were verified by reading the sources directly; this table has not been.

| Tool | Has | Lacks |
|---|---|---|
| **OpenFermion** | `FermionOperator`, `MajoranaOperator`, JW / BK / parity transforms, normal ordering, `MolecularData` | symbolic closed forms; mixed boson–fermion algebra; Python-only |
| **Qiskit Nature** | `FermionicOp`, mappers, electronic-structure drivers | symbolic algebra |
| **PennyLane** `qml.fermi` | `FermiWord` / `FermiSentence`, `jordan_wigner`, `bravyi_kitaev`, `parity_transform` | symbolic reduction beyond word arithmetic |
| **SymPy** `physics.secondquant` | `CreateFermion` / `AnnihilateFermion`, `NO`, `wicks`, `contraction` — **the real symbolic competitor** | Gröbner-based reduction; mixed statistics; numeric bridge |
| **QuTiP 5** | numerics only; JW by hand | any fermionic algebra module |
| **SecondQuantizedAlgebra.jl** | fermionic support claimed — **verify level** | — |

**Where QuantumFramework can be uniquely strong:** Gröbner-based reduction over *mixed* statistics, closed-form lifting of the resulting series (the `ClosedFormBCH` machinery from the bosonic project applies unchanged), and a symbolic → numeric → circuit path in one kernel. Where it will not win: raw speed on large electronic-structure Hamiltonians. Concede that up front.

---

## 5. Design decisions to lock before writing code

These are the choices that are expensive to reverse. A recommendation is given for each.

**D1 — Separate symbols, not a `"Statistics"` option.**
`AnnihilationOperator[size, order]` takes a *truncation dimension*; a fermionic operator takes a *mode count and index*. The arguments mean different things, and the Hilbert space is built differently (JW strings versus independent tensor factors). Overloading one symbol would make `AnnihilationOperator[4]` ambiguous.
→ **Recommend:** `FermionicAnnihilationOperator[k, n]`. Verbose, unambiguous, greppable.

**D2 — No `$FermionModes` global.**
`$FockSize` is safe because changing a truncation only degrades accuracy. Changing the mode count changes the JW strings and therefore the *operator itself*; a stale global would silently produce a wrong-dimension operator that still composes.
→ **Recommend:** mode count always explicit. If ergonomics demand a default, make omission an error rather than a guess.

**D3 — Statistics must be recoverable from the variable.**
`ExtractNCVars` in `Utils.m` collects any formal symbol. A mixed expression needs to know which of them are fermionic.
→ **Recommend:** key statistics on the **base symbol name** — `\[FormalA]`-derived is bosonic, `\[FormalC]`-derived is fermionic — with a `FermionicVarQ` predicate mirroring the existing `FormalSymbolQ`, and `FermionicFieldVariables` defaulting to `\[FormalC]`. This reuses the existing machinery, keeps expressions self-describing, and adds no registry state that can go stale.

**D4 — Fix the JW convention and document it in every usage message.**
Choose `c_k = (∏_{j<k} Z_j) ⊗ σ⁻_k`, with the string running over *lower* indices (OpenFermion-compatible). This is the convention verified in §2.
→ Cross-check numerically against OpenFermion output for a 4-mode Hamiltonian as an acceptance test.

**D5 — Fix the Majorana convention.**
`γ_{2k-1} = c_k + c_k†`, `γ_{2k} = i(c_k† - c_k)`, giving `{γ_a, γ_b} = 2δ_ab`. State the factor explicitly, as `CovarianceMatrix` already does with its `"QuadratureScaling"` option, and offer the same escape hatch.

**D6 — Reference pages: match precedent, or break with it?**
No `SecondQuantization` symbol currently has a reference page — the context is documented only by [Tutorials/SecondQuantization.nb](../../../QuantumFramework/Documentation/English/Tutorials/SecondQuantization.nb), while 46 symbol pages exist for the rest of the paclet.
→ **Recommend:** write reference pages for the fermionic symbols and treat that as the pilot for backfilling the bosonic ones. Decide before P3, since it roughly doubles the documentation effort.

---

## 6. Track A — kernel additions

New files under [QuantumFramework/Kernel/SecondQuantization/](../../../QuantumFramework/Kernel/SecondQuantization/), mirroring the existing split. Effort assumes the ~1 hr/day envelope; phases are dependency-ordered, with calendar dates to be slotted by the owner.

### P1 — Vertical slice (proves the design end to end) · ~5 days

| # | Symbol | File | Signature |
|---|---|---|---|
| 1 | `FermionicFieldVariables` | `FermionicUtils.m` | `[n]`, `[var, labels]` → `{c₁, c₁†, …}`, plus `FermionicVarQ` |
| 2 | `FermionicRelations` | `FermionicAlgebra.m` | `[vars]` → anticommutators plus nilpotency, mirroring `BosonicRelations` |
| 3 | `FermionicNormalOrder` | `FermionicAlgebra.m` | `[expr, vars, Method -> "GrobnerBasis", "Scalars" -> {}]` |
| 4 | `FermionicAnnihilationOperator` | `FermionicOperators.m` | `[k, n]` → `QuantumOperator`, memoized as `AnnihilationOperator` is |
| 5 | `JordanWignerTransform` | `FermionicMappings.m` | `[expr, vars]` → `QuantumOperator`; the `ToBosonicOperator` analogue |

**Exit criterion:** the six reductions and four CAR checks in §2 pass as `.wlt` tests, and `JordanWignerTransform[FermionicNormalOrder[…]]` agrees with the directly built numeric operator to machine precision.

### P2 — Algebra completion · ~6 days

| # | Symbol | Notes |
|---|---|---|
| 6 | `FermionicAntinormalOrder` | needs the fermionic analogue of `ladderSwap`; verify the automorphism argument survives nilpotency |
| 7 | `FermionicVEV` | Wick's theorem with signed pairings — cross-check against Gröbner reduction on random monomials |
| 8 | `FermionicMatrixElement` | `⟨m|…|n⟩` in the occupation basis |
| 9 | **`NormalOrder` (mixed)** ★ | dispatches on `FermionicVarQ`; bosons commute with fermions. **The lead capability.** Budget generously. |

### P3 — Numeric layer · ~6 days

| # | Symbol | Notes |
|---|---|---|
| 10 | `OccupationState` | `[{0, 1, 1, 0}]`; must respect the D4 sign convention |
| 11 | `MajoranaOperators` | `[n]` → `2n` operators; assert `{γₐ, γᵦ} = 2δ` in tests |
| 12 | `FermionicHamiltonian` | `[h]`, `[h, V]` from coefficient tensors |
| 13 | `FermionicThermalState` | Fermi–Dirac occupations |
| 14 | `SlaterDeterminantState` | `[U]` from an isometry |
| 15 | `FermionicCovarianceMatrix` | antisymmetric Majorana covariance; `"Scaling"` option per D5 |

### P4 — Mappings and bridges · ~5 days

| # | Symbol | Notes |
|---|---|---|
| 16 | `BravyiKitaevTransform` | acceptance test: agrees with JW up to unitary equivalence on spectra |
| 17 | `ParityTransform` | cheapest of the three |
| 18 | `FermionicBeamSplitter` | Givens rotation; mirrors `BeamSplitterOperator[…, Method -> …]` |
| 19 | Circuit bridge | JW output → `QuantumCircuitOperator`; check interop with `PauliStabilizer` |

### P5 — Should-have model constructors · ~4 days, drop without regret

`FermiHubbardHamiltonian`, `KitaevChainHamiltonian`, `TightBindingHamiltonian`. High demo and test value — each has an analytically known spectrum to validate against — but not required for the module to be coherent.

---

## 7. Track B — tests and validation

Mirror the existing layout: `Tests/SecondQuantization/Fermionic*.wlt`, registered from [Tests/SecondQuantization.wlt](../../../Tests/SecondQuantization.wlt).

| File | Covers |
|---|---|
| `FermionicAlgebra.wlt` | the §2 reductions verbatim; `[nᵢ, nⱼ] = 0`; random-monomial Gröbner-versus-Wick agreement |
| `FermionicOperators.wlt` | full CAR at `n = 2…5`; nilpotency; `γ` relations; the `AnnihilationOperator[2]` anti-pattern from §2 |
| `FermionicMappings.wlt` | JW versus BK spectra; JW versus the hand-built `Fold` construction; OpenFermion cross-check on one fixed 4-mode Hamiltonian |
| `FermionicStates.wlt` | occupation-basis signs; Slater determinant norm and antisymmetry; Fermi–Dirac trace |
| `MixedStatistics.wlt` | boson–fermion commuting; a Holstein-type term reduced two ways |

**Physics acceptance tests** — each has a closed-form answer, so a failure is unambiguous:

1. **Tight-binding chain** — numeric spectrum versus `-2t cos(kπ/(n+1))`.
2. **Hubbard dimer** — exact ground-state energy versus the analytic two-site formula.
3. **Kitaev chain** — Majorana zero modes appear in the topological phase and not in the trivial one.
4. **Free-fermion Wick check** — a four-point function from `FermionicVEV` versus the Pfaffian of the covariance matrix.
5. **`g²(0) = 0`** for any single-mode fermionic state — the sharpest statistics contrast with the bosonic layer.

---

## 8. Explicitly deferred (scope-creep guard)

| Deferred | Why |
|---|---|
| **Grassmann phase space** — fermionic Wigner / Husimi | Requires an anticommuting-number datatype; research-grade, and no clean truncated-matrix analogue. Would double the project. |
| **Fermionic displacement / squeeze operators** | Grassmann-valued parameters; belongs with the above. |
| **Electronic-structure integrals** (`MolecularData` equivalent) | Needs a quantum-chemistry backend. Accept hardcoded published integrals for the H₂/STO-3G cross-check instead. |
| **Anyons / parastatistics** | Interesting, unbounded. |
| **Fermionic tensor networks** | Depends on the TensorNetworks paclet's graded-tensor support. Revisit after P4. |
| **Performance work on large Hamiltonians** | Conceded to OpenFermion in §4; optimize only if a real user hits it. |

---

## 9. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| Gröbner reduction blows up on high-degree multi-mode monomials | **Medium–High** | `FermionicVEV` via Wick's theorem is an independent path that avoids Gröbner entirely; benchmark both in P2 and route `FermionicNormalOrder` by degree. Nilpotency *helps* here — the ideal is much smaller than the bosonic one. |
| Mixed-statistics reduction (item 9) does not close cleanly | Medium | Pre-decide: if it is not working in 3 days, ship pure-fermionic P1–P4 and move mixed statistics to a follow-up. The module is coherent without it. |
| JW sign convention mismatch discovered late | Medium | D4 fixes it in writing *and* in an OpenFermion cross-check test landing in P1, not P4. |
| Naming churn after release | Low–Medium | §5 settles it before code. Note that `AnnihilationOperator` is already exported and unqualified — the fermionic name cannot be retrofitted onto it later without a breaking change. |
| Documentation effort doubles under D6 | Medium | Decide D6 before P3. Tutorial-only is a legitimate answer that matches current precedent. |

---

## 10. Success criteria

**Must-have**

- All CAR relations verified numerically for `n = 2…5`, and symbolically via Gröbner reduction.
- `FermionicNormalOrder` and `FermionicVEV` agree on random monomials — two independent methods, one answer.
- `JordanWignerTransform` reproduces a published 4-mode Hamiltonian's Pauli decomposition exactly.
- The five physics acceptance tests in §7 pass against their closed forms.
- Every new symbol carries a usage message in the existing `SecondQuantization` style.
- A user who knows the bosonic layer can use the fermionic layer without reading new documentation.

**Should-have**

- Mixed boson–fermion normal ordering (item 9).
- `BravyiKitaevTransform` agreeing with JW up to unitary equivalence.
- H₂/STO-3G ground state from hardcoded integrals matching the literature value.
- Reference pages (D6).

---

## 11. One-line summary

Mirror the bosonic layer for fermions — the algebra engine and the tensor machinery both already work, as §2 verifies — then use Jordan–Wigner as the bridge that connects second quantization to the rest of QuantumFramework, and treat mixed boson–fermion normal ordering as the one capability no competing tool has.
