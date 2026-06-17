---
name: qf
description: >
  Deep working skill for the Wolfram QuantumFramework paclet. Trigger when the user writes,
  modifies, debugs, or reviews code using QuantumFramework symbols (QuantumState, QuantumOperator,
  QuantumBasis, QuantumCircuitOperator, QuantumChannel, QuantumMeasurement(Operator),
  QuantumPartialTrace, QuantumDistance, QuantumEntangledQ, QuantumEntanglementMonotone,
  QuantumEvolve, QuantumLinearSolve, FockState, AnnihilationOperator, CoherentState,
  GenerateParameters, ParametrizedLayer, ClassiqSetup, PauliStabilizer, StabilizerStateQ,
  StabilizerFrame, CliffordChannel, GraphState, LocalComplement, QuantumQASM, IBMJob,
  IBMJobSubmit, QiskitTarget, QuantumAdiabaticEvolve, the QuantumOperator KAK / Decompose
  properties, etc.), in any project. Encodes hard-won pitfalls from a complete kernel + doc +
  notebook audit. A bundled reference/ set (function-index, mistakes, idioms-and-performance,
  architecture, qf-tn-integration) carries the detail; reference/VERSION.md records the paclet
  version and commit the file:line cites were verified against. Prefer this over a generic
  quantum-computing skill whenever a QuantumFramework symbol is involved.
---

# QuantumFramework: Deep Working Skill

Anchored to a complete read-through of all kernel `.m` files, doc files, and supplementary
notebooks in the **QuantumFramework paclet** (`WolframResearch/QuantumFramework`). Captures
non-obvious facts you'd otherwise re-discover the painful way. Paclet version 2.0.0; see
[reference/VERSION.md](reference/VERSION.md) for the exact commit the `file:line` cites were
verified against and how to re-check them on your install.

## 0. The bundled reference set: read it, and read it first

The deep reference is the **active set** bundled with this skill under `reference/`. Read the
relevant file(s) before answering a non-trivial QF question (no sampling):

- [reference/function-index.md](reference/function-index.md): every public function with `file:line`
- [reference/mistakes.md](reference/mistakes.md): ~90 documented pitfalls with citations
- [reference/idioms-and-performance.md](reference/idioms-and-performance.md): fast paths, anti-patterns, notebook-confirmed workflows
- [reference/architecture.md](reference/architecture.md): object map, property dispatch, dimension semantics
- [reference/qf-tn-integration.md](reference/qf-tn-integration.md): the TensorNetworks bridge and the default circuit-application path
- [reference/VERSION.md](reference/VERSION.md): the paclet version + commit the cites were verified against, and how to re-check them

**The `file:line` cites are point-in-time** for the version in `reference/VERSION.md`. The
paclet on the user's machine may differ, so when an exact cite matters, re-verify it against
their own kernel (read the `.m` source, or test with `wolframscript`) before relying on it. If
the kernel has moved, the *behavior* notes remain the best guide; treat the line numbers as
approximate. The skill body below is self-contained, so it is usable even without opening the
reference files.

### Where signatures and behavior actually live

1. **`Kernel/<file>.m`**: the implementation, the only contract.
2. **Doc pages** (`Documentation/English/`): intent and examples; many symbol reference pages
   are 2022-era stubs with `XXXX` placeholders, but the Tutorials are substantive
   (`Tutorial.nb`, `Bellstheorem.nb`, `SecondQuantization.nb`, `QuantumOptimization.nb`,
   `TensorNetwork.nb`, `QPUServiceConnect.nb`).
3. **A live test** via `wolframscript -file` (never `Quiet`, never `Print` in committed code).

**Never read `Kernel/Usage.m` for what a function does.** It is auto-generated from the doc
pages (serialized `BoxData` usage boxes, not human-authored source). It is not evidence of
behavior, signatures, options, or return shapes, and hand edits to it are overwritten on
regeneration. Usage-text fixes go in the doc page's Usage section. This repeats the strict rule
in the repo's `CLAUDE.md` because the instinct to grep it dies hard.

## 1. Loading

```wolfram
Needs["Wolfram`QuantumFramework`"]
```

Sub-contexts (loaded by `PacletInfo.wl`, but auto-completion only kicks in after `Needs`):

- `Wolfram`QuantumFramework`SecondQuantization``: Fock-space, bosonic quantum optics
- `Wolfram`QuantumFramework`QuantumOptimization``: VQE, QAOA, QNGD, VQLS, `QuantumLinearSolve`, `QuantumAdiabaticEvolve`
- `Wolfram`QuantumFramework`Experimental``: `QuantumMPS`, `QuantumWignerTransform` (EXPERIMENTAL API)
- `Wolfram`QuantumFramework`DiagramPlot``: circuit diagram rendering
- `Wolfram`QuantumFramework`Gates``: gate library
- `Wolfram`QuantumFramework`Init``: boot

A bare `Needs` leaves the optimization symbols (`ParametrizedLayer`,
`QuantumNaturalGradientDescent`, parameter-shift, `QuantumLinearSolve`) inert in `Global``
unless the sub-package context is on the path; load the sub-package before concluding QF
"lacks" a feature.

Note the separately-installed paclet can lag the repo HEAD. To test against a local checkout of
the working tree, point `PacletDirectoryLoad` at the inner `QuantumFramework` paclet directory:

```wolfram
PacletDirectoryLoad["/path/to/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]
```

## 2. The mistakes that bite first

If you do nothing else, internalize these. (Line numbers are point-in-time, see
[reference/VERSION.md](reference/VERSION.md); ~90 more entries in
[reference/mistakes.md](reference/mistakes.md).)

### 2.1 Pre-2.0 list-form gate specs are dead, and they fail *silently*
- `QuantumOperator[{"R", θ, "XY"}]` does not error: it misparses the list as a state vector and
  builds a $d = 3$ "operator" whose amplitudes are the symbols `{R, θ, XY}`. Verified live.
- `QuantumOperator[{"RX", θ}]` is equally dead.
- Always use the call-form: `"R"[θ, "XY"]`, `"RX"[θ]`, `"U3"[θ, φ, λ]`, `"C"["1"]`.
- Any pre-2.0 doc cell or notebook using list-forms must be rebased before trusting it.

### 2.2 `"SX"` is a 2-qubit operator ($S \otimes X$), not $\sqrt{X}$
- For $\sqrt{X}$ use `"V"` (arity 1). `QuantumOperator["SX"]["Arity"]` is 2. Verified live.
- Related, from the named-operator audit: `"RX"`/`"RY"`/`"RZ"` with qudit dimension $d > 2$ are
  non-unitary; lowercase gate aliases are dead; `"Z"[d]` is the inverse clock convention
  relative to `"HeisenbergWeyl"`; `"Deutsch"` uses the half-angle convention.

### 2.3 `qs["VonNeumannEntropy"]` and `qs["Purity"]` silently return `Indeterminate` for $\ge 10$ qubits
- Hardcoded `ConfirmAssert[qs["Dimension"] < 2^10]` at `QuantumState/Properties.m:461, :481`.
- For larger systems: bipartition first via `QuantumPartialTrace`, then compute on the
  reduction; or use `qs["TraceNorm"]` (no dimension limit); or route through `QuantumMPS`.

### 2.4 `QuantumDistance` with no measure returns $1 - F$, not the fidelity $F$
- `QuantumDistance.m:16`: default is $1 - \mathrm{Re}\,\mathrm{Tr}\,(\rho_1 \rho_2)^{1/2}$.
- For raw fidelity call `QuantumDistance[qs1, qs2, "Fidelity"]`, or use `QuantumSimilarity`
  (exported, defined at `QuantumDistance.m:63`).

### 2.5 `QuantumPartialTrace` densifies sparse density matrices
- `MatrixPartialTrace` (`Utilities.m:135`) uses `ArrayReshape` + `TensorContract`, which always
  materializes the dense $2^N \times 2^N$ tensor.
- For large $N$: route through `QuantumMPS` (Experimental) or compute the reduced state from a
  bipartition directly.

### 2.6 `eigensystem` is generic: no Hermitian path
- QF wraps `Eigensystem` with `Simplify` and `Chop` (`Utilities.m:154`), no
  `Method -> "Hermitian"`.
- For pure-numeric Hermitian matrices, drop to
  `Eigensystem[qs["DensityMatrix"], Method -> "Hermitian"]` directly.

### 2.7 `Inverse` silently falls back to `PseudoInverse`
- `MatrixInverse` (`Utilities.m:309-312`) wraps in `Quiet` + `Check`, then `PseudoInverse`.
  **No warning.**
- For potentially singular matrices, check `MatrixRank` first or call
  `Inverse[m, ZeroTest -> ...]` explicitly.

### 2.8 The Magic circuit injects an $i$ global phase on 2 of 4 mappings
- `QuantumCircuitOperator["Magic"][QuantumState["10"]]` returns $i\,|\psi^+\rangle$, not
  $|\psi^+\rangle$. Verified live.
- Mapping table: $|00\rangle \to |\Phi^+\rangle$, $|10\rangle \to i|\psi^+\rangle$,
  $|01\rangle \to i|\Phi^-\rangle$, $|11\rangle \to |\psi^-\rangle$.

## 3. Construction patterns

### 3.1 States

```wolfram
QuantumState[{α, β}]                         (* 1 qubit, computational *)
QuantumState["PhiPlus"]                       (* Bell state *)
QuantumState["GHZ"[5]]                        (* 5-qubit GHZ *)
QuantumState["RandomPure"[3]]                 (* 3-qubit random pure *)
QuantumState["RandomMixed", {3, 2}]           (* random mixed bipartite 3⊗2 *)
QuantumState["Werner"[p, 2]]                  (* Werner state *)
QuantumState["Dicke"[{2, 1, 2}]]              (* Dicke generalized *)
QuantumState["BlochVector"[{rx, ry, rz}]]     (* qubit from Bloch vector *)
QuantumState["Graph"[g]]                      (* graph/cluster state *)
QuantumState["010"]                           (* bit-string *)
QuantumState["+0L"]                           (* combined: + ⊗ 0 ⊗ L *)
```

### 3.2 Operators (108 named operators in `$QuantumOperatorNames`)

```wolfram
QuantumOperator["X"[3]]                       (* generalized Pauli-X in d=3 *)
QuantumOperator["U3"[θ, φ, λ]]                (* IBM U3 gate *)
QuantumOperator["R"[θ, "XY"]]                 (* generic R(θ, X⊗Y) two-qubit *)
QuantumOperator["RandomUnitary"[{3, 2}]]      (* random unitary 3⊗2 *)
QuantumOperator["Hamiltonian"[H, {L1, L2}, {γ1, γ2}]]   (* Lindblad form *)
QuantumOperator["Trace"[d]]                   (* partial-trace operator *)
QuantumOperator["ZSpider" -> {1,2} -> {1,2,3}]  (* spider with explicit i/o legs *)
```

Order: `{output, input}` (rows-then-columns) OR rule `{input} -> {output}`. The two are
equivalent. Call-form arguments only; see mistake 2.1.

### 3.3 Circuits (44 named circuits in `$QuantumCircuitOperatorNames`)

```wolfram
QuantumCircuitOperator[{
  "X" -> 2, "Barrier"[;;3], "Y" -> 1, "Z" -> 3,
  "CNOT" -> {1, 2},
  "R"[θ, "ZZ"] -> {2, 3},
  "BitFlip"[p] -> {1},
  {1, 2}                           (* sequential measurements on qubits 1, 2 *)
}, "Parameters" -> {θ, p}]
```

Named circuits include `"Bell"`, `"Toffoli"`, `"Magic"`, `"Cup"`, `"Cap"`, `"Graph"[g]`,
`"Grover"[bool]`, `"Fourier"[n]`, `"PhaseEstimation"[u, n]`, `"BernsteinVazirani"[s]`,
`"Simon"[s]`, `"DeutschJozsa"[k, n]`, `"Trotterization"[ops, order, steps, time]`,
`"PhaseNumber"[n, m, boo]`, `"CHSH"`, `"GrayOracle"[θs, prec, target]`,
`"Multiplexer"[op1, op2, ...]`, `"QuantumState"[state]`, `"StatePreparation"[state]`; the full
catalog with per-circuit signatures is in [reference/function-index.md](reference/function-index.md).

Apply: `qc[QuantumState["00"]]`, or `qc[]` to apply to the default register state
(`QuantumCircuitOperator["Bell"][]` returns the $|\Phi^+\rangle$ `QuantumState`).

### 3.4 Stabilizers

`PauliStabilizer` accepts the same named call-forms (`PauliStabilizer["GHZ"[5]]`,
`PauliStabilizer["Random"[n]]`; 9 names in `$PauliStabilizerNames`). The fast path for long
Clifford circuits is `ps["ApplyCircuit", qc]` / `PauliStabilizerApply` (compiled generator-major
gate fold; orders of magnitude faster than gate-by-gate application).

## 4. Circuit application paths: TensorNetwork is the default

Circuit application routes through a tensor-network contraction by default
(`TensorNetwork.m`); `Method -> "Schrodinger"` forces the dense statevector path. Details and
the bridge to the TensorNetworks paclet: [reference/qf-tn-integration.md](reference/qf-tn-integration.md).

- Contraction-path control: `qc[qs, Method -> {"TensorNetwork", "Path" -> spec}]`
  (`TensorNetwork.m:182, 215`).
- **Trap**: for QPE-style circuits with many controlled powers (e.g. HHL at precision $\ge 5$),
  the default TN path bloats the output qudit count and can exhaust memory; the workaround is
  `Method -> "Schrodinger"`.
- `Method -> "Schrodinger"` is also the right choice for small circuits and fidelity
  benchmarking, where TN overhead dominates.

## 5. Performance idioms

### 5.1 Caching

QF has **two independent global caches**:

| Cache | Where | What | Bypass |
|---|---|---|---|
| `$QuantumFrameworkPropCache` | `QuantumFramework.m:15` | Property lookups | Set `False`. State cache also bypassed when `ParameterArity > 0`. |
| `$QuditBasisCache` | `QuditBasis/NamedBases.m:27` | Named-basis construction | Set `$QuditBasisCaching = False`. Built-in bypass for `Random*MIC`. |

Setting one to `False` does NOT clear the other. For long sessions:

```wolfram
Wolfram`QuantumFramework`PackageScope`$QuditBasisCache = <||>
```

### 5.2 Parametric workflows

Always declare `"Parameters" -> {...}` at construction time. The framework propagates
parameters through compile / measurement / evolution. Prefer:

```wolfram
qc = QuantumCircuitOperator[{...}, "Parameters" -> {θ, φ}];
qc[<|θ -> π/3, φ -> π/4|>]                    (* substitute symbolically *)
qc[θval, φval]                                (* substitute positionally *)
```

over rebuilding the circuit each iteration. The cache is bypassed when `ParameterArity > 0`,
which costs less than re-instantiation.

### 5.3 VQE / QAOA cost-function pattern

```wolfram
ansatz = QuantumCircuitOperator[{...}, "Parameters" -> {θ1, θ2, ...}];
|φ⟩ = ansatz[];
⟨φ| = |φ⟩^†;
cost[θ1_, θ2_, ...] = ("Scalar" //. ⟨φ|@H@|φ⟩) // ComplexExpand // FullSimplify;
NMinimize[cost[θ1, θ2, ...], {θ1, θ2, ...}]
```

**Use `NMinimize`, not `FindMinimum`**: VQLS-like cost surfaces have wide barren plateaus that
break local-derivative methods ([reference/mistakes.md](reference/mistakes.md) has the demo).

### 5.4 Quantum Natural Gradient Descent

When the cost surface has near-singular regions (barren plateaus near $\theta = 0$), QNGD with
the Fubini-Study metric tensor converges in ~2 steps where vanilla gradient descent needs 50+:

```wolfram
metric = FubiniStudyMetricTensor[ψ];   (* use the variational state ψ, NOT the ansatz *)
qopt = QuantumNaturalGradientDescent[cost, metric, "InitialPoint" -> init,
                                     "LearningRate" -> 0.05, "MaxIterations" -> 50];
```

### 5.5 Solving linear systems

For `m . x = b` with `Length[b]` a power of 2:

```wolfram
QuantumLinearSolve[m, b]                     (* result *)
QuantumLinearSolve[m, b, All]                (* full Association: Result, Ansatz, CircuitOperator, GlobalPhase, OptimizedParameters, Parameters *)
QuantumLinearSolve[m, b, Method -> "DifferentialEvolution"]   (* alt heuristic *)
```

Internally uses `QuantumOperator[m]["PauliDecompose"]` + `"Multiplexer"`. The cost function
constrains only $|\langle b|\psi\rangle|^2 = 1$, so the result is correct **up to a global
phase**; the function post-processes via `globalphase = ψopt / (b/Norm[b])` to recover
amplitudes.

## 6. Phase space (Wigner): parity matters

`QuantumWignerTransform` is EXPERIMENTAL (per its doc page). For shipping code, use the
equivalent stable form:

```wolfram
QuantumState[ρ["Double"], QuantumBasis["Wigner"[d]]]
```

Phase-space dimensions depend on parity:
- **Odd $d$**: `qs["PhaseSpace"]` is $d \times d$, equal to `qs["QuasiProbability"]`.
- **Even $d$**: `qs["PhaseSpace"]` is $2d \times 2d$. Only the upper-left
  `[[1;;d, 1;;d]]` quadrant is the actual Wigner state vector. To recover position
  probabilities from the full matrix: `Total[qp][[;;;;2]]`.

Multi-qudit: `qs["PhaseSpace"]` is a tensor with paired dims $\{2 d_i, 2 d_i\}$ (even) or
$\{d_i, d_i\}$ (odd).

## 7. Hardware backends and interop

```wolfram
ibmq = ServiceConnect["IBMQ"]
ibmq["BackendQueue"]
ibmq["Backend", "ID" -> ibm_brisbane, "Property" -> "Configuration"][basis_gates]

(* Synchronous (kernel blocked): *)
qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane}]

(* Asynchronous: *)
LocalSubmit[
  Needs["Wolfram`QuantumFramework`"];
  qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane}],
  HandlerFunctions -> <|"TaskFinished" -> ((result = #EvaluationResult)&)|>
]
```

> **WARNING**: For $\ge 10$ qubits using `Method -> "Qiskit"`, the return is a raw counts
> `Association[bitstring -> count, ...]`, **not** a `QuantumMeasurement`. Test with
> `Head[result]` before further processing.

Interop hub facts (synced 2026-06-11):
- `QuantumQASM` is the native-WL OpenQASM import/export (no Python). The emitter rebases
  non-computational (X/Y) measurement bases; the importer handles physical qubits (`$N`),
  controlled `gphase`, scientific-notation angles, and drops `reset`.
- `IBMJob` splits Properties from Actions: `"ExecutedCircuit"`, `"Refresh"`, `"Cancel"` are
  Actions. `"ExecutedCircuit"` cannot work for QF-submitted jobs: IBM stores no
  transpiled-circuit object for ISA + V2 sampler submissions, so the fetch 404s by design.
- `IBMJobSubmit` pre-validates primitive options (`ibmValidatePrimitiveOptions`) and decodes
  `PrimitiveResult` natively.

## 8. Documented-but-nonexistent symbols

The Guide (`Documentation/English/Guides/WolframQuantumComputationFramework.nb`) lists three
symbols that do not exist in any kernel file:

| Documented (broken) | Real equivalent |
|---|---|
| `QuantumSchmidtDecomposition` | `qs["SchmidtDecompose"]` / `qs["SchmidtBasis"]` |
| `QuantumSpectralDecomposition` | `qs["SpectralBasis"]` / `qo["Diagonalize"]` |
| `QuantumCompile` | `qco["Compile"]` / `qco["QuantumOperator"]` |

Don't recommend the broken symbols: they 404 in user docs.

## 9. Known kernel filename typo

`Kernel/QuantumHamiltonianOperator/QuantumHamiltonialOperator.m` is missing an `n`
("Hamiltonial"). Still present at HEAD `6f953ad5`. The previously catalogued typos
(`QuantumOperartorQ`, `__U2Gate`, the `RZXGate` YX body, `QuatumPartialTranspose.m`) are all
fixed; don't cite them.

## 10. When in doubt

1. `?Wolfram`QuantumFramework`Symbol*` prints all definitions on the symbol.
2. `Symbol // InputForm` shows the internal Association structure.
3. `Symbol[arg]["Diagram"]` is a visual sanity check.
4. `obj["Properties"]` is the head's own property list; `obj["AllProperties"]` is the full
   union including delegated ones (e.g. `QuantumOperator`: 67 vs 209). Pre-2.0.0 code that
   expected the union under `"Properties"` will under-count.
5. [reference/function-index.md](reference/function-index.md) has the `file:line` for every implementation;
   [reference/mistakes.md](reference/mistakes.md) covers most footguns; [reference/idioms-and-performance.md](reference/idioms-and-performance.md) has
   notebook-confirmed patterns (Lindblad forms, expectation values four ways, hardware chains).
6. When uncertain about behavior, run a tiny example via `wolframscript -file` in Bash (load
   the working tree with `PacletDirectoryLoad` if the change is unreleased). Never assert a
   signature from memory, and never from `Usage.m` (§0).

## 11. Error handling and silent fallbacks

QF declares ~50 `Head::tag` messages across the kernel. Tags form a small reused vocabulary;
every public head with a `["Property"]` API has `undefprop` and `failprop`. The full strategy
(layer selection, body grammar, naming) lives in `wl-style-guide`; this section catalogs what
QF actually does so you can match the local style when editing.

### Shared tag vocabulary

| Tag | Heads using it |
|---|---|
| `invalidName` | `QuantumState`, `QuantumOperator`, `QuantumCircuitOperator`, `QuantumChannel`, `QuantumMeasurementOperator`, `QuantumHamiltonianOperator`, `QuditBasis` |
| `invalidArgs` | catch-all on the named heads (added in the 2.0 refactor) |
| `undefprop` + `failprop` | the seven above plus `QuantumBasis`, `QuantumMeasurement` |
| `dim` / `dimmismatch` | `QuantumCircuitOperator`, `CliffordChannel`, `ParametrizedLayer`, `QuantumLinearSolve` |
| `inconsistent*`, `wrongData`, `dependentElements`, `zeroDimension` | `QuantumBasis`, `QuditBasis` (constructor invariants) |
| `invalidorder` | `DisplacementOperator`, `SqueezeOperator` |
| `notpure` / `nonclifford` / `nonpaulibasis` / `nongraph` | regime mismatches in `QuantumState`, `PauliStabilizer`, `GraphState` |
| `pkg` / `link` / `PythonEvaluators` | external dependencies (`QuEST`, `ClassiqSetup`, `PythonTools`) |

Most files use a top-of-file message block (e.g. `QuantumBasis.m`,
`Stabilizer/PauliStabilizer.m`). Older files (`SecondQuantization.m`) co-locate; both styles
work but the block scales.

### Property dispatch (`QuantumState/Properties.m`, near the top)

```wolfram
QuantumState::undefprop = "property `` is undefined for this state"
QuantumState::failprop  = "property `` failed with ``"

(qs_QuantumState[prop_ ? propQ, args___]) /; QuantumStateQ[Unevaluated[qs]] :=
  With[{result = QuantumStateProp[qs, prop, args]},
    ...cache logic...
    /;
      (!FailureQ[Unevaluated @ result] || Message[QuantumState::failprop, prop, result]) &&
      (!MatchQ[Unevaluated @ result, _QuantumStateProp] || Message[QuantumState::undefprop, prop])
  ]
```

Two short-circuit checks emit a tagged message and reject the rule when the property table
didn't match. The same shape is repeated in `QuantumOperator/Properties.m`,
`QuantumBasis/Properties.m`, `QuantumCircuitOperator/Properties.m`, etc.; grep for `undefprop`
to find all copies.

### Silent fallbacks to watch for

QF deliberately stays quiet in places where you'd want a message (see also §2.7):

- `MatrixInverse` (`Utilities.m:309-312`): `Quiet[Check[Inverse[m], PseudoInverse[m]]]`.
  Singular matrices silently become pseudo-inverse.
- `eigensystem` (`Utilities.m:154`): generic `Simplify` + `Chop` wrapper, no
  `Method -> "Hermitian"` path.
- `quantumStateQ[___] := False` and friends: no `badform` message on a malformed
  `QuantumState[...]`; the only check is construction-time `ConfirmBy`.
- Many `_/; predicate` overloads reject silently to fall through to the next overload.
- A `CompiledFunction` (e.g. the stabilizer gate fold) silently drops to interpreted
  evaluation when handed unpacked arrays; benchmark with packed inputs.

When adding new public functions, follow `wl-style-guide`'s rules, and don't wrap a noisy
`Confirm` in `Quiet` to silence it.

## 12. Common idiomatic patterns (notebook-confirmed)

### Composing circuits

```
qc1 /* qc2 /* qc3                          (* left-to-right composition *)
QuantumMeasurementOperator[{1,2}] @* qc1   (* right-to-left, append measurement *)
```

### Expectation value (4 ways)

```
Tr[op["Matrix"] . qs["DensityMatrix"]]                           (* trace *)
QuantumPartialTrace[op[qs["Operator"]]]["Scalar"]                (* partial trace *)
QuantumMeasurementOperator[op][qs]["Mean"]                       (* measurement *)
(qs^† @ op @ qs)["Scalar"]                                       (* bra-op-ket *)
```

### Inner product / fidelity

```
(qc /* {ψ^†})[]["Norm"]^2                                        (* |⟨ψ|ψf⟩|² *)
QuantumDistance[qs1, qs2, "Fidelity"]                            (* same as ⟨ψ1|ψ2⟩ for pure *)
```

### Lindblad master equation

```
(* Use any of these: all equivalent *)
QuantumEvolve[H, {L1, L2}, {γ1, γ2}, ρ0, {t, 0, T}]            (* rates separate *)
QuantumEvolve[H, {Sqrt[γ1] L1, Sqrt[γ2] L2}, ρ0, {t, 0, T}]    (* rates folded *)
ℒ = QuantumOperator["Liouvillian"[H, {Ls}, {γs}]];
ρt = Exp[ℒ t][ρ0]                                              (* time-independent shortcut *)
```

### Hardware run patterns

```
(* IBMQ via QPY *)
qpy = BaseEncode @ qc["Qiskit"]["QPY", "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane];
ibmq["RunCircuit", {"QPY" -> qpy, "Backend" -> ibm_brisbane}]

(* IBMQ async (kernel free) *)
LocalSubmit[Needs["Wolfram`QuantumFramework`"];
  qc[Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> ibm_brisbane}],
  HandlerFunctions -> <|"TaskFinished" -> ((result = #EvaluationResult) &)|>]

(* AWS Braket *)
braket["CreateQuantumTask",
  "DeviceArn" -> deviceArn,
  "Action" -> <|braketSchemaHeader -> <|name -> "braket.ir.openqasm.program", version -> "1"|>,
                source -> qc["Qiskit"]["QASM", "Provider" -> "AWSBraket"]|>,
  "Shots" -> 100, "OutputS3Bucket" -> bucket]
```

### Comparing exact vs hardware results

```
quantum_pred = qc[]["ProbabilitiesList"];
qpu_pred = N @ Normalize[Values @ Normal @ counts, Total];
BarChart[Transpose[{quantum_pred, qpu_pred}]]
```

### Programmatic Hamiltonian for Max-Cut

```
index = MapApply[List, EdgeRules[g]] /. x_Integer :> {x, x};
op[I] = StringRepeat["I", VertexCount[g]];
sum = Total[(1/2)(op[I] - StringReplacePart[op[I], "Z", #]) & /@ index];
Hc = QuantumOperator[sum]
```

### Counts simulation (mimics QPU stats)

```
counts = Counts[exp["SimulatedMeasurement", 1024]]   (* exact dist sampled, not noise *)
```

### Other commonly used non-doc-listed idioms

- `qc[<|θ -> π/3|>]` substitutes parameters via `Association` (preferred over positional when
  there are many parameters).
- `qs["Number"]` extracts the complex scalar from a 0-qudit state object (e.g., output of an
  inner product).
- `qc[Method -> "Schrodinger"]` forces the dense statevector path (vs default tensor network);
  see §4 for when this is mandatory.
- `qc[]["Reverse"]` reverses qudit order in a measurement output (matches IBM's little-endian).
- `qc[[2;;-2]]` slices a `QuantumCircuitOperator` like a list: drops state-prep + measurement
  bookends.
- For analog/Rydberg: use `AquilaHamiltonian[sites, {Ω, Δ, φ}]` (custom function from
  `Notebooks/QuEra_WhitePaper.nb`).
