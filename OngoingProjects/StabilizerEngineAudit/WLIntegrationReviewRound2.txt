# ROUND: 2/3 — independent stabilizer WL correctness and quality review

## VERDICT

The cumulative optimization is a coherent use of Wolfram Language and QuantumFramework's existing architecture. The storage, dispatch and execution layers are meaningfully separated; the integer kernel is justified, and the new interpreted code uses native expression/array operations. Exact Clifford and inverse behavior, symbolic binary signs, storage equivalence, matrix-honest circuit plans, and the optimized long-circuit scale were reproduced.

One inherited canonical-association validation gap remains in the newly claimed checked-execution boundary. It needs a small stabilizer-only fix before claiming that every malformed public representation is stopped before the native fold. This finding predates the optimization: I reproduced the same boundary escape on the original snapshot. It is counted because the changed validator and checked boundary still rely on the inconsistent public property.

No target sources or tests were edited. All review evidence was written under `/private/tmp`. Main checkout: `/Users/mohammadb/Documents/GitHub/QuantumFramework`; HEAD `708dc533489f0b7d7a0f92b412c380ce137ff25d`. Reviewed bytes are recorded in `stabilizer-style-r2-manifest.json`.

## REPRODUCED

Runtime: Wolfram Language 15.0.1, MacOSX-ARM64; working-tree QuantumFramework 2.1.3, loaded from main's `QuantumFramework/Kernel/QuantumFrameworkMain.m`. The old audit anchor is `f77a6640`; implementation/signature claims below use current source and live output, not stale audit line references.

### Claims and evidence

1. **Full exact Clifford X/Z/sign action and inverse.** Independently enumerated the 24 one-qubit Clifford actions from exact H/S matrices, keeping the corresponding words. Every forward and inverse X/Z image, both composition identities, and inverse involution passed across canonical, chunk and generator storage. Separately tested 12 random entangling 27-gate circuits on 2, 3 and 4 qubits against independently constructed System matrix embeddings, with all forward/inverse generator images and exact unitarity. The source implementation of the inverse residual phase is `GateUpdates.m:152–180`, with phase composition in `Compose.m:29–53`.
2. **Physical V/V-dagger.** Exact matrix conjugation agrees with scalar V and its adjoint. The H/S/H and H/S-dagger/H definitions (`GateUpdates.m:125–126`) match the compiled encoder and the exact matrices in the frame dependency; no unintended register growth occurred.
3. **Resident storage and property equivalence.** Independent scalar AG recurrences and split evolution agree with bulk evolution at n=1,2,5,31,32,61,62,63,124,125,1024. These recurrences use ordinary bit arrays and arithmetic modulo two, not production packed kernels. Canonical/chunk/generator tableaux, signs, phase, matrix, stabilizer/destabilizer slices, Pauli symbols, entropy, expectation and measurement keys agree. Materialized state vectors normalize exactly. `CompiledEvolution.wlt` additionally reproduces the 1,024-qubit 1,000-gate scalar comparison.
4. **Symbolic binary phase semantics.** Independent exact symbolic phases `1 - 2 phaseBit[j]` were retained through a nontrivial two-qubit circuit and inverse. All 16 binary assignments agree with concrete evolution/inversion. Symbolic signs remain outside the machine-integer path. The focused suite also covers symbolic measurement/sampling and the non-Clifford frame boundary.
5. **Immutable validity and ordinary malformed inputs.** Mutation of a valid expression does not preserve validity incorrectly. Missing packed fields, malformed dimensions, illegal word values and sign lengths are rejected in ordinary packed layouts. Canonical signed tableaux and symbolic signs keep the canonical fallback. There is one exception involving inconsistent stored canonical `"Qubits"`, detailed below.
6. **Wire boundaries and exact integer encoding.** Focused tests cover scalar and bulk lower/upper bounds, equal two-qubit wires, and ±2^63 wire values. Extreme encoder outputs stay exact; out-of-register values do not enter the native fold in those ordinary layouts. Empty evolution preserves the object. Additional live reversed-wire QCO application matches direct CNOT action.
7. **Circuit plans and encoding caches.** All focused bound, eviction, clearing and injected hash-bucket-collision checks passed. The original circuit/specification is retained and compared on a hit (`GateMatrix.m:503–513`, `Compiled.m:217–230`), and limits are 16 entries. Matrix-honesty, wire placement and parameter-substitution tests passed. An additional live computational H followed by H stored in a PauliX basis correctly rejected the latter with `PauliStabilizer::nonclifford`; a previous computational cache entry did not change that result. Existing lower-level gate caches retain their existing design; the collision claim here concerns the two added plan/encoding caches.
8. **Current signatures and compile targets.** Live contexts for `PauliStabilizerApply`, `psFromGeneratorPacked`, `applyCompiledFold`, `stabilizerExecutionPlan` are the QF PackageScope context. `Options[PauliStabilizerApply]` is `{"Compress" -> None}`. Object properties and return heads match the source. C was used for the primary run. In a separate fresh kernel, forced WVM compilation passed a 63-qubit 300-gate circuit with V, adjoints and entangling gates against scalar application; the resulting object was a CompiledFunction.
9. **Performance retention.** Five clean serialized trials per workload, same driver and prepared inputs for frozen pre-style and current source; startup, circuit preparation and first compilation excluded. `ClearSystemCache[]` precedes every trial; cold cases additionally clear the stabilizer caches. Every result head is PauliStabilizer, all trial message lists are empty, and result digests match between versions. This supports retention of the previous optimized scale, not uniform improvement or a new speedup claim.

### Independently executed focused suite

- `CompiledEvolution.wlt`: 65/65 success.
- `InversePhase.wlt`: 330/330 success.
- `ExecutionContracts.wlt`: 65/65 success.
- Independent grouped checks: 23/23 success.
- Independent WVM fallback check: success.

The parent-run full-suite JSON was also inspected directly: 1,592/1,592 tests across 15 files, no failures and `{}` outside-message list. That full-suite execution is parent-run evidence; the independent execution above is mine.

### Independent timing results

Values are milliseconds: median [minimum, maximum], five trials each.

| Workload | Frozen pre-style | Current |
|---|---:|---:|
| Cold encoding, n=128, 10,000 gates | 9.174 [8.816, 9.285] | 8.961 [8.781, 9.107] |
| Warm direct, n=1024, 10,000 gates | 3.625 [3.413, 4.562] | 4.105 [3.530, 5.120] |
| Ten public calls, n=1024, 1,000 gates/call | 15.277 [14.646, 15.898] | 14.260 [13.871, 16.035] |
| V-containing direct, n=128, 10,001 gates | 2.853 [2.833, 3.029] | 2.804 [2.740, 2.844] |
| Cold encoding, n=1024, 10,000 gates | 9.778 [9.494, 9.846] | 9.860 [9.691, 10.304] |

Warm direct's current median is 13.2% higher in this sample, with overlapping ranges; the other representative medians range from 6.7% lower to 0.8% higher. Do not characterize these measurements as improvement on every path. The long-circuit execution remains in the same optimized millisecond regime.

### A–G quality grading

**A. Documented/current symbol use — PASS.** Verified the current QF carrier forms, property signatures, options and contexts in `PauliStabilizer.m`, `Properties.m`, `Conversions.m`, current State/Operator/Circuit/Basis definitions, and live signatures in the independent log. No new external dependency or unsupported signature was introduced. Packed array conversion and both compilation targets were exercised in the real kernel.

**B. Idiomatic functional WL — PASS.** `Compiled.m:207–215` uses one Dispatch/Replace pass, Catenate, packed array conversion and Transpose for column offsets; `Packed.m:175–185` constructs the zero register with SparseArray, Thread and a map over X/Z index sets. Pattern dispatch separates scalar gates from helpers; Fold remains the circuit abstraction. The retained Do at `Compiled.m:124` is inside Compile, where a sequential integer machine loop is appropriate. The unchanged greedy rank-selection Do at `Constructors.m:32` is a genuinely sequential rank-growth algorithm using an Internal Bag, not repeated AppendTo. There is no introduced interpreted loop or per-gate full-tableau rebuild in the bulk path. No concrete idiom rewrite is needed.

**C. Terse physics reasoning — PASS.** The inverse phase comment and composition express the residual-phase cancellation directly (`GateUpdates.m:152–174`). V uses the common primitive decomposition. Packed layout/conversion/validation now lives in Packed.m; bulk encoding/execution lives in Compiled.m; matrix classification/plan construction lives in GateMatrix.m. Trusted constructors are shared rather than duplicated. There is no extra competing model of a QF state and no dead scaffolding introduced in the changed paths.

**D. Physics and representation invariants — FINDING (same single issue below).** Exact unitarity, Pauli conjugation, symplectic tests, normalization, inverse identities and storage equivalence pass. Canonical register-size metadata can nevertheless disagree with the tableau and defeat the native wire bound. A physical register dimension is an input invariant and must not be read from inconsistent metadata.

**E. Nontrivial examples — PASS.** The evidence is not restricted to one qubit or an all-positive stabilizer state: all 24 sign-sensitive one-qubit Clifford actions, entangling 2–4-qubit dense references, all 64 sign assignments for the three-qubit inverse regression in the focused suite, large 1,024-qubit evolution, symbolic measurement signs, and malformed/overflow boundaries were executed.

**F. Symbolic-first derivations — PASS.** Changed execution uses exact binary/62-bit arithmetic and exact sign information, with no decimal physics approximation. Symbolic binary phases are retained and checked under every assignment in the independent two-qubit example. Dense reference matrices use exact rationals, radicals and I; numerics are not a substitute for the derivation. Machine-word packing is a storage representation of exact finite arithmetic.

**G. Required regimes — FINDING (same single issue below).** General arbitrary-register formulas and symbolic phases, all24 exactly solvable actions, entangling dense/scalar references, n→large storage behavior through1,024, boundary word sizes, and ordinary malformed/failure cases are covered. A continuous weak/strong-coupling or classical limit is inapplicable to this finite exact Clifford action; large register size is the relevant scaling regime. One malformed canonical metadata case still crosses the checked execution boundary.

## FINDINGS

### [P2, inherited] Reject inconsistent canonical register metadata before trusting native execution

Source: `QuantumFramework/Kernel/Stabilizer/PauliStabilizer.m:65–67` validates canonical tableau/sign shapes but does not check a stored `"Qubits"` key. The direct stored-key accessor at line132 permits that key to override the actual tableau dimension. `Properties.m:11–12` and `Compiled.m:244–256` then trust the resulting register size.

Safe reproduction (no malformed native execution):

```wolfram
base = PauliStabilizer[2];
p = PauliStabilizer[<|
    "Signs" -> base["Signs"],
    "Tableau" -> base["Tableau"],
    "Qubits" -> 99
|>];
{PauliStabilizerQ[p], psConcreteFastQ[p], p["Qubits"], Dimensions[p["Tableau"]]}
(* {True, True, 99, {2, 2, 4}} *)

calls = 0;
Catch[
    Block[{compiledCliffordFold = Function[Null,
        calls++; Throw["ReachedNativeBoundary", "NativeBoundary"]]},
        p["ApplyCircuit", {"H" -> 50}]
    ],
    "NativeBoundary"
]
(* "ReachedNativeBoundary"; calls == 1 *)
```

The native fold would be indexed by wire50 although there are only two X/Z qubit rows. The new wire-range guard sees99 and passes. I intercepted the call before the native kernel; no claim that an actual process crash was observed is made.

The identical safe reproduction was run against `/private/tmp/qf-stabilizer-before`: it also reports Valid=True, Fast=True, Qubits=99, dimensions={2,2,4}, NativeCalls=1. This is an inherited validation hole, not a regression created by the style work. It remains relevant to the new boundary guarantee.

## RECOMMENDED FIXES

In the canonical branch of `pauliStabilizerQ`, either reject a stored `"Qubits"` value that differs from `Dimensions[a["Tableau"]][[2]]`, or consistently derive canonical qubit count from the tableau regardless of a stored key. The first is a small extension to the representation validator and matches the packed validators' dimension consistency checks. Keep ordinary canonical metadata and symbolic signs supported.

Add a regression that constructs this inconsistent canonical association and mutates an already valid canonical association's stored `"Qubits"`; require invalidity/no fast path, and intercept the native fold with a counter to require zero calls. An equal stored qubit count should retain the chosen documented policy. No broad refactoring is needed.

## Executable evidence and commands

- `/private/tmp/stabilizer-style-r2-independent.wls`, `.log`, `.json`: independent exact/scalar/symbolic/storage/signature probes and all460 focused tests.
- `/private/tmp/stabilizer-style-r2-wvm.wls`, `.log`, `.json`: WVM fallback check.
- `/private/tmp/stabilizer-style-r2-boundaries.wls`, `.log`, `.json`: supplemental QCO/basis and metadata exploration. The exploratory `$Failed` stub intentionally triggered downstream unpack diagnostics; use the clean tagged-Throw reproductions below as the defect evidence.
- `/private/tmp/stabilizer-style-r2-metadata-baseline.wls`: common safe reproduction, run with `QF_REVIEW_REPO` set to original snapshot and current main. Clean logs: `stabilizer-style-r2-metadata-baseline.log`, `stabilizer-style-r2-metadata-current.log`.
- `/private/tmp/stabilizer-style-r2-benchmark.wls`: common five-workload/five-trial benchmark. `stabilizer-style-r2-timing-{before,after}.{log,json}` contain all individual trials, paths, compiler and matching output digests.
- `/private/tmp/stabilizer-style-r2-manifest.json`: reviewed source/test hashes.

Kernel command throughout:

```sh
'/Applications/Wolfram 15.01.app/Contents/MacOS/WolframKernel' -noinit -script /private/tmp/stabilizer-style-r2-independent.wls
```

Timing commands were serialized with `QF_REVIEW_REPO=/private/tmp/qf-stabilizer-pre-style-review QF_REVIEW_LABEL=before` and then `QF_REVIEW_REPO=/Users/mohammadb/Documents/GitHub/QuantumFramework QF_REVIEW_LABEL=after`, using the common benchmark script. No Quiet or message suppression was added by this reviewer. MUnit's expected-message assertions were used as supplied in the target tests.

OPEN ISSUES: 1
