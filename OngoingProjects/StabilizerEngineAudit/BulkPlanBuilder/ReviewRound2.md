# Independent bulk-plan review, round 2 of 3

## VERDICT

Correctness, all seven quality axes, and the scoped plan-construction performance claim pass on the reviewed production bytes. No unresolved defect found.

Reviewed production: `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/Stabilizer/GateMatrix.m`, SHA256 `e63d64ef3487badac05f03c99838f18d21572ea4e48e318cb8632cd2efc863d2`.

Reviewed new suite: `/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer/BulkPlanBuilder.wlt`, SHA256 `2749a6093f54bf6338c1f20b90ece00c0b2468b1f0c3419a5ff7846491210a4d`.

Kernel: Wolfram Language 15.0.1 for Mac OS X ARM (64-bit), July 2, 2026. The loaded main file was the requested main checkout's `QuantumFramework/Kernel/QuantumFrameworkMain.m`. No production or test file was edited by this reviewer. The probes used real C and fresh WVM compilation, without `Quiet` or message suppression.

## REPRODUCED

Probe source: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-probes.wls`. Logs: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-C.log` and `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-WVM.log`.

The C run passed all 18 physical, fallback, representation, and suite assertions. One initial instrumentation assertion incorrectly assumed that reversed relative wire order must bypass the general classifier. That assumption is stronger than the implementation contract: the primitive Dispatch matches canonical stored local data; a different relative order correctly retains the matrix classifier. The corrected canonical primitive assertion passed in the fresh WVM run, which passed 19/19 probe assertions, and subsequently passed independently in C with zero classifier calls (`/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-counts.log`). The new suite was independently rerun in each target: 56 succeeded, 0 failed. Neither run produced unexpected kernel messages.

The independent probes reproduced these claims:

- All 12 canonical exact primitive matrices, including identity, Y, V and its adjoint, controlled gates and SWAP, bypass the general classifier. Approximate matrices and exact global-phase variants enter the general classifier.
- Physical primitive and entangling action equals independently constructed dense `U P U†`, including signed Y images, on a three-qubit register.
- Four combinations of independently ordered input/output wires for a nonsymmetric two-factor operator agree with a dense reference. A changed picture produces distinct physical local data and agrees with the per-operator fallback.
- Four random four-qubit custom Clifford matrices carrying a nontrivial exact global phase, on a nonmonotone wire order, agree with dense conjugation. Three five-qubit random circuits of 35 gates also agree for every X/Z generator.
- Seven-qubit symbolic initial signs stay on the symbolic path and agree with independently folded scalar gates under 40 random binary assignments. The suite additionally exhausts all 64 assignments for its three-qubit sign example.
- A product of exact symbolic phases `Exp[I alpha]` and `Exp[I beta]`, reused on distinct wires with Hadamards and CNOT, yields a state projector agreeing with the dense circuit at nine nontrivial phase assignments.
- Symbolic tokens are reconstructed per occurrence. Mislabelled composite tokens, moved wires, a noncomputational basis, a four-level wire and a nonunitary matrix are refused with the expected message.
- Nested circuits, barriers and empty blocks preserve sequence. Empty packed evolution is identity.
- At 31, 32, 61, 62, 63, 124, 125 and 1024 qubits, 120 random gates agree exactly with scalar evolution from an explicitly canonical tableau. Reusing a 300-gate, 128-qubit circuit 40 times agrees with a single 12,000-gate direct fold.

### Benchmarked-control audit

`/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/PairedBulkPlanBenchmark.wls` has a valid paired control for its stated flat-operator QCO blocks:

- Lines 9–11 reconstruct the baseline flat `Catenate[element /@ operators]` semantics. The baseline GateMatrix snapshot and baseline paclet's GateMatrix file were byte-identical.
- Lines 27–32 switch only builder definitions, clear all six stabilizer engine caches, alternate mode order over seven trials, and call `ClearSystemCache[]` before evaluation.
- Lines 47–48 initialize compiled execution and ordinary primitive dispatch before ordinary rows. Lines 30–31 restore the original delayed OwnValue in the dedicated FirstPlan rows, so those timed After calls include rebuilding the primitive table.
- Lines 35–36 preserve tableau/sign signatures, result heads and message lists. Lines 74–75 reject a changed signature, an unexpected message or a result of the wrong head.
- Cold means stabilizer engine caches cold in an already-loaded runtime with existing QF property caches. FirstPlan measures dispatch-table construction in that same warm runtime. It does not measure fresh process startup. Warm rows mainly check that unchanged plan reuse does not regress. These distinctions should accompany performance claims.

### Independent performance reproduction

Command: `wolframscript -file /Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-performance.wls`; exit code 0. Log and raw trials: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-performance.log` and `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r2-performance.json`. This uses a different deterministic gate recipe with 800 gates per circuit, eight applications, C execution, five alternating paired trials for public cases and 31 paired trials for each direct control. Engine caches are cleared between modes and `ClearSystemCache[]` is called before each timed evaluation. All ten rows have one shared physical signature across modes/trials, and no messages. Startup is outside the timings.

| Qubits | Regime | Before ms, median [min, max] | After ms, median [min, max] | Before/After median |
|---:|---|---:|---:|---:|
| 128 | Cold | 130.692 [130.300, 436.507] | 40.596 [40.364, 40.796] | 3.219 |
| 128 | Warm | 4.948 [4.690, 5.123] | 4.644 [4.507, 4.813] | 1.065 |
| 128 | FirstInitialization | 129.400 [129.071, 130.909] | 61.231 [61.181, 62.751] | 2.113 |
| 128 | DirectCold | 6.491 [6.307, 14.711] | 6.475 [6.247, 9.847] | 1.002 |
| 128 | DirectWarm | 2.145 [2.071, 3.166] | 2.129 [2.084, 2.590] | 1.008 |
| 1024 | Cold | 137.906 [137.530, 138.603] | 48.881 [46.704, 51.361] | 2.821 |
| 1024 | Warm | 13.399 [10.292, 17.869] | 12.169 [11.711, 33.734] | 1.101 |
| 1024 | FirstInitialization | 226.817 [179.410, 300.100] | 152.183 [90.122, 212.413] | 1.490 |
| 1024 | DirectCold | 10.578 [7.825, 35.772] | 10.815 [7.322, 28.835] | 0.978 |
| 1024 | DirectWarm | 10.104 [3.127, 29.144] | 11.734 [3.273, 32.250] | 0.861 |

Cold construction improves by 3.22× at 128 qubits and 2.82× at 1024 qubits in this independent recipe. Including primitive-table construction still improves the paired medians by 2.11× and 1.49×. Every paired cold and first-initialization trial favors After. The parent's 22-row file was independently inspected: all 14 trials per row had matching signatures, correct result heads and empty message lists.

The direct control must not be interpreted from unpaired medians alone. At 1024 qubits, the direct-warm medians are 10.104 and 11.734 ms, but paired wins split 15/16, the median paired After/Before ratio is 1.0233, and the ranges overlap widely. Direct-cold medians are 10.578 and 10.815 ms, yet 21/31 pairs favor After; its median paired After/Before ratio is 0.9137. At 128 qubits both direct paired-median ratios are within 0.2% of 1. An observed process-state check found several unrelated background Wolfram kernels using CPU during the independent run. Absolute times and small warm/control differences therefore need the stated ranges; the evidence does not support a universal warm speedup or a causal direct-code regression.

`Compiled.m` and `PauliStabilizer.m` are byte-identical to the full baseline paclet (SHA256 `c5f4e17678b953788e0aa03f195f9642582ea10c8512d4ae80a8ab098144b3da` and `36e15a279123848412c3e3250785fe6a0a2f9ee2a082079743c97017743dc334`). Recursive comparison of the stabilizer directory finds only GateMatrix.m different. The public warm path can still use fewer encoded gates: for the parent's 1000-gate recipe the baseline emits 2001 encoded gates, while After emits 1290 at 128 qubits and 1286 at 1024 qubits. This was checked in the separate successful C structural-control kernel. Thus an improvement in a public warm plan does not imply a change to the compiled kernel or unchanged direct-list implementation.

## FINDINGS

No production correctness defect found in the checked scope. The narrower primitive-bypass qualification above is a reporting precision point, not a physical bug.

### A. Documented symbols — PASS

The changed builder uses core `With`, `PositionIndex`, `Values`, `SortBy`, `Catenate`, `MapThread`, `Thread`, `Transpose`, `Part`, `Ordering` and `Dispatch` with valid shapes. Package object structures were checked against `QuantumBasis/QuantumBasis.m` and `QuantumOperator/QuantumOperator.m`, and circuit operator extraction against `QuantumCircuitOperator/Properties.m:54`. The atomic QuantumBasis extraction at GateMatrix lines 507–510 works in both kernels. Existing classification and execution contracts at lines 439–493 and `Stabilizer/PauliStabilizer.m:215–273` determine refusal and application. No recalled-only package signature is needed.

### B. Idiomatic functional WL — PASS

GateMatrix lines 527–532 build the immutable reference Dispatch with `KeyValueMap`. Lines 537–539 gather all placement columns in bulk. Lines 552–569 restore occurrences and sequence with `MapThread`, `PositionIndex`, `SortBy`, and `Catenate`. No new `Do`, `For`, `While`, `AppendTo`, repeated growing-array copy or scalar loop appears in the production change. The remaining compiled loop is unchanged and appropriate to the compiled mutation kernel.

### C. Terse physics reasoning — PASS

Lines 505–515 explain what physical data defines a group and how exact matrix content certifies primitives. Lines 535–556 identify local-to-physical wire placement and the occurrence-specific token exception. The implementation follows that argument directly: collect physical data, classify each group, remap each occurrence, restore sequence. The existing physical refusal comments and classifier remain authoritative.

### D. Physics invariants — PASS

Exact signed Pauli conjugation agrees with dense `U P U†` through five-qubit interacting examples; this checks phase/sign conventions, generator ordering and the induced canonical commutation structure. The symbolic phase projector agrees with the dense unitary reference. Approximate/global-phase inputs are kept on the general certificate path, and nonunitarity, basis mismatch, dimension mismatch and moved wires retain refusals. C and WVM agree with scalar canonical execution across packed-word boundaries.

### E. Nontrivial examples — PASS

The reviewed artifact is not limited to a one-qubit or separable base case. The suite contains two-/three-/four-qubit dense checks, large packed registers and 10,000 gates; the independent probes add nonsymmetric input/output wire permutations, four-qubit custom entanglers with exact global phase, five-qubit dense generator references, symbolic phase projectors, 1024-qubit scalar checks and 12,000 reused gates.

### F. Symbolic-first derivations — PASS

The new production path contains no numerical approximation. Primitive matrices use exact rational, algebraic and imaginary entries at lines 516–525. General symbolic phases and signs remain symbolic through the existing classifier and fallback. Numerical phase substitutions and approximate H inputs are verification regimes, not substitutions for the symbolic route.

### G. Regimes and limits — PASS

The source did not provide a separate regime table; the smallest relevant table is reconstructed here from executable evidence before grading coverage.

| Regime | Executable evidence |
|---|---|
| General symbolic/arbitrary size | Seven-qubit symbolic signs, parameterized diagonal phases, group algorithm independent of register size; suite symbolic tokens/signs |
| Exactly solvable special case | Exact reference matrices for all 12 primitives, signed dense Pauli images |
| Large-size/limiting case | Packed boundaries 31–125, 1024-qubit comparison, 12,000 gates and reused plan; empty circuit identity |
| Independent numerical reference | Symbolic phase state projectors at nine exact phase assignments evaluated numerically; dense exact 4-/5-qubit conjugation is an even stronger reference |
| Failure/edge case | Moved wires, mixed input/output permutations, picture/basis distinctions, ququart, nonunitary matrices, dishonest labels, empty and nested blocks |

## RECOMMENDED FIXES

No production fix required. When reporting performance, distinguish engine-cold plan creation from warm plan reuse, include the one-time dispatch-initialization case and measured ranges, and qualify primitive bypass as matching exact canonical stored local data. Avoid a universal warm-performance claim. All reviewer kernels have exited.

OPEN ISSUES: 0
