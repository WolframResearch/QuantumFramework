# ROUND: 1/3 — independent correctness and quantum-quality review

## VERDICT

One optimization-contract defect remains. Physical execution checks pass, but the new physical grouping helper is not defined on the operators it receives. It leaves `stabilizerOperatorData[op]` unevaluated, so the group key still contains the complete operator, including its label and absolute wire numbers. Primitive dispatch therefore misses ordinary named/relabelled operators except an operator structurally identical to its default reference object.

Production and repository tests were not edited by this reviewer. All correctness kernels launched by this reviewer have exited. No timing was collected.

Reviewed source SHA1: `GateMatrix.m` = `3cb134873e48a03e12049fa3af636b85ac7c6686`; `BulkPlanBuilder.wlt` = `411bdcacee228b47235f467af896d880f9161ea6`. Baseline SHA1 = `ef2db2e29b493db532cdbac90a19e2ba5ff07724`.

## REPRODUCED

Environment: Wolfram 15.0.1 for Mac OS X ARM (64-bit), July 2, 2026; repository paclet loaded from `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/Kernel/QuantumFrameworkMain.m`; installed Arrays 1.4.1 and TensorNetworks 1.2.0.

Established test command, with output redirected to the named evidence logs:

```text
'/Applications/Wolfram 15.01.app/Contents/MacOS/wolframscript' -file Tests/RunTests.wls BulkPlanBuilder
'/Applications/Wolfram 15.01.app/Contents/MacOS/wolframscript' -file Tests/RunTests.wls CompiledEvolution
'/Applications/Wolfram 15.01.app/Contents/MacOS/wolframscript' -file Tests/RunTests.wls ExecutionContracts
```

Results: BulkPlanBuilder **53/53 GREEN**, CompiledEvolution **65/65 GREEN**, ExecutionContracts **84/84 GREEN**. Separately, a fresh kernel set `$stabilizerCompileTarget = "WVM"` before any native evolution and reproduced all **53/53** new tests. No Quiet or message suppression was used.

Claims independently tested:

- Matrix-honest physical action, order, signs, relabelled primitives, reverse wires, and general exact/custom/approximate/global-phase gates: the established 53 tests pass. My independent reference builds a full matrix by enumerating basis-index pairs and embeds each local operator on its sorted wires. It compares `U P U†` directly with signed Pauli matrices assembled from returned stabilizer strings, without circuit compilation, `stabilizerPauliFromMatrix`, or bulk helper internals. Four custom cases pass: asymmetric input/output permutations; repeated exact custom entanglers at different wire placements; algebraic global phases; mixed dense/sparse/machine/50-digit storage. Eight additional random 32-gate four-qubit exact circuits pass.
- Packed boundaries and nontrivial register sizes: independent scalar application agrees with public circuit application at n = 31, 32, 61, 62, 63, 124, 125, 257, 1024. A 1,000-gate block applied five times agrees with scalar application of all 5,000 gates. The established test separately checks a reused plan for 10,000 applied gates.
- Symbolic signs: ten symbolic sign bits agree with scalar application under every one of the 1,024 valid binary assignments. The emitted BitXor forms need not be syntactically identical before substitution.
- Symbolic phase: the circuit H1, P[a]1, CNOT1→2, H2, P[-a]2 agrees with the legacy builder. Its density matrix equals that of `{1, Exp[-I a], Exp[I a], -1}/2` for real a. This checks normalization and relative phase without mistaking an irrelevant global phase for a defect.
- Twelve edge tests pass: empty zero-qubit register, nested empty/mixed circuits, measurement refusal, two moved-wire cases, wire zero refusal, existing public register extension, noncomputational basis refusal, four-level-wire refusal, nonunitary diagonal refusal, symbolic missing-label refusal, and dishonest composite-label refusal. Expected refusal messages are recorded in the log.
- General custom classification happens once for sixty occurrences on distinct wires in the interception probe because the existing gate memo still reuses classifications. This does not establish bulk grouping: the direct grouping defect remains independently proven by the unevaluated key and named-H key inequality.
- Existing compiled target code is unchanged by this patch. C-mode regression suites and fresh WVM new-suite execution pass.

Evidence: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit/BulkPlanBuilder/qf-bulk-plan-review-r1-suite.log`, `-compiled.log`, `-contracts.log`, `-wvm.log`, `-probes.wls`, `-probes.log`, `-debug.wls`, `-debug.log`, `-structure.wls`, `-structure.log`, `-dispatch-debug.wls`, `-dispatch-debug.log`, `-dispatch-cause.wls`, `-dispatch-cause.log` (all share that prefix).

Probe corrections are retained transparently: an initial temporary script had a bracket syntax error; the corrected script runs to completion. Initial symbolic comparisons incorrectly used SameQ on sparse zero arrays and unresolved symbolic BitXor forms; Normal and exhaustive valid-bit assignments resolve those. An initial expected refusal on an undersized supplied register was incorrect: direct legacy comparison confirms both public circuit paths extend that register. The malformed-constructor checks in `-structure.log` assumed constructor contracts that QF does not promise; they are not production findings. `-dispatch-cause.log` also contains shadow warnings induced only by its temporary `$ContextPath = {}` formatting diagnostic; this is not a clean-load warning and was not suppressed.

## FINDINGS

### P1 — The physical grouping key helper has an evaluated pattern on its left-hand side

`QuantumFramework/Kernel/Stabilizer/GateMatrix.m:509–510` defines a deeply nested pattern using active QF constructors:

```wl
stabilizerOperatorData[QuantumOperator[QuantumState[amplitudes_, basis_QuantumBasis], order_]] := ...
```

The live DownValues diagnostic shows its left-hand-side constructor evaluation has become `stabilizerOperatorData[Failure["ConfirmationFailed", ... QuantumBasisQ ...]]`. Normal operators therefore do not match it. This is an evaluation trap, not a numerical classifier issue.

Concrete untimed commands in `-dispatch-debug.wls` assign `a = QuantumOperator["H", {2}, "Label" -> "T"]` and `b = QuantumOperator[Hmatrix]`. Calling the helper on a and b yields two **unevaluated helper calls**, rather than data triples. The two keys are unequal, despite identical stored amplitudes, physical basis, and relative order. Replacing those keys by the primitive dispatch yields `$Failed` for a and `{ "H" -> {1} }` for b, because b happens to equal the exact reference expression. Disabling the general classifier in `-structure.wls` records **12 classifier calls and refusal placeholders for twelve exact relabelled primitives**; the primitive table should certify them without that classifier.

Impact: lines 564–567 group complete operator expressions instead of the intended amplitude/basis/order fingerprint. This violates the main requested bulk contract and causes label/wire variation to split physical groups. Fallback preserves physical results, so the 53 behavioral tests do not expose it.

Smallest fix: make the helper boundary `op_QuantumOperator` and obtain its state/amplitudes, basis, and orders through Part. Keep the separate `stabilizerBasisData` dispatch boundary for the atomic basis. Add direct fingerprint equality and shape assertions for operators differing only by label/absolute wires, plus primitive-recognition coverage that intercepts the general classifier. Retain inequalities for different stored bases, dimensions, and relative order patterns.

## QUALITY AXES A–G

| Axis | Result | Evidence and smallest fix |
|---|---|---|
| A. Documented symbols/signatures | FINDING | Live implementation and DownValues show the active constructor pattern at line509 evaluates during definition, leaving no usable physical-key definition. Use an outer QuantumOperator pattern and Part extraction. The actual getters, order shape, matrix semantics, atomic QuantumBasis boundary, PauliStabilizer APIs and compiled target are supported by local source and real-kernel executions. |
| B. Idiomatic functional WL | PASS | Lines498–569 use With, Map, PositionIndex, KeyValueMap, MapThread, SortBy and Catenate. There is no imperative loop or repeated list append in the added builder. The single `Append` constructs a bounded twelve-entry immutable dispatch table. |
| C. Terse physics reasoning | PASS | Physical data, relative order, local specs, absolute wire remapping, and restoration of occurrence order are explicit. Labels are intentionally excluded from classification while symbolic/composite token fallback remains per occurrence. The one active evaluation defect is counted once under A. |
| D. Physics invariants | PASS | Independent exact/numeric Pauli conjugation checks all 2n generators, including their signs, for entangling four-qubit circuits. Symbolic density matrix matches an exact normalized reference. Dimension/basis/moved-wire refusals and compiled symplectic tests pass. |
| E. Nontrivial examples | PASS | Entangling custom matrices, asymmetric input/output permutations, random noncommuting gates, symbolic frame phases, up to1,024 qubits, packed boundaries and5,000/10,000 applied gates exceed the single-qubit/separable base case. |
| F. Symbolic-first derivations | PASS | Primitive tables use exact integers, I and radicals. Symbolic signs and a symbolic real phase survive; approximate matrices are checked separately. No premature numerical replacement was introduced. |
| G. Regimes/limits | PASS | Reconstructed regime table below gives executable evidence for every applicable row. |

| Regime | Evidence |
|---|---|
| General symbolic / arbitrary size | Ten symbolic sign bits under all binary assignments; symbolic real phase; arbitrary runtime n exercised through1,024. |
| Exactly solvable special case | Twelve primitives and independent four-qubit exact conjugations; closed-form symbolic two-qubit density matrix. |
| Asymptotic / limiting | Large-register packed-word boundaries and repeated5,000/10,000-gate application. Empty/identity is the zero-action limit. No continuous physical approximation is introduced by a discrete plan builder. |
| Numeric reference | Independent machine and50-digit matrix conjugation in the mixed-representation case. |
| Failure / edge | Empty/mixed/nested circuits, moved wires, nonunitarity, noncomputational basis, qudit dimension, dishonest labels, wire zero. |

## RECOMMENDED FIXES

1. Correct the physical fingerprint helper and add direct structural contract tests before timing or claiming the intended grouping improvement.
2. Have a fresh round2 reviewer rerun the revised bytes, especially primitive classifier interception, relabelled/wire-shifted key equality, both compile targets and physical references.
3. The paired benchmark control correctly restores the original flat-operator mapper and alternates modes while clearing plan caches. Assert all per-row output signatures agree before treating a row as valid; the parent reports adding that assertion. First-plan primitive-table initialization should be reported separately from the already-built immutable-table case. No performance claim has been evaluated in this round.

OPEN ISSUES: 1
