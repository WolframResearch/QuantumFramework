# ROUND 1/3 — independent verifier and WL quality critic

## VERDICT

PASS for the named-helper extraction in
`QuantumFramework/Kernel/Stabilizer/GateMatrix.m:541–562`: correctness, all seven WL quality
axes, and retained cold-planning gains. No consistent extraction regression appeared in the
cold representative checks. Warm plan hits provably bypass the helper. Warm timing variance,
especially at 1024 qubits, prevents a precise measured bound or a universal warm-speedup claim.

The reviewed source is SHA256
`b5cbb28228a7e7f3e49d205bc20abc5a6368d8ce83617b929536d5340754f2e1`, on main HEAD
`890517010618b3a2bb018d90758c3d96df3ed4c2`. The diff against the supplied before snapshot
contains only the extraction and the equivalent output-order Part selection. No production
file was edited by this reviewer.

## REPRODUCED

Kernel command:

```sh
'/Applications/Wolfram 15.01.app/Contents/MacOS/WolframKernel' -noinit -script /private/tmp/qf-group-helper-review-r1/Correctness.wls > /private/tmp/qf-group-helper-review-r1/correctness.log 2>&1
```

The command exited 0. The loaded paclet came from the requested main checkout, not the
installed paclet. Kernel: Wolfram 15.0.1, ARM64 macOS; QF 2.1.3. No messages were suppressed;
the final message list was `{}`.

- `Tests/Stabilizer/BulkPlanBuilder.wlt`: 56/56 successes, zero failures.
- Eight independently assembled before/after group calls were `SameQ`: identities;
  misleadingly relabelled H; reverse-wire CNOT; three-wire H decomposing into three local
  specs; exact global-phase H; symbolic P with one valid and one dishonest label; composite
  T.H with one valid and one dishonest label; and a moved-wire later occurrence.
- The primitive helper returned `{"H" -> {1}}` while making zero gate-key and zero general
  classifier calls. The fallback RHS remains lazy.
- A valid first occurrence plus a later mismatched input/output set returned `$Failed`
  without evaluating primitive dispatch or general classification.
- A nonprimitive globally phased H on two different absolute wires invoked the general
  classifier exactly once. The second helper call used the gate cache.
- Direct multi-spec remapping across three arbitrary sorted wire sets exactly matched an
  independent occurrence-by-occurrence `Replace`; empty specs returned three empty blocks,
  and the empty-spec/empty-occurrence case returned `{}`.
- The suite additionally exercised exact primitive dense U.P.U-dagger conjugation, random
  2–4-qubit circuits, all 64 binary assignments to six symbolic signs, machine and
  40-digit custom gates, packed word boundaries through 1024 qubits, a reused 10,000-gate
  plan, misleading labels, basis and qudit refusal, nested circuits and barriers,
  nonunitary refusal, the non-Clifford frame boundary, and all storage layouts.

Logs and machine-readable results: `correctness.log`, `correctness.json`. The exact
before definition was loaded from the supplied snapshot into the review kernel only,
compared, and restored in that kernel. Source hashes match the starting manifest.

### Independent performance reproduction

After the parent's timed kernel exited, the following separate kernel ran a different
deterministic eleven-gate mixed sequence, prebuilt as 1000 operators and applied ten times.
Construction and startup were outside the clock. Each row used seven trials per mode,
randomized mode order, `stabilizerClearCaches[]` and `ClearSystemCache[]`. Warm rows first
seeded the plan once. The old flat per-operator path and the exact before-helper definition
were substituted only inside this review kernel.

```sh
'/Applications/Wolfram 15.01.app/Contents/MacOS/WolframKernel' -noinit -script /private/tmp/qf-group-helper-review-r1/IndependentTiming.wls > /private/tmp/qf-group-helper-review-r1/independent-timing.log 2>&1
```

The command exited 0, using the C compilation target. All 84 outputs had the same exact
tableau/sign hash within each workload, every output was a `PauliStabilizer`, and all trial
and outside message lists were `{}`.

| Qubits | Plan | Before helper median (ms) | After helper median (ms) | Pre-bulk median (ms) | Median paired After/Before |
|---|---|---:|---:|---:|---:|
| 128 | Cold | 60.513 | 61.011 | 196.755 | 1.0072 |
| 128 | Warm | 6.357 | 6.597 | 6.989 | 1.0019 |
| 1024 | Cold | 118.472 | 87.212 | 208.110 | 0.8428 |
| 1024 | Warm | 24.869 | 75.605 | 80.869 | 2.0882 |

The cold after-helper path retained 3.22x and 2.39x median gains over the old flat planner.
The extraction did not show a consistent cold slowdown. These timings do not establish a
precise overhead bound: cold 1024 ranges were 74.469–160.242 ms before and 79.500–223.495 ms
after. Warm 1024 was particularly inconclusive, spanning 20.600–152.001 ms before,
18.649–128.547 ms after, and 21.917–143.952 ms pre-bulk. Its raw medians must not be presented
as a measured helper regression or as proof of universal warm gains.

A separate `WarmCacheProof.wls` kernel then warmed a circuit and replaced both the circuit-spec
builder and the extracted classification helper with counting poison functions. Ten subsequent
applications on each of 128 and 1024 qubits returned valid tableaux with **zero calls to either
function** and no messages. The warm plan-hit execution therefore does not execute the changed
helper. Command shape was identical to the two commands above; its output is saved in
`warm-cache-proof.log`, and it exited 0.

Supplemental parent evidence, distinguished from this independent run: the complete 1667-test
stabilizer suite passed; the parent's ten workload rows retained exact signatures across 270
outputs. Its mixed cold medians before/after were 75.857/75.419 ms at 128 and 61.905/62.718 ms
at 1024; after-helper gains over pre-bulk were 2.55x/2.87x. Its random warm 1024 row was slower
than pre-bulk with broad overlapping ranges, so there is no universal warm-speedup claim.

## FINDINGS

No correctness, source-quality, or reproducible extraction performance defect found in the
reviewed delta. The explicit warm measurement limitation above is retained rather than treated
as a precise no-regression result. No required universal warm-speedup claim is made or approved.

### A. Documented symbols — PASS

No external API is introduced. The new private helper's signature is explicitly
`operator_QuantumOperator, data_, orders_List` (`GateMatrix.m:541`), matching its sole
caller at lines 551–552. The live working-tree implementation and kernel probes confirm
its returned local-spec list or `$Failed`/Missing refusal behavior. The `QuantumOperator`
input/output property definitions confirm order storage is `{output,input}`; selecting
`operators[[All,2,2]]` therefore selects the same input orders as the prior intermediate
`orders[[All,2]]`. No shape or signature is inferred solely from audit prose.

### B. Idiomatic functional WL — PASS

The delta adds no `For`, `Do`, `While`, `AppendTo`, indexed mutation, or unnecessary
`Module`. `AllTrue` validates occurrence sets (`:542`); delayed `Replace` and `Lookup`
retain the existing conditional and cache semantics (`:543–546`); `With` binds the
constant local classification; and `MapThread` pairs results with original positions
(`:551–561`). The `If` expresses a required lazy boundary rather than a finite enum
dispatch. Existing compiled loops outside the target are unchanged.

### C. Terse physics reasoning — PASS

The new helper's name identifies classification; `stabilizerOperatorGroupSpecs` now
reads as classify, remap local wires or recover each occurrence, pair with positions.
`operator`, `orders`, and `local` name their actual domain objects. The occurrence-local
label warning is retained immediately above the fallback (`:557–559`). No dead state,
duplicated classification body, or numeric scaffolding was added.

### D. Physics invariants — PASS

The independently rerun 56-test suite certifies primitive and random exact Pauli
conjugations against dense matrices, including signs, reversed control/target placement,
and misleading display labels. The exact global-phase cases exercise the invariant
that conjugation ignores a unit-modulus scalar. Invalid unitary/basis/dimension/wire
cases retain their refusal boundary. Before/after `SameQ` on fallback and multi-spec
results guards against reordering a noncommuting gate sequence at the extraction.

### E. Nontrivial examples — PASS

Coverage extends beyond a single real one-qubit gate: CNOT/CZ/SWAP with reversed wires,
random entangling 2–4-qubit circuits, three local factor specs, symbolic signs, and
128/1024-qubit repeated circuits. Test IDs include `BPB-RandomDense-*`,
`BPB-PackedBoundary-*`, and `BPB-TenThousandGates-ReusedPlan`.

### F. Symbolic-first derivations — PASS

The helper only reorganizes exact classification and wire mapping. It adds no
`N`, tolerance, decimal constant, or numerical solver. Symbolic `P[ArcCos[theta]]`
recovery/refusal and all binary assignments to symbolic tableau signs passed.
Machine/40-digit input is confined to numerical cross-check regimes.

### G. Regimes and limits — PASS

| Regime | Executable evidence |
|---|---|
| General symbolic | `BPB-SymbolicSigns-AllBinaryAssignments`, `BPB-SymbolicTokensRemappedPerOccurrence`; before/after symbolic-token group equivalence |
| Exactly solvable | All primitive exact dense conjugation tests; random 2–4-qubit exact dense reference |
| Large/limiting | Packed boundaries 31,32,61,62,63,124,125,128,1024; identity and empty circuits; direct empty specs; reused 10,000 gates |
| Numeric reference | Machine and 40-digit custom H compared with an independently assembled exact dense action |
| Failure/edge | Later moved occurrence, dishonest symbolic/composite labels, noncomputational basis, qudit, nonunitary, nested/barrier path, non-Clifford frame boundary |

No regime is omitted as inapplicable.

## RECOMMENDED FIXES

None. Keep the performance statement scoped to retained cold-planning gains and the helper's
unchanged warm bypass; retain the warm-1024 timing limitation when reporting results.

OPEN ISSUES: 0
