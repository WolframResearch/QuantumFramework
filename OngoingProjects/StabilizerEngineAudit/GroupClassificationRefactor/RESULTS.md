# Stabilizer group classification refactor

`stabilizerOperatorGroupClassification` owns the wire-set check, exact primitive dispatch and lazy cached general classification. `stabilizerOperatorGroupSpecs` only remaps the local specs (or maps the existing occurrence-specific fallback) and pairs them with circuit positions. Both definitions use functional WL, lexical constants and list operations, with no procedural loops.

Only production change for this request: `QuantumFramework/Kernel/Stabilizer/GateMatrix.m`. Main HEAD remains `89051701`. No staging or commit.

## Correctness

All 1,667 tests in 16 stabilizer test files passed in a fresh Wolfram 15.0.1 kernel, with no outside messages. This includes the 56 bulk-plan tests and 330 inverse-phase tests. See `suite-after.json` for per-file counts and timings.

## Paired timings

Nine alternating trials per mode; medians in milliseconds. Each workload applies a prebuilt flat 1,000-gate Clifford circuit ten times (10,000 gates). QCO construction, native compilation and kernel startup are outside the timed region. Cold clears all six stabilizer caches; warm primes the same workload. FirstPlan additionally resets the immutable primitive dispatch initializer. Runtime and QF property caches are already initialized. `ClearSystemCache[]` runs before every timed sample.

PreBulk restores flat per-operator mapping, equivalent to the old planner for these flat QCO workloads. BeforeHelper installs the exact saved group definition; AfterHelper uses the extracted helper. The new helper adds no memoization. All 270 results have identical tableau/sign signatures within their workload and no messages.

| Workload | PreBulk ms | Before helper ms | After helper ms | PreBulk / After | Paired after/before median |
|---|---:|---:|---:|---:|---:|
| MixedCold/n=128/g=10000 | 192.174 | 75.857 | 75.419 | 2.55x | 0.988 |
| MixedFirstPlan/n=128/g=10000 | 187.265 | 95.190 | 97.973 | 1.91x | 1.000 |
| MixedWarm/n=128/g=10000 | 6.920 | 5.521 | 5.814 | 1.19x | 1.060 |
| RandomCold/n=128/g=10000 | 203.245 | 61.153 | 63.624 | 3.19x | 1.015 |
| RandomWarm/n=128/g=10000 | 7.267 | 6.671 | 6.370 | 1.14x | 0.955 |
| MixedCold/n=1024/g=10000 | 179.872 | 61.905 | 62.718 | 2.87x | 1.035 |
| MixedFirstPlan/n=1024/g=10000 | 174.909 | 85.951 | 96.303 | 1.82x | 1.043 |
| MixedWarm/n=1024/g=10000 | 18.750 | 18.371 | 16.830 | 1.11x | 0.937 |
| RandomCold/n=1024/g=10000 | 206.168 | 72.593 | 73.676 | 2.80x | 1.015 |
| RandomWarm/n=1024/g=10000 | 24.033 | 39.909 | 30.311 | 0.79x | 0.894 |

Background Wolfram jobs were active throughout. Full trial ranges and pairing are retained in `paired.json`; `background-load.txt` records the load. Small differences between helper versions should not be read as a causal change to native evolution, whose source is untouched. Warm plan hits bypass both group implementations. In particular, the randomized 1024-qubit warm row has broad overlapping ranges: its medians do not establish a gain over PreBulk, though the extracted helper version is faster than BeforeHelper in that row. This report does not claim a warm gain for every workload.

## Source identity

- Before extraction GateMatrix SHA256: `e63d64ef3487badac05f03c99838f18d21572ea4e48e318cb8632cd2efc863d2`
- After extraction GateMatrix SHA256: `b5cbb28228a7e7f3e49d205bc20abc5a6368d8ce83617b929536d5340754f2e1`

Independent correctness, functional-style and representative performance reproduction: see `ReviewRound1.md` and the reviewer scripts/logs.

## Independent reproduction

Round 1 converged with zero open issues. The independent reviewer reran all 56 bulk-plan tests, eight exact before/after group comparisons, cache/laziness/remapping probes, and seven paired trials per mode on a different 10,000-gate sequence. All 84 timed outputs matched exact tableau/sign signatures, without messages. Independent cold gains over PreBulk were 3.22x at 128 qubits and 2.39x at 1024 qubits.

Independent warm timings at 1024 qubits were inconclusive: before-helper median 24.869 ms (range 20.600–152.001), after-helper 75.605 ms (18.649–128.547), PreBulk 80.869 ms (21.917–143.952). These measurements do not establish a precise warm regression bound or a universal warm speedup. A separate live check applied warmed plans ten times on each of 128 and 1024 qubits with zero calls to both circuit-spec construction and the extracted helper. The changed code is bypassed on warm plan hits. All seven WL quality axes passed.
