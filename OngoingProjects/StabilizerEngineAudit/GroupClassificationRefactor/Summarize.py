from pathlib import Path
import json, hashlib
p=Path(__file__).parent
bench=json.loads((p/'paired.json').read_text())
suite=json.loads((p/'suite-after.json').read_text())
assert len(bench['Rows']) == 10
all_trials=[t for row in bench['Rows'] for t in row['Trials']]
assert len(all_trials)==270
assert all(t['Head']=='PauliStabilizer' and t['Messages']=='{}' for t in all_trials)
assert all(len({t['Signature'] for t in row['Trials']})==1 for row in bench['Rows'])
assert bench['OutsideMessages']=='{}'
assert sum(row['Count'] for row in suite['Rows'])==1667
assert all(not row['Failures'] for row in suite['Rows']) and suite['Messages']=='{}'
text='''# Stabilizer group classification refactor

`stabilizerOperatorGroupClassification` owns the wire-set check, exact primitive dispatch and lazy cached general classification. `stabilizerOperatorGroupSpecs` only remaps the local specs (or maps the existing occurrence-specific fallback) and pairs them with circuit positions. Both definitions use functional WL, lexical constants and list operations, with no procedural loops.

Only production change for this request: `QuantumFramework/Kernel/Stabilizer/GateMatrix.m`. Main HEAD remains `89051701`. No staging or commit.

## Correctness

All 1,667 tests in 16 stabilizer test files passed in a fresh Wolfram 15.0.1 kernel, with no outside messages. This includes the 56 bulk-plan tests and 330 inverse-phase tests. See `suite-after.json` for per-file counts and timings.

## Paired timings

Nine alternating trials per mode; medians in milliseconds. Each workload applies a prebuilt flat 1,000-gate Clifford circuit ten times (10,000 gates). QCO construction, native compilation and kernel startup are outside the timed region. Cold clears all six stabilizer caches; warm primes the same workload. FirstPlan additionally resets the immutable primitive dispatch initializer. Runtime and QF property caches are already initialized. `ClearSystemCache[]` runs before every timed sample.

PreBulk restores flat per-operator mapping, equivalent to the old planner for these flat QCO workloads. BeforeHelper installs the exact saved group definition; AfterHelper uses the extracted helper. The new helper adds no memoization. All 270 results have identical tableau/sign signatures within their workload and no messages.

| Workload | PreBulk ms | Before helper ms | After helper ms | PreBulk / After | Paired after/before median |
|---|---:|---:|---:|---:|---:|
'''
for row in bench['Rows']:
    m=row['MediansSeconds']
    text+=f"| {row['Label']} | {1000*m['PreBulk']:.3f} | {1000*m['BeforeHelper']:.3f} | {1000*m['AfterHelper']:.3f} | {m['PreBulk']/m['AfterHelper']:.2f}x | {row['MedianPairedAfterBefore']:.3f} |\n"
text+='''
Background Wolfram jobs were active throughout. Full trial ranges and pairing are retained in `paired.json`; `background-load.txt` records the load. Small differences between helper versions should not be read as a causal change to native evolution, whose source is untouched. Warm plan hits bypass both group implementations. In particular, the randomized 1024-qubit warm row has broad overlapping ranges: its medians do not establish a gain over PreBulk, though the extracted helper version is faster than BeforeHelper in that row. This report does not claim a warm gain for every workload.

## Source identity

'''
start=json.loads((p/'starting-manifest.json').read_text())
text+=f"- Before extraction GateMatrix SHA256: `{start['BeforeSourceSHA256']}`\n- After extraction GateMatrix SHA256: `{start['AfterSourceSHA256']}`\n"
text+='\nIndependent correctness, functional-style and representative performance reproduction: see `ReviewRound1.md` and the reviewer scripts/logs.\n'
(p/'RESULTS.md').write_text(text)
print(text)
