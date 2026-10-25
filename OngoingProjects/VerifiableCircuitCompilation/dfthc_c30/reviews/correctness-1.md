# C30 independent correctness verification — round 1/3

## VERDICT

PASS for the loader, its WL contracts, the released-data evidence, resource counts and combined local estimate. The complete standalone family was independently rebuilt and simulated in an isolated temporary tree. No target artifact or official competition code was changed. No prior review was read. No approval barrier, missing dependency, timeout or hidden WL message reduced coverage.

Python 3.13.1, NumPy 2.5.4, kickmix 0.1; real Wolfram kernel 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. No `Quiet`, WL `Print` or message suppression was used. Timings are run provenance, not a speed claim; no performance claim is made by the artifacts.

## REPRODUCED

### Target and dependency coverage

Read all 16 requested public target files: README.md, DOCUMENTATION.md, RESULTS.md, PLAN.md, NOTICE.md, LICENSE, loader.py, generator.py, verify.py, probe_loader.py, LoaderIdentities.wl, run.wls, build.json, verified.json, probes.json, wl-checks.json. No applicable AGENTS.md was present. Inspected the directly required sibling input pins, release/provenance records, DFTHC and DQI baseline/mixing/current cost records, and inverse-square-root verified comparison. Inspected the actual official verification, circuit simulation, metric, cost and fixed-shift implementations and the released gauge/ANF source sections used by this port. Installed BinaryReadList and PolynomialMod documentation plus live usage/options/attributes were checked.

### Explicit claims, examples and expected outputs

1. **Fixed inputs and provenance:** the version-2 public factor SHA256 is `9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193`; W has shape (65,60,30), H2_I has shape (65,30), the word table has shape (3976,75), U has shape (65,60,76), and the tree has shape (75,2). Rebuilt W/H2_I equal the released arrays exactly; H1_I/Hsym_lower/tree equal the baseline exactly; H1_I/Hsym_lower equal the official fixed shifts. All relevant official code hashes match the earlier pinned baseline. The source port contains no downloaded notebook-module import or factor fitting path. Its attribution agrees with the retained Apache-2.0 headers and source links.

2. **Only row signs change physically:** independent decoding gave 2011 negative row signs and maximum sign-adjusted vector drift `7.771561172376096e-16`. The stated norm calculation reproduces `6.661338147750939e-16` (an independent einsum summation gives `7.771561172376096e-16`, the expected rounding-order variation). The generated U is exactly the decoded two-body table. The `1e-12` zero-branch threshold and `1e-11` rejection cap agree with the released implementation. WL independently confirms symbolic sign-projector equality and Givens norm.

3. **Exact four-address algebra:** arbitrary GF(2) coefficients reproduce f0,f1,f2,f3 at (0,0),(1,0),(0,1),(1,1). Every one of the 512 encoded values in [-256,255] modulo 1024 has bits 9 and 8 equal, and compact/expand is exact. The interval endpoints, all 512 constant tables, 200 additional random four-by-seven tables and every released four-row block reconstruct exactly. Three compact differences use 81 bank wires; actual released circuits allocate 12 selection + 30 public + 81 bank + 5 scratch = 128 qubits.

4. **Dirty targets and phases:** symbolic integer measurement-phase cancellation, ANF distribution, dirty-target fanout and XOR involution all pass. Additional exhaustive Boolean fanout, 1024-by-512 word XOR-involution cases, and integer phase cases pass. The emitted circuits pass arbitrary dirty-target XOR, two applications restoring the input, clean workspace, unchanged selection and unflipped phase.

5. **Sequential deltas:** independent NumPy accumulation recovers all 25 consecutive three-word batches and reverse traversal clears the final batch. WL does the same over all 298200 words, and confirms signed difference bits.

6. **Deterministic standalone build:** an independent generator run recreated every fixed-array hash, all 25 profile entries, all 25 circuit bytes/hashes, the complete factorization NPZ hash `340efa960c4128a67f302cf38a3d24405c7d38c258ce570c263a6899a13c765c`, and word hash `8903db9b02b1c314cb3b6c8b760738791e42d4bbce7a7738170217b95ac84761`. Every actual circuit has 128 qubits and 1129 static Toffolis. Actual operation counts and reaction depths match every profile; total operations are 3241737, total reaction depth 250. Rebuilt build.json differs only in elapsed_seconds.

7. **Full official verification:** reran the unmodified verifier, default 64 shots on each of 3976 addresses, all 25 circuits forward and reverse. PASS; verified_shots 254464 (shots, not circuit applications), submission hash `4ced52799bc913e2cd92ff5724378838d461f268544690ee7c463d7eb098f20c`. Inspected the official loops: every circuit checks phase, selection, workspace and words; reverse restores the preceding batch and finally zero. This remains sampled measurement-branch evidence, as the documentation correctly states.

8. **Actual costs and metric gates:** all reported metrics and resource fields reproduced exactly: score `2051023203.3721037`, full Toffolis 65244, full qubits 661, ANGLE maximum 128, one-body error 13.911731985012415 mHa, two-body error 19.98067964067221 mHa, lambda_eff 47.55854233515931. The official score formula lambda_eff × Toffolis × qubits gives the same number. State preparation and shared temporary registers were checked against the live cost implementation. Remaining two-body headroom is 0.019320359327789305 mHa.

9. **Comparisons and combined estimate:** recomputed from sibling verified costs: inverse-square-root ratio 3.9385957568764147, DQI 1.2581323878369532, DFTHC baseline 1.6360612511835888; prior DFTHC improvement 1.6154555940004633 and reduction 38.09795801792165%. Their geometric mean is 2.0088882102784944; replacing the current DFTHC cost by the prior incumbent gives 1.7120801622097384. The estimate is below 8.43991. No leaderboard result is asserted or obtained by these commands.

10. **Supplementary published probes:** independently reran probe_loader.py: exactly 795200 released circuit applications, 4704 synthetic applications, all 14 malformed cases rejected; zero phase, workspace, selection and value errors. Synthetic zero/constant/nonconstant tables at 16,20,32,128 rows all pass. Added randomized dirty-target simulations at 24,28,36,60,68,124,132,252,260,508,516,1020,1028,2044,2052,4092,4096 rows with distinct seeds, all passing and at most Q128; 14 additional malformed build inputs also reject. Ten additional Möbius involution tests (truth tables of sizes 2 through 1024) pass. Batch padding shapes were checked at several widths and batch sizes.

11. **WL exact and full-data checks:** fresh `run.wls` returned success and reproduced all 17 named symbolic/edge checks, all six full-data checks, and word_count 298200. Independent tests checked no-argument and excess-argument calls, ragged/empty/wrong-row-count matrices, floats, unsigned-gap words, negatives, overflow, symbolic entries, booleans, complex values, infinity and indeterminate; FourAddress returns Failure for all. Nonempty valid four-row matrices return the original matrix exactly.

12. **Live WL/API contracts:** PolynomialMod integer reduction, Boolean polynomial moduli, BitGet zero-based indexing, Listable attributes of BitGet/BitXor/BitAnd/BitOr, and BinaryReadList's ByteOrdering option were checked in the real kernel. A constructed byte file reads `{256,255,768,1023}` as little-endian UnsignedInteger16, and the filename stream is closed afterward. Python inspect.signature confirms the live official verify/evaluate_words/load_inputs/resource_costs/passes_accuracy call forms used here. No QuantumFramework API is used. The recorded preceding commit exists and is contained in local origin/main; the new loader remains untracked/local.

### Commands and evidence

Commands were run in an isolated copy with environment-specific locations supplied as variables:

```sh
"$PYTHON" generator.py --data "$DATA" --public-factors "$PUBLIC_FACTORS"
"$PYTHON" verify.py --data "$DATA"
"$PYTHON" probe_loader.py
"$KERNEL" -noinit -noprompt -script run.wls
"$PYTHON" "$AUDIT/independent.py"
"$KERNEL" -noinit -noprompt -script "$AUDIT/independent.wls"
```

The audit directory contains generated reports/circuits, original-evidence copies, official-run.log, probe-run.log, wl-run.log, independent-python.json/log, independent-wl.json/log and comparisons.json. All rerun report fields match the originals after excluding elapsed_seconds and hashes that necessarily change when the fresh build's elapsed_seconds changes. Original executable target hashes remained unchanged throughout this verification.

## FINDINGS

No unresolved correctness, API-contract, example, numerical-metric or resource finding. No unsupported speed claim requires a benchmark.

The pipeline's final review manifest is written after the independent reviews; its completion and final source-hash binding should be checked by the parent before delivering the final artifact set. This is a remaining publishing-stage action rather than a loader/WL correctness defect. No prior review was read to assess its contents.

## RECOMMENDED FIXES

None required for the verified code or technical claims. Finish the planned review manifest before presenting review evidence as complete.

OPEN ISSUES: 0
