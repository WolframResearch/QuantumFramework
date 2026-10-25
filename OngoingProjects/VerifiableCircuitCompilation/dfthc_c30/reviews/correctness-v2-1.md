# Independent correctness review

ROUND: 1/3; correctness visit 2.

## VERDICT

PASS. The fixed-factor C30 loader, its public reports, its Wolfram Language identities, and the claimed full official verification were independently reproduced. No unresolved correctness or reproducibility defect was found.

The scope is exactly the sixteen public target files listed in `target-sha256.json`. Reviews and the subsequent publishing-stage evidence manifest were excluded. Original targets were never edited. All builds, report regeneration, circuit simulation and additional tests ran in this review's isolated temporary project copy. Target hashes were checked again after testing and remained unchanged.

## REPRODUCED

### Commands and environment

Executed in the isolated `project/dfthc_c30` directory, with the prescribed competition bundle and public factor file supplied read-only:

```sh
PYTHON generator.py --data DATA --public-factors PUBLIC_FACTORS
PYTHON verify.py --data DATA
PYTHON probe_loader.py
KERNEL -noinit -noprompt -script run.wls
```

Also executed the independently authored `independent.py` and `independent.wls` from this review directory. The runtime was Python 3.13.1, NumPy 2.5.4, kickmix 0.1, and Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. Every command completed successfully. Both Wolfram logs were empty: no kernel messages were suppressed or emitted. No `Quiet` or WL `Print` was used.

### Claim inventory and results

| Claim / target | Independent result |
|---|---|
| README / generator: standalone port without downloaded notebook-module imports | Read all target code; imports are standard library, NumPy, kickmix and the official verifier package. The local loader includes its required helper implementations. Fresh build succeeds. |
| NOTICE: released factor and code attribution | Public factor SHA256 is `9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193`. Rechecked source-module hashes against the sibling release audit. Downloaded notebook source blob matches its pinned SHA256, kernel ID 140416287 and version 8. The source files contain the declared Apache-2.0 and upstream construction attribution; version-2 factor metadata declares Apache 2.0. |
| Fixed factor dimensions `(65,60,30)` | Fresh W shape is `(65,60,30)`, U shape `(65,60,76)`, and angle words shape `(3976,75)`. |
| W/H2_I, official shifts and tree fixed | Compared W and H2_I element-for-element with the released factors. Compared H1_I and Hsym_lower against both baseline and packaged fixed shifts. Compared the tree against the supplied baseline. All exact matches. |
| Only a row-sign Givens gauge; branch-zero threshold 1e-12 | Audited the reverse-tree gauge against the hash-pinned released implementation. Independently decoded original and transformed words without calling the generator's gauge or official decoding function. Maximum component drift after row signs was `7.771561172376096e-16`; maximum row-projector drift was `8.881784197001252e-16`. Recomputed the generator's stated norm statistic through the fresh build as `6.661338147750939e-16`. A separate einsum summation gives `7.771561172376096e-16`, consistent with floating-point summation order and the 1e-11 cap. |
| Four-address arbitrary-coefficient GF(2) interpolation | Ran all 17 published exact/edge checks and separately checked additional tables with 1, 3 and 7 columns. Verified symbolic coefficients are retained and integer PolynomialMod reduces coefficients modulo 2. |
| Signed interval endpoints and 9-to-10-bit expansion | Verified all 512 encoded signed words, both endpoints 768/255, compact/expand identity, and all 262144 pairwise XORs remaining in the sign-extended subset. |
| Arbitrary target XOR; phase repair; clean workspace; selection preservation | Audited fanout-sandwich algebra, lookup ANF, AND construction/uncomputation and HMR corrections against live kickmix API documentation. Ran published dirty-target/involution probes plus independent boundary probes. Exact WL truth tables and integer measurement-phase cancellation also passed. |
| Sequential chaining and reverse clearing | Independently rebuilt batch deltas, checked all forward batches and final reverse cancellation, plus nondivisible word counts and padded batch sizes. The official simulation also exercised all forward/reverse circuit stages. |
| Q128 and circuit resources | All 25 rebuilt circuit byte hashes match `build.json`; all widths are 128 and static Toffoli counts are 1129. Fresh profiles, operation counts and reaction depths match completely. The build report differs only in elapsed provenance. |
| Full official simulation at default 64 shots | Executed the unmodified official verifier. Its live signature defaults to 64, and its source checks every address and every circuit in forward/reverse passes, phases, selection, workspace and loaded words. Result: PASS, 3976 addresses, 25 circuits, 254464 shots. This corresponds to 12723200 circuit applications. |
| Exact cost and full block-encoding resources | Reproduced score `2051023203.3721037`, Toffoli count 65244, qubits 661, ANGLE width 128, effective normalization `47.55854233515931`, one-body error `13.911731985012415` mHa and two-body error `19.98067964067221` mHa. Score also independently satisfies lambda_eff × T × Q. |
| Improvement claims and combined estimate | Recomputed from the required baseline/incumbent/DQI/inverse-square-root evidence: DFTHC baseline ratio `1.6360612511835888`, incumbent ratio `1.6154555940004633`, reduction `38.09795801792165%`, and geometric mean `2.0088882102784944`. Previous combined estimate `1.7120801622097384` is also reproduced. Inverse-square-root ratio was recalculated from its recorded baseline cost/current score rather than trusting the hardcoded ratio. |
| Published supplementary coverage | Reproduced exactly 795200 released and 4704 synthetic circuit applications with zero observed value, phase, workspace or selection failures, and 14 documented invalid cases rejected. `probes.json` is byte-for-byte reproducible. |
| 17 symbolic/edge and six complete-data WL checks | Real-kernel rerun passes all 23 checks, including all 298200 released words, forward/reverse batches, Q128 and fixed-factor flag. Report outcomes match; only regenerated build-hash provenance differs. |
| WL FourAddress contract and API claims | Live attributes confirm listable bit operations and zero-indexed BitGet. Live Options and independently exercised BinaryReadList confirm little-endian unsigned-16 import and automatic file closure for filename input. Installed reference pages were read for BinaryReadList and PolynomialMod. FourAddress rejects empty/ragged/wrong-row, float, symbolic, Boolean, unsigned-gap, negative, oversized and wrong-arity cases with Failure. |
| No factor training, Kaggle submission, commit or push | The reviewed scripts do not train factors or submit to Kaggle. This review made no such changes. Existing preceding-audit commit is present and contained in the saved origin/main reference; the new loader is local. Ignored local data/cache paths are correctly ignored by Git. |

### Additional adversarial cases

`independent.wls` passes 15 independently assembled checks, including exact truth tables for dirty fanout and involution, arbitrary-integer phase cancellation, all pairwise signed-word XORs, randomized FourAddress matrices, malformed inputs, file byte ordering/closure and live symbol/options contracts.

`independent.py` passes 31 checks. It adds 46368 dirty-target circuit applications at 16, 24, 36, 256, 1024, 2048, 4092 and 4096 rows, including the largest supported address size and incomplete address-space prefixes. All preserve target XOR/involution, selection, clean workspace and phase. The resulting widths are 120, 121, 122, 124, 126, 127, 128 and 128. It also checks unsigned/bool/complex/oversized malformed tables, nondivisible batch padding and Möbius involution through 1024 entries.

### Evidence integrity

The target hash bindings are in `target-sha256.json`. `final-bindings.json` records successful final comparison of unchanged original targets, all fresh official outcomes, all supplementary outcomes, Wolfram outcomes, source release hashes, notebook source metadata and empty Wolfram message logs. Original reports' source/build hashes bind to the actual reviewed targets. Official dependency hashes match both the feasibility input pins and the fuller baseline dependency pins, including `circuits.py` and `reference.py`.

Fresh generation independently reproduced the exact factorization archive, word cache and all 25 circuit bytes. The full official submission digest is `4ced52799bc913e2cd92ff5724378838d461f268544690ee7c463d7eb098f20c`, matching published evidence. The original `build.json` digest is `3d911ccb95018b64a9389ced0e1400d8f6f3cd3a61caf67c6fb9107456057a54`. Regenerated build digests differ because elapsed build provenance changes; all deterministic contents match.

Detailed evidence files: `generator.log`, `official.log`, `probes.log`, `wolfram.log`, `independent-python.json`, `independent-wl.json`, `kickmix-api.json`, `final-bindings.json`, the two independent test scripts, and regenerated isolated reports.

## FINDINGS

None. No severity-ranked correctness defects or unverifiable required claims remain.

The published limitation about unsampled measurement branches is accurate: the official and supplementary simulations sample measurement outcomes. Their success is not represented here as exhaustive certification of every branch. The exact algebraic checks and source/API audit supplement that sampling. Elapsed fields are provenance only; there is no runtime-performance claim requiring repeated timing trials. The excluded final publishing manifest is not treated as a missing verification dependency.

## RECOMMENDED FIXES

None required for the reviewed target hashes.

OPEN ISSUES: 0
