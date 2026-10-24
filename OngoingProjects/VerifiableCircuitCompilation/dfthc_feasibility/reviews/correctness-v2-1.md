# VERDICT

FAIL: the main feasibility and released-source audit reproduce, but one public API boundary case silently loses nonzero spectral modes. This is correctness visit 2, round 1/3. No target artifact was edited; all reproduction writes used a temporary copy, including its ignored caches and sibling incumbent summary. No prior review was read.

# REPRODUCED

Commands completed with exit status 0 using the supplied Python environment and Wolfram kernel:

```text
python <temporary-copy>/prepare.py --data <official-data> --released <public-release-audit>
WolframKernel -noinit -noprompt -script <temporary-copy>/run.wls
python <temporary-copy>/audit_loader.py --data <official-data> --source <pinned-module-source> --words <public-word-table>
WolframKernel -noinit -noprompt -script <independent-adversarial-tests>.wls
WolframKernel -noinit -noprompt -script <independent-real-numeric-edge-tests>.wls
```

Versions: Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026; NumPy 2.5.4; SciPy 1.18.1. Diagnostics used Echo; no Quiet, Print, message suppression, or performance benchmark was used. No runtime warning or kernel message occurred. The test scripts and complete command stdout are retained beside this report in temporary storage.

Claims and independently reproduced evidence:

- Full pair packing preserves the official N^4 Frobenius residual: full 0.018639570630973906 Ha, packed 0.0186395706309741 Ha. Exact pair and exchange symmetry errors are 0. The regenerated 2926-by-2926 little-endian matrix has SHA256 4c7e70360e876c9a13a7ffc3482e5cb7979b18c9cef5fa597923465563ae1efd. An additional exact, generic 3-orbital pair embedding satisfies E^T E=I and preserves every entry of an arbitrary six-by-six coefficient matrix in its full tensor Frobenius norm.
- The PSD/projector reasoning and cap Min[d,R Min[B,C]] follow from the official reconstruct_x/reconstruct_h2 definitions. Additional exact rank realization and cap saturation/asymmetry cases pass.
- All 23 stored symbolic/edge checks pass, including exact radical/constant ordering, four-address Boolean interpolation, all 512 signed words, invalid input, scratch failure and propagation.
- The full fresh WL eigensystem reproduces rank 290 at 20 mHa. Rank 289 lower bound is 0.020187308109326476 Ha; rank 290 is 0.019858425813203292 Ha. Negative floor is 0.0002976497842924418 Ha. Eigenpair Frobenius residual 3.624148053680244e-14; orthogonality Frobenius residual 2.5268021680905354e-12; independent SciPy spectrum difference 1.3048146537550791e-14.
- Shape (15,57,19) has rank cap 285 and bound 0.021548084729907786 Ha, so its REJECT classification is separated from 0.020 by more than the stated 1e-8 guard. The other eleven table classifications and all twelve README table rows agree with regenerated JSON.
- All 482 resource rows cover exactly the twelve shapes, five precisions and every positive legal batch width. They reproduce the official state-preparation/Givens model and full reference loader cost. An independent sweep of 576 more cases, including B=1, powers of two and adjacent B values, matches all seven returned PreparationCost fields. The separate ANGLE-cap boundary probe passes. Live official source confirms the separate Toffoli/scratch exponents and special power-of-two handling.
- Public C30/10-bit factor metrics reproduce exactly: error_1=0.01391173198501434, error_2=0.019980679640672218, Lambda=290.7542054993138, E_gap=3.9159356215095613, lambda_eff=47.55854233515931; official gates pass. Two-body margin is 0.01932035932778256 mHa; strict gap-floor margin is 0.0010376215095613261 Ha. The research baseline's five parameter arrays equal the sample arrays.
- All 25 pinned released circuits build within Q128 and pass their emitted static table/ANF assertions. Gauge drift is 7.771561172376096e-16 up to row signs. Both-pass static ANGLE count 56450, full Toffolis 65244, block-encoding Q661, conditional cost 2051023203.3721037 and gain 1.6154555940004633 reproduce exactly.
- Complete structural equality holds between each regenerated JSON and its original: inputs.json, comparison.json, released-loader.json. All thirteen target SHA256s remained unchanged through review.
- Official/starter-kit inputs and all nine inspected starter files match their recorded SHA256s; public dataset metadata pins starter version 4 and public words version 2. Public notebook metadata pins version 8 and kernel ID 140416287; the JSON source string SHA256 equals e4793278e7fd7f6828cb38a9d0983d2061c91a101e68a9fbb19c3ec44ae962a4. Historical release Git blob hashes reproduce after removing the one extra final newline present in the provided fetched-text audit copies.
- Installed current documentation Details/Options and live kernel behavior agree for Eigensystem eigenvectors-as-rows, eigenvalue magnitude ordering, explicit matrix Frobenius Norm, BinaryReadList Real64/ByteOrdering and file closing, and PolynomialMod coefficient reduction.

Explicit scope exclusions and untested claims: historical published sampling counts and certified scores are attributed literature/source statements, not local score validations. Their source contracts/count descriptions were checked against the supplied pinned release text, but their large experiments were not rerun. Full official forward/reverse circuit shot verification, dirty-target/measurement-phase correctness, new factor training and new verified competition scores are explicitly not claimed. The numerical guard remains a diagnostic guard, not an interval proof. These exclusions are consistent with the documented audit-only scope and are not open findings.

# FINDINGS

## [P2] Accepted zero-imaginary complex representations are silently discarded by PSDTail

Evidence: RankCost.wl:11,13-15. realNumberQ accepts a machine complex expression with zero imaginary part because NumericQ[x] is True and Im[x]==0. However Positive and Negative both return False for this representation, so both Select calls omit it.

Fresh-kernel reproduction:

```wolfram
NumericQ[1. + 0. I]                 (* True *)
Im[1. + 0. I] == 0                 (* True *)
Positive[1. + 0. I]                (* False *)
Negative[1. + 0. I]                (* False *)
DFTHCFeasibility`PSDTail[{1. + 0. I}, 0]   (* 0 *)
DFTHCFeasibility`PSDTail[{-1. + 0. I}, 0]  (* 0 *)
```

Both mathematical inputs are real unit modes, and both rank-zero Frobenius errors must be 1. The negative mode must also remain as a unit floor at every permitted rank. Returning Failure for the unsupported representation would also honor a restricted contract; returning zero silently does not. Current official real-matrix numerical results are unaffected, but the exposed public numeric contract is incorrect at this boundary.

# RECOMMENDED FIXES

Normalize validated real-valued spectra with Re before applying Positive/Negative and sorting, or reject Complex-head values explicitly in the validator. Add regression cases for positive and negative machine and arbitrary-precision zero-imaginary representations, including rank zero and a large rank. Re-run correctness after the executable change.

## SHA256 binding of reviewed target bytes

| File | SHA256 |
|---|---|
| RankCost.wl | 9a478208770a281ab2d79014371d48d8f107cce14bab2f20e2016f0018cee1e4 |
| run.wls | a9e295e743bd84abf44a642376706f2a542d823bceb0e89a195c97810b922757 |
| prepare.py | 646b4d94e58fec92bd6bd0173676f4a94dce047625750fd98e2f42960a780d50 |
| audit_loader.py | b7947facc9d2eda9e93038002156fa598eda35ed63b5a467907aaf3bdc843184 |
| DOCUMENTATION.md | 58a1a5d9fd581588a9d12b6b89a399ddda4451e5286c43be2c678075a9ead530 |
| README.md | 7be43f16dc13c409ee2abda9351fa83cfd765e36f7c8c80c44a39af73aacc21b |
| AUDIT.md | 77baa2a09dc91e4d3f010e82d2f694a344b049b045876f019f2e413d5a378195 |
| PLAN.md | 8a9e8571aa4cb99c0c45efbd6bd5534555dda04e1219f4ead4b091249b001bbf |
| inputs.json | e6d37f4bbf685562b9e837f9abd48de91e3ca5c2666b57e9140e99aa41e6906f |
| comparison.json | 87872061a3c9b575775c7a2bf7202f8e8a2092030fc077afe932b90313b6fac3 |
| released-loader.json | a045dcccfecfe6bd58746aebeb2ac0a6cffb0444045d52e2016dcd2c9dbbef0f |
| released-sources.json | 7463fd40fa460055507186c7e989e63ffe9ebc424c9bf5e8a99eda0ef09ccca9 |
| public-release-provenance.json | 34eb48a8887bceda16847b21ff4186dce2b3a5cc5ee64306129375748aa6bde2 |

OPEN ISSUES: 1
