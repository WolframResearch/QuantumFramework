# DFTHC rank and cost feasibility

Audit and necessary-condition screening only. No new verified score or submission.

The 2926-dimensional shifted pair matrix needs at least **290** unrestricted PSD modes at 20 mHa. Rank 289 has minimum error **20.1873081 mHa**; rank 290 has **19.8584258 mHa**. The negative-spectrum floor is **0.297649784 mHa**. This is a numerically validated bound with an explicit guard, not an interval-arithmetic certificate.

All **482** cost rows across twelve shapes, five angle precisions and legal batch widths agree between WL and the official Python cost model. Eigenpair residual 3.63e-14, orthogonality residual 2.53e-12, independent SciPy spectrum difference 1.31e-14. All 32 symbolic/edge checks pass, including explicit acceptance at Q128, rejection above the ANGLE cap, and normalization-domain checks.

These are conditional ordinary-unary-loader resource projections at fixed lambda_eff=48.68043460482008. Every row below selects beta=10, J=10; those precisions have not been shown accurate for hypothetical factors. Gains compare with the incumbent DFTHC cost 3313336907.312215. NOT_RULED_OUT is only a necessary spectral screen, not a feasible fit.

| R,B,C | Rank cap | Error lower bound (mHa) | Screen | T total | Q total | Conditional gain |
|---|---:|---:|---|---:|---:|---:|
| 65,60,60 | 2926 | 0.297650 | NOT_RULED_OUT | 76741 | 734 | 1.208x |
| 65,60,30 | 1950 | 0.297650 | NOT_RULED_OUT | 72378 | 731 | 1.286x |
| 15,57,19 | 285 | 21.548085 | REJECT | 19525 | 722 | 4.828x |
| 16,57,19 | 304 | 15.542820 | NOT_RULED_OUT | 20489 | 722 | 4.601x |
| 20,48,24 | 480 | 1.629534 | NOT_RULED_OUT | 21724 | 725 | 4.321x |
| 24,48,24 | 576 | 0.640356 | NOT_RULED_OUT | 25038 | 726 | 3.744x |
| 32,32,16 | 512 | 1.139588 | NOT_RULED_OUT | 22732 | 693 | 4.321x |
| 32,48,24 | 768 | 0.309572 | NOT_RULED_OUT | 31644 | 726 | 2.963x |
| 40,48,32 | 1280 | 0.297654 | NOT_RULED_OUT | 39037 | 730 | 2.388x |
| 48,48,32 | 1536 | 0.297650 | NOT_RULED_OUT | 45797 | 730 | 2.036x |
| 56,48,40 | 2240 | 0.297650 | NOT_RULED_OUT | 53564 | 733 | 1.734x |
| 65,48,30 | 1950 | 0.297650 | NOT_RULED_OUT | 59802 | 731 | 1.557x |

No smaller-rank fit is produced by this work. The reference unary loader is distinct from the released four-address loader audited in `released-loader.json`; see [AUDIT.md](AUDIT.md) for its directly reproduced metric gates and static projection.

Read [DOCUMENTATION.md](DOCUMENTATION.md) for the spectral argument, API contracts, numerical limitations, and regimes. Read [PLAN.md](PLAN.md) for the bounded workflow. Independent review results are saved under `reviews/`.

## Reproduction

Use a Python environment with kickmix, NumPy and SciPy, and a Wolfram 15.0 kernel. Supply the official competition download externally. Download the public release NPZ (version 2) and the official FeMoCo research baseline (starter-kit version 4) into an audit directory. `prepare.py --data DATA --released AUDIT_DIR` writes ignored matrix/spectrum caches and `inputs.json`. Then run `WolframKernel -noinit -noprompt -script run.wls` from this directory or by absolute path. The script locates its own files.

For the loader profile, extract the three hash-pinned Apache-2.0 modules from public notebook version 8, preserving its credits, into a separate audit directory. Run `audit_loader.py --data DATA --source SOURCE_DIR --words PUBLIC_WORDS_NPZ`. This builds the released circuits in memory, writes the aggregate profile, and does not export a submission. Full official circuit verification is a later required step.

No raw integral matrix, factor array, full eigenvector/eigenvalue file, downloaded source, personal path or credential belongs in the public artifact set.
