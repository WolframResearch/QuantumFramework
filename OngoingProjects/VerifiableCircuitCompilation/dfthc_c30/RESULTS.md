# Fixed-factor C30 results

The standalone C30 four-address family passed the unmodified official
DFTHC verifier: **254464 shots**, 3976 legal addresses, 64 measurement shots
per address, and 25 circuits in both forward and reverse directions.
No sampled value, phase, address or scratch violation was accepted. See
[verified.json](verified.json) for the full official breakdown and input hashes.

| Quantity | Verified result |
|---|---:|
| DFTHC cost (lower is better) | 2051023203.3721037 |
| Toffolis, full block encoding | 65244 |
| Qubits, full block encoding | 661 |
| Maximum ANGLE width | 128 |
| One-body error, mHa | 13.911731985012415 |
| Two-body error, mHa | 19.98067964067221 |
| Effective normalization | 47.55854233515931 |
| Improvement over original DFTHC baseline | 1.6360612511835888x |
| Improvement over prior verified DFTHC incumbent | 1.6154555940004633x |
| Cost reduction from prior incumbent | 38.09795801792165% |

The public weights and offsets are exactly fixed, as are the official shifts
and Givens tree. Only the released row-sign gauge is applied. Maximum vector
drift after accounting for signs is 7.771561172376096e-16; the unit-norm error
is 6.661338147750939e-16. No factor training was performed. The two-body
error has only 0.0193203593 mHa of headroom, so this result is bound to the
specified factor file, transformed words and circuit hashes.

The combined local estimate is now **2.0088882102784944**, up from
1.7120801622097384. It combines verified baseline cost ratios:
inverse square root 3.9385957568764147x, DQI 1.2581323878369532x,
and DFTHC 1.6360612511835888x, using their geometric mean. This remains
well below 8.43991. No new Kaggle submission or leaderboard result is claimed.

Supplementary checks passed 795200 released and 4704 synthetic circuit
applications, using dirty target words and checking involution, phases,
selection and clean workspace. WL passed 17 exact symbolic/edge checks and
six full-data checks, including all 298200 angle words. Measurement sampling
does not certify every possible branch.

The preceding audit, DQI integration and feasibility folders were pushed to
main as [df53dba](https://github.com/WolframResearch/QuantumFramework/commit/df53dba1d898082fcd2600a46fed1d505b00c0e3).
This new loader and its verification results are local. Review evidence is
recorded in `reviews/` and the final `evidence-manifest.json`.
