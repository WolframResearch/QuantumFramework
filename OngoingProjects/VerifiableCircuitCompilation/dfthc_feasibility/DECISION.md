# Decision after released-design audit and WL feasibility comparison

Choose the public **(65,60,30), ten-bit signed-Givens four-address loader**
as the next standalone implementation. Its fixed factors already pass the
official metric gates, and the released loader builds at the legal Q128
boundary. Attribution, source versions and hashes are in [AUDIT.md](AUDIT.md)
and `public-release-provenance.json`. This is external public work audited
here; it is not a newly discovered factorization.

| Quantity | Current full-verifier incumbent | Released design: static projection |
|---|---:|---:|
| Effective normalization | 48.68043460482008 | 47.55854233515931 |
| Total Toffolis | 93237 | 65244 |
| Total block-encoding qubits | 730 | 661 |
| Cost | 3313336907.312215 | 2051023203.3721037 |
| Gain over incumbent | 1 | 1.615455594 |

The projection is a 38.097958% cost reduction. It has **not** passed full
circuit verification. The metric-verified public two-body error is
19.980679640672218 mHa against a 20 mHa limit, so preserve its factors,
words, and tree while establishing that circuit baseline. A future standalone
port must preserve the Apache-2.0 credits and test arbitrary target words,
measurement phase repair, batch chaining, and all official forward/reverse
shots before claiming a verified improvement or publishing a submission.

Our combined local estimate remains **1.712080162**. If the projected cost
survives full circuit verification, the combined estimate becomes
**2.008888210**, still far below the requested 8.43991 benchmark. No Kaggle
submission or leaderboard score is claimed by this audit.

## What the WL comparison changes

The fixed-shift target needs at least **290** unrestricted positive modes at
20 mHa. The paper's (15,57,19) shape supports at most 285 and is numerically
rejected: even the relaxed minimum error is 21.548084730 mHa. Do not train or
port that shape under this competition contract.

The other eleven shapes are only **not ruled out** by the necessary bound.
For larger-gain research, prioritize bounded pilots of **(20,48,24)** and
**(32,32,16)** after establishing the released-loader baseline. Their relaxed
error floors are 1.629534 and 1.139588 mHa. At beta=12 their ordinary-unary
reference resource projections are approximately 3.568x and 3.563x cheaper
than the incumbent, conditional on retaining its normalization. No accurate
fit or new normalization is established for either shape. Include nearby
address-width transitions in a later search rather than treating these
sampled ranks as globally optimal.

Future fitting should use the derived projector Gram identity
`G[b,d] = (u_b . u_d)^2`, with the official fixed shifts, quantized residuals,
gap, normalization and full T*Q cost all evaluated together. The promising
public BM decoder port already handles the competition field and variable
mask; audit that port before rebuilding it from the older fixed-mask source.
Its full correctness and cost have not been reproduced here.

## Evidence and scope

All **32** built-in WL checks and all **482** full-cost comparisons pass.
Fresh correctness reviews converged on the final source; the final fresh
wl-quality round passes **A–G**, with **OPEN ISSUES: 0**. See
[review history](reviews/PIPELINE.md), [comparison](README.md),
[mathematics and regimes](DOCUMENTATION.md), and `evidence-manifest.json`.

This stage produced no factor training, new submission files, commit or push.
Existing inverse-square-root, DQI and DFTHC verification evidence is unchanged.
Raw matrices, spectra and external factor/source data remain outside the
public artifact set. The public files were checked for personal paths and
credential markers.
