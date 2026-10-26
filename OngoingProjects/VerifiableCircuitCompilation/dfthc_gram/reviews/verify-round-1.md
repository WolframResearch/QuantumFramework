# VERDICT

ROUND: 1/3. The mathematical feasibility results, frozen-array certificates, all eight candidate metrics, numerical screens, and conditional DFTHC resource counts reproduce. One contextual three-task estimate lacks independently reproducible provenance. No mathematical or numerical defect was found in the studied DFTHC results.

All work ran in an isolated temporary replica containing only the named artifacts and permitted dependencies. No circuit implementation, circuit verifier, optimization-budget extension, commit, or push was performed. The target artifacts were not edited.

# REPRODUCED

## Claim inventory and coverage

- **Core representation and invariants:** the augmented core blocks, two-body reconstruction, offset contraction/norm, PSD quadratic form, orthogonal-copy invariance, exact rank example, and the same-core/different-LCU-normalization counterexample. The official math source confirms the native-offset conversion, one-body absorption, SOS energy, normalization, and effective-normalization formulas.
- **Schur and projection statements:** singular range condition `q2 == 0 && tau >= q1^2`, rejection outside that range, normal equations, nonnegative Schur projection residual, the exact `7/5` epigraph optimum, and its independent numerical SDP value `1.399999996255239`.
- **Dictionary calculus and obstruction:** projector Gram identity, Euclidean gradient, eliminated-weight gradient, left inverse, rank-obstruction minor/support, noncommuting real projectors, identity-rotation limit, and vanishing-error limit. The full-rank/common-support argument supports the stated per-block obstruction; it does not exclude approximate or globally redistributed representations.
- **Offset compression and Schur rescaling:** all eight compression identities and all six rescaling/convexity checks, including the coupled residual/gradient, offset/norm preservation, copy-rank bound, PSD Hessian, and nonnegative boundary optimum. Together with the 23 core checks, all 37 exact checks pass.
- **Frozen data:** all 65 orbital frames have rational modular rank 60, all 65 augmented arrays have rank 30, and every block has a fully supported copy. All recorded singular values and reconstruction residuals reproduce. These certificates concern the frozen Float64 arrays interpreted as binary rationals, as stated.
- **Eight pilots:** every official metric (`error_1`, `error_2`, `Lambda`, `E_SOS`, `E_gap`, `lambda_eff`, and domain flag), gate flag, fit history/cap, quantization/refit flag, and conditional cost row reproduces. The B48/B32 histories each contain 50 monotone steps. The offset pilot reaches 200 iterations without claiming optimizer convergence or global optimality.
- **Bounds:** all 13 global PSD rank floors and their guarded decisions reproduce. The quantized B32 fixed-dictionary bound reproduces with the stated global-target/fixed-direction scope. Neither screen establishes impossibility for reoptimized dictionaries or for all shapes that pass the relaxation.
- **Resources/conclusions:** all 21 ordinary-unary rows are the minimum modeled cost over legal batches under Q128, and the unchanged-loader row matches the pinned static loader/preparation model. The incumbent and conditional DFTHC objective products, 156-Toffoli saving, `0.1519406375%` modeled worsening, zero passing candidates, and absence of a new score reproduce. The three-task combined estimate is the sole unverified numerical claim, detailed below.
- **Public function signatures and boundaries:** source and live-kernel checks cover `AugmentedCore[w,v]`, `ProjectorGram[u,v]`, `ProjectorGradient[u,w,u0,w0]`, and the three zero-argument check functions. Exact/radical, zero, rectangular, empty, ragged, complex, nonfinite, symbolic-input, and incompatible-shape cases behave consistently with their numeric-real domains. The tau-zero Schur boundary forces the offset contraction to vanish. No notebooks or runtime benchmark claims are present.

## Reproduction commands and outputs

The following commands ran from the replica's `dfthc_gram` directory. `PY`, `WL`, and `DATA` denote the supplied Python executable, Wolfram kernel, and official bundle respectively; private absolute paths are intentionally omitted.

```sh
"$PY" study.py prepare --data "$DATA"
"$WL" -noinit -noprompt -script run.wls
"$PY" -m unittest test_study -v
"$PY" study.py run --data "$DATA"
"$WL" -noinit -noprompt -script verify-data.wls
```

All five commands exit 0. Preparation reports `PINNED_INPUTS`; WL reports `PASS`, four passing all-data checks, and the SDP value above; all five Python tests pass; the numerical study reports eight candidates and zero metric passes; the final WL cost check reports `PASS`, 21 cost rows and one unchanged-loader row. No Wolfram messages were emitted. The regenerated `inputs.json`, `symbolic.json`, `results.json`, and `wl-data-checks.json` are byte-identical to the targets.

Versions: Wolfram `15.0.0 for Mac OS X ARM (64-bit) (February 17, 2026)`, Python `3.13.1`, NumPy `2.5.4`, SciPy `1.18.1`.

Additional verifier-only checks ran in fresh kernel/Python processes:

- A direct symbolic derivative of explicitly flattened rectangular rank-one factor matrices equals the implemented gradient. Live probes confirm `MatrixRank`'s `Modulus` option, `CharacteristicPolynomial`'s sign convention, exact binary rationalization, the prime modulus, and the boundary cases above. One initial verifier-only `SameQ` probe compared factored and expanded polynomials structurally; replacing it with exact polynomial equality resolves that harness error.
- Independent Python modular elimination, using `float.as_integer_ratio()` and modular inverses rather than WL rank calls, reproduces all 130 rank certificates. The smallest released offset norm is `1.1089580842005626`, confirming the nonzero-offset regime.
- Direct symmetric-pair reconstruction reproduces all eight official two-body errors within `4e-15` Ha. The offset candidate keeps every incumbent angle word identical, with coefficient drift `8.881784197001252e-16` and norm drift `2.842170943040401e-14`.
- An independent eigenvalue driver reproduces all 13 PSD tails within `1e-10` Ha. Pivoted QR, changing the dictionary pivot order, gives the B32 bound `0.13469057200509504` Ha, minimum absolute QR diagonal `0.020112784339302938`, and orthogonality residual `1.3654382401560638e-13`.
- The offset candidate reproduces error `0.020271423640563186` Ha and normalization `47.744962460733106`; its conditional unchanged-loader totals are `65088` Toffolis, `661` qubits, and cost `2054139541.1018138`. Its NNLS KKT residual is `3.469446951953614e-17`.

# FINDINGS

1. **[P2] The displayed three-task estimate is an unpinned constant.** `OngoingProjects/VerifiableCircuitCompilation/dfthc_gram/study.py:347` inserts `current_combined_local_estimate = 2.035845221582375` without loading component costs or an aggregation formula. `RESULTS.md:3` presents its rounded value as the current three-task local estimate. The permitted provenance covers the DFTHC incumbent, official math inputs, and shifted matrix only. The permitted historical C30 report contains a different combined estimate, `2.0088882102784944`; that does not prove the new constant wrong, but it cannot substantiate it. Re-running the script reproduces the literal rather than independently verifying its value. All DFTHC-only counts and conclusions remain reproducible.

# RECOMMENDED FIXES

Remove the unrelated three-task estimate from `study.py`, `results.json`, and `RESULTS.md`, or pin all three task costs and the exact aggregation formula and derive the value from them. Regenerate the hash-bound reports after a source change. No additional fitting or optimization is needed to resolve this evidence issue.

OPEN ISSUES: 1
