# VERDICT

PASS for the stated mathematical-feasibility scope. The entire README reproduction succeeded in an isolated temporary replica using only the named artifacts and direct dependencies. The regenerated `inputs.json`, `symbolic.json`, `results.json`, and `wl-data-checks.json` are byte-for-byte identical to the target files. No target source or evidence file was changed.

# REPRODUCED

All paths below are relative to `OngoingProjects/VerifiableCircuitCompilation/`.

**Commands and environment.** From the replica's `dfthc_gram/`, using the supplied Python environment and installed Wolfram kernel:

```text
python study.py prepare --data DATA
WolframKernel -noinit -noprompt -script run.wls
python -m unittest test_study -v
python study.py run --data DATA
WolframKernel -noinit -noprompt -script verify-data.wls
```

All five commands exited 0. Salient outputs: pinned-input preparation succeeded; WL returned `PASS` with all four data checks true and SDP value `1.399999996255239`; all five Python tests passed; numerical execution returned eight candidates and zero metric passes; final WL verification returned ten true checks, 21 reference-cost rows and one unchanged-loader row. No WL messages occurred. Versions: Wolfram 15.0.0 for Mac OS X ARM (February 17, 2026); Python 3.13.1; NumPy 2.5.4; SciPy 1.18.1. There are no runtime benchmark claims or notebooks in scope.

**Claim inventory and exact mathematics.** Read all twelve target artifacts completely. Independently checked the formulas and reproduction results for:

1. Augmented-core blocks, pair-tensor reconstruction, offset contraction, offset norm, PSD quadratic form, rank/minimum-factor interpretation, and orthogonal-copy gauge invariance.
2. The singular Schur range condition, infeasibility outside the range, and the nonsingular Schur epigraph optimum `7/5`.
3. Projector Gram contractions, eliminated-weight normal equations, Schur projection residual/nonnegativity, Euclidean projector gradient, and the envelope-gradient use after least-squares elimination.
4. The full-support factor rank obstruction using a left inverse, its support/minor identities, and the limited applicability to exact unchanged individual blocks.
5. Same-core/different-normalization example (`1/2` versus `1`), noncommuting projectors, identity-rotation limit, and vanishing-projection-error limit.
6. Offset-orthogonal removal: preserved contraction/norm, PSD projector, one-copy rank reduction, rank-one core subtraction, coupled quartic objective, and both gradients.
7. Positive Schur rescaling: core identity, preserved offsets/norm, convex nonnegative scale objective, PSD Gram Hessian, and nonnegative-fit boundary. The stated `tau>0` condition and separate `tau=0` caveat are correct.
8. Invalid core/projector/gradient shapes and the declared numerical-domain guards.

The 23 `GramSymbolicChecks` flags, eight `OffsetCompressionChecks` flags, and six `SchurRescaleChecks` flags all reproduced exactly: **37/37**.

**Live WL signatures and boundary cases.** Inspected the three public functions' live downvalues, `MatrixRank` and `SingularValueList` options, and exercised empty/ragged/mismatched input, complex/nonfinite/symbolic rejection, exact radicals, zero cores and gradients, and rectangular projector output. All passed without messages. Independently confirmed prime modulus `2147483647`, the exact binary-rational interpretation of `SetPrecision[0.1,Infinity]`, the `Det[m-t I]` characteristic-polynomial convention, and the zero-scale PSD boundary.

**Frozen-array claims.** Reproduced all 65 complete records, including modular ranks `rank(U)=60`, `rank([w;v])=30`, full-support-copy flags, singular values and four reconstruction diagnostics. An additional exact modular check found `rank(w)=30` in every block. Every offset is nonzero (minimum squared norm `1.2297880325137824`); maximum unit-row norm-squared error is `1.11e-15`. These substantiate the frozen-binary-rational statements only; no exact unrounded-trigonometric claim is made.

**Derivative and fit claims.** All five small-array tests passed, including independent finite differences for variable projection and coupled removal, greedy-projection monotonicity, and WL/Python gradient agreement. An additional fresh-seed variable-projection finite-difference check had maximum absolute error `3.83e-8`; a singular dictionary was rejected explicitly. Dictionary normal-equation residuals after quantization were between `2.98e-15` and `4.67e-15`. Both 50-step histories, all accepted-step counts, and the 200-iteration cap with unsuccessful optimizer status reproduced exactly. No optimization budget was extended.

**Official decoded-word metrics.** Checked the native-offset conversion and every used metric/cost signature against official implementation source. A second evaluation of every saved candidate through the official metric entry point matched every reported metric exactly. All word values are integers in `[0,1023]`; saved two-body vectors exactly match their decoded words. Reproduced two-body errors and effective normalizations:

| Candidate | Error (mHa) | Effective normalization | Metric pass |
|---|---:|---:|---|
| dictionary-B59 | 72.09889194326 | 59.38268566617 | False |
| dictionary-B56 | 214.04513341035 | 95.48011194989 | False |
| dictionary-B48 | 527.37632731318 | 163.37160966515 | False |
| dictionary-B32 | 1604.34061118401 | 304.67853117365 | False |
| copies-C29 | 20.33738341688 | 51.60183484209 | False |
| copies-C24 | 24.58318657568 | 51.56620827737 | False |
| copies-C16 | 62.35583636933 | 53.24262751023 | False |
| offset-preserving-C29 | 20.27142364056 | 47.74496246073 | False |

The incumbent's seven metric fields, including `19.98067964067221` mHa error, also reproduced. Offset-preserving C29 retains every incumbent word, misses the error gate by `0.27142364056` mHa, has KKT residual `3.469446951953614e-17`, contraction drift `8.88e-16` and norm drift `2.842170943040401e-14`. Its predicted and official physical errors agree within `1.5e-16` Ha. All remaining metric fields, fit fields and flags in `results.json` matched, not just this table.

**Necessary bounds and scope.** All thirteen PSD rank caps/floors and their screens reproduced exactly. A separate full-spectrum calculation using the `evd` LAPACK driver differed by at most `1.89e-17` Ha. The `(8,32,32)` floor is `33.76848165684` mHa, so its rejection is well beyond the declared `1e-8` Ha guard. The full-size rank floor is only `0.29764979020` mHa. An independent pivoted-QR fixed-direction calculation gives `134.69057200502` mHa, versus recorded `134.69057200510` mHa; span dimension is 2080 and orthogonality residual is `1.37e-13`. The stated restrictions to numerical screening, fixed proposed directions, and absence of a theorem excluding reoptimized/global representations are appropriate.

**Resources and recorded context.** Reproduced every field of all 21 conditional loader rows, the legal-batch minimizations, complete official preparation/Givens costs, shared temporary-space maxima, and the ANGLE Q128 cap. Static loader profiles, incumbent/build hashes and all six recorded PASS-report identities agree. The unchanged-loader C29 estimate is `65088` Toffolis, `661` total qubits, width `128`, and cost `2054139541.1018138`; it saves 156 Toffolis but is `0.15194063746%` costlier. The incumbent cost is `2051023203.3721037`; the ordinary-unary C30 estimate is `72378` Toffolis, `731` qubits and `2516242481.485071`. Shape-only normalization assumptions and rejected-candidate estimate labels are explicit and correct.

The six report costs give component ratios `4.099287418831162`, `1.2581323878369532`, and `1.6360612511835888`; their geometric mean is `2.035845221582375`. This is independently reproduced in WL from the recorded reports. Historical circuit outcomes are accepted only as identified recorded evidence, as requested; those circuit runs were not repeated. No circuit implementation or new competition score was produced.

# FINDINGS

None within the declared scope. The claims that no feasible improvement was found among these eight bounded pilots, that the incumbent and recorded combined context remain unchanged, and that global approximate refactorization remains open are supported. Recommendations about future global factorization are research judgments, not claims of demonstrated feasibility or optimality.

# RECOMMENDED FIXES

None required.

OPEN ISSUES: 0
