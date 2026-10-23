VERDICT

PASS for the stated feasibility scope. No new submission or full repeated-shot circuit verification is inferred. ROUND 1/3.

REPRODUCED

- Read the complete feasibility source, prose and structured artifacts, the official arithmetic generator/verifier, matching immutable baseline/incumbent reports, and installed documentation reading copies for the exercised signatures. Original artifact hashes agree with their recorded provenance.
- Ran copied `run-feasibility.wls` and `verify-feasibility.wls` with `WolframKernel -noinit -noprompt -script`, from outside the source directory. Runtime: Wolfram Language 15.0.0, Mac OS X ARM, February 17, 2026. No WL messages occurred. Every regenerated design/phase JSON equals the retained original semantically; 25/25 WL checks passed.
- All 348 intervals independently regenerated, including all 1,044 exact relative-error, derivative and Horner-register certificates. Exact coefficient rounding, interval coverage, Horner and input-truncation budgets, Newton identity/limit/uniform relative bound, rational square-bit drop bound and its attainment, and row-versus-whole-product witness reproduce.
- Independently derived the arithmetic perturbation with ell=sqrt(m): `(1+e)(ell^2 ds+dm)/2 - ell dp/2`. Exact Resolve independently proves each exported A arithmetic allowance for `1<=ell<=2`, the certified seed-error range, square omission `0<=ds<=SquareDropBound`, `0<=dm<=20*2^-b`, and `0<=dp<=w*2^-b`. Therefore the reported Newton bounds remain conservative for their specified arithmetic model.
- Live MiniMax return shape, its WorkingPrecision option, exact Resolve domains, BooleanConvert ANF/ESOP, BooleanMinimize full-adder identity, PolynomialReduce over rational coefficients, integer-domain LinearOptimization syntax and result shape all reproduce. Tested malformed/inexact/symbolic/nonpolynomial certificate inputs, invalid interval/eta boundaries, zero polynomial, no omitted square bits and no Horner rows.
- Both B relative bounds and margins reproduce: cubic `32285/1099511627776 = 2.9363036447e-8`, quadratic `40489/1099511627776 = 3.6824530980e-8`; each is below 2^-24. Their remaining fractions are approximately 50.7369995% and 38.2186890%.
- General parity identities for sizes 1 through 6, two-/three-bit identities, arbitrary real angle diagonal unitarity, distinct Boolean/phase coefficient domains, modulo-one reductions and the exact zero-state phase reproduce. An independent rational monomial expansion checks every exported term in all 12 phase blocks against the retained partial products, modulo integer turns. Additionally, the actual C all-exponents circuit matches the incumbent integer model exactly on 38 boundary/exponent inputs with two independent HMR seeds; phase, restored x and clean scratch all agree.
- Ran the original `resource_probe.py --data "$QF_INVSQRT_DATA" --output review-workdir/resource-results.json` with the prescribed Python environment. Its entire result equals the retained result after removing elapsed_seconds only. This reproduces compiler operation counts, internal ancilla charging, padded joint coefficient lookups, disjoint incumbent stage attribution, zero classically conditioned priced operations, sampled executed/static equality, synthesis prices and all resource estimates.

| Probe | Q | Toffolis | Fine rotations | Exhaustive model |
|---|---:|---:|---:|---|
| A quadratic nr36 | 300 | 10497 | 44 | Pass |
| A quadratic nr28 | 260 | 8361 | 36 | Pass |
| A linear nr32 | 256 | 8256 | 40 | Pass |
| A linear nr28 | 236 | 7180 | 36 | Exactly one failing x: 1048575 |
| B cubic | 256 | 6638 | 43 | Pass |
| B quadratic | 220 | 4449 | 43 | Pass |
| C fixed s=0 diagnostic | 217 | 5081 | 403 | No full-task eligibility claim |
| C all exponents | 217 | 5099 | 8875 | No full circuit verification claim |

- All six A/B variants cover all 3,133,452 inputs with zero register-capacity violations. The compiler/model samples reproduce the reported 81 distinct quadratic-seed inputs and 93 for each other design, two seeded shots each. Predicted objectives and improvement ratios reproduce, including B quadratic objective 1,043,220.2106 and C hypothetical half-margin objective 16,399,035.5275.
- Prior baseline/incumbent full-verifier PASS records were audited against current circuit/verifier hashes; their historical 200,540,928-shot executions were not rerun because they are explicitly immutable inputs to this comparison. New full repeated-shot verification is expressly unclaimed. No runtime/performance speedup claim requires benchmarking.

FINDINGS

None within the stated feasibility scope. The distinction between exact arithmetic/model evidence, sampled circuit evidence, conditional costs and prior full-verifier evidence is accurate. Architecture choice remains deferred.

RECOMMENDED FIXES

None required for correctness in this round. Independent scripts and copied/regenerated evidence are retained at `review-workdir/`.


Publication note: local paths have been replaced with portable references.
The findings and quantitative review evidence are unchanged; supplemental
temporary reviewer scripts are not included in the public folder.

OPEN ISSUES: 0
