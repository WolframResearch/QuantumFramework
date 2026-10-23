# Independent wl-quality review, round 1/3

Verdict: PASS for the stated architectural-feasibility scope. No source artifact was changed. No new circuit is certified for the official repeated-shot task by this review.

## Scope and authority

Read the complete `Feasibility.wl` (88 lines), both WL runners, `resource_probe.py` (590 lines), `COMPARISON.md`, `DOCUMENTATION.md`, and `PLAN.md`. Loaded all of `designs.json`, `phase-terms.json`, `wl-checks.json`, and `resource-results.json`; inspected all rows structurally and rederived every interval certificate. Read the full directly used official `nr_reference_circuit.py` and `verify_submission_invsqrt.py`, the relevant installed documentation and live kickmix docstrings. Applicable repository instructions were checked. The reviewed code uses core WL and `FunctionApproximations`, with no new QuantumFramework API.

All source names below refer to `OngoingProjects/VerifiableCircuitCompilation/feasibility`. Supplemental review evidence is in `review-workdir`.

Runtime: Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. Used the real kernel with `-noinit -noprompt -script`; messages were unsuppressed. No `Quiet` or `Print` was used in the WL review scripts. Original scripts/JSONs were copied to the review folder before execution.

Reproduction:

```bash
WolframKernel -noinit -noprompt -script review-workdir/run-feasibility.wls
WolframKernel -noinit -noprompt -script review-workdir/verify-feasibility.wls
WolframKernel -noinit -noprompt -script review-workdir/independent-quality.wls
python3 OngoingProjects/VerifiableCircuitCompilation/feasibility/resource_probe.py --data "$QF_INVSQRT_DATA" --output review-workdir/resources-reproduced.json
python3 review-workdir/phase-circuit-samples.py
```

Full design and phase-term regeneration matched the submitted bytes exactly. Reproduced resource evidence matched every field after excluding elapsed times. All stored WL, official-source, and resource-probe hashes matched current bytes. The supplied WL check runner exited 0 with `passed -> 25, failed -> {}`. All supplemental independent checks in `independent-quality.json` are True, including a separately reconstructed exact polynomial test for every one of 348 intervals.

## A. Documented symbols — PASS

Evidence: `Feasibility.wl:31-39` uses `Needs["FunctionApproximations`"]`, `MiniMaxApproximation[expr,{t,{0,h},degree,0},WorkingPrecision->50]`, and extracts `approx[[2,1]]` and `approx[[2,2]]`. The installed `MiniMaxApproximation` reference page, read via `installed-documentation/MiniMaxApproximation.md:3-28`, explicitly documents this argument form, relative error objective, return `{abscissa,{approximation,maxerror}}`, and WorkingPrecision. Live kernel probes confirm the `{List,List}` outer heads, degree 2, error precision 50, and WorkingPrecision option. The package implementation `FunctionApproximations/Approximations.m` also confirms relative error and iterative discovery behavior.

`CertifyPolynomial` at `Feasibility.wl:19-26` follows the installed Resolve authority: quantified polynomial inequalities over Reals, exact rationals, and no finite WorkingPrecision. The installed Resolve page warns that finite WorkingPrecision may yield an incorrect result (`installed-documentation/Resolve.md:734-751`); the source avoids that option. The independent squared-certificate equivalence returned True for all `0<eta<1`, `p>0`, `m>=1`.

`verify-feasibility.wls:61-65,76` agrees with the installed BooleanMinimize ANF form, PolynomialReduce rational coefficient/remainder shape, and LinearOptimization integer-domain syntax. Supplemental kernel checks reproduced equivalence, zero rational remainder and integer optimum 4. LinearOptimization remains an API exercise, never an exact architecture-optimality certificate. Live kickmix docstrings confirm clean-target joint lookup and deletion, little-endian order, modulo-register addition, Z_POW exponent conventions, INFO_TAG marks, and static operation statistics; the complete resource runner executes those APIs successfully.

Defect: none. Minimal fix: none.

## B. Idiomatic functional WL — PASS

Evidence: `Feasibility.wl:14-17,32-48,59,67-68,73-85` uses direct list arithmetic, Map, Select, FoldList, Lookup, Catenate, GatherBy, and declarative construction. Horner evaluation is represented by Fold in `verify-feasibility.wls:15-17`, directly expressing its recurrence. The source search found no Do, For, While, AppendTo, Quiet, or Print in the WL artifacts. The remaining Table constructions at `Feasibility.wl:54,79-80` and verifier lines 43-48,59-60 enumerate width grids, dependent retained-bit bounds, or integer parameter regimes; they do not build arrays by mutation. The two Modules at lines 29-32 genuinely allocate local symbolic `t` and intermediate mutable state; the constant phase block and assertion scopes use With.

The fixed 8-by-8 diagonal matrix is an exact three-bit special-case invariant, not an unnecessarily dense representation of a scalable state-space simulation. No sparse-data issue affects this example. Applied the wl-native-style and wl-style-guide rubrics without treating every Table or necessary Module as a defect.

Defect: none. Minimal fix: none.

## C. Terse physics reasoning — PASS

Evidence: `COMPARISON.md:44-88` gives the normalized inverse-square-root target, exact polynomial certificate, Horner row-error bound, scale invariance, and Newton relative-error map in the order needed to assess the model. `Feasibility.wl:12-17` defines the three mathematical error objects, lines 19-26 define certificate domain and proof, lines 28-68 build the quantified architectures, and lines 71-85 derive and collect rational parity phases. Names such as seedBound, hornerBound, deltaBound, SquareDropBound and NewtonRelativeError denote mathematical quantities. The Python code separates circuit generation, exact integer arithmetic, sampled circuit checks, and pricing instead of making proof claims from a resource count.

The shared B design records retain unused Newton metadata, but it does not obscure the arithmetic argument, alter certification, or justify an additional refactoring requirement for this bounded feasibility artifact. The report labels all estimates, samples, failures and selection limits explicitly (`COMPARISON.md:75-78,95-101,137-152,165-169,196-199`).

Defect: none. Minimal fix: none.

## D. Physics invariants — PASS within the stated evidence level

Evidence: `verify-feasibility.wls:40-42,66-75` asserts normalization-scale invariance, exact diagonal unitarity for arbitrary real theta, absolute zero-state phase, retained-row phase equality and modulo-one phase. These all reproduce True. Supplemental kernel checks prove the squared relative-error certificate equivalent to `Abs[p Sqrt[m]-1]<=eta`, symbolic phase norm preservation on all eight three-bit basis values, and the arithmetic Newton residual

`e_out = -3 e^2/2-e^3/2+(1+e)(m ds+dm)/2-Sqrt[m] dp/2`.

The associated exact coefficient bounds reproduce True for `|e|<=h<=1`, `1<=m<=4`, and `1<=Sqrt[m]<=2`. The source's Newton truncation allowance is conservative relative to this derivation. The one-step Horner error inequality also returned True, supporting geometric accumulation of row drops rather than whole-product rounding.

`resource_probe.py:257-320` records finite-register capacity before modular reduction; the full-domain rerun checks 3,133,452 inputs for each of six A/B variants with zero capacity violations. `resource_probe.py:366-407` directly asserts restored x, clean scratch and exact phase/model agreement on 162 quadratic-seed shots and 186 shots for each other A/B variant. Every assertion reproduces. The independent supplemental C circuit sample used 63 distinct inputs, 504 shots, and all exponents 0 through 10; output was `x_restored=true`, `scratch_clean=true`, `phase_matches_independent_official_model=true`, `phase_mismatch_inputs=[]`.

The relevant Hilbert-space contract is a diagonal phase on the preserved 22-bit x register with clean ancillas, not a density-operator or Hamiltonian calculation; trace/positivity and Hamiltonian Hermiticity are inapplicable. The reviewed feasibility comparison explicitly withholds an all-input/all-measurement-branch circuit proof, so absence of the excluded full verifier is not an unfulfilled invariant claim.

Defect: none. Minimal fix: none. Limit retained: sampled circuit evidence does not certify all HMR branches or the full official task.

## E. Nontrivial examples — PASS

The base case could have been a scalar zero-error Newton map or a two-bit phase identity. The artifacts go beyond both: four nonconstant piecewise-polynomial architectures spanning 12, 96, 48 and 192 intervals; independent quantized row truncation; six finite-register domain scans; eleven coherent normalization exponents; and a genuine failing width boundary. This structure is visible in `run-feasibility.wls:3-7`, `Feasibility.wl:54-68,78-85`, and `resource_probe.py:542-563`.

The supplied parity identity is tested symbolically for sizes 1 through 6 (`verify-feasibility.wls:59-60`); the independent review extends that polynomial identity through size 8. The 28-bit linear seed fails for exactly x=1,048,575, with reproduced error/tolerance 1.062493771317591 and minimum margin -67,102,208 units. Thus the examples expose a real precision failure rather than hiding behind zero angles, a single input, or an exact seed.

Defect: none. Minimal fix: none.

## F. Symbolic-first derivations — PASS

Evidence: `Feasibility.wl:12-26` retains arbitrary relative error and exact rational symbolic certificates; only coefficient discovery at lines 33-38 uses 50-digit minimax. The resulting coefficients are integers over a power of two, eta is rational over 2^40, and all 348 intervals are subsequently checked by exact Resolve. Supplemental direct quantifier elimination reconfirmed every positivity, squared-error, derivative and Horner-capacity inequality.

Newton's arbitrary-e identity, zero-error special case, limit and scale invariance are exact (`verify-feasibility.wls:35-42`), with numerics used after the derivation. Parity supports and turns in `Feasibility.wl:71-85` use rational arithmetic and reduce modulo one, keeping GF(2) Boolean coefficient arithmetic separate from arbitrary-angle phase arithmetic. Fixed widths and exponent ranges represent this specific benchmark instance; they are not premature substitution of a symbolic physical parameter. Numeric minimax is a discovery tool whose output is certified independently; no minimax/global synthesis optimum is inferred from it.

Defect: none. Minimal fix: none.

## G. Regimes and limits — PASS

The required regimes table exists at `COMPARISON.md:176-184` and is supported as follows:

| Regime | Executed/source evidence |
|---|---|
| General symbolic | Arbitrary e and normalization scale (`verify-feasibility.wls:35,38-42`); arbitrary real theta (`66-69`); polynomial identity sizes 1..6 (`59-60`), independently extended to 1..8. Supplemental rational certificate equivalence is symbolic in p,m,eta. |
| Exact special | Zero seed error (`36`), exact bit products (`57-60`), exact diagonal unitarity and zero-state phase (`66-70`). All True. |
| Limit | `Limit[NewtonRelativeError[e]/e^2,e->0]===-3/2` (`37`), independently reproduced True. |
| Numerical reference | High-precision polynomial samples (`Feasibility.wl:47`); exhaustive bit-exact scans (`resource_probe.py:323-363`), compiler samples (`366-407`), and independently sampled C against the separate official model. All regenerated evidence matches stored results. |
| Failure/edge | Inaccurate, inexact and empty certificate inputs (`verify-feasibility.wls:54-56`); independent invalid-a rejection; row-vs-whole-product witness (`49`); exactly one failing linear width input; normalization boundaries, modulo-one phases and absolute zero phase (`70,75`). |

Boundary assignment and both polynomial endpoints are documented (`COMPARISON.md:111-115`); the exact certificates include endpoints. No row is waived because it is difficult or deferred.

Defect: none. Minimal fix: none.

## Reproduced quantitative evidence

| Probe | Q | Toffolis | Fine rotations | Model minimum absolute margin (2^-64 turns) | Predicted Q*C |
|---|---:|---:|---:|---:|---:|
| A quadratic, NR36 | 300 | 10497 | 44 | 620330513 | 3237996.82728617 |
| A quadratic, NR28 | 260 | 8361 | 36 | 7078021 | 2244985.109086959 |
| A linear, NR32 | 256 | 8256 | 40 | 598402767 | 2182381.513220222 |
| A linear, NR28 | 236 | 7180 | 36 | -67102208 | ineligible |
| B cubic | 256 | 6638 | 43 | 537782754 | 1773730.398377167 |
| B quadratic | 220 | 4449 | 43 | 413820069 | 1043220.210631028 |
| C all exponents | 217 | 5099 | 8875 | hypothetical half-tolerance | 16399035.52753608 |

The two B bounds independently regenerate as `32285/1099511627776` and `40489/1099511627776`; both are below 2^-24. Stored incumbent/baseline verification reports agree with the supplied immutable costs and 200,540,928-shot counts. Those prior full runs were not rerun; all new results remain feasibility/model/sample evidence. No new full-circuit PASS, runtime speedup, global synthesis optimum, or final architecture selection is established.


Publication note: local paths have been replaced with portable references.
The findings and quantitative review evidence are unchanged; supplemental
temporary reviewer scripts are not included in the public folder.

OPEN ISSUES: 0
