# Correctness review — round 3/3

## VERDICT

All stated mathematical and aggregate-data results reproduced. No artifact correctness defect found. The fresh-context requirement is qualified by accidental exposure to prior-review excerpts, recorded below as one procedural issue. This report is not an unqualified fresh-context sign-off.

Targets are bound by the SHA-256 inventory below; ordered `filename hash\n` bundle: `1934259ceb2960aa49d1bd4a3f4a3080ef9df6ed58183da61b3ee373cd534ada`.

| File | SHA-256 |
|---|---|
| `Representations.wl` | `85817e6845e6dfaf53673e9773d2b5ca168ebbb2be0c56e611e081e45ed50a19` |
| `ProjectiveDecoder.wl` | `4b15a77c2cda15f57cd30e93cf08313bc90c1201e8e172bc8cb0f434fe35fd16` |
| `run.wls` | `2d85571bb326105edf2a47c2969a78cbcafe3c69487156cde8bd2c34a6cec205` |
| `curvature_screen.py` | `821a6119382eaa0d7329688ae833444fc5c39a6196bf3ace84ec50ec6a91f0f3` |
| `affine_rounding.py` | `8bd41280d90b06e502101977b3e088b1c66af8612b8d9eff012e3de8933c35de` |
| `table_screen.py` | `60b730ca9286d2e2e6935a1308ec618ca123e5d7321d76103cb77004f7e1d051` |
| `RESULTS.md` | `9498310af24638addfd28cf2a12a5ffc0c5dc563c4ff23d36b3a6d2c6e6cc328` |
| `PLAN.md` | `4d2c0833bd60bd73d921aa54dfe2e83151447a2c7332a1ccf4e0084dfda6a745` |
| `symbolic-checks.json` | `e8cdfee7631ff37b982ced8d7d4b95e47abb3b6f9f1ce7b2726d24824cd42622` |
| `curvature-screen.json` | `9e8e79ecfbc3390a26595ac9168062139f3a28e97c86ae3bd9b981f21dd3174d` |
| `affine-rounding.json` | `baa761187fb066ee3e9b3266b41b7225847f04340d5636c9dbfac2dbd722402b` |
| `table-screen.json` | `a2ba5abfbe00d5e4b0ee8e8d1b10325eb304c9f0eafaf7ad1c0071ed2813cda5` |

## REPRODUCED

Environment: Wolfram **15.0.0**, Mac OS X ARM 64-bit (February 17, 2026); Python **3.13.1**, NumPy **2.5.4**, kickmix **0.1**. Native field modulus **4179**, primitive element **2**. Applicable ancestor `AGENTS.md` was empty. No notebook is among the targets.

The target directory and directly required generator/coefficient dependencies were copied to an isolated directory, denoted by `SCRATCH`. `DATA` denotes the supplied official verifier bundle and `WORDS` the fixed C30 `words.u16`; no incumbent was edited. Commands, with `KERNEL` and `PYTHON` denoting the supplied runtimes:

```sh
"$KERNEL" -noinit -noprompt -script "$SCRATCH/representation_search/run.wls"
"$PYTHON" "$SCRATCH/representation_search/curvature_screen.py" --data "$DATA"
"$PYTHON" "$SCRATCH/representation_search/affine_rounding.py" --data "$DATA"
"$PYTHON" "$SCRATCH/representation_search/table_screen.py" --words "$WORDS"
"$KERNEL" -noinit -noprompt -script "$SCRATCH/extra.wls"
"$PYTHON" "$SCRATCH/independent_models.py"
"$PYTHON" "$SCRATCH/native_compare.py"
```

All final commands exited zero. The packaged Wolfram runner returned **PASS**, all **36 Boolean checks true**, with two scope strings. All four regenerated JSON reports matched the target JSON objects exactly, including hashes. The derived candidate NPZ also matched byte-for-byte. Corrected kernel runs emitted no messages; no `Quiet` was used. Preliminary verifier-only probes emitted `Function::slot1` and `AssociateTo::argrx` from temporary test typos; these were corrected and rerun. Python signature introspection initially failed on the pybind11 builtin and a nonexistent `inv` method; live `mul`/`pow` docstrings resolved the actual API, and final probes use `pow(x,-1)`.

The explicit claim/example inventory and independently reproduced results are:

| Claim / regime | Result and scope |
|---|---|
| Curvature identity, uniform bound, equioscillation, affine minimax certificate, Δ→0 (`RESULTS.md:10–19`, `Representations.wl:72–95`) | Exact identity, extrema `±h²/8`, finite-difference certificate `1/2`, and zero-change limit pass. The triangle inequality gives the stated `1/8` optimum and factor-eight envelope improvement. |
| Full inverse-square-root integer model and register bounds (`RESULTS.md:21–28`, `curvature_screen.py:19–48`) | All six candidate aggregates reproduce. Naive f8: **1,590,097** failures; compensated f8: **385,356**; f9: **186**; f10/f11/f12: **0**. A second full-domain implementation, using bit-length normalization and integer division, independently checks all **3,133,452** inputs at f10: worst ratio **0.8505589342620424**, x=**1,064,928**, A range **536,871,167..1,073,741,797**, B range **525,831..4,194,000**. Initial positive coefficients and monotone subtractions prove every partial accumulator also fits. Candidate physical widths **(30,22,9)**; incumbent generator widths **(30,22,15)**. This is an integer model, with no qubit/circuit saving established. |
| Finite affine-rounding rejection (`RESULTS.md:30–36`, `affine_rounding.py:12–36`) | Exact report: repairs `(row,ΔA,ΔB)` = **(2,3,0), (7,−1,0), (158,−1,0), (183,0,1), (184,−2,0)**; row **170** remains unresolved. Independent scalar interval intersections cover its **5,461** legal inputs at all **65** B offsets; every interval is empty, minimum integer gap **1**. Only the specified fixed-C finite family is rejected. |
| Cached quotient signature/domain/ideal certificate (`Representations.wl:6–32,96–105`) | Live `PolynomialExtendedGCD` shape `{g,{s,t}}` and Bezout identity pass over GF(16); live `PolynomialReduce` shape `{quotients,remainder}` and reconstruction pass. `CachedQuotient` returns exactly `denominator,inverse,quotient`. Degree-35 h34, degree-one h0, degree-five h71 and ideal remainder checks reproduce. Repeated roots, forbidden zero root, constant/zero/nonpolynomial locators, nonpolynomial evaluator, malformed arity/variable and negative/noninteger/machine-real/symbolic exponent are rejected. A simple zero root at h0 and a zero evaluator work. |
| Native DQI convention and phases (`RESULTS.md:48–65`, `Representations.wl:38–65`) | Released scanner source gives n=70, m=73, t=35, evaluator slots **37..71** with exponents **34..0**, and odd locator slots equal **zσ′**. A separate kickmix field implementation evolves these actual slots over all **4,095** primitive-power positions. At all **35** roots, the cached h1 quotient agrees with the native numerator/denominator; all **4,096** polynomial-basis masks per root agree after explicit trace-dual conversion: **143,360** phase comparisons. The other **4,060** conditioned phases are +1. Trace pairing is bijective over all masks; packaged degree 1/2/17/35 phase cases and the incorrect-h34 counterexample also pass. Fixed-point evaluation is linear in cache coefficients; root flag times a linear bit gives the checked unit-modulus CZ sign pattern. |
| Inversionless classical recurrence/scale/Forney (`ProjectiveDecoder.wl:7–38,44–50`, `RESULTS.md:89–106`) | Packaged 35-, one-, and zero-error cases pass. An independently implemented normalized BM recurrence agrees with every normalized projective locator at every one of **71** states in **11** further cases: GF(16) errors 1..7; native-field errors 1,2,17,35 with deterministic randomized nonzero weights/distinct locations. Final locators, Forney values, degree bounds and nonzero scales all pass. Source recurrence contains no division; testing normalization/Forney does use divisions. Zero scale is outside the valid representation. |
| Actual shifted matrices/formal inverse/truncation (`ProjectiveDecoder.wl:39–43`, `RESULTS.md:111–120`) | Direct matrix action matches both branches. Determinants are **δ z^m** and **γ**. Formal rational inverse multiplies back to identity. For the displayed four-coefficient truncated map, nonzero vector **(0,0,0,1)** maps to zero: the old template's linear coefficient is lost. Thus rational-function invertibility does not prove fixed-register reversibility; tails, branch history and inverse-update cost remain explicitly unresolved. |
| Conditional static budgets (`RESULTS.md:71–83`) | Inspected native scan and GF8 backend: five 42-operation products plus 12 root-masked ANDs per position give **909,090**. Existing manifest supplies **4,503,030** static operations and **1,903** qubits. Equal `T Q²` width is **2130.128471604313**; an added 420-qubit block gives objective ratio **1.1892874201342571** even with free setup. A degree-<d product has d² field products and monic reduction at most d(d−1), hence **2d²−d**: **2,415** at d35, four expressions cost **405,720** at 42 each. All exclusions in the text are necessary and present. |
| Fixed DFTHC table/support/ranks (`table_screen.py:47–73`, `RESULTS.md:125–143`) | Pinned word hash matches. All **25** batches reproduce 994 legal rows, **983** nonlinear support, degree histogram **[1,10,45,120,210,252,206,110,35,5,0]**, output rank **108**, difference rank **81**. Direct sign-bit checks and the existing loader's independent Möbius transform reproduce completion, unchanged prefix coefficients under randomized invalid suffixes, and support preservation under an invertible output shear. For every prefix, subset index s⊆m implies s≤m; the general proof is valid. No general AND/XOR lower bound follows. |
| Provenance/scope (`RESULTS.md:66–69,92–94,145–176`, `PLAN.md`) | The [published DQI design](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html) explicitly describes derivative-inverse caching; authoritative indexed text of the [iBM paper, §2.2.2 and Algorithm 2](https://eprint.iacr.org/2022/469.pdf) confirms scalar rescaling. Direct full-PDF retrieval returned HTTP 403, but the referenced passage was available in the primary-source search extract. Existing verification metadata contains combined local estimate **2.0088882102784944**; historical incumbent verification was inspected, not rerun. Git confirms derived `local/` tables are ignored. No speed benchmark, compiled candidate, workspace cleanup, native BM adapter, reversible decoder, or new score is asserted; those are not untested success claims. |

## FINDINGS

**P2 — Procedural independence limitation; no target-code defect.** A broad dependency `rg` inadvertently traversed `representation_search/reviews/` before an explicit exclusion was added. The tool result exposed the following prior-review snippets. The review files were not subsequently opened, and these excerpts were not used as evidence for the reproductions above.

| Prior file | Exact exposed lines | Excerpt subjects |
|---|---|---|
| `reviews/correctness-2.md` | **27,29,30,31,43,77,82,83,92** | Prior quotient/API/native-mask checks; static budgets; table support/ranks; scope/provenance; h71 test status; claim-bundle and source/report hashes. |
| `reviews/correctness-1.md` | **5,7,50,53,56,70,78,83,90,108** | Prior source hashes; abstract quotient cases; budgets; historical score/ignored tables; native-slot mismatch discussion, numerator formula, corrected h1 convention, counterexample setup and repair recommendation. |

This exposure cannot be undone in this reviewer context. It limits the required pristine fresh-context independence, despite the separate fresh-kernel, external-field and full-domain reproductions. The procedural issue is the only open finding; mathematical/code issues found: **0**. No missing runtime or data dependency prevented any core check.

## RECOMMENDED FIXES

No target-code changes are recommended from this round. Retain the explicit scope qualifications and report the independence limitation with this evidence. Do not describe this round as an unqualified pristine fresh-context review or restart a fourth round to evade the three-round cap. Any later independently commissioned review should exclude all prior-review paths from discovery before reading dependency-search results.

OPEN ISSUES: 1
