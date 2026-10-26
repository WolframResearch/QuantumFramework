# Independent WL quality review, round 2/3

The feasibility conclusion is supported: the eight bounded pilots find no
accurate improvement. All physics, symbolic, numerical and resource checks
below pass. One small WL idiom issue remains under the explicit Listable-function
criterion. It does not change the physical conclusion.

## Scope and execution

Read every required target artifact in full: `GramCores.wl`, `run.wls`,
`verify-data.wls`, `study.py`, `test_study.py`, `README.md`, `RESULTS.md`,
`PLAN.md`, `inputs.json`, `symbolic.json`, `results.json`, and
`wl-data-checks.json`, all under
`OngoingProjects/VerifiableCircuitCompilation/dfthc_gram/`.
Also read the permitted `dfthc_feasibility/RankCost.wl` dependency and relevant
official `verify.py` and `scoreinfra/dfthc_math.py` definitions.

Checked exact applicable ancestor `AGENTS.md` paths. No applicable nonempty
repository instructions were found. Did not read other review reports or
representation-search material. Used isolated copies for all scripts producing
reports; did not modify target implementation or evidence. No circuits, full
circuit verifiers, fitting extensions, commits or pushes were run.

Executed:

- `WolframKernel -noinit -noprompt -script run.wls`: exit 0, PASS; all 37 exact
  checks, all 65 block rank/reconstruction checks, and numeric SDP fixture pass.
- `python -m unittest test_study -v`: five tests, all pass.
- `WolframKernel -noinit -noprompt -script verify-data.wls`: exit 0, PASS; all
  21 ordinary-loader cost rows, one unchanged-loader row, source hashes, recorded
  report identities and combined context pass.
- Independent read-only Python checks of the eight SHA-bound saved candidates,
  recomputing their official decoded-word metrics and physical invariants.
- Independent PSD-rank-tail recomputation and an SVD-based fixed-dictionary
  relaxation, rather than the implementation's QR route.

The kernel identifies itself as Wolfram Language 15.0.0 for Mac OS X ARM,
February 17, 2026. No messages were suppressed.

## A. Documented symbols — PASS

The six exported `DFTHCGram` symbols have definitions and usage text in
`GramCores.wl:2-7`; the three array APIs' numeric domains, shapes, return formulas
and failure branches are explicit at lines 9-21. The exact-check entry points
are defined at lines 23-117 and executed in the isolated kernel. The independent
gradient tests pass, including the direct WL/Python fixture comparison.

The only other non-core WL APIs used are `DFTHCFeasibility`'s `PreparationCost`
and `LookupCost`. Their signatures and association return fields are defined in
`OngoingProjects/VerifiableCircuitCompilation/dfthc_feasibility/RankCost.wl:4-5,26-62`.
The cost rerun confirms the assumed return fields and all 21 rows against the
recorded official Python model outputs.

Live `Information[SemidefiniteOptimization,"Usage"]` confirms the symbolic
objective/constraints/variables form and optional solution-property argument.
`Options[SemidefiniteOptimization]` returns `MaxIterations`, `Method`,
`PerformanceGoal` and `Tolerance`. The actual property request at `run.wls:39-40`
returns scalar `1.399999996255239`, matching the independently obtained exact
minimum `7/5`. Live `Options[MatrixRank]` includes `Modulus` and `Tolerance`;
the exact modular-rank calls run successfully. The binary-rational rank scope is
explicit at `README.md:70-73` and `run.wls:48`. No unsupported shape, option,
package symbol or memory-only API dependency was found.

## B. Idiomatic functional WL — FINDING (one issue)

`GramCores.wl:51` uses

```wolfram
Map[Factor, rotated.Transpose[rotated]-k, {2}]
```

although `Factor` is Listable. Installed-kernel evidence:

```wolfram
Attributes[Factor]
(* {Listable, Protected} *)
```

For the full symbolic 3-by-2 augmented array and rational orthogonal rotation
used by the artifact, the independent kernel check

```wolfram
Factor[residual] === Map[Factor, residual, {2}]
```

returns `True`. The smallest fix is `Factor[rotated.Transpose[rotated]-k]`.
There is no semantic change for this nested matrix. This is a minor idiom
finding, not a correctness or physics failure.

Otherwise the WL implementation uses direct matrix operations, `MapThread`,
`With`, `Select`, `AllTrue` and declarative exact algebra. No `Do`, `For`,
`While` or `AppendTo` occurs in the target WL. The small `Table` expressions
enumerate principal-minor cardinalities and are appropriate. The test
procedures' `Module` locals are genuinely assigned symbolic expressions.

## C. Terse physics reasoning — PASS

`README.md:10-36,40-68,77-102` develops the physical objects, eliminated weights,
rank obstruction and offset-preserving compression with explicit assumptions.
The implementation distinguishes factor residual from the global physical
two-body residual (`study.py:88-96,205-267`) and uses coupled global residuals
for copy removal. `RESULTS.md:27-31,43-59` gives the iteration cap, accuracy
failure, normalization penalty and limitations of each rejection. The scalar
normalization is kept factor-dependent. No dead executable calculation or
verbosity that obscures a substantive physics step warrants a finding.

## D. Physics invariants — PASS

Exact evidence from `GramCores.wl:47-69,76-98,102-117` covers augmented-core
reconstruction, real PSD quadratic forms, orthogonal gauge invariance, Schur
range restrictions, one-body contraction, offset norm and rank-one removal.
The all-data kernel rerun confirms the 65-block records rather than relying
only on saved PASS flags.

Independent saved-candidate checks use official `reconstruct_x`, which
normalizes the real projector rows (`scoreinfra/dfthc_math.py:76-85`). Across all
eight candidates:

- maximum unit-row squared-norm error is `7.771561172376096e-16`;
- maximum Hermiticity/symmetry error of reconstructed factors is
  `4.163336342344337e-17`;
- minimum computed augmented-core eigenvalue is
  `-1.2410167430624446e-14`, consistent with floating-point evaluation of
  the exact PSD construction `A.Transpose[A]`;
- every official metric reproduces the corresponding `results.json` value
  with maximum difference `0.0`, and all eight metric gates reject.

For the offset-preserving C29 candidate, direct contractions reproduce maximum
one-body coefficient drift `8.881784197001252e-16` and offset-norm drift
`2.842170943040401e-14`. These agree with `results.json` and the assertions at
`study.py:256-258`. Official metric code derives the one-body operator from this
contraction and SOS energy from its eigenvalues and offset norm
(`scoreinfra/dfthc_math.py:148-164,328-359`). No circuit-evolution invariant is
needed for the authorized feasibility scope.

## E. Nontrivial examples — PASS

The work could hide behind commuting projectors, a scalar core or zero offsets,
but does not. `GramCores.wl:63-64` tests both normalization change under a
same-core rotation and noncommuting rank-one real projectors. The 3-by-3
singular-Schur example has a genuine out-of-range obstruction; the nonsingular
SDP has off-diagonal coupling and optimum `7/5` (`run.wls:35-41`). The coupled
removal test has two outer blocks, three copy coordinates and five pair-space
coordinates (`test_study.py:26-42`), exposing cross-block interactions. The
full 65-block numerical instance is checked. Complexity is tied to the claimed
mechanisms rather than added independently of them.

## F. Symbolic-first — PASS

Core identities, rational rotations, projection residual, gradients,
compression and convex-Hessian identities are derived with exact symbolic
parameters (`GramCores.wl:23-117`). The exact `Minimize` optimum precedes its
numerical SDP cross-check (`run.wls:35-41`). Numerical work is appropriately
confined to frozen decoded data, bounded local pilots and explicitly guarded
spectral screens. The exact binary rational interpretation of the frozen
Float64 rank certificate is stated without extending it to unrounded
trigonometric data. No premature numerical substitution masks a symbolic
derivation claimed by the report.

## G. Regimes — PASS

The requested regimes table is present at `README.md:137-143`.

- General symbolic: the 37 exact core, Schur, gradient, compression and
  rescaling checks rerun successfully.
- Exact special: noncommutation, same-core/different-normalization and the SDP
  optimum `7/5` are executable and pass.
- Limits: identity rotation and vanishing projection residual are explicit
  at `GramCores.wl:65-66`. The rescaling identity is polynomial in its scale,
  hence includes zero scale; an independent exact fixture with unit
  `o={3/5,4/5,0}` and `v={-4/5,3/5,1}` additionally gives `True` for
  zero-scale core equality, PSD, unchanged offset contraction and unchanged
  offset norm. `README.md:96-97` correctly excludes division by zero offsets.
- Numerical reference: all eight official metric rows reproduce exactly;
  the Python derivative tests and all-data WL checks pass. Independent
  recomputation matches all 13 PSD-rank tails within `1e-10` Ha. Replacing QR
  by SVD gives fixed B32 span dimension 2080, minimum span singular value
  `0.010091641114343308`, and relaxed error `0.1346905720044907` Ha, matching
  the recorded `0.13469057200509507` Ha and lying far above the 0.020 Ha gate.
- Failure/boundary: the exact out-of-range Schur constraint, malformed API
  calls and nonnegative-fit boundary pass. The zero-removal-direction Python
  test passes; a direct singular-dictionary call raises
  `ValueError('Singular or unresolved dictionary Gram.')`. Actual metric
  rejection and the independent rank rejection are confirmed.

The report consistently treats passing relaxations as NOT RULED OUT and the
failed local pilots as evidence about these bounded attempts. No new score or
global impossibility claim is inferred from their failures.

OPEN ISSUES: 1
