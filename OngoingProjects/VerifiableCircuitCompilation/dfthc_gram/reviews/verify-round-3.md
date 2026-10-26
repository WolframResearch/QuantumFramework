# ROUND: 3/3 — independent WL verification

## VERDICT

PASS within the stated mathematical-feasibility scope. The complete README reproduction succeeded in an isolated temporary replica containing only the named artifacts and dependencies. The original target files were preserved. No circuits, historical circuit verification, new competition score, optimization-budget extension, incumbent changes, commits, or pushes were performed.

Paths below are relative to `OngoingProjects/VerifiableCircuitCompilation/`.

## REPRODUCED

**Claim inventory and coverage.** Read the complete twelve target artifacts, the required `dfthc_feasibility/RankCost.wl`, pinned input fields, recorded build/resources/cost reports, and the official metric/preparation implementations needed by the study. Checked applicable exact `AGENTS.md` locations; no additional instructions applied. Claims covered: augmented-core identities and rank obstructions; singular/nonsingular Schur conditions; normalization noninvariance; projector and coupled-removal gradients; offset-preserving compression and convex rescaling; frozen-array rank/reconstruction evidence; eight official decoded-word metric pilots; numerical PSD/fixed-dictionary bounds; conditional whole-task resources; six-report score-context arithmetic; and the stated negative feasibility conclusion. There are no notebooks or runtime benchmarks in this target.

**Commands and environment.** With `PY`, `WL`, and `DATA` denoting the supplied Python environment, Wolfram kernel, and official bundle, ran the following commands from the replica of `dfthc_gram/`, in README order:

```text
$PY study.py prepare --data $DATA
$WL -noinit -noprompt -script run.wls
$PY -m unittest test_study -v
$PY study.py run --data $DATA
$WL -noinit -noprompt -script verify-data.wls
```

All five commands exited 0. Versions: Python 3.13.1; NumPy 2.5.4; SciPy 1.18.1; Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. The production reproduction emitted no warnings or kernel messages. No `Quiet` or message suppression was used. WL diagnostics used `Echo`.

- **Exact mathematics and all-data evidence:** `run.wls` returned `PASS`; all 37 symbolic checks passed, all 65 orbital-frame modular ranks were 60, and all 65 augmented-array modular ranks were 30. The prime was independently checked. Live `MatrixRank` options include `Modulus`; exact binary rationalization was checked directly. All full-support-copy and reconstruction checks passed. Across the stored blocks, minimum singular values were approximately 0.000526403 and 0.0000541554; the largest recorded reconstruction residual was 9.60e-15. These certify only the frozen binary-rational arrays, as the reports explicitly state.
- **Schur fixture:** the independently evaluated quadratic form with the inverse leading block is exactly 7/5. `Minimize` and `SemidefiniteOptimization` reproduced 7/5 and 1.399999996255239 respectively. The singular range condition, invalid out-of-range vector, zero-offset boundary, and rank-drop characteristic-polynomial convention were checked.
- **Derivatives and boundaries:** all five packaged Python tests passed. Additional fresh-kernel tests checked an independently differentiated rectangular factor residual against `ProjectorGradient`, exact rectangular core/Gram outputs, zero cores, empty/incompatible inputs, and nonfinite/complex/symbolic rejections. The direct dependency's boundary suite also passed. Additional independent Python directional differences over five random fixtures gave maximum absolute derivative discrepancies of 1.42e-8 for variable projection and 7.18e-8 for coupled removal. Singular/zero dictionaries, duplicate dictionary spans, and zero removal directions were rejected. Small exact-spectrum fixtures independently checked the fixed-dictionary bound with rank 0, positive-rank truncation, a negative eigenvalue, and target support outside the dictionary.
- **Numerical pilots:** all eight official decoded-word metric rows reproduced exactly, with zero metric passes. Two-body errors in mHa were 72.098892, 214.045133, 527.376327, 1604.340611, 20.337383, 24.583187, 62.355836, and 20.271424 in the report's order. B48/B32 retained their 50-step monotone histories. Offset-preserving C29 used exactly 200 direction iterations, reported optimizer failure at its cap, and preserved all angle words. Its Schur KKT residual was 3.46945e-17; coefficient and offset-norm drifts were 8.88178e-16 and 2.84217e-14. The reported predicted physical residual agreed with the official metric. All offset norms were nonzero; the minimum was 1.108958.
- **Bounds:** reproduced every shape row and the fixed quantized B32 dictionary floor 0.13469057200509507 Ha, with orthogonality residual 1.37768e-13 and the explicit 1e-8 guard. The (8,32,32) unrestricted PSD floor was 0.03376848165684026 Ha, correctly rejected. The other displayed smaller shapes remain only `NOT_RULED_OUT`. Source inspection and independent small fixtures support the stated relaxation scope; these are numerical screens, not interval certificates or bounds on reoptimized dictionaries.
- **Resources and recorded context:** `verify-data.wls` returned `PASS` for all ten checks, covering 21 ordinary-unary resource rows, their optimal legal batches and Q128 caps, and one unchanged-loader row. Static loader counts and the official preparation model reproduce C29's 65,088 Toffolis, 661 qubits, and conditional cost 2,054,139,541.1018138: 0.1519406% above the incumbent. The six SHA-bound reports have recorded `PASS` status; their ratios reproduce 4.099287418831162, 1.2581323878369532, and 1.6360612511835888, with geometric mean 2.035845221582375. No new competition score is asserted.
- **Artifact consistency:** complete parsed `inputs.json`, `symbolic.json`, `results.json`, and `wl-data-checks.json` from the replica equal their target counterparts. Source hashes and all eight stored candidate hashes match. No comparison depended solely on ZIP hashes.

**Evidence limits.** Historical full-circuit success remains recorded evidence: only its hashes, statuses, costs, static counts, and arithmetic were checked, as requested. Circuit execution, unrounded-trigonometric exact ranks, rigorous spectral intervals, and global optimization are explicitly outside the study and were not treated as reproduced claims. No required dependency, kernel, or approval barrier reduced the requested coverage.

## FINDINGS

None. The conclusion that these eight bounded pilots found no feasible improvement is supported; it does not establish impossibility for other factorizations or dictionaries. Conditional resource estimates and unchanged historical context are labeled appropriately.

## RECOMMENDED FIXES

None required.

OPEN ISSUES: 0
