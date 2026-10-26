# Augmented Gram and projector-dictionary feasibility results

**No feasible improvement was found in the eight bounded pilots.** No circuits were implemented or verified, and nothing was committed or pushed. The verified DFTHC incumbent remains cost **2,051,023,203.37**, with the three-task local estimate **2.0358452216**. This study produces no new competition score.

## What WL established

- All 37 exact symbolic checks pass: augmented-core reconstruction, singular Schur range condition, copy-rotation invariants and noninvariant LCU normalization, projector gradients, offset-preserving rank reduction, coupled global residual, and convex Schur rescaling.
- In all 65 frozen blocks, the augmented coefficient array has rank 30 and the orbital frame has rank 60 over a prime field, interpreting the stored Float64 values as exact binary rationals. Every block has a factor with all 60 nonzero weights. The unchanged exact block therefore needs at least 30 copies and 60 shared projectors. This does not rule out approximate compression or redistribution across blocks.
- The singular Schur condition needs both q in Range[M] and tau >= q^T M^+ q. The exact convex fixture has optimum 7/5; the independent numeric SDP gives 1.39999999626.
- An orthogonal copy rotation can preserve the whole augmented Gram core while changing LCU normalization. Rank reduction alone cannot predict the competition objective.

## Actual decoded-word metric results

All candidates use 10-bit angle words and the official metric routines. The two-body error gate is 20 mHa; the incumbent error is 19.980680 mHa. Smaller dictionary candidates also recompute the one-body eigenvectors and angle words. No candidate passes all metric gates.

| Pilot | Two-body error (mHa) | Effective normalization | Metric gates |
|---|---:|---:|---|
| dictionary-B59 | 72.098892 | 59.382686 | FAIL |
| dictionary-B56 | 214.045133 | 95.480112 | FAIL |
| dictionary-B48 | 527.376327 | 163.371610 | FAIL |
| dictionary-B32 | 1604.340611 | 304.678531 | FAIL |
| copies-C29 | 20.337383 | 51.601835 | FAIL |
| copies-C24 | 24.583187 | 51.566208 | FAIL |
| copies-C16 | 62.355836 | 53.242628 | FAIL |
| offset-preserving-C29 | 20.271424 | 47.744962 | FAIL |

The B=48 and B=32 searches each take 50 monotone variable-projection steps; they minimize factor residual locally, not the global physical error. B=59 and B=56 use a greedy projector selection and least-squares refit. These failures do not prove that every dictionary of those sizes fails.

The offset-preserving C29 candidate uses 200 coupled direction steps followed by one convex nonnegative Schur-rescaling fit. The direction optimizer hits its iteration cap; no global optimality is claimed. The scalar fit satisfies its recorded KKT residual below 4e-17. It reduces the error from 20.337570 to 20.271424 mHa, still **0.271424 mHa above the gate**. Its one-body contraction and offset norm drift by less than 3e-14, and every incumbent angle word remains unchanged.

For the specific quantized B32 directions, even a global PSD fit with all block restrictions relaxed has error at least 134.690572 mHa (numerical screen, 1e-8 Ha guard). That fixed dictionary is rejected. Other reoptimized dictionaries remain open.

## Whole-task cost estimates

Costs include forward/reverse angle loading, both state-preparation stages, Givens rotations, and shared temporary space. They are the official objective product lambda_eff * Toffolis * qubits. Every loader estimate meets the 128-qubit ANGLE cap. These are estimates; no new circuit family is built.

| Representation | Total Toffolis | Total qubits | Effective normalization | Whole DFTHC cost | Status |
|---|---:|---:|---:|---:|---|
| Verified incumbent C30 | 65,244 | 661 | 47.558542 | 2,051,023,203 | Previously fully verified |
| Offset-preserving C29, existing loader | 65,088 | 661 | 47.744962 | 2,054,139,541 | Accuracy FAIL; estimate only |
| C30, ordinary-unary reference | 72,378 | 731 | 47.558542 | 2,516,242,481 | Analytic reference only |

C29 saves 156 Toffolis, but the normalization increase makes the modeled cost **0.152% worse** than the incumbent. Merely shrinking C does not remove the large R*B angle dictionary. The ordinary-unary reference is costlier than our existing optimized loader; it is included to make comparisons consistent across shapes.

The following shapes have **no accurate fitted candidate**. Their error floors relax the shared-projector and block restrictions; passing a floor only means NOT RULED OUT. Costs assume the incumbent normalization 47.558542, 10-bit words, and the ordinary-unary reference. A new fit could change the normalization, gap, error and optimal loader substantially.

| Shape (R,B,C) | Relaxed error floor (mHa) | Total Toffolis | Total qubits | Conditional whole cost | Rank screen |
|---|---:|---:|---:|---:|---|
| (32, 32, 16) | 1.139588 | 22,732 | 693 | 749,202,844 | NOT_RULED_OUT |
| (20, 48, 24) | 1.629534 | 21,724 | 725 | 749,042,286 | NOT_RULED_OUT |
| (16, 32, 20) | 11.743128 | 14,037 | 690 | 460,629,689 | NOT_RULED_OUT |
| (10, 32, 30) | 16.697068 | 10,126 | 700 | 337,104,460 | NOT_RULED_OUT |
| (8, 32, 32) | 33.768482 | 9,115 | 443 | 192,038,778 | REJECT |

The 0.34–0.75 billion rows identify the scale of a useful mathematical target, not achieved performance. The (8,32,32) row is ruled out by the necessary PSD rank screen despite its attractive resource count. The full size (65,32,30) shape has a weak rank floor of 0.298 mHa, whereas our fixed fitted directions have a much stronger 134.691 mHa floor; rank alone misses the central difficulty.

## Consequence for the representation search

The next substantive mathematical problem is a joint factorization of the **global** pair tensor, with augmented offset constraints, into fewer outer blocks and shared rank-one projectors. A per-block exact change of coordinates cannot shrink these full-rank blocks. A useful objective must control the global physical residual and the factor-dependent LCU one-norm together; the bounded factor-residual pilots do not do that. No such implementation is authorized or included here.

## Evidence and scope

`symbolic.json` contains the exact identities, all-data rank/reconstruction checks and SDP fixture. `results.json` contains every metric result, numerical bound, fit cap and resource row. `wl-data-checks.json` independently prices the complete rows in WL against the official Python cost model, including reuse of the static verified loader counts. Five Python tests check the Gram contractions, gradients, projection monotonicity and malformed inputs. Raw arrays and candidate files remain in ignored `local/`.

The unchanged combined context is calculated from six SHA-bound recorded PASS reports. Their baseline/current ratios are 4.099287418831162 (inverse square root), 1.2581323878369532 (DQI), and 1.6360612511835888 (DFTHC). The geometric mean is 2.035845221582375. WL independently reads the component reports and checks that arithmetic; none of those historical circuit runs is repeated here.

Fresh wl-verify and wl-quality outcomes are recorded in the separate review reports. Full circuit verification, exact unrounded-trigonometric rank claims, interval-certified spectral bounds and a globally optimal representation are outside this study.

## Independent review outcome

Fresh wl-verify passes on round 3 with zero open issues. Fresh wl-quality reaches its round-3 cap with one open symbolic-depth item; it is **not a quality pass**. Axes A-E and G pass. The saved exact projector-gradient check fixes numeric weight matrices whose cross-weight diagonal vanishes, so that check does not establish the formula for arbitrary symbolic weights. The reviewer independently differentiated a fully symbolic rectangular case and confirmed the formula, but the saved derivation should retain symbolic weight entries. This is a coverage gap, not a numerical or physical error found in the results.

The remaining source change and its fresh reviews are outside the exhausted review budget. The implementation and evidence remain at the final independently verified bytes. No circuit implementation is selected. See `reviews/verify-round-3.md` and `reviews/quality-round-3.md` for the final verdicts.
