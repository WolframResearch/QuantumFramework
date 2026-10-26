# DFTHC representation research plan

Saved October 26, 2026. Research is paused; these are the next proposed mathematical experiments.

The goal is to discover better mathematical representations using WL's symbolic tools. The [global refit](dfthc_global_16/RESULTS.md) produced no accurate candidate: its best two-body error is 280.233521 mHa against the 20 mHa limit. Its offset parameterization imposes an extra restriction, so the failed bounded searches leave other representations open. The verified C30 implementation remains the incumbent.

## Diagnose the source of the error

Compute relaxed bounds for the learned projector dictionary. Compare unrestricted PSD cores with the current block, rank and offset restrictions to distinguish insufficient projector directions from restrictive offsets and incomplete local optimization. Label numerical screens separately from certified bounds. Use this diagnosis to choose the next representation before another large nonlinear fit.

## Derive general augmented Gram cores in WL

Use WL to derive the Schur-complement constraints, including singular and zero-offset cases, and eliminate dependent variables symbolically. Permit offset arrangements excluded by the tight-frame construction. Check the general symbolic case, exact examples, limits, independent numerical references and failure cases before fitting the full data.

## Compress gradually from the verified C30 solution

Start from the verified factors and reduce blocks and shared projectors gradually, retaining the best accurate intermediate solutions. Compare (20,48,24) with (16,32,20); the smaller shape is a candidate rather than a fixed choice. Track both official reconstruction errors, the SOS gap, effective normalization and the complete task cost at each stage.

## Select and verify a feasible improvement

Require floating accuracy comfortably below 20 mHa before testing angle precision. Evaluate the actual decoded words, recompute normalization and the one-body rows, and price the whole task with legal loader widths and shared workspace. Select only a candidate that passes both accuracy gates and improves the incumbent cost.

Run fresh wl-verify and wl-quality on any new mathematical work. Circuit implementation and full circuit verification follow only after a feasible mathematical candidate is selected and that further work is authorized. Preserve the existing verified implementations and keep raw data, factors, words and logs outside publication.
