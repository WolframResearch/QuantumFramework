# Compensated quadratic implementation

Exact WL minimax derivation and independent coefficient correspondence ->
compiled resource comparison -> exhaustive integer model and boundary/cleanup
checks -> unmodified full official verifier -> fresh wl-verify -> fresh
wl-quality -> hash-bound actual cost and combined-score report.

The bounded schedule comparison uses padded arithmetic with and without
coefficient recycling, direct controlled arithmetic, and one smaller adder
budget. The approved representation and precisions remain (30,23,10).
No further coefficient or architecture search belongs to this implementation.

Quality code fixes return to a fresh correctness visit. Each invocation has
at most three fresh rounds; at most three correctness visits are allowed.
Notebook generation is inapplicable: this is a code deliverable.

All incumbent sources, reports, frozen coefficients and other tasks remain
intact. Raw binaries stay in ignored local/. No commit, branch, push or Kaggle
submission is authorized in this stage.
