# Augmented Gram cores and smaller projector dictionaries

Deliverable: mathematical feasibility and conditional full-resource estimates.
No circuits, push, commit or changes to the verified incumbents.

Graph: pinned-input audit -> exact WL identities/Schur conditions/gradients ->
all-data core ranks and reconstruction checks -> bounded dictionary/copy
compression pilots -> official decoded-word metric gates -> full-task resource
estimates -> fresh wl-verify -> fresh wl-quality -> bound report.

The exact representation is K = [w;v^T][w;v^T]^T per outer block. Its blocks
control two-body coefficients, the one-body offset contraction and offset
norm. The LCU one-norm also depends on the chosen factor of K. It is never
inferred from K alone.

Bounded pilots: B = 59,56,48,32 at R=65,C=30; physical-Gram copy truncation
C = 29,24,16 at R=65,B=60. B=48,32 receive at most 50 monotone variable-
projection steps, with weights eliminated by the projector Gram normal
equations. Test exact quantized vectors and recompute one-body rows, offsets,
gap and both accuracy gates. No claim of global optimization.

The first copy-truncation pilot misses the two-body gate narrowly. One bounded
extension removes a direction orthogonal to the native identity offsets in
each copy space, preserving the one-body contraction and offset norm. Optimize
the coupled full-target rank-one-removal objective for at most 200 L-BFGS
iterations. This changes only the mathematical factorization; no circuits.
After those fixed direction steps, one convex nonnegative least-squares solve
rescales the retained Schur-complement part of each core. The offset column
and norm stay fixed. No additional direction-optimization rounds are allowed.

Relaxations: unrestricted PSD pair-rank bounds; fixed-dictionary span bounds
with block/rank restrictions relaxed. A passing relaxation is only NOT RULED
OUT. Local block approximation errors are not summed into a false global
impossibility theorem. Additional rank/cost shape comparisons are conditional,
not fitted candidates.

Price the ordinary-unary loader reference plus the official full preparation,
Givens and shared-qubit model, enforce the ANGLE Q128 cap, and state the
assumed normalization. No compiled-circuit score will be claimed.

Correctness and quality have at most three fresh rounds each. Executable
quality fixes return to correctness, with at most three correctness visits.
No notebook is required. Raw arrays/candidates stay in ignored local/; only
aggregate evidence may be shareable.
