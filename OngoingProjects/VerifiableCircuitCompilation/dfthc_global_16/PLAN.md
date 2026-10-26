# Global (16,32,20) mathematical refit

New authorized study, with fresh correctness and quality invocations. Previous
study directories and verified incumbents are historical inputs. No circuits,
full circuit verifier, commit, branch, or push. Raw arrays stay in ignored local/.

Graph: pin data -> fix the fully symbolic rectangular gradient -> derive global
Gram gradients and the offset-preserving tight-frame/polar pullback in WL ->
independent finite differences and exact fixtures -> two bounded global fits ->
quantized official metrics and complete resource pricing -> fresh wl-verify ->
fresh wl-quality -> evidence/privacy manifest.

Representation: the incumbent global offset Z is negative definite. For
L^T L=-Z and a 512-by-76 frame F with F^T F=I, put S=F L,
u_i=s_i/||s_i||. Assign the first copy in each block weights
w_{r,b,0}=-||s_{r,b}||^2/alpha_r and offset v_{r,0}=alpha_r.
All remaining offsets vanish; the other 19 weight columns are free.
Positive alpha are parameterized with sum alpha_r^2 equal to the incumbent
offset norm. Thus sum v X=Z and sum v^2=tau before quantization.
All 16 block covariances and all projector directions fit the global target.
Offset-allocation logits are bounded to [-4,4] for numerical conditioning.

Bounded numerical budget: two deterministic seeds (released-direction polar
frame and spectral-direction polar frame), at most 1200 L-BFGS iterations per
seed, at most 20000 function evaluations each. Initialize free covariance
cores by at most four full global residual coordinate sweeps; each fixed-frame
core update is the best PSD rank-19 fit in its whitened projector span.
No claim of global optimization. Record early stops and iteration limits.

Quantization budget: beta=10,12,14,16. Take the lowest-error floating fit;
refit free weights and offset allocations against each decoded dictionary for
at most 300 L-BFGS iterations. The floating tight-frame equality need not
survive angle quantization: report its measured drift, recompute one-body
eigenvectors, gap, both official errors and normalization. No frozen-gap claim
is made for quantized candidates. At most one 150-iteration copy-gauge fit for
the best accuracy candidate may reduce LCU normalization while preserving its
two-body tensor, native offset contraction and norm.

Price every complete candidate with its actual lambda_eff and the official
preparation/Givens/shared-workspace costs plus all legal ordinary-unary loader
batches under ANGLE Q128. No candidate is a new circuit-verified score.

Independent review budget for this new invocation: at most three fresh
wl-verify rounds and three fresh wl-quality rounds, with at most three
correctness visits. Executable quality fixes return to correctness. Report
remaining issues if a cap is reached; do not silently relax the criteria.
