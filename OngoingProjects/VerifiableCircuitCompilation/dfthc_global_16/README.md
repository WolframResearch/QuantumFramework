# Global DFTHC augmented-core refit

This is a mathematical (R,B,C)=(16,32,20) refit of the full released shifted
Hamiltonian. Accuracy and normalization are measured by the official decoded-
word routines. All matrices and coefficients in this representation are real.
Whole-task resource counts are analytical estimates. The prior study and
verified incumbents remain historical inputs. No quantum circuits,
circuit simulations, commits or pushes are part of this study. Raw factors,
angle words and reference arrays stay in ignored local/.

## Symbolic gradient repair

GramCores.wl is the corrected version of the earlier study's mathematical
module. Its gradient derivation uses independent symbolic U of size 2 by 3,
U0 of size 3 by 3, W of size 2 by 2 and W0 of size 3 by 2. No weight entries
are frozen and the diagonal cross-weight terms are nonzero symbolic
polynomials. Direct differentiation agrees with the general formula. The
previous study is an unchanged historical snapshot; this corrected module
is the source used here.

## Offset-preserving mathematical representation

Let T be the full 2926-by-2926 symmetric-pair target, with orthonormal svec
packing (off-diagonals multiplied by sqrt(2)). Let Z=sum_rc v_rc X_rc and
tau=sum_rc v_rc^2 be the incumbent one-body offset matrix and offset norm.
The native v is stored H2_I minus the sum of the lowercase weights.

On these data Z is negative definite. Choose L with L^T L=-Z, and a tall
512-by-76 frame F satisfying F^T F=I. Put Y=F L. For each row y_i define

    u_i = y_i/||y_i||,       a_i = -||y_i||^2.

Partition the rows into 16 groups of 32. For positive alpha with
sum_r alpha_r^2=tau, reserve one of the 20 copies per group:

    w_(r,b,0) = a_(r,b)/alpha_r,    v_(r,0)=alpha_r,
    v_(r,c)=0 for c>0.

The other 19 weight columns are free. Then

    sum_rc v_rc X_rc = sum_i a_i u_i u_i^T
                    = -Y^T Y = -L^T F^T F L = Z,
    sum_rc v_rc^2 = tau.

The one-body operator and SOS gap are therefore fixed before quantization,
while all 512 projector directions and all free covariance cores can move.
Each augmented core K_r=[w_r;v_r^T][w_r;v_r^T]^T remains PSD, with rank at
most 20. The offset-carrying factor in each block is negative semidefinite;
this is a restricted family of global refactorizations. A failed bounded fit
does not exclude every (16,32,20) representation.

The route requires Z negative definite, tau>0, at least 76 frame rows, and
nonzero y_i. Positive or indefinite Z is outside this parameterization.
The released offset spectrum and all domain checks are recorded. Allocation
uses alpha_r=sqrt(tau softmax(gamma)_r), with gamma in [-4,4].

## Global objective and exact gradients

Write P_r for the matrix of svec(u_b u_b^T) columns and X for the concatenated
P_r w_r columns. The optimizer minimizes the complete physical residual

    f = ||X X^T - T||_F^2,       d f/d X = 4 (X X^T-T) X.

Every block interacts with every other block in this objective. It is not a
sum of separate factor approximation errors. GlobalRefit.wl differentiates
the full multi-block polynomial with all weights and directions symbolic,
checking the weight and projector chain rules including cross-block terms.

Use a full-column-rank raw matrix A and its thin polar decomposition A=F H
to parameterize F. For an upstream derivative G, set

    H K + K H = F^T G - G^T F,
    d f/d A = (I-F F^T) G H^-1 + F K.

H is positive definite and K is skew. The Sylvester trace identity proves
the adjoint independently of dimension; the diagonal-H eigenbasis check is
fully symbolic. Python evaluates the formula in the singular-vector basis.
The chain through Y, its row norms, the constrained first weight column and
the alpha allocation is included. Independent finite differences test the
complete implemented map, not only isolated matrix formulas.

For fixed projector directions, let G_r=P_r^T P_r and
Q_r=P_r G_r^-1/2. For the global residual R after removing the other blocks
and the fixed offset factors,

    ||Q_r M Q_r^T-R||_F^2
      = ||M-Q_r^T R Q_r||_F^2 + ||R||_F^2-||Q_r^T R Q_r||_F^2.

Thus the best PSD rank-19 core retains the 19 largest positive eigenvalues
of Q_r^T R Q_r. WL proves the orthogonal-span identity and an exact mixed-
sign rank boundary. Four coordinate sweeps initialize each fit; the
subsequent polar-frame/weight search fits all variables jointly.

## Bounded fits, quantization and normalization

Two deterministic initializations use released directions selected across
the 65 original blocks, and eigenvectors of the leading global pair modes.
Each receives at most 1200 L-BFGS iterations. PLAN.md fixes all budgets;
iteration limits are reported without global-optimality claims.

The lower-error floating fit is encoded at 10,12,14,16 bits. For each decoded
dictionary, free weights and offset allocations receive at most 300 global
refit iterations. The offset coefficients a_i stay fixed. The offset norm
stays tau, but the exact floating Z equality need not survive angle rounding.
Measured drift, newly computed one-body eigenvectors, both official errors,
SOS gap and normalization are reported for the actual words. Already encoded
two-body words are reused exactly; vector signs do not define a different
rank-one projector.

LCU normalization depends on factors, not only their Gram cores. For the
lowest-error decoded candidate, one bounded orthogonal gauge fit rotates the
19 zero-offset columns inside each block. This preserves the represented
pair tensor, offset contraction and offset norm. The piecewise LCU gradient
is derived symbolically with arbitrary sign-chamber parameters and tested by
finite differences. The nonsmooth optimizer is not claimed globally optimal.

## Whole-task resource model

For A=76+16*32=588 and s=ceil(log2 A), a reference ordinary-unary loader has
ANGLE width 2s-1+J beta and forward/reverse Toffoli count
2 ceil(75/J)(A-2). Every legal J is priced. Include the official outer and
inner preparation plus Givens resources; add persistent workspaces and
share temporary space with max. Cost=lambda_eff * total Toffolis * total
qubits, using the candidate's measured normalization. verify-data.wls
independently checks these complete costs against RankCost.wl. All retained
loader rows obey ANGLE Q128. These are resource estimates, not compiled scores.

The unchanged three-task context comes from six SHA-bound recorded PASS
reports. Its geometric-mean arithmetic is independently recomputed in WL;
those historical circuit runs are not repeated.

## Required regimes and reproduction

| Regime | Evidence |
|---|---|
| General symbolic | All gradient coefficients, global cross-block residual, frame identity, polar/Sylvester adjoint, allocation and LCU chamber gradients |
| Exact special | Noncommuting projectors, mixed-sign PSD rank truncation, same-core/different-normalization example |
| Limit | Identity copy rotation, vanishing projection residual, zero covariance scale, polar radial derivative zero |
| Numerical reference | Independent full-map finite differences, actual 76-orbital fits, official decoded-word metrics and complete cost model |
| Failure/boundary | Rank-deficient polar input, unsupported offset sign, zero projector rows, singular dictionaries, iteration limits and actual accuracy rejection |

Relevant invariants are real Hermiticity, unit projector rows, PSD augmented
cores, tensor reconstruction, global offset contraction and offset norm.
Circuit evolution and cleanup are outside this mathematical study.

From this directory, with the existing Python environment and Wolfram 15:

    python refit.py prepare --data DATA
    WolframKernel -noinit -noprompt -script symbolic.wls
    python -m unittest test_refit -v
    python refit.py run --data DATA
    WolframKernel -noinit -noprompt -script verify-data.wls

Primary WL references consulted: [D](https://reference.wolfram.com/language/ref/D.html),
[Tr](https://reference.wolfram.com/language/ref/Tr.html),
[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html),
[GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html),
[SingularValueDecomposition](https://reference.wolfram.com/language/ref/SingularValueDecomposition.html),
and [LeastSquares](https://reference.wolfram.com/language/ref/LeastSquares.html).
The supplied package definitions and live kernel establish the actual
signatures and return shapes used in the study.
