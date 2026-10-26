# DFTHC augmented Gram cores: mathematical feasibility

This study changes the mathematical factor representation and measures its
feasibility. It builds no quantum circuits, modifies no incumbent, and makes
no push or new competition-score claim. Raw arrays and numerical candidates
are in ignored `local/`. Only aggregate reports are shareable.

## Core representation

For one outer block let the columns of P be svec(u_b u_b^T), with real unit
vectors u_b and the orthonormal symmetric-matrix packing. Let w have shape
B by C and let v be the C native identity offsets. The native offsets are
`H2_I - Sum_b w_b`, not the stored `H2_I` array itself. Define

    A = [w; v^T],       K = A A^T = [M, q; q^T, tau].
    M = w w^T,         q = w v,          tau = v^T v.

The shifted two-body pair matrix is Sum_r P_r M_r P_r^T. The one-body
correction uses Sum_r P_r q_r. The offset contribution to the SOS energy
uses Sum_r tau_r. Holding these objects fixed preserves the one-body
operator and its spectrum, and hence the SOS gap. K is positive semidefinite
and its rank is the minimum number of columns needed to factor that exact K.

The LCU normalization is **not** determined by K. Its two-body term is

    Lambda_2 = (1/4) Sum_c (Sum_b |w_bc| + |v_c + Sum_b w_bc|)^2.

Every candidate therefore retains explicit factors and is evaluated with the
official metric routine. An exact two-projector example in GramCores.wl has
the same K under an orthogonal copy rotation but Lambda_2 changes from 1/2
to 1. Reducing factor redundancy is not by itself a cost improvement.

For singular M, PSD requires q in Range[M] and tau >= q^T M^+ q. Omitting
the range condition admits impossible cores. WL Reduce proves the singular
fixture M=diag(1,0): q_2=0 and tau>=q_1^2. A separate nonsingular Schur
epigraph has exact optimum 7/5 and is checked by SemidefiniteOptimization.

## Projector dictionary and eliminated weights

The projector Gram matrix is G_bd=(u_b.u_d)^2. For fixed proposed directions
U and target factor matrices P_0 w_0, the least-squares weights solve

    G w = G_cross w_0.

This removes the linear weights from the nonlinear dictionary search. The
factor residual and its exact Euclidean gradient use only Gram contractions:

    L = Tr(w^T G w) - 2 Tr(w^T G_cross w_0) + constant,
    dL/dU = 4 [((w w^T) * (U U^T)) U
               - ((w w_0^T) * (U U_0^T)) U_0].

Here * is elementwise multiplication of equal-shape matrices. The envelope
theorem permits this gradient after solving the weight normal equations.
The numerical pilot takes tangent steps on unit-vector rows and refits the
weights after quantizing the final directions. It is a local optimization
of factor-matrix residual, not a global optimizer of the two-body error or
of the full competition objective. Both official reconstruction errors,
the SOS gap and normalization are recomputed afterwards, including new
one-body eigenvectors and their quantized angle words.

There is an exact obstruction to simply shrinking the unchanged blocks.
If U has B independent rows and some copy has every w_b nonzero, then
X=U^T Diag(w) U has rank B. A left inverse of U^T exposes Diag(w), giving
both rank bounds. Any equal block Gram tensor has the same factor span;
all factors from B' shared projectors have a common orbital support of
dimension at most B'. The factor span contains this rank-B X, so B'>=B.
This applies to equality of each individual block, not to a refactorization
that redistributes contributions between outer blocks or permits error.

The all-data modular rank certificates interpret the frozen Float64 values
as exact binary rationals via SetPrecision[...,Infinity]. They do not claim
exact ranks for unrounded trigonometric vectors. The numerical singular
values and reconstruction residuals are also recorded.

## Offset-preserving approximate core compression

An approximate C->C-1 reduction removes a unit copy-space direction o with
o.v=0. Put Q=I-o o^T. Then

    M' = w Q w^T,       q' = q,       tau' = tau.

The removed pair-space factor is y=P w o. If E is the current global
two-body residual and rows of Y are the removed y vectors, then

    ||E-Y^T Y||_F^2 = ||E||_F^2 + ||Y Y^T||_F^2 - 2 Tr(Y E Y^T).

WL derives the coupled quartic and its gradient. The direction pilot is
bounded at 200 L-BFGS iterations on offset-orthogonal spheres; reaching the
iteration cap is reported and never called a global optimum.

For tau>0 write the core as M=q q^T/tau + S, S>=0. After removal, keep
q,tau fixed and rescale each retained S by t_r>=0. This preserves PSD, the
copy-rank bound, the one-body correction and the SOS energy. The represented
two-body tensor is affine in t, so its global least-squares fit is a convex
nonnegative least-squares problem. WL checks the core identity, the unchanged
offsets and norm, and the PSD Hessian. For tau=0, PSD forces q=0; the division
form is inapplicable. The released data have nonzero offsets throughout.

This candidate keeps every existing angle word unchanged. Whole-task pricing
can consequently use the existing verified ANGLE family's static gate counts
and width with the new official preparation costs. This remains an estimate
until that factorization is fully verified; no new circuit run is claimed.

## Feasibility bounds and resource estimates

The unrestricted PSD rank tail of the fixed shifted pair matrix is a necessary
bound. Passing it is only NOT RULED OUT: it omits shared rank-one-projector,
block, quantization, normalization and one-body restrictions.

For a fixed proposed dictionary, relax all block restrictions and project the
entire target into the span of its projectors. Retain the best positive
eigenmodes up to the allowed total pair rank. Its residual is a stronger
necessary bound for those fixed directions. It does not exclude a different
optimized dictionary. These full-size numerical screens use an explicit
1e-8 guard and diagnostics; they are not rigorous interval certificates.
Per-block approximation errors are never summed into a global rejection.

The ordinary-unary reference loader is priced analytically, without building
it, for all legal batch sizes. With A=76+RB and s=ceil(log2 A), its angle
width is 2s-1+J beta and its forward/reverse Toffoli count is
2 ceil(75/J)(A-2). The official preparation and Givens costs are included.
Temporary space is shared with max, as required by the official model. All
retained rows meet the ANGLE Q128 cap. The WL implementation in the earlier
RankCost.wl independently checks every estimate against the official Python
cost routine. Shape-only rows hold lambda_eff at the C30 value explicitly;
candidate rows use their measured decoded-word normalization. Rejected
candidates have no valid competition score, whatever their resource estimate.

The unchanged three-task context is derived from six recorded PASS reports:
baseline and current costs for each task. Their SHA-256 hashes and component
ratios are stored in results.json. The combination is the geometric mean of
the three baseline/current ratios; the WL data check independently reproduces
it from those files. Those historical circuit runs are not repeated here.

## Regimes and reproduction

| Regime | Evidence |
|---|---|
| General symbolic | Core identities, singular Schur condition, projection residual, projector gradients and coupled-removal/Schur-rescaling identities |
| Exact special | Noncommuting real projectors; same-core/different-normalization example; exact SDP optimum 7/5 |
| Limit | Identity copy rotation; vanishing Schur projection residual; zero scale on a retained covariance component |
| Numerical reference | All 65 released blocks, all official decoded angle words, independent derivative tests, full shifted target and complete cost rows |
| Failure / boundary | Singular dictionary, invalid shapes, out-of-range Schur vector, zero removal direction, nonnegative fit boundary, rejected metric and rank screens |

The relevant invariants are real symmetry/Hermiticity of factor matrices,
unit projector vectors, PSD Gram cores, the reconstructed coefficient tensor,
one-body correction and offset energy. Quantum state evolution, circuit
unitarity, workspace cleanup and measurement branches are outside this study.

Use the existing project Python environment and a Wolfram 15 kernel. DATA is
the user-supplied official competition bundle. From this directory:

    python study.py prepare --data DATA
    WolframKernel -noinit -noprompt -script run.wls
    python -m unittest test_study -v
    python study.py run --data DATA
    WolframKernel -noinit -noprompt -script verify-data.wls

Official documentation consulted for [MatrixRank](https://reference.wolfram.com/language/ref/MatrixRank.html),
[SetPrecision](https://reference.wolfram.com/language/ref/SetPrecision.html),
[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html),
[CharacteristicPolynomial](https://reference.wolfram.com/language/ref/CharacteristicPolynomial.html),
and [SemidefiniteOptimization](https://reference.wolfram.com/language/ref/SemidefiniteOptimization.html).
CharacteristicPolynomial uses Det[m-t I]; that sign convention is included
in the rank-drop check. Exact binary rationalization is deliberately used
instead of precision-equivalent continued-fraction rationalization.

Results and independent review outcomes are recorded separately. This source
does not imply that any candidate passes the full circuit verifier.
