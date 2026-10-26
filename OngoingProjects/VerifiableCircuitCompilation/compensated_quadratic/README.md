# Compensated quadratic: full circuit verification

The compensated representation passes the unmodified full official inverse-
square-root verifier: **200,540,928 shots**, 64 for every legal integer input.
Measured cost is **506618.35132005374**, **3.9199901% below** the previous
verified reduced quadratic's 527287.988387595. This is a modest gain.
The combined local estimate rises from **2.0088882102784944** to
**2.035845221582375**, using the existing verified DQI and C30 DFTHC results.
The other two tasks were not rerun or changed. No Kaggle submission or push
was made.

## Mathematical representation

The frozen 36-bit coefficient source is unchanged. For P(t)=A-Bt+Ct² on
0≤t≤h, h=1/64, choose a coarser C+Δ and transform the affine terms together:

    A_new = A + Δh²/8,    B_new = B + Δh.
    P_new-P = Δ(t²-ht+h²/8).

WL `Solve` derives the alternating-extrema solution for s²+u s+v as
u=-1, v=1/8. WL `Reduce` eliminates the affine parameters from the full
quantified uniform-error constraint and returns r≥1/8. Thus the curvature
perturbation bound |Δ|h²/8 is eight times tighter than changing C alone.
This is a standard minimax approximation argument, not a novelty claim.
It is the representation change that permits the selected fractional widths
(30,23,10), with physical coefficient widths (30,22,9), versus (30,22,15).
Python exact Fractions reproduce all 192 WL-derived coefficient rows.

Two nested Horner stages still round each binary partial-product row. WL
keeps the signed error e0-t e1 and accounts for rows that are exactly aligned.
The conservative continuous envelope proves 191 intervals. For interval 0,
`grid.wls` checks all 16,384 normalized binary offsets using exact integers
after clearing denominators in the squared relative-error inequalities.
Together these prove the 2^-24 relative bound on the complete normalized
20-bit grid in [1,4). This is not a continuous certificate for interval 0.
The official circuit verifier separately checks its rounded uint64 targets,
input restoration, clean scratch and the remaining rotation-synthesis margin.

## Measured resources and cost

| Quantity | Previous verified quadratic | Compensated quadratic |
|---|---:|---:|
| Physical coefficient widths | 30,22,15 | 30,22,9 |
| Peak circuit qubits | 142 | 151 |
| Toffoli count | 3460 | 3097 |
| Fine rotations | 37 | 37 |
| Toffoli equivalent, including rotation synthesis | 3713.295692870388 | 3355.0884193381044 |
| Official cost | 527287.988387595 | 506618.35132005374 |

The smaller curvature word changes the best memory schedule. Retaining it
avoids 279 Toffolis of release/reload lookup work, while using nine more peak
qubits than the previous recycling schedule. The full cost includes that
tradeoff. Six fewer coefficient bits do not mean six fewer peak qubits.
The same compensated table with recycling costs 516040.55554601084 by the
exact-model margin, so it loses to the retained schedule. The four bounded
schedule comparisons and their boundary/cleanup checks are in screening.json;
only the selected schedule has the new full official verification result.

The full-run worst error/tolerance is 0.8505589342620424 at x=1064928.
Its cost agrees with the independent exact integer model's synthesis margin.
The gain over the original inverse-square-root baseline is 4.099287418831162.
The combined score is the geometric mean of the three baseline/current cost
ratios: 4.099287418831162, 1.2581323878369532, and 1.6360612511835888.
It remains far below the reference 8.43991. This change does not solve the
larger representation problems in DQI or DFTHC.

## Invariants and regimes

| Regime | Evidence |
|---|---|
| General symbolic | Quantified minimax solution, scaled compensation identity, signed Horner error and normalization identity |
| Exactly solvable | Three alternating residual extrema; exact ties-to-even rounding |
| Limit | Δ→0 restores the original polynomial |
| Numerical/exact reference | 192 coefficient rows; continuous/grid certificates; all 3,133,452 integer inputs; all 200,540,928 official circuit shots |
| Failure/boundary | Invalid tables/widths; noncertified continuous interval remains explicit; every reachable interval/normalization boundary tested in four compiled schedules |

The circuit implements a diagonal phase on the promised input domain.
Phase modulus is one. Positive coefficients and monotone partial subtractions
fit the certified registers. The compiled boundary probes match the integer
phase exactly, restore x, clear scratch and reconcile priced gate counts.
The full verifier repeats restoration and cleanup checks for every shot.
No Hamiltonian dynamics or many-body state simulation is claimed.

## Reproduction and provenance

From this directory, with the parent project dependencies and kickmix:

    WolframKernel -noinit -noprompt -script run.wls
    WolframKernel -noinit -noprompt -script grid.wls
    python generator.py
    python screen.py --data DATA
    python verify.py --data DATA --workers 4

The builder imports the existing normalization and Horner primitives; it does
not modify those files or import competition code. It checks the frozen
source hash, independently derived table, WL source bindings and capacities.
The standalone default builds the selected retained-curvature circuit;
screening checks that its canonical binary matches the selected schedule.
The verification runner pins the official verifier, the full specification,
the candidate binary and all executable inputs before and after the run.
All original comparison reports stay intact. Binary output is ignored local/.

Documented APIs checked in WL 15 and against official references:
[Solve](https://reference.wolfram.com/language/ref/Solve.html),
[Reduce](https://reference.wolfram.com/language/ref/Reduce.html),
[Resolve](https://reference.wolfram.com/language/ref/Resolve.html).
The source package contains the exact executable argument. Reviews and their
hash bindings are recorded separately; this README does not imply their status.
