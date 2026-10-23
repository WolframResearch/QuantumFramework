# Reduced-width quadratic: full verification

The optimized quadratic circuit **passed all 200,540,928 shots** of the
unmodified official local verifier: 64 shots for each integer input from
1 through 3,133,452. Input restoration, scratch cleanup and phase accuracy
passed with the default problem specification and the same declared salt
as the baseline and earlier quadratic implementation.

| Quantity | Supplied baseline | Earlier quadratic | Reduced quadratic |
|---|---:|---:|---:|
| Qubits | 260 | 220 | **142** |
| Raw Toffolis | 7,735 | 4,449 | **3,460** |
| Fine rotations | 37 | 43 | **37** |
| Toffoli-equivalent C | 7,987.593207 | 4,741.910048 | **3,713.295693** |
| Actual objective Q*C | 2,076,774.233715 | 1,043,220.210631 | **527,287.988388** |

This reduces cost **49.4557%** from the previous quadratic and **74.6102%**
from the supplied baseline: improvement factors **1.978464** and **3.938596**,
respectively. Rotation synthesis is priced using the actual full-verifier
margin, with rotation epsilon `1.3917528078052015e-12`. The result equals
the exhaustive-model cost estimate exactly.

The full run's worst error/tolerance ratio is **0.7178548753050588** at
`x = 1,052,634`, matching the integer model. The observed run took 198.03
seconds with four workers; this is a measured duration on the exercised
environment, not a general speed guarantee.

## Why it is smaller

The original hash-pinned 36-bit coefficients remain unchanged. Their rounded
representations use fractional widths `(30,23,16)` and physical register widths
`(30,22,15)`. WL independently certifies these representations on all 192
complete intervals. It keeps the signed Horner truncation error `e0-t*e1`,
which is tighter than adding both error magnitudes. Four more aggressive width
choices fail that certificate and are excluded.

Product scratch is allocated per stage, and phase extension bits are live only
during phase extraction. The quadratic coefficient is measurement-uncomputed
after its first use and reloaded before inverse arithmetic. The cost includes
both added lookup operations. Compiled comparisons also considered retaining
the coefficient, direct controlled addition and smaller adder-workspace budgets.
The selected design wins the bounded comparison; global optimality is unclaimed.

## Evidence

- [Full official result, input hashes and comparisons](results/optimized-quadratic-verified.json)
- [Generator build manifest](results/optimized-quadratic.build.json)
- [Symbolic argument, regimes, documentation and reproduction](optimization/README.md)
- [All eight width certificates](optimization/width-certificates.json)
- [Exhaustive integer screening and resource comparisons](optimization/screening.json)
- [Independent correctness review](optimization/reviews/correctness-1.md)
- [Independent quality review](optimization/reviews/quality-1.md)

The fresh correctness reviewer reproduced all 21 WL checks, every saved
certificate and exact integer error/capacity identities across all inputs.
The Python suite passes 21 checks across the initial, quadratic and optimized
implementations, including nine optimization checks and 16,380 compiled
boundary shots. The full run verifies the serialized binary whose SHA-256 is
`3e8da3474cb20e7db77bb51a86514c3ef422c9a5673b0955a8106d76969396b4`.
Its source and certificate hashes remained unchanged during verification.

This is a full local verifier PASS with the declared salt. It is not a
Kaggle submission or an exhaustive proof over all measurement outcomes.
With the other two tracks at baseline, this inverse-square-root improvement
would yield aggregate factor **1.579236**. The earlier generators and their
verification records are preserved. Generated binaries remain ignored by Git,
and downloaded competition sources remain external.
