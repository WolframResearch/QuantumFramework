# Correctness review, round 1

This review is bound to the original artifact set:
`Representations.wl` SHA-256
`166a96be01652c277e9f3f3aa7b121e2eccfeea94ea8a31940374c385465acf2`
and `RESULTS.md` SHA-256
`a4abe8b97977ee41753359a8e207b483ae7409f96372220e6c48e2c6826d174e`.
Later revisions are outside this verdict.

## VERDICT

Changes required: the general cached-quotient identity is correct, but the claimed
mapping to the released DQI scanner uses the wrong evaluator ordering and
prefactor. One native-field counterexample changes an observable phase.
The compensated-curvature and fixed-table screens reproduce.

All eleven requested source, claims, and report files were read, together with
their directly required generator, verifier, loader, and existing resource
evidence. No applicable AGENTS.md imposed additional rules. Target source,
reports, local tables, circuits, and factors were left unchanged; all regenerated
reports were written to a temporary directory. Original-file hashes remained
unchanged before this review file was added.

## REPRODUCED

Environment: Wolfram **15.0.0 for Mac OS X ARM (64-bit), February 17, 2026**;
the supplied Python environment with NumPy and kickmix; native field modulus
**4179**, primitive element **2**. Commands, with executable and external-data
locations supplied by the caller:

```text
$WOLFRAMKERNEL -noinit -noprompt -script TMP/representation_search/run.wls
$PYTHON TMP/representation_search/curvature_screen.py --data DATA
$PYTHON TMP/representation_search/affine_rounding.py --data DATA
$PYTHON TMP/representation_search/table_screen.py --words ../dfthc_c30/local/words.u16
$WOLFRAMKERNEL -noinit -noprompt -script TMP/independent-native-check.wls
```

Each command exited successfully without WL messages. All four regenerated
JSON reports equal their supplied counterparts as parsed objects, including
source and input hashes. The claim inventory and reproduction are:

| Claims or examples | Independent result |
|---|---|
| Curvature identity, uniform bound, alternating extrema, affine minimax lower certificate, zero-change limit | Exact WL checks pass; the second difference is 1/2, yielding the claimed lower bound 1/8. |
| Integer-intercept interval equivalence | Exact WL check passes. |
| All six exhaustive nested integer screens and their register bounds | All 3,133,452 inputs rerun; every aggregate matches. The f10 compensated model has widths (30,22,9), zero phase/capacity violations, and worst ratio 0.8505589342620424 at x=1064928. |
| Too-coarse curvature rejection | Reproduced failure counts: f8 naive 1,590,097; f8 compensated 385,356; f9 compensated 186. f11/f12 both pass with the reported ratios. |
| Finite affine search, B offsets -32..32 | Reproduces exactly five changed rows and unresolved row 170. Its stated finite-search limitation is correct. |
| General quotient ideal identity; GF(16) degree-one h0, degree-five h71; GF(4096) degree-35 h34 | All pass as abstract algebraic identities for the chosen exponent. They do not establish the native scanner mapping; see finding 1. |
| Repeated-root, zero-root h>0, constant-locator, malformed-exponent failures | Reproduced. Additional malformed polynomial/exponent cases return Failure; a zero root with h0 is correctly accepted. A symbolic rational-coefficient quotient certificate also passes. |
| Native trace Gram nondegeneracy; all 4096 masks; quadratic phase signs | Gram determinant is nonzero; its linear-algebra consequence covers all masks. Phase signs are {1,1,1,-1}, with unit squares. This verifies only the phase identity, as explicitly scoped. |
| Static removal ceiling, width ceiling, polynomial multiplication budget | 4095(5*42+12)=909090; 1903 sqrt(4503030/(4503030-909090))=2130.1284716043; 2*35^2-35=2415; four products at 42 cost 405720. The 1903/4503030 incumbent figures agree with existing evidence. These remain conditional arithmetic budgets. |
| DFTHC exact support and ranks, all 25 batches | Every batch reproduces legal prefix 994, nonlinear support 983, output rank 108, difference rank 81, and the complete degree histogram. The loader uses the same packing and signed-word convention. |
| Prefix freedom and invertible output changes | WL triangular/involution checks pass. The subset argument proves arbitrary-prefix invariance, and an invertible linear map preserves zero/nonzero coefficient vectors. This gives no general AND/XOR lower bound, as stated. |
| Ignored derived tables, historical combined score | Git confirms derived tables are ignored. Existing verified cost ratios reproduce the stated combined arithmetic 2.0088882102784944. Historical circuit verification was inspected rather than rerun; no new circuit/score claim is made here. |

The documented three-argument PolynomialExtendedGCD return shape and finite-field
support agree with the [official reference](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html).
PolynomialReduce's quotient/remainder shape agrees with its
[reference](https://reference.wolfram.com/language/ref/PolynomialReduce.html).
The [cited DQI design](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html)
does describe caching the derivative inverse. No speed benchmark or notebook
behavior is asserted by these targets, so neither requires a performance rerun
or notebook inspection. No dependency or kernel limitation reduced these checks.

## FINDINGS

1. **[P1] Native DQI mapping does not match the released scanner.**
   `RESULTS.md:48-51,64-65` claims that relabeling evaluator slots 37..71 gives
   h=34, and `Representations.wl:65` labels its generic test as native evidence.
   The released scanner's actual multipliers and summation are in
   `verifiers/dqi/dqi_reference_circuit.py:545-560`. Put
   **beta = primitive^j**, let **S(z)=sum(b_k z^k,k=0..35)** from locator slots,
   and let **w_k** be evaluator slots after prealignment. The scanner evaluates

   ```text
   N(beta) = sum(w_k beta^(71-k), k=37..71)
   D(beta) = sum(b_k beta^k, k odd) = beta S'(beta).
   ```

   Thus a directly matching convention is
   **W(z)=sum(w_k z^(71-k),k=37..71)** with **h=1** and locator **S(z)**.
   Relabeling slots in ascending order produces a different polynomial.
   A reciprocal-locator convention must also transform the derivative and
   account for the locator degree; it cannot justify the stated h34 merely
   from the index range.

   Exact counterexample over the native field: choose
   **S(z)=product(z-F[r],r=2..36)**, **w_37=1**, and all other evaluator slots zero.
   At the locator root **beta=F[2]**, the released scanner produces
   **beta^33/S'(beta)=F[3677]**. The stated ascending relabeling gives Omega=1,
   whose h34 cache produces **beta^(-34)/S'(beta)=F[2751]**. Classical mask 2
   gives opposite phase signs, +1 and -1. A real-kernel check returns
   `released_scan_vs_h34 -> False`; reversing the evaluator and using h1 gives
   `released_scan_vs_reversed_h1 -> True` at all 35 roots.

   The current test independently verifies the algebra for an assumed h; it
   never constructs the released numerator or its odd-index derivative sum.
   Therefore all supplied WL booleans can pass while this native claim fails.

## RECOMMENDED FIXES

Define the locator variable and the evaluator coefficient order explicitly.
Use the released scanner expressions above as an independent oracle, including
the odd-slot derivative accumulator. Test the proposed cached quotient against
that oracle over native-field locators of several degrees, including degree 35,
and retain a nontrivial classical-mask phase comparison. Correct the h34/h71
claims and regenerate the named native check and report. The generic
CachedQuotient implementation and the reproduced curvature/table screens need
no change for this finding.

OPEN ISSUES: 1
