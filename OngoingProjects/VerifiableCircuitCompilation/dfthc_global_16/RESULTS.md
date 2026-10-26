# Global (16,32,20) refit results

**No passing candidate was found.** The symbolic gradient gap is repaired, and the new global tight-frame representation preserves the incumbent offsets before quantization, but the bounded global fits remain well above the 20 mHa two-body gate. No circuits, full circuit verification, commit or push were performed.

## Actual accuracy, normalization and whole-task estimates

Each row below is evaluated by the official decoded-word routine with recomputed one-body rows. Both error gates are 20 mHa. All one-body errors and domain checks pass; every row fails the two-body gate. Costs use the measured effective normalization, not the incumbent value.

| Candidate | One-body error (mHa) | Two-body error (mHa) | Lambda | lambda_eff | Total Toffolis | Total qubits | Whole-task cost estimate |
|---|---:|---:|---:|---:|---:|---:|---:|
| Verified C30 incumbent | 13.911732 | 19.980680 | 290.754205 | 47.558542 | 65,244 | 661 | 2,051,023,203 |
| global-beta10 | 13.309520 | 285.336564 | 416.063936 | 58.448531 | 14,037 | 690 | 566,105,005 |
| global-beta12 | 3.233974 | 280.545638 | 415.761185 | 57.070049 | 15,809 | 698 | 629,749,840 |
| global-beta14 | 0.804187 | 280.242245 | 415.760982 | 56.940139 | 18,753 | 688 | 734,645,322 |
| global-beta16 | 0.203080 | 280.233521 | 415.742596 | 56.932929 | 21,697 | 686 | 847,397,803 |
| global-gauged-beta16 | 0.203080 | 280.233521 | 354.583899 | 52.557348 | 21,697 | 686 | 782,271,025 |

The lowest two-body error is **280.233521 mHa**, about **14.01 times the gate**. Its gauged effective normalization is **52.557347620**, and its modeled whole-task cost is **782,271,025.12**. This is an accuracy-rejected resource estimate, not a valid competition score or an achieved improvement.

All resources include forward/reverse angle loading, outer and inner preparation, Givens rotations, persistent space and the shared maximum temporary space. The reference ordinary-unary loader is optimized over every legal batch size. The ANGLE widths are 119, 127, 117, 115, 115 qubits, all within Q128. No new loader was compiled.

## Floating fit and constraints

| Initialization | Global error (mHa) | Iterations | Offset-matrix drift (Ha) | Offset-norm drift |
|---|---:|---:|---:|---:|
| released | 280.329457 | 1200 | 2.283e-13 | 2.274e-13 |
| spectral | 316.179505 | 1200 | 2.263e-13 | 0.000e+00 |

Both floating fits stop at the 1,200-iteration cap with optimizer_success=False. The four decoded dictionaries each receive 300 global weight/allocation iterations and also hit that cap. No convergence or global optimality is claimed. The best floating fit is already far above the gate, so angle precision is not the main limitation in these attempts.

The incumbent offset Z has eigenvalues in [-12.319681640, -6.100522580] and tau=1078.660024986199. The tight-frame map preserves both Z and tau before rounding; this keeps the floating SOS gap at 3.9159356215095613. It permits all 512 projector directions to move while fitting the full target jointly.

After quantization the offset matrix is recomputed rather than assumed fixed. Its Frobenius drift falls from 0.239765319 Ha at 10 bits to 0.003770371 Ha at 16 bits. The corresponding actual SOS gaps are 4.125872680 and 3.916725397. The offset norm remains fixed up to roundoff; all reported normalizations use these actual gaps.

## What the new mathematical representation establishes

- The saved factor-residual derivative now retains every weight entry symbolically in a rectangular example. Its diagonal cross-weight terms are nonzero symbolic polynomials. The earlier numeric-weight specialization is absent from this study.
- The global objective is the full shifted pair-tensor residual. WL derives its complete projector and weight gradients, with cross-block interactions and arbitrary symbolic coefficients.
- A transformed Parseval frame parameterizes an exact global one-body offset and offset norm. The polar pullback, allocation gradient and PSD core update are derived symbolically; independent finite differences check the full numerical chain.
- Fixed-frame core initialization has a closed-form best PSD rank-19 update against the current global residual. The nonlinear frame/weight fit remains a local search.
- The offset-carrying factor in each block is restricted to be negative semidefinite. This covers a specific offset-preserving family, not every possible (16,32,20) factorization. Two failed bounded searches are not an impossibility theorem.

## Normalization gauge

The one 150-iteration copy-gauge fit reduces Lambda_H2 from 402.440273578 to 341.281577140. Lambda_eff falls from 56.932929206 to 52.557347620. The pair-tensor drift is 1.001e-14 Ha; the native offset contraction and all two-body words remain unchanged. Both official reconstruction errors are unchanged. The gauge optimizer hits its cap and is not claimed globally optimal.

## Verification evidence and score status

The local kernel passes 53 exact checks. Nine Python tests check the complete polar/global/fixed-dictionary gradients, copy-gauge gradient, frame/offset invariants, direct global residual, core-update monotonicity and domain failures. verify-data.wls independently checks all five full-resource estimates, legal batches, actual normalizations, provenance and recorded score arithmetic. Fresh independent review verdicts are recorded separately.

The verified incumbent cost remains 2,051,023,203.37; the combined local estimate remains 2.035845221582. No new circuit-verified score is claimed. Raw inputs, floating candidates and five quantized factor/word tables remain in ignored local/. Old study directories and the verified incumbents have not changed.

Full methods and reproduction commands are in README.md. Machine-readable values, iteration histories, caps, source hashes and candidate hashes are in results.json, symbolic.json, inputs.json and wl-data-checks.json.
