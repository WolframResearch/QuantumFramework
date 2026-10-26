# VERDICT

PASS for the authorized mathematical global (16,32,20) study. No unresolved in-scope correctness or reproducibility issue was found. No candidate passes the official accuracy gates, and the artifacts correctly make no new circuit-verified score claim.

# REPRODUCED

All thirteen named target artifacts were read completely; numeric histories were parsed and checked in full. Directly required helper, resource, official mathematical/preparation/reference definitions and the six named recorded reports were inspected. Only named dependencies were copied into an isolated temporary replica; large raw inputs were linked read-only. Original sources/results were checked against the frozen copies after verification and remained identical.

The README workflow was executed in that replica with the supplied Python runtime and Wolfram kernel, using `OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1` for Python:

```sh
"$PYTHON" refit.py prepare --data "$DATA"
"$KERNEL" -noinit -noprompt -script symbolic.wls
"$PYTHON" -m unittest test_refit -v
"$PYTHON" refit.py run --data "$DATA"
"$KERNEL" -noinit -noprompt -script verify-data.wls
```

Here `DATA` denotes the supplied official bundle; runtime aliases denote the supplied executables. Versions were Python 3.13.1, NumPy 2.5.4, SciPy 1.18.1 and Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026. Every reproduction command exited 0. Both WL scripts emitted PASS, with no kernel messages. The nine numerical tests all passed.

Claim coverage and independent cross-checks:

- **Rectangular factor derivative and augmented-core identities:** all 53 saved exact checks reproduced, including every symbolic coefficient, nonzero diagonal cross-weight polynomials, singular Schur range condition, PSD/rank/gauge identities, projection errors, compression/rescaling identities and their limiting/failure fixtures. An additional independently formed 3-row versus 4-row, 3-dimensional factor residual differentiated exactly to the stated general rectangular formula.
- **Global objective and gradients:** the generic symbolic cross-block projector/weight checks reproduced. Additional full-map finite differences covered shapes (3,3,3) with four orbitals and (2,4,4) with three orbitals, for both moving and fixed dictionaries. Twelve independent directions per map, central step 1e-5, gave maximum normalized discrepancy 4.12e-9. The full shifted pair target was separately rebuilt from the official reference and fixed shifts; maximum entry discrepancy from the pinned pair matrix was 1.39e-17.
- **Tight frame, polar/Sylvester, allocation and LCU chamber:** the exact frame/offset, tangent, diagonal-basis pullback, dimension-independent Sylvester trace, positive allocation and arbitrary sign-chamber checks passed. An additional rational orthonormal frame with a non-diagonal positive-definite H verified the polar adjoint exactly. Numerical polar checks also passed at repeated singular values, for a tall frame and a square reflection. Independent finite differences covered the implemented gauge derivative through the polar map.
- **Fixed-frame PSD update:** the orthogonal-span identity and mixed-sign rank-one boundary reproduced exactly. A separate indefinite-residual fixture compared the implemented rank-two coordinate update against an independently whitened and spectrally truncated optimum; matrices and direct residual norms agreed within 1e-12. Initialization histories were monotone, including all four recorded sweeps for both full-size seeds.
- **Signatures, return shapes and domains:** live supplied definitions were checked, rather than inferred from documentation. AugmentedCore returned the expected square augmented matrix, ProjectorGram a rectangular matrix, and ProjectorGradient the original direction-array shape; malformed, complex and empty WL arguments returned Failure. Relevant Python failures rejected rank-deficient/nonfinite polar inputs, unsupported offset signs, nonpositive tau, zero projector rows, singular dictionaries and wrong parameter lengths. Empty zero-orbital Python evaluation also rejects, by IndexError; no general validation/error-type contract for zero-dimensional models is claimed or relied upon by this fixed 76-orbital study. The live RankCost boundary, degeneracy and malformed-input checks passed.
- **Full numerical reproduction:** all 3,869 numeric leaves and 97 other leaves of results.json matched exactly after excluding elapsed_seconds. Prepared arrays, both floating fit arrays and all five decoded candidates were bitwise equal to their frozen counterparts; their recorded hashes also matched. Both floating fits reached 1,200 iterations; every decoded refit reached 300; the sole gauge fit reached 150. No optimization budget was increased.
- **Independent raw-candidate reconstruction:** 134 separate checks passed. These used saved factors/words, independent pair and offset contractions, native offset conversion, official loading/evaluation, separately recomputed one-body eigenvalue matching, LCU normalization, SOS energy/gap and effective normalization. Realness, row normalization, augmented-core PSD/rank, offset norm and gauge invariants were checked. Both floating offsets remained within about 2.3e-13 Ha of the incumbent and their SOS gap remained 3.9159356215095613.

The independently recovered final candidate results were:

| Candidate | One-body error (mHa) | Two-body error (mHa) | lambda_eff | Toffolis / qubits / ANGLE width / batch | Modeled cost |
|---|---:|---:|---:|---|---:|
| global-beta10 | 13.309519895 | 285.336563705 | 58.448531437 | 14037 / 690 / 119 / 10 | 566105004.688 |
| global-beta12 | 3.233973954 | 280.545637746 | 57.070048760 | 15809 / 698 / 127 / 9 | 629749839.796 |
| global-beta14 | 0.804186922 | 280.242245082 | 56.940139336 | 18753 / 688 / 117 / 7 | 734645321.880 |
| global-beta16 | 0.203080104 | 280.233521245 | 56.932929206 | 21697 / 686 / 115 / 6 | 847397802.785 |
| global-gauged-beta16 | 0.203080104 | 280.233521245 | 52.557347620 | 21697 / 686 / 115 / 6 | 782271025.115 |

All five pass the one-body gate and official normalization/gap domain floors, and all five fail the 20 mHa two-body gate. Independent resource enumeration over batches 1 through 75 recovered every selected legal minimum, using the actual candidate lambda_eff, official preparation/Givens costs, persistent space and maximum shared temporary space. The independent WL data script also passed all eleven checks and priced all five rows.

The gauge preserved the angle words and augmented cores; independent pair-tensor drift was 1.01e-14 Ha. Native offsets were unchanged. Lambda_H2 decreased from 402.440273578 to 341.281577140 while both official errors remained unchanged. Quantized offset drift and recomputed gaps matched the saved values, including 0.239765319 Ha at beta=10 and 0.003770371 Ha at beta=16.

The six historical report hashes, PASS statuses, individual score products and three baseline/current ratios matched. Their geometric mean is 2.035845221582375 and the recorded incumbent cost remains 2051023203.3721037. These historical circuit runs were not rerun, as required. Raw local artifacts are ignored by Git.

# FINDINGS

None requiring a fix within the authorized scope. The evidence supports the representation, gradients, core update and quantization statements in `OngoingProjects/VerifiableCircuitCompilation/dfthc_global_16/README.md`, together with the numerical, resource and scope statements in `OngoingProjects/VerifiableCircuitCompilation/dfthc_global_16/RESULTS.md`.

The scope limits are accurate and material: the offset-carrying factor is negative semidefinite in each block; logits are bounded; all nonlinear fits hit their caps without a global-optimality claim. Therefore the failures exclude these bounded attempts, not every (16,32,20) representation. The best error is approximately 14.01 times the gate. Resource rows are analytical estimates for accuracy-rejected word tables, not achieved competition improvements. Circuit construction, simulation, cleanup verification and historical circuit revalidation are outside this study. elapsed_seconds is run metadata, not a performance-improvement benchmark; no speed claim was evaluated.

# RECOMMENDED FIXES

None for the mathematical study or its recorded evidence. Preserve the present numerical-failure and score-status qualifications.

OPEN ISSUES: 0
