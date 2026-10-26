# Correctness review, round 2/3

## VERDICT

PASS for the stated mathematical, fixed-table and conditional resource claims. The complete target sources and claim reports were read; no prior reviews were consulted. All reproductions used temporary copies, and the target files were unchanged. No circuit, reversible decoder, native BM integration, performance improvement or new score is certified by this review.

## REPRODUCED

Environment: Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026; Python 3.13.1; NumPy 2.5.4; kickmix 0.1. The native field is GF(2^12), modulus `0x1053`, primitive element 2.

Commands, with machine-specific executable/data locations supplied externally and the working tree copied into a temporary directory:

```text
KERNEL -noinit -noprompt -script TEMP/representation_search/run.wls
PYTHON TEMP/representation_search/curvature_screen.py --data DATA
PYTHON TEMP/representation_search/affine_rounding.py --data DATA
PYTHON TEMP/representation_search/table_screen.py --words WORDS
KERNEL -noinit -noprompt -script TEMP/independent.wls
```

All four project commands exited 0. Their regenerated JSON structures exactly match the four supplied reports, including source hashes. The project Wolfram run emitted all 34 Boolean checks as `True`, status `PASS`, with no kernel messages. Independent verification used `Echo`, no `Quiet`, and compared every state against a separately coded normalized BM recurrence.

Claim inventory and independent checks:

1. **Curvature** (`RESULTS.md:10–36`): expansion gives the asserted residual; its three alternating extrema and exact bound are correct. The second-difference certificate gives `4E >= 1/2`, hence the minimax bound `E >= 1/8`. The zero-change limit is zero. All 3,133,452 legal integers reproduce every candidate report: failure counts are 1,590,097, 385,356, 186, 0, 0, 0 for the listed candidates. The selected `(30,23,10)` model has physical widths `(30,22,9)`, maximum ratio `0.8505589342620424` at x=1,064,928, and no nonpositive or overflowing stage values. Positive decrement terms make each stage monotone, so stage-end positivity also bounds intermediate steps. Independent exact scalar nested arithmetic agrees on 20,139 boundary/random inputs. Current frozen quantization has physical widths `(30,22,15)`, confirming six fewer coefficient bits in the model.
2. **Affine rounding** (`affine_rounding.py:12–33`, `RESULTS.md:30–35`): independently derived ceil/floor bounds match the implemented interval formula, including x=1's phase wrap. Exhaustive regeneration reproduces five repaired rows and row 170 as the sole unresolved interval over all 65 searched B offsets. The restricted-family rejection is correctly qualified.
3. **Cached quotient and native adapter** (`Representations.wl:6–65`, `RESULTS.md:40–85`): live `PolynomialExtendedGCD` returns `{g,{s,t}}` in the expected shape, and its Bezout identity holds over GF(16). Live `PolynomialReduce` returns `{quotients,remainder}`; the ideal remainder is zero. The finite-field certificates, degree-1 and degree-5 cases and rejected domains reproduce. Additional probes reject zero/constant/nonpolynomial locators, nonpolynomial evaluator, negative/noninteger/machine-real h, invalid variable and arity; a simple zero root is valid at h=0, and a zero evaluator returns zero. The released scanner source confirms n=70, m=73, t=35, evaluator slots 37..71 with exponents 34..0, and odd locator slots summing `z sigma'`. Independently evaluated all 4095 actual primitive-power scan points using kickmix field operations. All 35 roots agree with the Wolfram cache; all 4096 actual polynomial-basis masks at each root agree after explicit trace-dual conversion: 143,360 mask comparisons. Nonroots have the same conditioned phase +1. The incorrect h34 adapter's counterexample is detected. The degree-1/2/17/35 packaged phase cases also pass.
4. **Projective recurrence** (`ProjectiveDecoder.wl:7–48`, `RESULTS.md:89–115`): the recurrence itself contains no division. Independently compared every projective state, after normalization, with a separately implemented ordinary BM recurrence for degrees 0,1,2,3,17,34,35. All locators, Forney values, nonzero constant scales and degree bounds agree; the degree-2 case used equal weights to force an initial zero discrepancy. The general update and scale-cancellation identities and both determinants reproduce. The nonzero-scale/domain restrictions and absence of a reversible/native-syndrome implementation are explicit.
5. **Conditional budgets** (`RESULTS.md:71–83`): inspected the released scanner and current GF(8)-tower inverse chain. Five 42-operation products plus 12 root-masked ANDs per point give `4095*222 = 909090`. The current manifest supplies 4,503,030 operations and 1903 qubits. Equal-objective width is `1903 sqrt(4503030/3593940) = 2130.128471604313`; adding a full 420-qubit block would multiply the objective by about 1.18929 even with free setup. For degree-<d residues and a monic degree-d modulus, d² product operations plus d(d−1) reduction operations give `2d²-d`; d=35 yields 2415 and four expressions at 42 cost 405,720. These are expression budgets with the stated exclusions.
6. **DFTHC support and rank** (`table_screen.py:47–73`, `RESULTS.md:119–137`): the pinned word hash matches. All 25 batches reproduce the full degree histogram, 994 legal rows, 983 nonlinear support, total output rank 108 and difference-bank rank 81. Independent use of the existing loader's transform reproduces prefix completion, unchanged prefix coefficients under random invalid suffixes, and support preservation under an invertible output shear. The original word sign-bit redundancy was checked directly. The subset argument proves the arbitrary-prefix statement; no claim about general AND/XOR lower bounds follows.
7. **Scope/provenance**: incumbent score 2.0088882102784944 agrees with existing verification metadata; the incumbent circuits were not rerun because no new circuit or score is claimed. The published [DQI design](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html) explicitly caches the derivative inverse; the [iBM paper](https://eprint.iacr.org/2022/469) describes the scalar rescaling and Algorithm 2. Plans and statements of future work are scope declarations, not computed outputs. No performance benchmark is claimed, so repeated timing trials are inapplicable. No notebook is present.

Complete packaged Boolean-check enumeration:

- `curvature_identity`: reproduced `true`.
- `curvature_uniform_bound`: reproduced `true`.
- `curvature_equioscillation`: reproduced `true`.
- `affine_minimax_lower_certificate`: reproduced `true`.
- `integer_intercept_interval`: reproduced `true`.
- `curvature_zero_change_limit`: reproduced `true`.
- `GF12_degree35_general_quotient_h34`: reproduced `true`.
- `GF16_degree1_prefactor0`: reproduced `true`.
- `GF16_degree5_prefactor71`: reproduced `true`.
- `general_quotient_ideal_certificate`: reproduced `true`.
- `repeated_root_rejected`: reproduced `true`.
- `zero_root_with_prefactor_rejected`: reproduced `true`.
- `empty_locator_rejected`: reproduced `true`.
- `malformed_prefactor_rejected`: reproduced `true`.
- `trace_dual_basis_all_4096_masks`: reproduced `true`.
- `cached_root_phase_unitary`: reproduced `true`.
- `prefix_independent_of_invalid_suffix`: reproduced `true`.
- `mobius_involution_sizes_2_to_64`: reproduced `true`.
- `invertible_output_basis_preserves_support`: reproduced `true`.
- `native_scan_reversed_evaluator_h1`: reproduced `true`.
- `native_scan_masked_phases_degrees_1_2_17_35`: reproduced `true`.
- `native_derivative_accumulator_is_z_sigma_prime`: reproduced `true`.
- `incorrect_h34_adapter_counterexample_detected`: reproduced `true`.
- `projective_general_update_identity`: reproduced `true`.
- `projective_Forney_scale_cancels`: reproduced `true`.
- `major_update_matrix_determinant`: reproduced `true`.
- `minor_update_matrix_determinant`: reproduced `true`.
- `native_field_35_errors_locator`: reproduced `true`.
- `native_field_35_errors_Forney`: reproduced `true`.
- `native_field_35_errors_degree_bound`: reproduced `true`.
- `native_field_35_errors_nonzero_scale`: reproduced `true`.
- `projective_one_error_case`: reproduced `true`.
- `projective_zero_error_case`: reproduced `true`.
- `zero_projective_scale_is_invalid`: reproduced `true`.

## FINDINGS

None within the stated scope. The two string-valued checks explicitly limit phase/unitarity and classical-decoder conclusions. Target scripts emitted no messages. Verifier-only exploratory representation queries were corrected before the independent clean rerun; they did not expose a target defect.

Source and claim binding uses SHA-256 of each file's bytes. Bundle hashes are SHA-256 of the corresponding filename-to-hash JSON mapping, with sorted keys and compact separators. The source bundle contains the six executable files; the claim bundle contains the two Markdown and four report files.

- Source bundle: `c8d27b881fbb145f6b2c09df485218520263559c72a8b22f3d5b699519a39f17`
- Claim bundle: `9236d277815864cc2a9cae1437ca6396d48743d26ddaf13be65df77f3b6ca83b`

| File | SHA-256 |
|---|---|
| `Representations.wl` | `85817e6845e6dfaf53673e9773d2b5ca168ebbb2be0c56e611e081e45ed50a19` |
| `ProjectiveDecoder.wl` | `e0a03b7d7b1c394936d273747ce39c88c556426cf40b9d7525f22e69dff554bb` |
| `run.wls` | `2d85571bb326105edf2a47c2969a78cbcafe3c69487156cde8bd2c34a6cec205` |
| `curvature_screen.py` | `821a6119382eaa0d7329688ae833444fc5c39a6196bf3ace84ec50ec6a91f0f3` |
| `affine_rounding.py` | `8bd41280d90b06e502101977b3e088b1c66af8612b8d9eff012e3de8933c35de` |
| `table_screen.py` | `60b730ca9286d2e2e6935a1308ec618ca123e5d7321d76103cb77004f7e1d051` |
| `RESULTS.md` | `bd1062c8e69ef41b78bcf261cd36ad06ef203f5903970fa0acbf9c4873f9743d` |
| `PLAN.md` | `4d2c0833bd60bd73d921aa54dfe2e83151447a2c7332a1ccf4e0084dfda6a745` |
| `symbolic-checks.json` | `78035f544058ad635240b030fe4fb828ed1b9b00a2990df09bca87429a1d8273` |
| `curvature-screen.json` | `9e8e79ecfbc3390a26595ac9168062139f3a28e97c86ae3bd9b981f21dd3174d` |
| `affine-rounding.json` | `baa761187fb066ee3e9b3266b41b7225847f04340d5636c9dbfac2dbd722402b` |
| `table-screen.json` | `a2ba5abfbe00d5e4b0ee8e8d1b10325eb304c9f0eafaf7ad1c0071ed2813cda5` |

## RECOMMENDED FIXES

None required. Keep the existing limitations on reversible setup, branch history, workspace, rotation synthesis and uncompiled resource budgets when choosing an implementation.

OPEN ISSUES: 0
