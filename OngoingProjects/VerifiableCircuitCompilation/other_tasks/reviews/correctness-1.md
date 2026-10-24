# Correctness review — round 1 of 3

## VERDICT

PASS. The stated symbolic identities, compiled-fragment resource and simulation
claims, unchanged full baselines, bounded mixing result, evidence bindings and
test counts were independently reproduced. No target artifact was edited. All
generators ran in temporary mirrors; verifiers, circuit data and integrals stayed
external to the repository.

The review covers `OngoingProjects/VerifiableCircuitCompilation/other_tasks/`:
all four WL/WLS sources, all five Python sources, README, PLAN and every JSON
evidence file. The ignored candidate NPZ was read only to check its recorded
SHA-256. Relevant external verifier contracts, mathematical reconstruction,
native multiplication/inversion implementation and live method documentation
were also inspected. No applicable AGENTS.md was found.

## REPRODUCED

### Environment and commands actually run

- Real kernel: Wolfram Language 15.0.0 for Mac OS X ARM (64-bit), February 17,
  2026. The direct `WolframKernel -noinit -noprompt -script` invocation was used.
- Python: the pinned kickmix environment, kickmix 0.1, NumPy 2.5.4.
- The external kickmix source HEAD is
  `638777fac4379bc62e2a1e4a52cf1a1753895d48`.
- `$PYTHON`, `$KERNEL`, `$DATA` and `$MIRROR` below denote the exercised Python
  runtime, direct kernel executable, external official-data directory and
  temporary mirror; personal absolute paths are omitted.

From the mirror's parent project directory:

```sh
"$PYTHON" other_tasks/profile_baselines.py dqi --data "$DATA" --workers 1
"$PYTHON" other_tasks/profile_baselines.py dfthcsos --data "$DATA"
"$KERNEL" -noinit -noprompt -script other_tasks/run-symbolic.wls
"$KERNEL" -noinit -noprompt -script other_tasks/verify-symbolic.wls
"$PYTHON" other_tasks/field_probes.py
"$PYTHON" other_tasks/profile_stages.py --data "$DATA"
"$PYTHON" other_tasks/component_mixing.py --data "$DATA"
"$PYTHON" -m unittest other_tasks.test_feasibility -v
```

The baseline and generator stdout logs were parsed in their entirety: they
contain only the expected running/static records and JSON reports, with no
warning/error text. The final WL verification reports 25 succeeded and 0 failed;
the Python suite reports 8 tests, OK. Package WL generator and final verification
runs emitted no kernel messages. No Quiet or message suppression was added.

An initial mirror-only verification overlapped the DQI generator's RUNNING
report and consequently failed its two baseline-confidence comparisons. Running
the established verifier against a separate stable mirror produced 25/25 and
the exact recorded `wl-checks.json`; the eventual full DQI report independently
matched that stable baseline. This was a reviewer scheduling collision, not an
artifact finding.

### Full official baseline and candidate results

| Run | Full coverage | Reproduced result |
|---|---|---|
| Unchanged DQI | 102,400 shots | PASS; zero value and phase errors; Q=1,893; mean executed Toffolis=5,725,107.220107422; static Toffolis=5,725,461; score=24,439,458,623,294.45 |
| Unchanged DFTHC | 254,464 shots, 3,976 addresses × 64 | PASS; Q=730; Toffolis=93,237; score=3,355,599,588.3155365 |
| Pair-coordinate mixing, original angle circuits | 254,464 shots | PASS; score=3,313,336,907.312215; factor=1.0127553225601815; lambda_eff=48.68043460482008 |

Every non-duration field of all three regenerated reports equals its original
record, including all metrics, resources, operation counts, submission hashes,
search counts and candidate hash. The unchanged DFTHC breakdown includes 79,480
angle Toffolis, 13,757 preparation/Givens Toffolis and shared temporary peak 544.

The DQI confidence calculation independently gives epsilon-plus
0.0008660587301657926 and repetition factor 1.1912604648115832 at delta=2^-128,
despite zero observed errors. The score formulas match the official source and
the recorded arithmetic.

### Field arithmetic and compiled fragments

- Live native kickmix field is GF(4096), modulus `0x1053` or
  `x^12+x^6+x^4+x+1`, primitive element 2. WL's matching polynomial field has
  characteristic 2, extension degree 12, 4,096 elements and polynomial
  representation. All 4,096 coefficient/index round trips passed.
- The GF(64), GF(16), GF(8) plans regenerate exactly, with native subfield-root
  indices 270, 1971, 937, respectively. Root, basis, interpolation and complete
  1,728-entry coefficient-tensor certificates are true for each plan. Evaluation
  map dimensions are {3,6,12}, {5,4,12}, {7,3,12}; the subproduct counts are 3,5,7.
- The GF(8) six-product polynomial identity is exact modulo `x^3+x+1`.
- All 4,095 nonzero inverse identities and the explicitly extended zero-to-zero
  convention passed. Generator order was checked against all prime divisors of
  4,095. Unsupported/inexact degrees and wrong argument counts are rejected.
- `FiniteField`, element coefficient/index properties and finite-field matrix
  inversion were checked against installed documentation and a real kernel.
  Independent modular-binary and finite-field matrix inverse probes passed.

Each cell below is Q / Toffolis. All twelve generated fragments match
`field-results.json` byte-for-byte, including operation counts and checks.

| Route | Multiply | Controlled multiply | Inverse phase |
|---|---:|---:|---:|
| Native | 36 / 63 | 49 / 75 | 72 / 315 |
| GF(64) | 54 / 63 | 67 / 75 | 90 / 315 |
| GF(16) | 48 / 45 | 61 / 57 | 84 / 225 |
| GF(8) | 46 / 42 | 59 / 54 | 82 / 210 |

Each fragment was freshly simulated on 8,192 basis shots, covering every native
field input with independent partner words, dirty XOR targets, both control
values where applicable, restored inputs, clean scratch and correct phase.
The Python suite additionally reran GF(8) inverse phase with masks 0 and 0xA55.
Live builder documentation and native source confirm XOR-target semantics,
controlled-product scaffolding, zero-Toffoli measurement cleanup and the
five-multiplication inverse scaffold.

As a separate oracle, a hand-written polynomial-bit multiplier, independent of
both WL finite-field arithmetic and kickmix's multiplication implementation,
checked all 144 basis pairs (1,728 output coefficients) and 4,096 additional
random products for each of the three exported routes. All passed. It also
checked all 64 GF(8) input pairs against the exported six-product formula.

### Stage profile and conditional decoder arithmetic

The regenerated `dqi-profile.json` equals the original complete JSON. This
reproduces the rounded flame-chart stage profile, 58,211 static native
multiplications and fraction 0.6405236189714679. The model arithmetic gives
factors 0.9811897837870002, 1.2085537281672043 and 1.258063132670707 for the
three routes. The source and README explicitly retain the necessary assumptions
about static rather than executed counts, scratch peak and repetition factor.

### Component mixing

- The orthogonal, bilinear two-body tensor, native-offset contraction and
  squared-offset norm identities all reproduce as exact symbolic identities.
  Identity rotation, limiting sign flip and exact 2× two-body normalization
  example reproduce.
- The pair objective has the stated constant + A cos(2 theta) + B sin(2 theta)
  form on fixed sign regions. Source inspection confirms that the numerical
  search evaluates sign boundaries and interior stationary candidates over a
  pi/2 period. Dense-grid reference tests pass, including zero and degenerate
  pairs. The code does not claim a global mixing optimum.
- Two rounds of 60 prioritized pairs over 65 groups give 7,800 attempts and
  5,944 accepted rotations. Raw Lambda, effective normalization, original angle
  words, dimensions, fixed shifts, float32 conventions and full official
  acceptance reproduce. Relative augmented-Gram drift is
  3.2140392044921713e-16. The spectral control worsens cost to
  3,496,650,989.2667217, as reported.
- Additional live WL boundary tests covered ten real exact/numeric parameters,
  all yielding 2×2 orthogonal matrices to the stated numerical tolerance; eight
  complex, malformed or nonfinite rotation inputs and eleven invalid field
  degrees were rejected. Both public functions reject wrong argument counts.

### Evidence integrity and timing treatment

All recorded official-input hashes, baseline generator hashes and twelve
source/dependency evidence bindings match the frozen source files. The ignored
candidate NPZ matches its recorded SHA-256, and `git check-ignore` confirms it is
excluded. The public artifact set contains source and reports, with no vendored
official circuits, integrals or verifier source.

`symbolic.json`, `field-results.json`, `wl-checks.json` and `dqi-profile.json`
reproduce exactly. The full baseline and mixing reports reproduce every field
except `elapsed_seconds`. Historical durations are observations, not portable
speed claims. Reviewer durations were 498.746403541998 seconds for DQI,
75.04663916601567 seconds for DFTHC and 110.2936766660423 seconds for mixing;
these runs shared host resources. No speed benchmark is claimed or certified.

## FINDINGS

None. There are no unresolved correctness findings, unavailable dependencies,
failed target checks or unverifiable required claims.

The documented limits remain material: no production DQI decoder or whole-decoder
substitution was implemented or verified; conditional decoder factors are not
candidate scores; floating pair search is not a certified global optimizer; the
measurement-based fragment simulations sample branches rather than exhaust
them. There is no notebook deliverable or notebook-behavior claim. No commit,
push, external submission or production change was performed by this review.

## RECOMMENDED FIXES

None required for the stated baseline-profiling and feasibility scope.

OPEN ISSUES: 0
