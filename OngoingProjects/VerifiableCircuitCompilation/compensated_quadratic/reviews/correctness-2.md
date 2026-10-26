# VERDICT

PASS — independent correctness round 2/3. The complete official verification and all declared mathematical, table, capacity, grid, circuit, and cost claims reproduce. The 13 target artifacts and candidate binary were bound before and after this review and remained unchanged. No earlier review, review summary, neighboring review file, or representation-search material was inspected. No applicable ancestor AGENTS.md was present.

# REPRODUCED

All target sources and complete JSON structures were read. Numeric tables were parsed and checked across every row: 192 frozen rows, 192 compensated rows, 192 source interval records, and 192 derived certificate records. The bounded JSON parse covered 16 explicitly named dependency/target JSON files and 14,289 leaves; no table was assessed from a sample.

The executed commands below use environment placeholders to avoid personal absolute paths. Two isolated copies separated the immutable full-verifier inputs from regeneration. `KERNEL` was the installed Wolfram 15.0 kernel, `PY` the specified kickmix virtual environment, and `DATA` the supplied official download.

```sh
"$KERNEL" -noinit -noprompt -script "$WLCOPY/run.wls"
"$KERNEL" -noinit -noprompt -script "$WLCOPY/grid.wls"
"$PY" "$WLCOPY/generator.py"
"$PY" "$WLCOPY/screen.py" --data "$DATA"
"$PY" "$COPY/verify.py" --data "$DATA" --workers 4
```

Every command exited 0. Target WL runners emitted only their expected Echo reports, without kernel messages. Versions: Wolfram 15.0.0 for Mac OS X ARM (64-bit), February 17, 2026; Python 3.13.1; NumPy 2.5.4; kickmix 0.1. Independent temporary WL and Python checks were also executed, without changing the targets.

| Claim or expected result | Reproduced evidence |
|---|---|
| Global affine minimax radius and alternating extrema | Live `Solve` returns `{u -> -1, v -> 1/8, r -> 1/8}`. Independent quantified `Reduce` returns exactly `r >= 1/8` (`CompensatedQuadratic.wl:18–24`). |
| Eightfold curvature-error reduction | Independent exact `MaxValue` checks give radii 1 and 1/8 on the unit interval; their ratio is 8. Scaled uniform bound holds for every real `h>0` and `0<=t<=h`. |
| Scaled compensation, signed Horner error, limit, normalization and phase modulus | All 11 named derivation checks reproduce as True, including the polynomial identity, uniform bound, zero-change limit, `e0-t e1`, relative-error normalization, unit phase modulus, ties-to-even, and malformed-input checks (`CompensatedQuadratic.wl:23–34`). |
| Coefficients and precision | WL and independently written exact Python Fraction computations agree for every one of the 192 rows. Fractional widths are `(30,23,10)` and physical widths `(30,22,9)`. |
| Exact capacity and row-alignment error bounds | All 192 capacity predicates reproduce; a separate Fraction implementation also checks all row bounds and all 192 inner-error values. The outer bound is exactly `98305/8796093022208` (`CompensatedQuadratic.wl:36–54`). |
| Continuous and binary-grid coverage | Exactly 191 continuous envelopes are True; interval 0 remains explicitly unproved continuously. WL and a separately written Fraction/integer checker both pass all 16,384 interval-zero offsets, with positive accumulators. The independent cleared-integer squared inequalities have minimum slack 21,536,728,635,567,136. Together the certificates cover the full normalized 20-bit grid in `[1,4)` (`grid.wls:9–30`). |
| Exhaustive legal integer inputs | All 3,133,452 inputs pass; zero phase failures and zero partial-accumulator violations. Minimum uint64 phase margin is 159,224,414. Worst error/tolerance is 0.8505589342620424 at `x=1064928` (`screen.py:40–56`). |
| Four compiled schedules | All schedules below pass exact-model phase matching, input restoration, scratch cleanup, register fit and priced-count reconciliation on 4,095 distinct boundary inputs repeated four times, 16,380 shots each, seed 371811 (`screen.py:59–81`). |
| Selected standalone default | Both `generator.py` and `screen.py` produce the selected retained-curvature schedule and identical canonical binary. Its SHA-256 is `21828cd5ccc1238e688a5b06bdda36316daf65ee8d6b976cecafe8b4bd45c007` (`generator.py:60`, `screen.py:100–104`). |
| Full official verification | Unmodified pinned verifier SHA-256 `ddef520dec6b68d1cf0a6f328ddd36d79eb84570b9192cf3f7acaf073f06733a` passes **200,540,928 shots**, with the complete default specification and salt `qf-prototype-local-validation-v1`. No sampling substitute was used (`verify.py:90–114`). |
| Selected resources | 151 qubits, 3,097 Toffoli-like gates, 1 T-like gate, 37 fine rotations, reaction depth 1,393, 34,426 operations, synthesis epsilon `7.328888008135618e-13`, Toffoli equivalent `3355.0884193381044`. |
| Cost and retained-register tradeoff | Official cost is `506618.35132005374`, versus incumbent `527287.988387595`, a `3.9199901235655754%` reduction. Recycling adds 279 Toffolis and saves nine qubits; its exact-model-margin cost is `516040.55554601084`. |
| Combined comparison arithmetic | Independently recomputed ratios are `4.099287418831162`, `1.2581323878369532`, `1.6360612511835888`. Their geometric mean is `2.035845221582375`, versus previous `2.0088882102784944`. DQI/DFTHC reports were reused and remained unchanged (`verify.py:57–82`). |
| Reproducible reports and provenance | Regenerated `coefficients.json`, `grid-certificate.json`, `build.json`, `screening.json`, and binary are byte-identical. Every full rerun report field equals stored `verified.json` except `elapsed_seconds`. All bound source/comparison hashes match. The candidate binary is ignored by Git as documented. |

| Arithmetic | Recycle curvature | Adder budget | Qubits | Toffolis | Rotations | Cost from exact-model margin |
|---|---|---:|---:|---:|---:|---:|
| padded | false | default | 151 | 3097 | 37 | 506618.35132005374 |
| padded | true | default | 142 | 3376 | 37 | 516040.55554601084 |
| padded | true | 130 | 130 | 3812 | 37 | 529111.4945139536 |
| controlled | true | default | 115 | 4944 | 37 | 598240.168223882 |

Function signatures were checked against the complete implementation and live execution: `CompensatedRows[rows,widths]`, `CompensationDerivation[]`, `CoefficientCertificate[integerRows,widths]`, Python `quantized_rows(widths=(30,23,10))`, and keyword-only `build_circuit(arithmetic='padded', recycle_coefficient=False, adder_budget=None)`. Live `Options[Solve]`, `Options[Reduce]` and `Options[Resolve]` were inspected; the exact domain/quantifier forms used in the package execute successfully.

Additional regime checks passed: empty, short, nonpositive, approximate and symbolic WL tables; zero, excessive, short, approximate and symbolic widths; empty/nonpositive certificate inputs; widths 1 and 36 preserving the `(192,3)` return shape; and width 36 exactly restoring the absolute original coefficients. Five malformed Python width cases and four malformed model-input cases are rejected. The empty model input returns empty arrays.

The full verifier uses width 22, `max_x=3133452`, `beta=1`, precision 24, 64 shots per input, batch size 256, four workers, and the pinned salt. Its implementation checks restoration and scratch cleanliness on every shot and reserves the remaining phase margin for rotation synthesis. The fresh full run took **185.48657279100735 seconds** inside the runner. This is one observed duration, excluding Python startup; it is not a speed benchmark. No execution-speed claim was made, so repeated timing trials were inapplicable. No notebook deliverable exists. No novelty, global circuit optimality, remote submission, Hamiltonian simulation, or quality-pass conclusion is drawn.

Historical absence of prior remote actions cannot be established from these local artifacts; this review made no commit, branch, push, submission, or other remote action. The stated correctness result concerns the frozen local mathematical and circuit artifacts. Other task verification was intentionally not rerun.

Before/after SHA-256 bindings (each unchanged):

| Target | SHA-256 |
|---|---|
| `CompensatedQuadratic.wl` | `c00cb5caf066d09fa1f6f32ae138aa8a1d477574301cfbec5cd0a88a701b543e` |
| `run.wls` | `165d8bfb78005c0bf80a7a6e472cfb56fb9e909b1d78b77195c42a101cec2de7` |
| `grid.wls` | `8ed16cee7824137d763967984a20de7b0ac013dfdc0e4212ac0b3548f32c1a43` |
| `generator.py` | `35bf116960bc425f972fe734e78eea87738cce3ebb60741c1720f932eef9cacb` |
| `screen.py` | `50ef2c80bc96a83386a6ae965fb9cca721db83cf4481de2402e73d408b167d1c` |
| `verify.py` | `dc3283a50e4f34e0c736f5a2c077c4384b112d80d1b88db8734d0a6dda68b62f` |
| `README.md` | `8b88be9dd37cc2377d148d798a851a8472f4ec9e1a3b235b70556aa9f2162173` |
| `PLAN.md` | `8611c9f968b201f9a76819a4448ac1a6cf199141ac781a7cf8630be2811e91f8` |
| `coefficients.json` | `7f033cfb62fce13e12d40a63405a314e086bcdcd6ec67f7e213e23b431870be2` |
| `grid-certificate.json` | `5ce76b02813cf04ff3cb88bff88a885e3c70aa79ab75bed81ce52c17961cb19b` |
| `build.json` | `e9f7dbbdad8bcef71fe84e4b5c97a86bd39341659cc12843c3289473c60a0fc1` |
| `screening.json` | `0264fa11cb48bb745f5c9e66380a2269cee2722089149bc05e988a9b1cd9fa45` |
| `verified.json` | `b3f7d53e053d992f647572461922dc5514310e315ebe005c5658b95db2b8efb8` |
| `local/compensated-quadratic.kmb` | `21828cd5ccc1238e688a5b06bdda36316daf65ee8d6b976cecafe8b4bd45c007` |

# FINDINGS

No unresolved correctness or reproducibility findings in the declared scope. The continuous-envelope gap for interval 0 is documented accurately and closed on the full required binary grid. The three nonselected schedules have bounded screening evidence; only the selected circuit has the full new official verification result, as documented.

# RECOMMENDED FIXES

None required for correctness.

OPEN ISSUES: 0
