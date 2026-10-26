# Correctness review — round 1/3

## VERDICT

PASS. A fresh independent verifier reproduced the supplied implementation in an isolated copy. No correctness issue remains. The full unmodified official run completed with all **200,540,928 shots**; the resulting report equals the supplied `verified.json` in every field except nondeterministic elapsed time.

Read the complete target sources and claims, parsed every coefficient/certificate row, inspected the named arithmetic and official-verifier dependencies, and followed ancestor instructions. No prior review or `representation_search/` material was opened or searched. No subagents were used. All 32 snapshotted original target/dependency files remained byte-identical; report creation is the only repository write by this reviewer.

## REPRODUCED

Commands executed in the isolated `project/compensated_quadratic` copy, with the supplied Python environment, Wolfram runtime and official `DATA` directory:

```text
WolframKernel -noinit -noprompt -script run.wls
WolframKernel -noinit -noprompt -script grid.wls
python screen.py --data DATA
python verify.py --data DATA --workers 4
WolframKernel -noinit -noprompt -script independent.wls
python independent.py
```

All commands exited 0. No Wolfram messages, Python tracebacks, suppressed-message wrappers, missing dependencies or approval barriers occurred. WL outputs were captured through `Echo`; no WL `Print` or `Quiet` was added.

Runtime: Wolfram **15.0.0 for Mac OS X ARM (64-bit), February 17, 2026**; Python **3.13.1**, NumPy **2.5.4**, kickmix **0.1**. The full verifier used its unchanged default batch size 256, four workers, the supplied salt and the complete default specification `(width,max_x,precision_bits,shots_per_x)=(22,3133452,24,64)`.

Claim and output inventory:

| Claim / artifact lines | Independent result |
|---|---|
| Exact affine minimax derivation, `CompensatedQuadratic.wl:17–34`, `README.md:14–24` | Live `Solve` returned `{{u -> -1, v -> 1/8, r -> 1/8}}`; quantified `Reduce` returned `r >= 1/8`. Independent extrema were `{1/8,-1/8,1/8}`; a radius `1/9` was infeasible. Exact maximum-ratio calculation returned 8. |
| Scaled compensation identity, zero-change limit, signed Horner error, normalization, phase modulus and rounding, `CompensatedQuadratic.wl:25–34` | All 11 saved derivation checks reproduced as `True`, including `e0-t e1`, normalization with real `y,s` and positive `x`, unit phase modulus and ties-to-even. |
| 192-row coefficient correspondence and widths, `generator.py:31–58`, `README.md:25–27` | All rows matched exact Python `Fraction` derivation. Shape `{192,3}`, physical widths `{30,22,9}`; frozen source and source bindings matched. |
| Domain/failure signatures, `CompensatedQuadratic.wl:2–15,36–55` | Invalid/empty/wrong-shape, wrong-sign, real and symbolic tables were rejected; invalid width lengths, 0, 37, real and symbolic widths were rejected. Width endpoints 1 and 36 returned 192-by-3 tables. Invalid certificate inputs returned `Failure`. |
| Capacity and continuous envelope, `CompensatedQuadratic.wl:36–54`, `README.md:29–37` | Every capacity check passed. Outer error bound exactly `98305/8796093022208`. Exactly 191 continuous envelopes were `True`, interval 0 was `False`, and none timed out. Rebuilt `coefficients.json` was byte-identical. |
| Exact grid closure, `grid.wls:9–30` | Interval 0 passed all 16,384 offsets, both accumulators positive, exact squared relative bound satisfied. All coverage checks were true; rebuilt `grid-certificate.json` was byte-identical. A separate integer implementation additionally checked **all 3,145,728 normalized grid points** in `[1,4)`, with zero phase failures and zero partial-register violations. |
| Full integer reference, `screen.py:20–56`, `README.md:59–60,74` | All **3,133,452** legal integer inputs passed; zero phase/capacity failures; minimum remaining margin **159224414** uint64 phase units; worst ratio **0.8505589342620424** at **1064928**. Out-of-domain and non-1D model inputs were rejected; empty model input returned two empty vectors. |
| Four compiled schedule probes, `screen.py:59–109`, `README.md:50–57,75–81` | Each schedule checked **4,095 unique inputs / 16,380 shots**, seed 371811. Input restoration, scratch cleanup, exact integer-phase equality, register fit and priced static/measured counts all passed. Rebuilt `screening.json` was byte-identical. |
| Selected resources, `build.json`, `README.md:41–57` | Padded retained-C2 schedule: **151 qubits, 3097 Toffolis, 1 T, 37 fine rotations, reaction depth 1393, 34426 operations**. Recycling: 142 qubits / 3376 Toffolis / cost 516040.55554601084. The 279-Toffoli and nine-qubit schedule differences reproduced. Smaller-budget and controlled alternatives also reproduced their saved counts and costs. `build.json` and candidate binary were byte-identical. |
| Full official circuit verification, `verify.py:84–123`, `README.md:3–8,36–37,59–63,77–81` | Unmodified full verifier returned `PASS`, **200540928 shots**, **506618.35132005374** cost, **3355.0884193381044** Toffoli equivalent and rotation epsilon **7.328888008135618e-13**. Worst ratio/input matched the integer model. Restoration and scratch cleanup were checked on every shot. |
| Provenance/stability, `generator.py:45–58,128–138`, `verify.py:14–53,89–113`, `README.md:93–98` | Frozen source, WL derivation, Python generator, arithmetic dependencies, candidate, official verifier and full specification matched. Official and candidate hashes matched after execution; all recorded incumbent input hashes matched current files. |
| Cost/combined-score arithmetic, `verify.py:56–81`, `README.md:5–8,61–65` | Baseline/current ratios **4.099287418831162**, **1.2581323878369532**, **1.6360612511835888**; geometric mean **2.035845221582375**, previous **2.0088882102784944**. Incumbent cost reduction **3.9199901235655754%**, combined improvement **1.341887077934678%**. All saved comparison fields reproduced independently. |
| Delivery constraints and plan, `PLAN.md:8–19`, `README.md:9–10,64–65,82,98,105` | Four bounded schedules only; no additional precision/architecture search, notebook target, Hamiltonian claim or external submission was involved. Incumbent, DQI and DFTHC artifacts were unchanged. This review does not assert a quality-review outcome. |

`Solve[expr,vars,Reals]`, quantified `Reduce[expr,vars,Reals]` and `Resolve[expr,Reals]` signatures and return forms were checked against live kernel results and current primary documentation: [Solve](https://reference.wolfram.com/language/ref/Solve.html), [Reduce](https://reference.wolfram.com/language/ref/Reduce.html), [Resolve](https://reference.wolfram.com/language/ref/Resolve.html).

The independent arithmetic implementation reconstructed both truncated Horner stages using integer shifts, checked partial accumulators, and tested the exact integer inequality `((2^24-1)^2) 2^32 <= M_numerator outer^2 <= ((2^24+1)^2) 2^32`; it did not call the author's integer model for this grid proof. Temporary independent script SHA-256 values: WL `eadf5db9c8e4262298ba6ce8b1b84874f259bee9947fc461064f98ea37a377c4`; Python `38fb31c929c4215b5c23c137c7726efd91f208429f19515fc66a9615cb767e5c`.

The complete official rerun took **184.327376250003 seconds**. This is an observed run duration, not a speed benchmark; deterministic circuit resources and official cost, rather than wall-clock improvement, are the claims assessed here.

Scope limits: DQI, DFTHC and prior inverse-square-root reports were consumed as existing evidence and their relevant hash bindings checked; those prior benchmarks were not rerun, as requested. The combined estimate is local arithmetic, not a Kaggle score or rank. The constant 8.43991 was treated as the supplied comparison reference. No independent external-account audit of historical push/submission flags was attempted; this review performed neither action. No dependency or coverage limitation reduced the current candidate's full verification.

## FINDINGS

None. In particular, the failed conservative continuous envelope for interval 0 is explicitly reported and completely covered by the exact finite-grid proof; it is not presented as a continuous certificate.

## RECOMMENDED FIXES

None required for correctness. Preserve the current scope distinctions between continuous/grid proof, schedule screening, full official verification and local combined arithmetic.

## Hash binding

Paths below are relative to the target directory. Target hashes were checked before and after reproduction.

| Target | SHA-256 |
|---|---|
| `CompensatedQuadratic.wl` | `c00cb5caf066d09fa1f6f32ae138aa8a1d477574301cfbec5cd0a88a701b543e` |
| `run.wls` | `165d8bfb78005c0bf80a7a6e472cfb56fb9e909b1d78b77195c42a101cec2de7` |
| `grid.wls` | `8ed16cee7824137d763967984a20de7b0ac013dfdc0e4212ac0b3548f32c1a43` |
| `generator.py` | `9f53e80624fec06474fafae9482d44fc38334d41f6b46e92736a754890bb9cc2` |
| `screen.py` | `0091f3b2f9ec491726d9a9ca45676f0490100537eeb8ad3a2af14b07e49773d4` |
| `verify.py` | `ca51dc81ac99981e492512625f8a9fe23ab8a47d1d24380c6670643507aa2bde` |
| `PLAN.md` | `8611c9f968b201f9a76819a4448ac1a6cf199141ac781a7cf8630be2811e91f8` |
| `README.md` | `57f65df4e4726f78ad48fdcd261bea2e40a5da93333aff01646144bdaa66ecdd` |
| `coefficients.json` | `7f033cfb62fce13e12d40a63405a314e086bcdcd6ec67f7e213e23b431870be2` |
| `grid-certificate.json` | `5ce76b02813cf04ff3cb88bff88a885e3c70aa79ab75bed81ce52c17961cb19b` |
| `build.json` | `45b27e0f7892226481d4dd3b2ea1fac2f838817403bcc8c6fb0ff3893e52130c` |
| `screening.json` | `b383d6cbef40567dbb41b0fdd8ff03f6dfecddfef227082476411d272d4710eb` |
| `verified.json` | `62b8589ad016b44d880388a0f82458be0d94ea943eee585811e7d174f048d470` |

| Direct / hash-record dependency | SHA-256 |
|---|---|
| `../quadratic.py` | `62cfa195356f83c7e9cca17898a74b78185fe42e0a1520d0cac8ab1153c55ce2` |
| `../optimized_quadratic.py` | `9ca63c842c066d957aa8b7d175b93d22eb40d7dadb0ca9fdcbab46cf8a954f36` |
| `../verify_quadratic.py` | `a0ddece35aa17687614043c0746d198ff22bd4889302ee8e3ad567e197a7d933` |
| `../feasibility/B-quadratic-192.json` | `2b768e13c7054d0e58bc4edf7afc0ff64559685aa309ad250b275fd9fd6a4c8d` |
| `../results/optimized-quadratic-verified.json` | `8576ee713ef79db0ffcd628af1b67617b6127adbd21b7719bf76061ac41f4109` |
| `../results/baseline-verified.json` | `2eb7fd500e2cd3efaff5e97b02ee61027f36a37661924b7108ad727634d2ccbb` |
| `../results/candidate-verified.json` | `ff31729451dc5b1c575cdae8853201e82f64546233feef27492c7d03f23baabd` |
| `../dqi_gf8/verified.json` | `3bbe0bac3b0dde8d1b591afffd8e8bb04085047609e53ad8b6b8e1d9de97d5fc` |
| `../dfthc_c30/verified.json` | `8412005457b6ee96e17460309d3906b15d3da7843022f18a33a0ae6c787a2ef5` |
| `../other_tasks/dqi-baseline.json` | `86f3134485e054806407e51b1cd34fadd89ece5a09871d4fe370ea7f750915bc` |
| `../other_tasks/dfthcsos-baseline.json` | `c714c2c36862305e7aedb420def79361ac3bb767d924f148bcc83f14161205fc` |
| `../verify_optimized_quadratic.py` | `12282432b2bacb5f75844d2d119b342907185ae46f5808cf433ee87c5ad3859f` |
| `../optimization/WidthOptimization.wl` | `971d23be409a4b85a804dd7f65e85dc9f40144e7f77961c473616622b0c2ad3a` |
| `../optimization/run-widths.wls` | `ce7fa1c74f9081df97b9de8dea20a14230a93f919d57468c7f0735fff4fca6e6` |
| `../optimization/verify-widths.wls` | `5560d4b91368e1d8fdc32b08be98c596d995d9e8c8b644356691ee16b36327e7` |
| `../optimization/width-certificates.json` | `ba544c74e32b49ccf31653c180b0827f4f4cb6b5eaec194a87eba19e60dfaa2d` |
| `../optimization/wl-checks.json` | `73a38fca44eb3e99628d65ee4bc618ba246b0e1bd0079412c5c48d2540b53b70` |
| `../optimization/analyze.py` | `ef3f90e95213fdf1ff3fae96b8193651fb66587f9091a9a8f108e4f1e9b86538` |
| `../optimization/screening.json` | `f0a80fe803c41bb1cd5d919f27f02fd1bea77cb57848d11fe593e235399a1ce9` |

Official verifier SHA-256: `ddef520dec6b68d1cf0a6f328ddd36d79eb84570b9192cf3f7acaf073f06733a`.
Candidate canonical KMB SHA-256: `21828cd5ccc1238e688a5b06bdda36316daf65ee8d6b976cecafe8b4bd45c007`.
Supplied official baseline KMB SHA-256: `e2bf8eb64350a4460a271fc7f1d591ef125b1e276377e9a0f64b96ab3ca321e1`.

OPEN ISSUES: 0
