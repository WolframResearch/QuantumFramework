# Independent WL quantum quality review, round 1

Scope: the thirteen supplied compensated-quadratic artifacts and their named direct dependencies. Reproduction used an isolated temporary copy; original artifacts and the candidate binary were not edited. No prior reviews, research-search directories, thread history, or other reviewers were consulted. The scoped claim is a promised-domain diagonal inverse-square-root phase oracle with a fixed 192-row realization. No Hamiltonian or many-body claim is assessed.

Runtime: Wolfram Language 15.0.0, Python 3.13.1, NumPy 2.5.4, kickmix 0.1. The Wolfram kernel ran without initialization files. No messages were suppressed. Only the selected circuit received a new full official run; the four schedules received the specified compiled boundary probes.

Reproduction commands, executed in the isolated artifact directory:

```text
WolframKernel -noinit -noprompt -script run.wls
WolframKernel -noinit -noprompt -script grid.wls
python screen.py --data DATA
python verify.py --data DATA --workers 4
```

## A. Documented symbols — PASS

The WL package uses core System symbols and its own three exported functions; it has no QuantumFramework or other external WL dependency. Public argument/return contracts are given in `CompensatedQuadratic.wl:2–4`, definitions and failure fallbacks at lines 10–15, 17–34, and 36–55. Executed calls returned the stated integer table, derivation Association, certificate Association, and Failure objects.

The installed kernel's `Information[symbol,"Usage"]` was inspected for Solve, Reduce, Resolve, ForAll, Exists, Round, IntegerDigits, and FileHash. It confirms the explicit-domain solving and quantifier forms actually used in `CompensatedQuadratic.wl:18–19,26,48–50`, and the fixed-length digit form in `grid.wls:10`. The exact calls execute with the documented return shapes: Solve returns rules, Reduce returns `r >= 1/8`, Resolve returns booleans. There are no nondefault solver options to validate.

For the Python boundary, the installed kickmix 0.1 docstrings were inspected for the builder, allocation, lookup/unlookup, and/un-and, addition/subtraction, rotations, Z_POW, simulator I/O/phase, serialization, statistics, and reaction depth methods. They support the used signatures, little-endian registers, Fraction exponents, phase in half-turns, clean lookup targets and cleared unlookup results. `generator.py:99–100` includes the required factor of two between turns and Z_POW half-turn exponents. The checked official reference explicitly uses the same little-endian x register and half-turn conversion. No recalled API or unsupported return-shape defense is needed.

Smallest fix: none.

## B. Idiomatic functional WL — PASS

`CompensatedQuadratic.wl:10–14` maps the row transformation and uses dependent lexical constants with With. Lines 38–54 use column reduction, MapIndexed, list arithmetic, and AllTrue. `grid.wls:9–20` uses digit matrices, Dot, listable Floor, and mapping over shifts to construct the needed two-dimensional product array. The mapping on line 13 is across exponents; it is not an unnecessary map of Floor over scalar elements. Both Modules provide fresh symbolic variables for quantifiers/algebra, rather than mutable loop counters or constant-only local bindings.

The complete three WL files contain no Do, For, While, AppendTo, mutable indexed table construction, unnecessary dense quantum operator, or replacement of a simpler built-in by an imperative implementation. The dense coefficient/digit arrays are appropriate to their actual nonzero structure. Runner globals are script outputs and provenance records, outside a reusable numerical inner loop.

Smallest fix: none.

## C. Terse physics reasoning — PASS

The argument proceeds from exact affine minimax compensation (`CompensatedQuadratic.wl:18–26`) to signed dependent Horner error (line 28), the normalization identity (29–30), phase modulus (31), then actual register and phase envelopes (36–54). Named intermediates `coarse`, `delta`, `innerError`, `outerError`, `origin`, and `p` expose the mathematical roles. The README's derivation (`README.md:14–37`) explains the representation and separates the continuous envelope from interval 0's finite-grid closure. Its cost account (`README.md:41–64`) explicitly includes the memory tradeoff and the unchanged other-task reports.

An independent kernel calculation gives residual extrema `{h^2/8,-h^2/8,h^2/8}` at `{0,h/2,h}` and the critical point `t -> h/2`; C-only perturbation has maximum factor `h^2` for h>0. Thus the eightfold statement concerns curvature perturbation, as claimed, rather than the entire rounded circuit error. No dead definition or repeated physical derivation requires contraction.

Smallest fix: none.

## D. Physics invariants — PASS

The following independent kernel probes returned the displayed exact outputs:

```wolfram
Reduce[Exists[{u,v},ForAll[t,0<=t<=1,-r<=t^2+u t+v<=r]],r,Reals]
(* r >= 1/8 *)
Resolve[ForAll[{t,h,d},h>0 && Element[d,Reals] && 0<=t<=h,
  Abs[d (t^2-h t+h^2/8)]<=Abs[d] h^2/8],Reals]
(* True *)
Expand[(a-t (b-t c+e1)+e0)-(a-t (b-t c))]
(* e0 - e1 t *)
FullSimplify[y 2^(s-10) Sqrt[x]-y Sqrt[x 4^s/2^20],
  x>0 && Element[{s,y},Reals]]
(* 0 *)
FullSimplify[Abs[Exp[2 Pi I phi]]^2,Element[phi,Reals]]
(* 1 *)
```

All eleven shipped derivation checks return True. The fresh coefficient run proves all 192 capacity checks with physical widths `{30,22,9}`. The 191 true continuous certificates and the false certificate at index 0 reproduce exactly. The grid runner checks all 16,384 offsets of that remaining interval with positive accumulators and exact squared relative-error inequalities.

An additional independent scalar integer realization checks every one of the 3,145,728 normalized points in [1,4), using Python arbitrary-size integers. With `n=(64+i) 2^14+offset`, its stages are `b-sum(bit_m floor(c 2^(13-m)))` followed by `a-sum(bit_m floor(inner 2^(7-m)))`, for m=7..20. It checks

```text
(2^24-1)^2 2^32 <= n outer^2 <= (2^24+1)^2 2^32
```

and positive accumulators at every point. Output: zero violations; minimum inner accumulator 525831; minimum outer accumulator 536870975. Independently recomputed all 192 inner-error fractions and the outer bound `98305/8796093022208` agree with the saved certificate.

The exhaustive promised-input integer model returns zero phase failures and zero partial-accumulator violations, with minimum remaining margin 159224414 units of 2^-64 turns. Every compiled schedule restores input, clears scratch, matches the exact integer phase, and reconciles priced gate counts on 4,095 distinct boundary inputs / 16,380 seeded shots. These checks are asserted in `screen.py:20–81`, not merely printed.

The unmodified default official verifier was independently rerun through `verify.py --data DATA --workers 4` (`verify.py:90–120`): **PASS, 200,540,928 shots**, 64 seeded repetitions for each of 3,133,452 legal inputs. It reproduces cost **506618.35132005374**, Toffoli equivalent **3355.0884193381044**, 151 qubits, 3097 Toffolis and 37 fine rotations; worst error/tolerance **0.8505589342620424** at **x=1064928**. Input restoration and clean scratch are asserted for every simulated shot by the pinned official verifier. Every field of the reproduced `verified.json` equals the supplied report except `elapsed_seconds`. The recorded duration of this run is an observation, not a speed comparison.

Unit modulus gives the diagonal promised-domain oracle's norm preservation. Positive real inverse-square-root values, capacity, basis order and phase conventions are the relevant invariants. Hermiticity, energy conservation, density-operator positivity, and Hamiltonian spectra are inapplicable because no such observables or dynamics are claimed. The official run exercises its documented seeded measurement-uncomputation contract; it is not an enumeration of every conceivable measurement branch.

Smallest fix: none.

## E. Nontrivial examples — PASS

The evidence does not end at zero compensation, endpoints, or one Horner stage. Independently recomputed curvature changes are nonzero in every one of the 192 rows: 97 positive and 95 negative. Across the complete normalized grid, 3,054,976 inputs have nonzero inner-stage rounding, 3,142,886 have nonzero outer-stage rounding, and 3,053,774 have both. These counts were computed from exact Fraction differences against the two unrounded stage expressions; the signed error is therefore exercised in its dependent two-stage regime.

`screen.py:59–65` probes both sides of every reachable table/normalization boundary and both promised endpoints. The exhaustive model covers the full legal x range, while the default full official protocol repeats every legal x 64 times. All four bounded schedules in `screen.py:14–17` are represented; only the selected schedule claims new official verification.

Smallest fix: none.

## F. Symbolic-first derivation — PASS

The minimax Solve/Reduce calculation retains symbolic u,v,r,t. The uniform bound retains arbitrary positive h and real perturbation; signed Horner error retains symbolic a,b,c,e0,e1; normalization retains symbolic x,y,s; phase modulus retains real phi (`CompensatedQuadratic.wl:17–31`). The independent exact extrema and C-only maximum confirm the scaling rather than assuming it from a numerical fit. All coefficient quantization uses exact powers of two and Round; Python independently uses Fraction (`generator.py:31–41`). The continuous certificate and grid closure use exact integer/rational arithmetic. Numeric eps and cost are computed only after an exact error-margin calculation, at the required rotation-synthesis cost boundary. There is no premature N, NIntegrate, NDSolve, numerical diagonalization, or replacement of a derivation by fitted constants.

Smallest fix: none.

## G. Regimes — PASS

The explicit regimes table at `README.md:69–75` is supported by executable evidence:

| Required regime | Reproduced evidence |
|---|---|
| General symbolic | Quantified global radius r>=1/8; exact scaled absolute bound for h>0; general signed Horner, normalization and real-phase modulus identities |
| Exactly solvable | Alternating extrema at 0,h/2,h; critical point h/2; exact ties-to-even including negative and positive half-integers |
| Limit | `Limit[d (t^2-h t+h^2/8),d->0]` returns 0 |
| Independent reference | All 192 exact coefficient rows and per-row error fractions; 191 continuous certificates plus interval-0 closure; independent full normalized-grid integer check; exhaustive promised inputs and compiled schedules; full official selected-circuit run |
| Failure / boundary | WL rejects malformed table dimensions/signs/types and widths of wrong length, zero, >36, machine real or symbol; certificate rejects empty table and zero width; Python rejects zero/>36/float/bool/wrong-length widths, unknown schedule, impossible qubit budget, zero/above-max/rank-two input; compiled boundary probes cover legal endpoint and interval/normalization transitions |

The fresh WL invalid-width/table/certificate vectors are all True. The ten independent Python invalid-input calls all raise ValueError. The conservative index-0 continuous rejection remains visible and is closed only on the claimed finite input grid; no continuous proof is inferred there. General symbolic representation parameters are retained until the fixed released realization, as required; an arbitrary-size Hamiltonian example would not test this task's scoped physics.

Smallest fix: none.

## Numeric claims and provenance

All numeric/table rows in the thirteen target artifacts were parsed. All coefficient rows, interval fractions/booleans, grid rows, and four schedule resource/cost rows reproduce. The README resource table agrees with the pinned previous and selected reports: 142->151 qubits, 3460->3097 Toffolis, 37 fine rotations in both cases, 3713.295692870388->3355.0884193381044 Toffoli equivalents, and 527287.988387595->506618.35132005374 cost. The same-table recycling schedule adds 279 Toffolis and has estimated cost 516040.55554601084.

Recomputing baseline/current ratios gives `4.099287418831162`, `1.2581323878369532`, `1.6360612511835888`; their geometric mean is `2.035845221582375`, versus `2.0088882102784944` using the previous quadratic. Cost reduction versus that quadratic is `3.9199901235655754%`. Other-task reports were read and hash-bound as comparison inputs; their circuits were not rerun or modified. The combined value remains a local estimate, with no submission claim. No wall-clock speedup claim is made.

All thirteen original target SHA256 values and the candidate hash were bound before/after. All original target/candidate bytes remained unchanged; each row below has identical before and after SHA256. The isolated rerun also reproduces byte-for-byte `coefficients.json`, `grid-certificate.json`, `build.json`, `screening.json`, and the candidate binary. Only the isolated official report's elapsed duration changes.

| Artifact (relative to this deliverable) | Before = after SHA256 |
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


OPEN ISSUES: 0
