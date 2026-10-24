# Independent WL quality and physics-depth review — round 1/3

Verdict: A–G PASS. No unresolved or unverifiable quality finding within the supplied narrow arithmetic-adapter scope.

## Scope and independence

I read the complete wl-quality skill and critic mandate, wl-native-style and wl-style-guide rubrics, and all requested artifacts. I found no applicable AGENTS.md in the repository or its ancestor/project directories. No prior reviewer report was read. I did not modify repository artifacts, official inputs, native kickmix sources, the historical baseline, or prior feasibility evidence.

Artifact root: `OngoingProjects/VerifiableCircuitCompilation/dqi_gf8`. Complete files read: `__init__.py`, `gf8_backend.py`, `generator.py`, `verify.py`, `test_integration.py`, `run-certificates.wls`, `certificates.json`, `build.json`, `verified.json`, and `README.md`. Required sibling files read completely: `other_tasks/FieldArithmetic.wl`, `symbolic.json`, `field_probes.py`, and `dqi-baseline.json`. The unrelated mixing values embedded in symbolic.json were read as dependency bytes, not treated as new claims to certify here.

I also inspected installed WL 15.0 FiniteField/FiniteFieldElement documentation extracted in `FiniteField.md` and `FiniteFieldElement.md`; native kickmix arithmetic implementation and Python binding source; and the external official decoder/verifier source used by the relevant protocols. Native source HEAD was independently checked as `638777fac4379bc62e2a1e4a52cf1a1753895d48`.

Claimed mathematical object: an exact bilinear GF(4096) product in the fixed native polynomial basis, implemented by seven GF(8) products of bilinear rank six; linear Frobenius maps; and the native-compatible inverse scaffold. Relevant quantum object: reversible classical XOR arithmetic and its intended diagonal phase, with measurement-based cleanup and restored input/control/classical registers. No Hamiltonian, density operator, propagation approximation, continuum limit, or general builder replacement is claimed.

## What actually executed

`$KERNEL`, `$PYTHON`, `$DATA` and `$MIRROR` identify the exercised direct kernel, pinned runtime, external competition bundle and temporary mirror. Personal absolute paths are omitted.

All generators/tests ran in a unique mirror, `$MIRROR`. Only read-only binding audits touched the original circuit and sources.

1. `"$KERNEL" -noinit -noprompt -script dqi_gf8/run-certificates.wls`, from that mirror: exit 0, status PASS, all eleven checks True. The regenerated certificates.json was byte-identical to the original. Checked names: same_tower_plan, same_six_product, chain_composition, native_scaffold_exponents, inverse_exponent, all_nonzero_inverses, zero_extension, frobenius_invertible, frobenius_period, frobenius_inverse, and frobenius_composition.
2. `QF_INVSQRT_DATA="$DATA" "$PYTHON" -m unittest dqi_gf8.test_integration -v`: all six tests passed, `Ran 6 tests in 0.937s`, exit 0. This executed both multiplication control modes, all inverse inputs, both division condition branches, all six Frobenius maps on all field words, malformed inputs, and the unmodified reduced official decoder protocol with 2,048 shots and one worker.
3. Fresh independent WL probe script `quality-probes.wls` in the mirror: exit 0, eight exact checks True. It confirmed native-polynomial irreducibility, polynomial index/coefficients round-trip for all 4,096 words, all six Frobenius matrices by independent PolynomialRemainder rather than FiniteField exponentiation, the six-product identity, exact native scaffold exponents, the arbitrary positive-integer chain identity, zero/one boundary values, and invalid subfield-degree Failure. Kernel reported `15.0.0 for Mac OS X ARM (64-bit) (February 17, 2026)`, polynomial representation, degree 12, and `Options[FiniteField] == {}`. Evidence is saved in `quality-probes.json`.
4. Fresh compiled Python probe `scaffold-probe.py`: for every native input 0…4095, directly read all four generated scaffold registers and the inverse output. Eight assertions passed: source preserved; scaffold powers 3, 15, 1020, and 2046; inverse output; no extra phase; zero temporary scratch. Evidence is saved in `scaffold-probe.json`.
5. Read-only report/binding audit: twelve checks passed, saved in `report-audit.json`: build source hashes; verification source hashes; historical baseline hash; official input hashes; original binary SHA-256 and byte size; official cost formula; comparison ratio; static Toffoli counts; multiplication call count; total operation count; and mean Toffoli count reconstructed from executed shot operation totals.

I did not regenerate the 1.18 GB full circuit or rerun the expensive 102,400-shot full simulation during this quality review. The full candidate and historical baseline results below are report audits, not new full-run observations. This limit does not leave a requested invariant unexecuted: arithmetic invariants were rerun directly and the reduced official decoder was executed; the task expressly permits auditing the retained full report.

## A. Documented symbols — PASS

Evidence:

- FieldArithmetic.wl:21–31 and run-certificates.wls:7–11 use the documented `FiniteField[2, polynomialFunction]`, integer field indexing, and `FiniteFieldElement["Coefficients"]`/`["Index"]` accessors. Installed FiniteField documentation states polynomial indexing is reverse lexicographic in coefficients, and the FiniteFieldElement page explicitly documents these two properties and field-valued matrix `Inverse`. That supports the little-endian native polynomial basis and the interpolation inverse at FieldArithmetic.wl:36.
- The direct WL probe confirmed the representation and degree in the actual 15.0 kernel, and all 4,096 polynomial index round-trips returned True. Independent polynomial-remainder matrices agreed exactly with all consumed Frobenius matrices.
- The only package symbols in the WL certificate are ``OtherTasksField`TowerPlan[3]`` and ``OtherTasksField`SixProductProof[]``. Their full definitions, usages, accepted-degree guard, return Associations, and failure case are visible in FieldArithmetic.wl:3–4, 11–18, and 20–58. The caller consumes exactly those shapes.
- Native `gen_gf_inverse.cc:108–193` and Python binding `gf_arithmetic.pybind.cc:907–980` establish the ordered inverse scaffold contract. The generated four words were independently read and matched `{3,15,1020,2046}` on all 4,096 inputs. Native division binding at 1000 onward and `gen_gf_div.cc` establish the inverse-word-plus-chain layout returned by gf8_backend.py:169–182.

No unsupported signature, guessed return shape, stale package API, or handwritten substitute for a materially simpler documented field operation was found. The custom selected-subfield/native-basis interpolation maps are part of this implementation's certificate, rather than an attempt to duplicate generic finite-field arithmetic.

Defects: none. Required fix: none.

## B. Idiomatic functional WL — PASS

Evidence:

- A literal scan of both WL target files found no `Do[`, `For[`, `While[`, `AppendTo[`, `Quiet[`, `NIntegrate[`, or `NDSolve[`.
- FieldArithmetic.wl:11–18 uses With for immutable formula constants. TowerPlan:20–57 uses Module for multiple dependent exact objects and local helper functions whose definitions capture that invocation's basis; this is justified lexical isolation.
- FieldArithmetic.wl:26–46 uses SelectFirst, listable field powers, Catenate, Partition, Map, MapThread, NestList, and explicit rectangular Table constructions. The Table indices construct basis/coefficient tensors or select interpolation tensor entries; they do not emulate a loop appending to lists.
- run-certificates.wls:9–29 is a compact map/AssociationThread/exact-reduction certificate. Its small dense 12×12 binary maps and at most 11×11 interpolation matrices do not present a significant sparse linear-algebra problem.
- The two one-level `If[p === "infinity", ...]` cases encode the finite/infinite interpolation-point definition. They are not finite-choice control cascades or validation logic hidden in a numerical inner loop.

Python builder loops are circuit emission and binary row elimination, outside the WL functional-idiom mandate. They do not provide evidence of Python-shaped WL.

Defects: none. Required fix: none.

## C. Terse physics reasoning — PASS

Evidence:

- FieldArithmetic.wl:20–56 reads as field/subfield → native basis → coordinates/embedding → evaluation/interpolation → coefficient-tensor identity. Names describe the mathematical objects. SixProductProof:11–18 places the six products next to their exact output contraction and remainder identity.
- run-certificates.wls:12–16 gives `beta[k] = 2^k - 1`, then the symbolic composition identity, the exact scaffold powers, and the inverse power. These are the needed mathematical explanation for the compiler schedule.
- gf8_backend.py:96–123 keeps input evaluation, six-product emission, output accumulation, and native cleanup in one small computation. Lines 140–155 expose the native inverse schedule without duplicating native cleanup code. Lines 78–79 delegate other operations.
- README:43–68 states the adapter contract, exact native field, native-compatible scaffold, and retained phase-mask semantics. README:26–39 and 82–94 distinguish cost, Clifford/qubit tradeoffs, observed duration, and an actual submission score.

No dead derivation or repeated scaffolding was found in the supplied artifacts. The public prose does not pretend the specialized degree-12 adapter implements the entire builder API.

Defects: none. Required fix: none.

## D. Applicable quantum invariants — PASS

Evidence:

- The exact coefficient tensor in FieldArithmetic.wl:43–56 certifies all 12³ native multiplication coefficients, rather than only selected arithmetic examples. Basis and interpolation inverse identities are independently asserted at 52–55. The regenerated selected plan equals the stored plan exactly.
- run-certificates.wls:23–29 asserts all 4,095 nonzero inverse values, zero extension, binary Frobenius invertibility, period, inverse pairs, and composition. All were executed and returned True.
- test_integration.py:18–27, via field_probes.py:128–164, checks dirty XOR targets, preserved operands and quantum control, absence of unintended phase, and zero scratch across 8,192 shots per control mode. Nonlinear Toffoli counts are checked as 42 and 54.
- test_integration.py:29–54 checks all inverse inputs, preserved independently sampled classical masks, the diagonal inverse phase, and native cleanup's zero scratch. The new readout probe additionally checked every scaffold value directly.
- test_integration.py:56–90 checks every divisor on both classical-condition branches, restored input and mask/condition registers, quotient, product phase, and clean scratch. It asserts the native-compatible chain size 60 and static count 252.
- test_integration.py:92–104 checks all field words under all six generated Frobenius circuits against independent native field arithmetic. The reduced official protocol at 122–133 passed with zero value/phase errors.
- Native `gen_gf_mul.cc:246–279` explicitly measures computed product targets in the X basis and pays back their byproduct phase with the native product-phase circuit. The adapter retains this cleanup only after the exact corresponding product has been computed, at gf8_backend.py:120 and 137. Source compatibility plus exact tensor identity supports the sampled phase tests.

For an arbitrary superposition, the relevant specification is the corrected reversible map or intended diagonal phase, extended linearly. A Hermiticity, density positivity/trace, or Hamiltonian-energy assertion would be unrelated to this artifact. Measurement branches are sampled, as README:110 explicitly says, rather than falsely advertised as exhaustively enumerated. No statevector normalization benchmark or dense unitary construction is necessary to validate this finite-field arithmetic contract.

Defects: none. Required fix: none.

## E. Nontrivial examples — PASS

The base case this work could have hidden behind would be a clean-target GF(2) product, input 0/1, an unconditioned arithmetic fragment, or a phase-free test. It goes substantially beyond all of them:

- FieldArithmetic.wl:20–56 proves a native 12-bit bilinear tensor reconstructed through a quartic extension over GF(8), with seven subproducts, and every coefficient checked.
- All 4,096 inverse words and every divisor are exercised, including singular zero extension; all four nontrivial scaffold powers are now directly checked.
- Dirty XOR targets, quantum controls, independently sampled classical phase masks, disabled/enabled classical conditions, measurement cleanup, and a reduced official decoder are executed by the supplied six tests.
- The retained full default decoder report is bound to the actual circuit hash and current inputs/sources; it is not a count obtained only by extrapolating a trivial arithmetic fragment.

Defects: none. Required fix: none.

## F. Symbolic-first derivations — PASS

Evidence:

- SixProductProof uses six symbolic input coefficients, Expand, PolynomialRemainder, and PolynomialMod at FieldArithmetic.wl:11–18. No floating-point arithmetic substitutes for the GF(8) identity.
- TowerPlan uses exact finite-field elements, exact binary basis inversion and exact field-valued interpolation. Its complete bilinear coefficient-tensor identity is a symbolic proof for arbitrary field inputs in the fixed basis, even though a finite basis is used to represent it.
- run-certificates.wls:12–16 keeps `u,v` symbolic under positive-integer assumptions in the chain identity. The final stored exponents specialize only the declared native degree-12 schedule. The fresh FullSimplify probe reproduced the general identity as True.
- The finite degree, modulus, subfield, and inverse exponent are exact design parameters required for native compatibility. They are not physical parameters prematurely replaced with decimal numerical values. All floating-point values appear in observed sampling statistics/cost calculations or run durations after the exact construction.

Defects: none. Required fix: none.

## G. General/special/limit/numeric/failure regimes — PASS

The required table exists at README:98–104. Its rows match executable evidence:

| Regime | Evidence checked in this review |
|---|---|
| General symbolic | Symbolic six-product identity; complete native bilinear coefficient tensor; arbitrary positive-integer chain composition; regenerated exact selected plan. |
| Exactly solvable special | Native zero and one, every nonzero inverse, exact four scaffold powers, and all 4,096 word values under each Frobenius map. |
| Limit/decoupling | Literal `control=False` emits no circuit; randomly varied quantum control includes disabled/enabled branches; each divisor appears in both classical-condition branches; zero-inverse extension is asserted exactly. A continuum/asymptotic approximation is inapplicable to this fixed finite-field construction. No asymptotic optimality or speed claim is made. |
| Independent numeric/computational reference | Compiled arithmetic compared to native GF arithmetic; Frobenius matrices separately reconstructed by polynomial remainder; reduced unmodified official protocol rerun. Full candidate report and historical baseline audited rather than rerun. |
| Failure/edge | Wrong field, short register, invalid target string, singular binary map and invalid decoder dimensions rejected in supplied tests; invalid TowerPlan degree independently returned Failure. Register disjointness is an explicitly stated precondition, as in the native generator contract. |

The inverse mask test samples one independent mask per input, rather than all 4,096 masks per input. README:106–110 accurately describes this coverage and sampled measurement branches, so no exhaustive-mask claim needs repair. The source-level field/basis and phase contract supplies the exact structural argument underlying those sampled checks.

Defects: none. Required fix: none.

## Retained full-result audit and practical limits

All source/input/binary bindings matched the current files. The original binary size is 1,178,484,163 bytes and SHA-256 is `bc54c73b8725acbf6d246df5205f6a2dff2ea5783f335154833e57d9d93c8668`. The report fixes m=4095, n=70, degree=12, r=2016, 102,400 shots, and one worker. It records status PASS, zero value/phase errors, Q=1903, static Toffolis=4,503,030 and total operations=117,029,650.

Independent arithmetic audit confirmed:

- `32,408 + 888*6 + 4,095*5 == 58,211` multiplication computations.
- Static CCX + CSWAP totals equal the static Toffoli-like count.
- Executed CCX + CSWAP shot totals divided by 102,400 equal 4,502,782.126171875.
- Official source formula `R * Q^2 * T`, with R=1.1912604648115832, reproduces 19,425,188,366,155.996 exactly under the report's floating-point operations.
- The retained historical baseline 24,439,458,623,294.45 divided by the candidate cost equals 1.2581323878369532; the report's 20.517108559676135% reduction is consistent.

These are consistency and binding checks on a retained full run. They are not a new baseline run, a new full-simulation run in this review, a submission to Kaggle, a global-optimality proof, or a portable wall-clock speed benchmark. Nothing in the reviewed public text requires those stronger interpretations.

Unresolved concrete defects: 0. Unverifiable required findings: 0. No artifact changes are requested.

OPEN ISSUES: 0
