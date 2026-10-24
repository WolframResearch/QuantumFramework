# Independent WL quality review — round 1/3

Verdict: PASS on axes A–G; zero unresolved or unverifiable quality findings within the stated feasibility scope.

This review used the complete `wl-quality` skill and critic mandate, the `wl-native-style` and `wl-style-guide` rubrics, repository `CONTRIBUTING.md`, and the applicable Codex instructions (the ancestor `AGENTS.md` was empty; no subordinate `AGENTS.md` was present). No target artifacts were edited. All supplied executable generators/tests ran against a temporary copy.

## Artifact set and scope

Root reviewed: `OngoingProjects/VerifiableCircuitCompilation/other_tasks`.

Read in full: `FieldArithmetic.wl`, `ComponentMixing.wl`, `run-symbolic.wls`, `verify-symbolic.wls`, `symbolic.json`, `wl-checks.json`, `README.md`, `PLAN.md`, `field_probes.py`, `field-results.json`, `component_mixing.py`, `mixing-results.json`, `profile_baselines.py`, `dqi-baseline.json`, `dfthcsos-baseline.json`, `profile_stages.py`, `dqi-profile.json`, and `test_feasibility.py`.

The work claims exact finite-field multiplication identities and compiled reversible arithmetic fragments, together with real orthogonal DFTHC component transformations preserving the represented SOS Hamiltonian. It is feasibility research; a modified DQI decoder and global optimization are explicitly outside the claim (`README.md:3–5,63–68,88–93`; `PLAN.md:3–5`). No QF API is used, so the repository's QF audit anchor is not an API dependency of these artifacts.

The temporary review mirror is `$MIRROR/other_tasks`. Independent reviewer programs and results are its sibling files `independent.wls`, `independent-wl.json`, `independent.py`, `independent-python.json`; a fresh complete fragment run is in `fragments.log`.

## Executed evidence

`$KERNEL`, `$PYTHON` and `$MIRROR` denote the exercised direct kernel executable, pinned Python runtime and temporary review mirror; personal absolute paths are omitted.

1. Real WL command:

   ```sh
   "$KERNEL" -noinit -noprompt -script $MIRROR/other_tasks/verify-symbolic.wls
   ```

   Exit 0; kernel `15.0.0 for Mac OS X ARM (64-bit) (February 17, 2026)`; 25 tests, 25 succeeded, 0 failed. The reproduced symbolic SHA-256 is `150019e663e1fa3f0ca452d81250676ab6435dd947c22b2c6d2d1ad2d395ce78`. This executes all three exact tower constructions, complete tensor comparisons, every nonzero native inverse, primitive-generator checks, symbolic mixing identities, exact gain, limits, malformed inputs and confidence-bound calculations (`verify-symbolic.wls:5–37`).

2. Supplied Python command, run from the temporary mirror's parent:

   ```sh
   "$PYTHON" -m unittest other_tasks.test_feasibility -v
   ```

   Exit 0; `Ran 8 tests in 0.483s`, `OK`. Coverage includes dense independent pair-search references, synthetic Gram invariance, spectral-control Gram invariance, malformed vectors, source bindings, compiled controlled multiplication, inverse phase, zero and multibit masks.

3. Complete supplied fragment command, again in the mirror:

   ```sh
   "$PYTHON" other_tasks/field_probes.py
   ```

   Exit 0. The parsed regenerated report equals the original `field-results.json` exactly. All 12 fragments ran 8,192 shots each, totaling 98,304 fragment shots. Resource tuples reproduced: multiplication Q/T `36/63,54/63,48/45,46/42`; controlled Q/T `49/75,67/75,61/57,59/54`; inverse phase Q/T `72/315,90/315,84/225,82/210` for native, degree 6, degree 4, degree 3 respectively. All report clean scratch, restored inputs and correct phase.

4. Independent reviewer Python reference:

   ```sh
   "$PYTHON" $MIRROR/independent.py
   ```

   Exit 0. A separate carryless integer-polynomial multiplier, with explicit shift/XOR reduction, reconstructs every saved bilinear basis pair using the saved maps. All 1,728 bit coefficients per plan (5,184 total) agree with direct reduction modulo integer 4179; all 64 GF(8) input pairs agree with the six-product identity modulo 11. This reference does not call WL or kickmix multiplication for its arithmetic. Native `km.GF2Field(12).modulus` independently reports 4179. The best saved candidate retains every baseline array except `W` and `H2_I`; float64 augmented-Gram relative drift is `3.069217553325251e-16`. The small difference from the report's `3.2140392044921713e-16` is contraction-order roundoff. All source-hash checks attempted passed.

5. Independent reviewer WL extension:

   ```sh
   "$KERNEL" -noinit -noprompt -script $MIRROR/independent.wls
   ```

   Exit 0. The final corrected probe emits no diagnostic messages and records eight true checks: a noncommuting base, both component Hermiticities, mixed Hermiticity, exact shifted-SOS preservation, independent polynomial GF(8) arithmetic, identity limit, and rejection of the singular complex rotation parameter. The operator example is `H1=a X+b Z+c I`, `H2=d Y+e Z+f I`, with arbitrary real parameters and offsets `w,v`; the rotated shifted operators obey `Together[sum(Qnew.Qnew)-sum(Q.Q)] == 0` entrywise. Thus the Hamiltonian argument was checked beyond commuting scalar surrogates.

## A. Documented symbols — PASS

Evidence:

- `FieldArithmetic.wl:21,24–31,36–46` uses `FiniteField[2,polynomialFunction]`, field elements created by integer indices, coefficient/index properties, modular `Inverse`, and finite-field matrix `Inverse`.
- Installed WL 15 reference documentation was consulted via the supplied complete extraction files. `qf-next-FiniteField.md:3–5` documents `FiniteField[p,d]`, `FiniteField[p,f]` for an irreducible polynomial, and polynomial/exponential representation. `qf-next-FiniteFieldElement.md:3–5,11–21` documents element construction and the `"Index"`/`"Coefficients"` properties. Its line 29 explicitly supports `Inverse` and related linear algebra on finite-field matrices.
- The independent live kernel returns native irreducible polynomial `1+#1+#1^4+#1^6+#1^12&`, `field[2]["Coefficients"] == {0,1}`, and `field[2]["Index"] == 2`. `Options[FiniteField] == {}`; `Options[Inverse] == {Method->Automatic,Modulus->0,ZeroTest->Automatic}`. `Inverse[{{1,1},{0,1}},Modulus->2]` equals `{{1,1},{0,1}}`. For a nontrivial GF(8) matrix, mapping `"Index"` over `Inverse[m].m` returns the ordinary 2×2 identity. Saved return shapes are also checked by `verify-symbolic.wls:19`.
- No undocumented external WL package symbols are present. Public local symbols have explicit definitions and usage strings (`FieldArithmetic.wl:3–4,11–58`; `ComponentMixing.wl:3–8`). Their argument/error behavior is exercised by tests.
- Pinned kickmix source is clean at commit `638777fac4379bc62e2a1e4a52cf1a1753895d48`. Its `src/kickmix/py/circuit/circuit_builder.pybind.cc:2138–2215` documents the XOR/dirty-target multiplication contract and measurement/CZ cleanup; lines 2312–2345 and 2391–2422 document returned inversion scaffolds and zero-Toffoli deletion. Live method documentation confirms the keyword arguments actually used. The full fragment run confirms those contracts in context.
- The explicit basis/evaluation/output maps are the intended certificate artifacts. WL arithmetic and inversion use their documented built-ins rather than replacing them with an ad hoc field implementation.

Concrete defect: none. Smallest fix: none required.

## B. Idiomatic functional WL — PASS

Evidence:

- An `rg` scan of all target `.wl`/`.wls` files found no `Do`, `For`, `While`, or `AppendTo`.
- `FieldArithmetic.wl:26–39` uses `SelectFirst`, `Catenate`, `Map`, `MapThread`, `Partition`, and `NestList` for selection, embedding and interpolation maps; `:44–46` uses `MapThread` for matched output/evaluation contraction. `run-symbolic.wls:4,8–9` maps the plans and source hashes directly.
- The remaining `Table` uses construct genuine index grids (the product basis, polynomial coefficients, and the full 12×12×12 bilinear tensor), where explicit exponent/axis indices are appropriate. They do not mutate accumulators. The small `If` distinguishes a finite interpolation point from infinity (`:33–39`); it is not a dispatch cascade.
- `TowerPlan` needs lexical local functions and staged dependent values; its `Module` is justified. `MixingProofs` needs fresh arbitrary component symbols and local dependent expressions. `SixProductProof` uses `With` for literal constant coefficient arrays. Invalid inputs use pattern dispatch and structured `Failure`.
- The matrices are at most 12×12 and their nontrivial maps have substantial density. No materially sparse large operator is represented densely.

Concrete defect: none. Smallest fix: none required. No semantic refactor is justified solely to remove the readable tensor-index tables.

## C. Terse physics reasoning — PASS

Evidence:

- `FieldArithmetic.wl:21–46` follows the mathematical construction in order: native/subfield, compatible root, product basis, coordinate embedding, evaluation points, Vandermonde inverse, output maps, independent expected tensor and reconstructed tensor. Intermediate names correspond to these objects. The certificates follow at `:47–56`.
- `ComponentMixing.wl:9–29` states the rotation, transformed component/offset vectors, sign-region objective, invariant identities, then the exact normalization example. The scalar names in the bilinear proof are appropriate arbitrary coefficients rather than opaque temporary data.
- `README.md:31–68,70–108` connects each identity to its actual compilation or energy implication and states the measured/resource limitations. `profile_stages.py:43–59` cleanly separates conditional objective arithmetic from empirical baseline measurements. No unexplained production speedup or global-optimality conclusion is smuggled into the report.

Concrete defect: none. Smallest fix: none required.

## D. Physics invariants — PASS

Evidence:

- Arithmetic dimensions and basis ordering are checked exactly: `FieldArithmetic.wl:51–56` asserts the chosen root, binary basis inverse, interpolation inverse, and complete tensor equality. These execute in the 25/25 verifier, and all three native-basis tensors independently pass the polynomial reference described above.
- `field_probes.py:133–163` explicitly checks reversible XOR output with dirty targets, unchanged operands and controls, no spurious multiplication phase, the inverse's diagonal phase, restored source, and zero scratch. These assertions passed on the fresh 12-fragment run. `test_feasibility.py:79–84` additionally passed zero and multibit phase masks. Relevant quantum invariants are asserted rather than merely mentioned.
- `ComponentMixing.wl:16–23` asserts real orthogonality, arbitrary bilinear preservation, the identity-offset contraction and squared offset norm, identity and limiting sign flip. All pass the real kernel. These preserve the full SOS coefficient tensor, the one-body offset contribution and the SOS energy.
- Independent arbitrary-real Hermitian 2×2 operators with noncommuting Pauli components retain both Hermiticity and the full shifted SOS matrix exactly. The augmented-Gram check on the actual saved candidate and the unchanged metadata checks independently substantiate the implementation connection (`component_mixing.py:90–103`).
- Positivity/trace conditions of density matrices and state-evolution normalization are not physical dependencies here: no density operator or Hamiltonian time evolution is constructed. The arithmetic maps are reversible XOR/diagonal-phase operations; DFTHC preserves an existing real SOS Hamiltonian. This applicability distinction is already stated at `README.md:120–124`.

Concrete defect: none. Smallest fix: none required.

## E. Nontrivial examples — PASS

Potential base cases are scalar multiplication alone, zero rotation alone, or a normalization change without identity offsets.

The artifacts exceed them: complete native degree-12 bilinear tensors for three different subfields; dirty-target controlled reversible fragments; all 4,096 source words; nontrivial phase masks; arbitrary component and offset symbols; an exact two-component factor-of-two LCU-H2 example (`ComponentMixing.wl:26–28`); synthetic rectangular Gram tests (`test_feasibility.py:31–41`); and the full-dimensional real factorization with a deliberately adverse spectral control (`mixing-results.json:29–59`). The independent noncommuting Hermitian example also excludes dependence on commuting operators. The fresh execution results above support these claims.

Concrete defect: none. Smallest fix: none required.

## F. Symbolic-first derivations — PASS

Evidence:

- The field modulus, subfield arithmetic, inverse maps, GF(8) identity and coefficient tensors are exact (`FieldArithmetic.wl:11–56`). Fixed characteristic and degree are native problem data, not avoidable numerical substitutions for mathematical parameters.
- Pair-rotation parameter `t`, component/offset entries, and objective angle `theta` remain symbolic (`ComponentMixing.wl:6–28`). `TrigReduce`/`TrigExpand` derive the sign-region objective before its numerical search is used. The exact gain is 2 rather than an approximate decimal assertion.
- The independent matrix extension retains arbitrary real coefficients, offsets and rotation parameter and produces exact zero differences. It corroborates the symbolic argument without numerical substitution.
- `component_mixing.py:23–47` numerically enumerates sign boundaries and stationary points for finite submitted data. `README.md:88–93` correctly calls this a numerical instantiation and bounded coordinate search. It does not claim a floating-point optimality certificate or substitute it for the symbolic invariant proof.
- Numerical baseline/resource observations and confidence estimates are appropriate final evaluations. `verify-symbolic.wls:8–10` begins the confidence calculation with exact `delta=2^-128` and the actual fixed shot count before numerical evaluation.

Concrete defect: none. Smallest fix: none required.

## G. Regimes and limits — PASS

The supplied regimes table exists at `README.md:110–118`. Each required row has executable evidence:

| Regime | DQI evidence | DFTHC evidence |
|---|---|---|
| General symbolic | Exact full bilinear tensors for all three field plans (`FieldArithmetic.wl:43–56`); 25/25 WL reproduction; independent 5,184-bit tensor comparison | Arbitrary component/offset bilinear, contraction and norm identities (`ComponentMixing.wl:16–20`); independent arbitrary-real noncommuting SOS check |
| Exactly solvable special | Six-product GF(8) polynomial identity (`:11–18`), all native nonzero inverses (`verify-symbolic.wls:21`), independent complete 64-pair GF(8) reference | Exact 2× LCU-H2 example (`ComponentMixing.wl:26–28`, `verify-symbolic.wls:30`) |
| Limit/decoupling | Both controls in compiled XOR tests (`field_probes.py:147–156`); zero-mask identity phase (`test_feasibility.py:79–84`) | Identity rotation and `t->Infinity` sign flip (`ComponentMixing.wl:21–22`), reproduced by WL; independent identity check |
| Independent numeric | Fresh all-12 compiled simulation with dirty targets, input/scratch/phase checks; independent carryless-polynomial reference | Dense 4,097-angle pair references and independent synthetic Gram tests (`test_feasibility.py:18–41`); actual-candidate Gram and metadata check |
| Failure/edge | Unsupported/inexact degrees and zero extension of inversion (`verify-symbolic.wls:22,25–26`); zero does not claim a multiplicative inverse | Complex and malformed rotation rejection (`verify-symbolic.wls:33–34`); malformed vector tests (`test_feasibility.py:43–45`); adverse spectral-control outcome is clearly retained |

Concrete defect: none. Smallest fix: none required.

## Coverage limits retained in this verdict

The expensive full unchanged DQI baseline and complete 254,464-shot DFTHC verifier were not repeated in this quality review. Their actual full-run records were read completely, and their generator/source bindings were checked; no empirical full-run reproduction beyond that is claimed. The complete compiled fragment simulation, supplied WL verifier and supplied Python unit suite were rerun. The verifier's quantum measurement branches are sampled; the README correctly does not claim exhaustive branch enumeration. No modified decoder, peak-scratch integration test, new submission, implementation recommendation or global optimizer has been produced or inferred by this review.

Those limits agree with the artifact's declared feasibility scope and are not unexecuted claimed physics invariants. Every applicable invariant needed for the quality verdict has fresh executable evidence described above.

OPEN ISSUES: 0
