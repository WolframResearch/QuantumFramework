# Independent WL quality review — round 1/3

Reviewed all twelve named artifacts in full against the critic mandate, without reading prior reviews or editing targets. The applicable user-level AGENTS.md is empty; no nearer AGENTS.md was found. Reproduction used isolated temporary copies, Wolfram 15.0.0 ARM (February 17, 2026), and the supplied Python environment and DATA. No raw inputs or derived tables are included here.

## Bound inputs

| File | SHA-256 |
|---|---|
| Representations.wl | 85817e6845e6dfaf53673e9773d2b5ca168ebbb2be0c56e611e081e45ed50a19 |
| ProjectiveDecoder.wl | e0a03b7d7b1c394936d273747ce39c88c556426cf40b9d7525f22e69dff554bb |
| run.wls | 2d85571bb326105edf2a47c2969a78cbcafe3c69487156cde8bd2c34a6cec205 |
| curvature_screen.py | 821a6119382eaa0d7329688ae833444fc5c39a6196bf3ace84ec50ec6a91f0f3 |
| affine_rounding.py | 8bd41280d90b06e502101977b3e088b1c66af8612b8d9eff012e3de8933c35de |
| table_screen.py | 60b730ca9286d2e2e6935a1308ec618ca123e5d7321d76103cb77004f7e1d051 |
| RESULTS.md | bd1062c8e69ef41b78bcf261cd36ad06ef203f5903970fa0acbf9c4873f9743d |
| PLAN.md | 4d2c0833bd60bd73d921aa54dfe2e83151447a2c7332a1ccf4e0084dfda6a745 |
| symbolic-checks.json | 78035f544058ad635240b030fe4fb828ed1b9b00a2990df09bca87429a1d8273 |
| curvature-screen.json | 9e8e79ecfbc3390a26595ac9168062139f3a28e97c86ae3bd9b981f21dd3174d |
| affine-rounding.json | baa761187fb066ee3e9b3266b41b7225847f04340d5636c9dbfac2dbd722402b |
| table-screen.json | a2ba5abfbe00d5e4b0ee8e8d1b10325eb304c9f0eafaf7ad1c0071ed2813cda5 |

## Execution evidence

From the isolated reproduction directory:

```text
WolframKernel -noinit -noprompt -script run.wls
    status -> PASS; every Boolean check -> True; exit 0
python curvature_screen.py --data DATA
    f8 naive: 1,590,097 failed; f8 compensated: 385,356 failed
    f9 compensated: 186 failed
    f10 compensated: 0 failed; widths (30,22,9); worst ratio 0.8505589342620424
    f11/f12: 0 failed; worst ratio 0.7464338467124986
    all candidates: 0 nonpositive and 0 overflowing accumulators; exit 0
python affine_rounding.py --data DATA
    modified [(2,3,0),(7,-1,0),(158,-1,0),(183,0,1),(184,-2,0)]
    unresolved [170]; exit 0
python table_screen.py --words WORDS
    nonlinear_support [983]; all_output_rank [108]; difference_bank_rank [81]; exit 0
```

Parsed JSON equality with the four supplied reports returned no different fields. No messages occurred in the target run. Additional live probes below used no message suppression.

## A — PASS: documented signatures and return shapes

Representations.wl:6–18 correctly interprets the extended-GCD result as `{gcd,{s,t}}`, including the constant gcd scale, and reduces the inverse and quotient. Live `PolynomialExtendedGCD` over the native three-root field case returned outer/inner lengths `{2,2}`; `Expand[{s,t}.{sigma,z D[sigma,z]}-gcd]==field[0]` returned `True`. `CachedQuotient` returned an Association with keys `{"denominator","inverse","quotient"}`. `Head[field]` and `Head[field[2]]` returned `FiniteField` and `FiniteFieldElement`.

The live `::usage` results and `Options` matched the signatures used. Extended-GCD options were `{Method->Automatic,Modulus->0}`; PolynomialReduce options included exact rational functions, modulus 0 and tolerance 0. Current primary documentation confirms the custom-polynomial constructor and finite-field extended GCD: [FiniteField](https://reference.wolfram.com/language/ref/FiniteField.html), [PolynomialExtendedGCD](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html), [PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html). The quotient ideal check correctly extracts the remainder. There are no external WL package symbols in these targets.

## B — PASS: functional WL

The projective recurrence is expressed with `FoldList`/`Fold` (ProjectiveDecoder.wl:19,35), and the transform construction with `Nest` (Representations.wl:35). Root, coefficient and invariant checks use `Map`, `MapThread`, `Outer` and `AllTrue`. No `Do`, `For`, `While` or `AppendTo` occurs in the WL targets. `Table` at ProjectiveDecoder.wl:18 is a direct index-dependent syndrome sequence, not an imperative list builder. Modules require isolated polynomial symbols or named intermediate objects. The matrix sizes are bounded at 64; retaining a dense test matrix is reasonable here. No mandatory rewrite identified.

## C — PASS: terse physical reasoning

Names identify locators, evaluators, roots, discrepancies, quotients, masks and phases. RESULTS.md:10–19 gives the compensation identity followed immediately by its minimax proof; :40–64 connects the quotient certificate to the phase. The projective algebra and its explicitly unresolved reversible costs are separated at :96–115. The arithmetic ceilings reproduce as 909,090 removed operations, width ceiling 2130.12847, 2415 field products and 405,720 budgeted operations. These are explicitly budgets, not implementation evidence. The determinant defect is counted once under D.

## D — FINDING: major-update determinant omits its polynomial shift

**Q1.** ProjectiveDecoder.wl:11,13 actually applies

```text
{Cnew,Bnew} = {g C - d z^gap B, C}
```

but :39 tests `Det[{{g,-d},{1,0}}]==d`. RESULTS.md:111–113 then describes this as the determinant of the update and infers nonsingularity. At the already nontrivial special case `gap=1`, the actual polynomial matrix gives the live result

```wl
Det[{{g,-d z},{1,0}}]
(* d z *)
Det[{{g,-d z},{0,1}}]
(* g *)
```

Thus `d` is not the major-update determinant. Moreover, nonzero `d z^gap` establishes invertibility over rational functions, not a bijection of arbitrary fixed-width coefficient registers. The two-coefficient truncation of that same shifted update has

```wl
Det[{{g,0,0,0},{0,g,-d,0},{1,0,0,0},{0,1,0,0}}]
(* 0 *)
```

because a highest coefficient is shifted out. This does not refute the successful forward recurrence tests or imply a failed claimed circuit: no reversible decoder is claimed. It does invalidate the supplied nonsingularity justification.

**Smallest fix:** test and state the actual shifted matrix, giving determinant `d z^gap`; qualify the algebraic domain and say recovery requires division by `z^gap` on the actual image, with boundaries and register retention unresolved. Alternatively remove both determinant tests and the nonsingularity claim, retaining the existing explicit reversible-cost limitation. A scalar determinant is insufficient as coefficient-register reversibility evidence.

The stated physical invariant otherwise passes. Representations.wl:30–32 proves quotient equality at roots; :48–63 checks the released reversed evaluator and `z sigma'` conventions. The released source scan at lines 545–584 independently confirms those slot exponents and the root-masked product. The trace Gram determinant checks all binary linear functionals. A new native-field three-root check explicitly tested **all 4096 trace masks** against the direct Forney expression and returned `True`. `ConjugateTranspose[DiagonalMatrix[{1,1,1,-1}]].DiagonalMatrix[{1,1,1,-1}]==IdentityMatrix[4]` returned `True`. Hamiltonian, density-matrix and workspace invariants are inapplicable because this deliverable claims phase identities only (RESULTS.md:149–151).

## E — PASS: nontrivial examples

The base cases are degree one, one error and zero compensation. The targets additionally test the native field at locator degrees 1,2,17,35, arbitrary trace masks, degree-35 recovery from 70 syndromes and all 35 Forney values (Representations.wl:38–65; ProjectiveDecoder.wl:16–46). All reproduced `True`. Exhaustive nested integer evaluation covers all 3,133,452 legal inputs, not only polynomial extrema. All 25 fixed tables exhibit the complete nonlinear-support/rank obstruction. No extra unrelated quantum example is needed for the explicitly classical-algebra scope.

## F — PASS: symbolic first

Curvature, uniform bounds, equioscillation, minimax lower bound, quotient ideal membership and projective scaling remain exact symbolic expressions (Representations.wl:72–101; ProjectiveDecoder.wl:36–38). The integer screen starts with `Fraction` and uses integer shifts/subtractions; floating point enters only in the final displayed error/tolerance ratio. Fixed n=70, field degree 12 and table sizes identify the actual benchmark; they do not replace the general identities. A numeric check at `d=1/7,h=1/64` produced equal extrema `4.3596540178571428571*10^-6`.

## G — PASS: regimes and limits, within the stated scope

RESULTS.md:139–147 supplies the required table. Its rows have executable evidence: general exact identities; degree-one quotient and one-error recurrence; `d->0` and no-error locator; full-domain official integer references and native-field comparisons; coarse precision, finite-search rejection, repeated roots, zero roots, empty locator and malformed exponent. An additional live `Limit[d (t^2-h t+h^2/8)/.t->s h,h->0,Assumptions->Element[{d,s},Reals]]` returned `0`, checking shrinking intervals as well as vanishing curvature change. Projective zero scale is explicitly outside the representation; native adaptation and reversible implementation are expressly unimplemented, so their missing circuits are not counted as hidden successes. Q1 remains the one open issue.

OPEN ISSUES: 1
