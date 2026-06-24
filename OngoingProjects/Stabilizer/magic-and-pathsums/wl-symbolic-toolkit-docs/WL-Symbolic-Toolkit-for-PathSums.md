# Wolfram Language symbolic toolkit for the path-sum simulator

A mapping of 91 Wolfram Language function documentation pages to the machinery of `../Beyond-SymPhase-Symbolic-Algebra-For-Magic.md`. The pages were collected from the [FiniteFields](https://reference.wolfram.com/language/guide/FiniteFields.html) and [PolynomialAlgebra](https://reference.wolfram.com/language/guide/PolynomialAlgebra.html) guides, the `Hold*` evaluation-control family, and `NonThreadable`; each is saved in full in this folder as `<Name>.md`. Every function below is read in full and tagged **Core** / **Supporting** / **Tangential** / **Unrelated** relative to the path-sum pipeline.

> **Framing caveat (see `../WL-Symbolic-Novelty-Reexamination.md`).** Tagging each function "relative to the path-sum pipeline" answers *which primitive implements an existing step*. A later re-examination found that this framing under-recorded two findings that surface only under the opposite question, *what can these primitives do that the competing tools cannot*: exact symbolic-parameter path sums (a capability no competitor has) and the Gröbner/variety point count (a distinct, unbuilt tractability axis, the concrete form of an Amy-Stinchcombe open problem). Read the `GroebnerBasis` entry in Stage 4 and the exact-arithmetic entries in Stage 3 with that in mind: they are not only implementation upgrades.

The path-sum simulator carries a Clifford+$T$ amplitude as a Gauss sum $\langle y\lvert C\rvert 0\rangle = 2^{-h/2}\sum_v \omega^{\varphi(v)}$, with $\omega = e^{i\pi/4}$, qubit values multilinear over $\mathbb{F}_2$, and the phase $\varphi$ a degree-$\le 3$ pseudo-Boolean polynomial valued in $\mathbb{Z}_8$. It reduces the sum by a confluent rewrite system (the $[HH]$ Hadamard-elimination rule).

## Three headline connections (kernel-verified)

1. **The field trace is the additive character of the Gauss sum.** The additive character of $GF(2^k)$ is $\chi(x) = (-1)^{\mathrm{Tr}(x)}$, where $\mathrm{Tr}$ is the absolute field trace. The report's amplitude is exactly a sum of such characters. WL exposes this natively: `FiniteFieldElementTrace`. Verified: $\sum_{x\in GF(2^4)}(-1)^{\mathrm{Tr}(x)} = 0$ (the trace is balanced, $8/8$). This is the non-obvious bridge from WL's `FiniteField` machinery to the report's character sums, and the natural route to generalize the simulator from $\mathbb{F}_2$ to $GF(2^k)$.

2. **The "variety-valued" generalization the report names is `GroebnerBasis` over the Boolean ideal.** Adjoining the idempotency relations $x_i^2 - x_i$ and working `Modulus -> 2` makes the radical Boolean ideal, whose normal forms (via `PolynomialReduce`) are the multilinear monomials and whose finite quotient ring turns solution counting into linear algebra. For *linear* polynomials a Gröbner basis is exactly Gaussian elimination, so the all-Clifford fold is row reduction over $\mathbb{F}_2$. Verified: `PolynomialReduce[x^3 + x, {x^2 - x}, {x}, Modulus -> 2]` returns the idempotent normal form $0$.

3. **`Collect` performs the $[HH]$ split, and `Cyclotomic`/`Root`/`RootReduce` are the exact $\mathbb{Z}[\omega]$ amplitude arithmetic.** `Collect[φ, x_v]` writes $\varphi = \varphi_0 + x_v\,C_v$ in one call, exactly the cofactor split the elimination needs (verified). `Cyclotomic[8, x] = 1 + x^4` is $\omega$'s minimal polynomial (verified), and `RootReduce` canonicalizes amplitudes exactly in $\mathbb{Z}[\omega, 1/\sqrt2]$ with no floating point.

---

## Stage 1: the $\mathbb{F}_2$ / $GF(2^k)$ carrier (qubit values, branch variables)

- **PolynomialMod** [Core] `PolynomialMod[poly, m]` reduces coefficients mod `m`; the prototype's $\mathbb{F}_2$ (`,2`) and $\mathbb{Z}_8$ (`,8`) reducer. The list form `PolynomialMod[poly, {x^2-x, 2}]` imposes idempotency and the modulus together, giving the $\mathbb{F}_2[V]/(x^2-x)$ normal form directly. Does not divide (unlike `PolynomialRemainder`), so it never inverts the modulus.
- **Modulus** [Core] the option `Modulus -> p` that turns `Solve`/`Reduce`/`Factor`/`GroebnerBasis`/`MatrixRank`/`LinearSolve`/`Inverse`/`Det` into $\mathbb{F}_p$ operations. The single switch behind every "over $\mathbb{F}_2$" step.
- **FiniteField** [Core] `FiniteField[2, k]` is $GF(2^k)$ with native exact arithmetic; `FiniteField[2]` is $\mathbb{F}_2$ itself. The native substrate the prototype currently hand-rolls with `PolynomialMod[.,2]` plus manual power reduction.
- **FiniteFieldElement** [Core] `ff[{c0,...,c_{k-1}}]` is a length-$k$ bit vector, the natural carrier of a branch-variable assignment $v\in\mathbb{F}_2^k$; `"Coefficients"` extracts the coordinate vector. Note the index-based API: `ff[123]` means index 123, `ff[{...}]` means a coefficient vector.
- **Expand** [Core] flattens gate compositions into a monomial sum so `MonomialList`/`CoefficientRules` can read them; `Modulus -> p` reduces in place. Does *not* impose $x^2 = x$ by itself; idempotency is a separate rewrite (matching the prototype).
- **FromFiniteFieldIndex** [Supporting] `FromFiniteFieldIndex[Range[0, 2^h-1], FiniteField[2, h]]` (Listable) materializes every assignment, the natural "enumerate the domain of the Gauss sum" primitive.
- **FiniteFieldIndex** [Supporting] the inverse map, element to integer; for $\mathbb{F}_2$ the binary integer of the coefficient vector. Bookkeeping for enumeration.
- **ToFiniteField** / **FromFiniteField** [Supporting] inject integer/polynomial data into $GF(2^k)$ and read it back; the coercion entry/exit points. `FromFiniteField[FiniteFieldElementTrace[x], ff]` returns the $0/1$ trace bit (used in the character verification above).
- **FiniteFieldEmbedding** [Supporting] embeddings between subfield towers $GF(2^k)\hookrightarrow GF(2^{km})$; supports relative traces and summing over subfields (the tower side of the variety generalization).
- **IrreduciblePolynomialQ** [Supporting] `IrreduciblePolynomialQ[f, Modulus -> 2]` is the validity check to build $GF(2^k)$ from a chosen modulus $f$.
- **PrimitivePolynomialQ** [Supporting] needed only for the `"Exponential"` (cyclic-group) field representation; the polynomial (bit-vector) representation only needs irreducibility.
- **MinimalPolynomial** [Supporting] bridges an element to its $\mathbb{F}_2$-algebraic data; $\mathrm{Tr}(a)$ and $N(a)$ are read off its coefficients, and its degree tells you whether an element generates $GF(2^k)$.

## Stage 2: the phase polynomial $\varphi$ (degree-$\le 3$, $\mathbb{Z}_8$)

- **Collect** [Core] `Collect[φ, x_v]` is the $[HH]$ split $\varphi = \varphi_0 + x_v\,C_v$ in one step; `Modulus` keeps it in $\mathbb{Z}_8$ / $\mathbb{F}_2$.
- **Coefficient** [Core] `Coefficient[φ, x_v]` is the cofactor $C_v$ and `Coefficient[φ, x_v, 0]` is $\varphi_0$ (the multilinearity makes $\varphi$ linear in each $x_v$).
- **CoefficientRules** [Core] the sparse `{exponent-vector -> coefficient}` representation of $\varphi$: each rule is one monomial, the exponent vector says which branch variables it contains, the coefficient is its weight (the $4$ in $4x_ix_j$ from $CZ$, the linear $T$ term). Degree $\le 3$ is "max total exponent $\le 3$". `Modulus -> 8`/`2` keep the ring.
- **MonomialList** [Core] the term-by-term handle on $\varphi$ for the monomial-level rewrites; already used by the prototype.
- **Variables** [Core] enumerates the active branch variables (the summation index set of the Gauss sum) and fixes the variable order the other readers need; already used. Caveat: its output is not always sorted, so fix an explicit order.
- **CoefficientList** [Supporting] dense coefficient tensor; for a single variable the length-2 list is $\{\varphi_0, C_v\}$, but `CoefficientRules` is better for the full multilinear form.
- **Exponent** [Supporting] degree bookkeeping; `Exponent[φ, x_v]` should be $\le 1$ (idempotency invariant), and total degree checks the degree-$\le 3$ assumption.

## Stage 3: the amplitude (Gauss/character sum, exact in $\mathbb{Z}[\omega, 1/\sqrt2]$)

- **FiniteFieldElementTrace** [Core, keystone] `FiniteFieldElementTrace[a]` is the absolute trace $\mathrm{Tr}(a) = a + a^2 + \cdots + a^{2^{k-1}} \in \mathbb{F}_2$; `(-1)^FiniteFieldElementTrace[a]` is the additive character. A Gauss/Weil sum is `Total[(-1)^(FiniteFieldElementTrace /@ values)]`. Every $\mathbb{F}_2$-linear functional is a trace, so all linear/affine exponents in a character sum are traces. Relative-trace and embedding forms support subfield/tower sums.
- **FrobeniusAutomorphism** [Core] $x \mapsto x^2$ on $GF(2^k)$: the $\mathbb{F}_2$-linear field homomorphism that builds the trace ($\mathrm{Tr} = \sum_k \phi^k$) and the idempotency $x^2 = x$ on $\mathbb{F}_2$, and whose orbit is the conjugate/root set of `MinimalPolynomial`.
- **Cyclotomic** [Core] `Cyclotomic[8, x] = 1 + x^4` is $\omega$'s minimal polynomial; the defining relation the whole exact arithmetic reduces modulo. $\mathbb{Q}(\omega) = \mathbb{Q}[x]/(x^4+1)$.
- **Root** [Core] the exact algebraic carrier for an amplitude; $\omega$ is `Root[1 + #^4 &, ...]`, exactly comparable and arbitrary-precision-evaluable, no floating point.
- **RootReduce** [Core] canonicalizes any algebraic combination of amplitudes to a single `Root` (already used: `{1/Sqrt2, 1/2 + I/2}`); canonical form means two amplitudes are exactly equality-testable. Caveat: the reduced minimal polynomial's degree can blow up for tall towers.
- **Extension** [Core] the option (e.g. `Extension -> {(-1)^(1/4)}`) that puts `Factor`/`PolynomialGCD`/`Together` into the cyclotomic field $\mathbb{Q}(\omega)$ instead of treating $\omega$ as an inert symbol.
- **ToRadicals** [Supporting] renders a `Root`-form amplitude as an explicit radical ($1/\sqrt2$, $\tfrac12 + \tfrac{i}{2}$); a presentation layer (every $\mathbb{Z}[\omega]$ amplitude is degree $\le 4$, always solvable by radicals).
- **GaussianIntegers** [Supporting] `GaussianIntegers -> True` $\equiv$ `Extension -> I`: the $\mathbb{Z}[i] = \mathbb{Z}[\omega^2]$ Clifford-only sublevel. Cannot reach the full $8$th-root ring (needs $\sqrt2 = \omega + \omega^{-1}$).
- **RootSum** [Supporting, with caveat] `RootSum[f, form]` is $\sum_{f(x)=0}\mathrm{form}(x)$, a closed form for rational `form` *without finding the roots*; note `f` and `form` are both pure functions, so the cyclotomic case is `RootSum[1 + #^4 &, form]` (or `Cyclotomic[n, #] &`), not `RootSum[Cyclotomic[n, x], form]` (which errors, since `Cyclotomic[n, x]` is an expression, not a function). It then sums over the $n$th roots of unity, capturing the *linear* root-of-unity power sums and rational-of-roots sums. **Honest limit:** it does NOT give the quadratic Gauss sum $\sum_k e^{i\pi k^2/n}$ in closed form, because the $k \mapsto k^2$ exponent is not "one fixed function over the roots of one polynomial." Use the report's own closed-form quadratic-Gauss-sum identities for that.
- **FiniteFieldElementNorm** [Tangential] the multiplicative companion $N(a) = \prod_k a^{2^k}$; enters multiplicative (Jacobi-sum) characters and the determinant of "multiply by $a$", not the additive character on the critical path.

## Stage 4: the $[HH]$ elimination / rewrite engine

- **GroebnerBasis** [Core] the named "variety-valued" future direction: `GroebnerBasis[polys ∪ {x_i^2 - x_i}, vars, Modulus -> 2]` is the radical Boolean ideal; basis $\{1\}$ means no common zero (inconsistent constraints, amplitude $0$). For linear polynomials it *is* Gaussian elimination; the 3-argument form eliminates variables (sums one out).
- **PolynomialReduce** [Core] `PolynomialReduce[poly, gb, vars, Modulus -> 2]` gives the unique normal form against a Gröbner basis; remainder $0$ is the ideal-membership test (a Clifford constraint already implied). Verified idempotent reduction $x^3 + x \to 0$.
- **Reduce** [Core] `Reduce[constraints, vars, Modulus -> 2]` (or domain `FiniteField[2, k]`) returns the complete, finite $\mathbb{F}_2$ solution set, deciding consistency and enumerating the support for the sum; uses Gröbner bases internally.
- **Resolve** [Core] `Resolve[Exists[v, constraints], Booleans]` (or `FiniteField[2, 1]`) decides satisfiability of the support and projects out (sums out) a variable over $\mathbb{F}_2$, keeping inequations (unlike `Eliminate`).
- **Solve** / **SolveValues** [Supporting] `Solve[linearConstraints, vars, Modulus -> 2]` solves the determined $\mathbb{F}_2$ system; `SolveValues` returns bare value tuples to feed the sum. Generic solutions only, so prefer `Reduce` for a complete count.
- **FindInstance** [Supporting] a single witness assignment in the support (or a Boolean-SAT witness); does not count.
- **Eliminate** [Supporting] `Eliminate[eqns, vars]` (modular via `Modulus == 2`) projects away the summed variable, returning the implied equation-only constraints (Zariski closure).
- **Resultant** [Supporting] pairwise variable elimination between two polynomials over $\mathbb{F}_2$; cheaper than a full Gröbner basis for one variable and two polynomials.
- **PolynomialGCD** [Supporting] the univariate special case of `GroebnerBasis` (`Modulus -> 2`): a common-root relation between two single-variable constraints.
- **PolynomialQuotient** / **PolynomialRemainder** / **PolynomialQuotientRemainder** [Supporting] univariate $\mathbb{F}_2[x]$ division; subsumed by `PolynomialReduce` in the multivariate setting. `PolynomialRemainder` divides (can invert coefficients), so it is not interchangeable with `PolynomialMod`.
- **PolynomialExtendedGCD** [Tangential] Bézout cofactors / modular inverses in $\mathbb{F}_p[x]/(q)$; the simulator needs membership and counting, not Bézout certificates.
- **PolynomialLCM** [Tangential] no role in summing over a Boolean variety.

## Stage 5: exact linear algebra (cap the over-count, fold the Clifford part)

- **MatrixRank** [Core] `MatrixRank[m]` on the matrix of $2^{\#T}$ component state vectors with exact $\mathbb{Z}[\omega]$ entries reads off the true span dimension $\le 2^n$, capping the eager stabilizer-frame over-count; `MatrixRank[m, Modulus -> 2]` does the $\mathbb{F}_2$ rank. Caveat: it assumes free symbols are independent (moot for concrete cyclotomic entries).
- **RowReduce** [Core] the literal $\mathbb{F}_2$ Gaussian-elimination engine (`Modulus -> 2`) under the all-Clifford fold (the "linear Gröbner basis = Gaussian elimination" fact); the common substrate of rank/null-space/solve.
- **NullSpace** [Core] `NullSpace[m, Modulus -> 2]` is the full solution space of the homogeneous $\mathbb{F}_2$ constraints $Q = 0$; on the rank side, `Length[NullSpace] == cols - MatrixRank`.
- **LinearSolve** [Core] `LinearSolve[m, b, Modulus -> 2]` is the inhomogeneous $\mathbb{F}_2$ solver (already used); over $\mathbb{Z}[\omega]$ it expresses one component as an exact combination of others. Returns one particular solution; pair with `NullSpace` for completeness.
- **Dot** [Supporting] builds the component-vector matrix ($U.\psi$) and verifies solutions. Two caveats: no `Modulus` option (use an explicit `Mod[., 2]`), and it is non-Hermitian (use `Conjugate[a].b` for genuine state overlaps).
- **Det** / **Inverse** / **LUDecomposition** [Tangential] all carry `Modulus` and work on exact entries, but they assume square / full-rank / non-degenerate matrices, which is exactly the structure the report does not have (it measures rectangular rank-deficiency); at most cheap sanity checks.

## Stage 6: implementing the rewriter (evaluation control)

These are how to carry the *formal* sum and write the rewrite rules idiomatically; the `Hold*`, `NHold*`, `SequenceHold`, `NonThreadable` are **attributes** (set with `SetAttributes`), the rest are **wrapper functions**.

- **HoldPattern** [Core] wrap a rule's left-hand side that contains evaluating heads (`Plus`, `Power`, an inert `Sum`) so the rule matches the inert structure instead of collapsing at definition time. The load-bearing tool for writing the $[HH]$ / Amy-Stinchcombe rules.
- **Inactive** [Core] `Inactive[Sum][ω^φ[v], {v, ...}]` is an inert sum that pattern-matches and rewrites while still indexing like a sum; individual operators stay inert while others evaluate (matches peeling operators one at a time). Caveat: `Inactive[Sum]` lacks `Sum`'s `HoldAll`, so guard the bound index.
- **Inactivate** [Core] `Inactivate[expr, Sum | Plus | Power]` freezes exactly the operators the rewriter manipulates in one pass, leaving `List`/`Rule` structure live.
- **Activate** [Core] the inverse: re-enable chosen operators after reduction to collapse the residual sum to the closed-form amplitude; selective `Activate[., patt]` releases one operator class at a time.
- **HoldAll** [Core] `SetAttributes[PathSum, HoldAll]` so a custom path-sum head does not pre-evaluate its symbolic phase into a number when built.
- **Hold** + **ReleaseHold** [Core] the hold-then-release lifecycle: keep the sum unevaluated while rules fire (substitution and upvalues still act inside `Hold`), then `ReleaseHold` to compute.
- **HoldComplete** / **HoldAllComplete** [Supporting] seal a term against upvalues and `Sequence` leakage when `Hold`/`HoldAll` is too permissive.
- **HoldFirst** / **HoldRest** [Supporting] finer-grained holding (e.g. hold the phase slot, not a numeric qubit count).
- **Evaluate** / **Unevaluated** [Supporting] per-call escape hatches: force one slot of a `HoldAll` head (`Evaluate`), or pass one argument unevaluated through a non-holding function (`Unevaluated`).
- **SequenceHold** [Supporting] stop `Sequence` objects in the argument list from auto-splicing (a robust carrier whose argument boundaries are preserved).
- **NonThreadable** [Supporting] mark the path-sum head non-scalar so it does not thread element-wise over a list payload (its own coefficient vector, a batch of amplitudes) and shatter; keeps the object a single algebraic atom under `+`, `*`.
- **HoldForm** [Supporting] display the formal sum literally (invisible wrapper) in report cells and traces.
- **HoldCompleteForm** [Tangential] display-only fully-shielded hold (sequences/upvalues visibly preserved); rarely needed.
- **NHoldAll** / **NHoldFirst** / **NHoldRest** [Tangential] keep exact data ($\omega$, integer exponents mod 8, $\sqrt2$ prefactors) from being floated by `N`; protect exactness, do nothing for the symbolic rewriting.

## Tangential and unrelated (honest)

- **Factor**, **FactorList**, **FactorSquareFree**, **FactorSquareFreeList** [Tangential] multiplicative factorization / root multiplicities; the phase is manipulated additively (by monomials and the $\varphi_0 + x_v C_v$ split), never as a product of irreducibles. At most incidental simplification of a closed-form coefficient.
- **Together**, **Cancel** [Unrelated] rational-function normalization (common denominator, gcd cancellation); the path-sum forms have no denominators.
- **Decompose** [Unrelated] functional composition $P_1 \circ P_2$ of univariate polynomials; the report decomposes $\varphi$ additively, not as a composition.
- **Discriminant** [Tangential] `Discriminant[Cyclotomic[n, x], x]` relates to the $\sqrt{n}$ magnitude of a quadratic Gauss sum (a number-theoretic cross-check), but computes a single number, not the sum.
- **SymmetricPolynomial**, **SymmetricReduction** [Tangential] Vieta / Newton-Girard themes (symmetric functions of roots of unity) already covered by `RootSum`/`RootReduce`; not pipeline components.
- **MultiplicativeOrder**, **FiniteFieldElementPrimitiveQ** [Tangential] multiplicative-group / discrete-log properties; the additive character does not need a primitive element.
- **CharacteristicPolynomial** [Tangential] $\det(m - xI)$; the report's linear algebra is rank/null-space over $\mathbb{Z}[\omega]$, not eigenvalues.
- **CylindricalDecomposition**, **CylindricalDecompositionFunction** [Unrelated] real semialgebraic quantifier elimination; the report's domain is the complex cyclotomic ring and $\mathbb{F}_2$, not real inequality regions.
- **AsymptoticSolve** [Unrelated] perturbative series solutions; the report needs exact finite algebra.
- **NSolveValues** [Unrelated] numerical root-finding; the floating point the exact-cyclotomic approach is built to avoid (useful only as verification scaffolding).
- **PolynomialSumOfSquaresList** [Unrelated] real sum-of-squares positivity certificates; no contact with the complex cyclotomic amplitude arithmetic.

## What this suggests for the prototype (the payoff)

The current prototype hand-rolls $\mathbb{F}_2$ arithmetic with `PolynomialMod[., 2]` plus a manual $x^2 \to x$ rewrite, and brute-forces the residual sum. WL's native toolkit offers four concrete upgrades, all grounded in the functions above:

1. **Replace the hand-rolled Boolean reduction** with `GroebnerBasis`/`PolynomialReduce` against the ideal $(x_i^2 - x_i)$ at `Modulus -> 2` (or with `FiniteField[2]` arithmetic). This is also the literal realization of the report's named "sum over varieties / ideals" generalization.
2. **Do the $[HH]$ split with `Collect`/`Coefficient`** and represent $\varphi$ by its `CoefficientRules`, removing the manual monomial bookkeeping.
3. **Generalize the character sum to $GF(2^k)$** via `FiniteFieldElementTrace` ($\chi = (-1)^{\mathrm{Tr}}$), the native form of the Gauss sum.
4. **Carry the formal sum with `Inactive[Sum]` + `HoldPattern` rules + `Activate`**, marking the path-sum head `HoldAll` and `NonThreadable`; do the over-count cap with `MatrixRank` over `RootReduce`d $\mathbb{Z}[\omega]$ entries, and the all-Clifford fold with `RowReduce`/`NullSpace`/`LinearSolve` at `Modulus -> 2`.

The one honest non-result: `RootSum` is not a closed-form evaluator for the quadratic Gauss sum, so that part stays with the report's explicit identities.

---

*Source pages saved in this folder (`<Name>.md`), one per function, fetched from the installed Wolfram release via the documentation resolver. Headline behavioral claims kernel-verified.*
