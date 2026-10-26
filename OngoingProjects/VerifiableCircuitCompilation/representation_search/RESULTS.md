# Representation search: mathematical results, no new circuits

The useful change is to optimize the represented function before synthesizing
its gates. Three candidates survive algebraic tests. None has a compiled cost
or a new competition score. Existing circuits, coefficient files, factors and
verification evidence are unchanged. Nothing from this search has been pushed.

## 1. Compensated curvature: a passing smaller integer model

For an interval of width h, write P(t)=A-Bt+Ct², 0≤t≤h. Changing C to C+Δ
while also changing A to A+Δh²/8 and B to B+Δh gives

    P_new(t)-P(t) = Δ(t²-ht+h²/8).

Its maximum absolute value is |Δ|h²/8, eight times smaller than the bound
|Δ|h² from changing C alone. This is the minimax affine compensation: for
e(s)=s²+u s+v, e(0)-2e(1/2)+e(1)=1/2. If |e|≤E, the triangle inequality
gives E≥1/8. The proposed residual attains this bound with alternating signs.
This is a standard minimax approximation argument, not a novelty claim.

With h=1/64, the original frozen 36-bit table generates a new candidate table.
The integer model reproduces the two nested, individually rounded Horner
stages and the phase scaling, rather than evaluating a real polynomial.
All 3,133,452 legal x pass with fractional widths (30,23,10), physical widths
(30,22,9), and worst error/tolerance 0.8505589342620424. Both intermediate
accumulators remain positive and fit their registers. The current verified
design uses (30,22,15): this candidate removes six coefficient bits from the
model. It does not establish a six-qubit reduction in peak circuit width.

The 8-bit and 9-bit fractional curvature candidates fail. A second finite
search chooses the integer intercept from the intersection of the exact
per-input admissible intervals, at B offsets -32 through +32, with C fixed
at its compensated 9-bit value. It repairs five rows but cannot satisfy row
170. That finite family is rejected; it is not an impossibility proof for
all coefficients. New tables stay in ignored local/. A future circuit must
also budget rotation-synthesis error from the remaining tolerance.

## 2. DQI: cache a quotient polynomial, then kick back a quadratic phase

Let σ be a nonconstant squarefree locator with nonzero roots. Put

    M = z^h σ' mod σ,
    Ψ = M^(-1) mod σ,
    Π = Ω Ψ mod σ.

Then Π(β)=β^(-h)Ω(β)/σ'(β) at every locator root β. PolynomialExtendedGCD
constructs Ψ exactly. This is a polynomial certificate, not a fit. In the
released decoder n=70, m=n+3=73, t=35; evaluator coefficients occupy indices
37 through 71. With z=primitive^j, the locator is Σ b_k z^k, the numerator is
Ω(z)=Σ w_k z^(71-k), and the odd locator slots accumulate zσ'(z). Thus the
direct native adapter reverses the evaluator slots and uses h=1. A direct
test compares the cache with both of these released scan sums. The first
independent review caught and rejected an incorrect ascending-slot/h34
adapter; its counterexample is retained. The abstract h34 identity by itself
had not tested the native convention.

For each classical measurement mask μ_j, evaluate Π at the fixed field point
using linear maps. The phase is (-1)^(χ_j L_μj(Π(β_j))), where χ_j says σ(β_j)=0
and L is the required binary linear functional. After caching, this is a
quadratic Boolean phase, implementable with classically selected CZ gates.
It avoids the per-position field inverse and the 12 root-masked evaluator ANDs.
The trace Gram matrix is invertible in the native GF(2^12), confirming that
all 4096 masks admit the trace-functional representation. The locator tests
still cost work; the cached polynomial must be created and erased reversibly.

Caching σ'^(-1) modulo σ already appears in the
[published DQI design](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html).
The candidate tested here caches the complete product, including the native
position factor. No claim is made that this is a new decoding identity.

The static scan removes at most 4095(5·42+12)=909,090 Toffoli-like operations
from the present 4,503,030. This ceiling excludes setup, erasure and all new
workspace. At unchanged width, setup plus erasure must cost less than that.
Under a T Q² comparison, even with free setup, total width must stay below
about 2130 qubits, versus the current 1903. Simply charging an additional
35·12=420-qubit polynomial block would lose. In-place conversion or proven
reuse of live registers is consequently a prerequisite, not a detail.

For monic degree d, ordinary product plus reduction uses at most 2d²-d
field products in an arithmetic expression. At d=35 this is 2415; four such
expressions at 42 operations/product give 405,720 before the inverse setup,
normalization, control and reversible scratch. This is an operation budget,
not a constructed quantum algorithm or a certified resource estimate.
The no-error constant locator needs a separate zero-phase branch. Repeated
roots and zero roots with h>0 are outside the inverse's domain and rejected.

## 3. DQI: projective coefficients instead of normalized coefficients

A more structural second candidate is to represent the locator only up to
nonzero scale. Scaling (σ,Ω) to (aσ,aΩ) leaves roots and Ω/σ' unchanged.
Normalization divisions can consequently be removed from the classical
locator recurrence. This is the established inversionless Berlekamp–Massey
idea, not a new algorithm; see the scalar-rescaling discussion and Algorithm 2
in [the authors' paper](https://eprint.iacr.org/2022/469).

If the normalized update is C_new=C-(δ/b) z^m B, keep scaled C~=aC and B~=cB
with γ=cb and δ~=aδ. Then

    C~_new = γ C~ - δ~ z^m B~ = γ a C_new.

The major update B~_new=C~, γ_new=δ~ preserves the same relation to the next
normalized step; the minor update retains B~ and γ. The recurrence in
ProjectiveDecoder.wl has no division. Exact tests in the native GF(4096)
recover a degree-35 locator from 70 synthetic syndromes, its 35 Forney values,
all intermediate degree/nonzero-scale invariants, plus one-error and zero-error
cases. A general symbolic identity proves the rescaled update and phase ratio.

This replaces the mathematical state representation of the Euclidean decoder.
It is not a reversible decoder implementation: the actual native syndrome
adapter, quantum control history, live register schedule, and reverse update
costs are unresolved. The major polynomial update matrix is
{{γ,-δ~ z^m},{1,0}}, with determinant δ~ z^m; the minor matrix has determinant γ.
The major update has a formal inverse over rational functions, requiring
division by δ~ z^m. This is not a fixed-register reversibility proof. Truncating
both polynomials to two coefficients with m=1 already gives a singular
four-coefficient map; a correct implementation must preserve the needed tails
or prove their promised zeros and retain the necessary branch information.
This qualification and counterexample were added after the first quality
review. Implementing inverse updates may also require field inversions.
Classical absence of division does not imply a quantum Toffoli saving. No
projected score is assigned to this candidate.

## 4. DFTHC: two tempting basis changes are exhausted

Every one of the 25 fixed C30 batches has 994 legal high-address rows and
983 nonzero nonlinear address-monomial coefficients. The 111 output bits
have rank 108, exactly reflecting three redundant sign bits; the 81
difference-bank bits have full rank 81 over GF(2).

For an ANF coefficient indexed by m, a_m is the XOR of truth-table rows
s⊆m. If m<994, all such s also satisfy s<994. Thus changing invalid rows
994..1023 cannot change these coefficients. The present completion already
zeros coefficients outside that prefix. An invertible output basis maps
nonzero coefficient vectors to nonzero vectors. These transformations cannot
reduce the support. WL checks the triangular Möbius identity and Python
measures the actual tables. BooleanMinimize expression length does not price
a reversible shared-AND circuit.

This rejects two approaches within the address-monomial representation.
It does not prove that the table is incompressible or that 983 AND gates are
necessary for a general factored Boolean circuit. A tensor decomposition of
the Boolean map, or a new physical factorization with fewer distinct vectors,
remains an open search. This stage has not found a better DFTHC representation.

## Regimes and checks

| Regime | Evidence |
|---|---|
| General symbolic | Curvature identity and minimax certificate; quotient ideal identity; projective update/phase identities; subset proof for arbitrary prefixes |
| Exactly solvable | Degree-one GF(16) quotient and one-error projective recurrence; three alternating curvature extrema |
| Limit | Δ→0 restores P; no-error constant locator is rejected by the quotient API and has constant locator in the projective recurrence |
| Numerical/reference | Exhaustive official inverse-square-root reference; native-field degree-35 quotient versus released scanner; synthetic 35-error projective decoder; all 25 fixed DFTHC batches |
| Failure/edge | Too-coarse curvature; finite affine search rejection; repeated/zero roots; malformed exponent; empty locator; invalid zero projective scale; singular truncated major update |

The physical invariant is equality of the kicked-back phase, including arbitrary
classical masks. Diagonal phase signs have unit modulus. No compiled circuit,
workspace-cleanup, Hamiltonian-fit, speed benchmark or new score is asserted.

## Reproduce

Use a Wolfram 15 kernel to run `run.wls`. Built-ins and signatures were checked
against [PolynomialExtendedGCD](https://reference.wolfram.com/language/ref/PolynomialExtendedGCD.html),
including its FiniteField example, and
[PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html).
Python requires NumPy, the existing project generator dependencies, and the
official competition verifier bundle supplied as DATA. From this directory:

    python curvature_screen.py --data DATA
    python affine_rounding.py --data DATA
    python table_screen.py --words ../dfthc_c30/local/words.u16

Only aggregate reports are public candidates. Raw derived tables and factor
inputs remain ignored. Current verified combined local score: 2.0088882102784944,
unchanged. The next implementation candidate supported by this search is the
compensated quadratic. Both DQI representations need a viable reversible
schedule before implementation is chosen. A larger DFTHC gain is unresolved.
