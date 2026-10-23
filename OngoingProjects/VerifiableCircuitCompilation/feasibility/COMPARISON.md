# Inverse-square-root phasing: architectural feasibility

This is a comparison before production implementation. Core WL discovers and
certifies alternative arithmetic; the official kickmix compiler supplies resource
counts. Exhaustive integer-model checks and small circuit samples are separate
from the official full repeated-shot verifier. No new submission is claimed.

## Contract and incumbent

For every integer `1 <= x <= 3,133,452`, restore the 22-bit input and every clean
scratch qubit while applying `exp(2 pi i / sqrt(x))`. The phase error modulo one
turn must not exceed `2^-24 / sqrt(x)`. The objective is
`Q * (Toffolis + (base_T + synthesized_T)/4)`. The synthesis price depends on
the *minimum absolute* remaining phase-error margin, not the maximum relative
error alone. Source: the downloaded official verifier, whose hash is recorded in
resource-results.json. No other competition task is optimized here.

The supplied baseline costs 2,076,774.23. The incumbent costs 1,941,280.56 and
has already passed the full 200,540,928-shot local verifier. The new designs below
have not run that verifier.

The incumbent's actual compiled Toffoli counts are:

| Stage | Compute | Uncompute | Total |
|---|---:|---:|---:|
| Normalize | 120 | 110 | 230 |
| Interval selection | 3 | 0 | 3 |
| Coefficient load | 0 | 0 | 0 |
| Cubic Horner | 954 | 954 | 1,908 |
| Square seed | 442 | 442 | 884 |
| Multiply normalized input | 1,001 | 1,001 | 2,002 |
| Form Newton correction | 27 | 27 | 54 |
| Final multiply | 833 | 833 | 1,666 |
| Denormalize | 290 | 290 | 580 |
| Total | | | 7,327 |

Thus the Newton square and two multiplications account for 4,552 Toffolis,
62.1% of the incumbent. The existing compiler already uses in-place Horner,
column truncation and measurement uncomputation. The profile is disjoint and
reconciles with official operation_statistics. There are no classically
conditioned Toffoli or fine-rotation counts in these probes; the sampled
executed counts agree. HMR correction Cliffords are conditioned.

## Shared symbolic argument

Normalize as `x = 4^k m`, `1 <= m < 4`. Then the phase is `2^-k / sqrt(m)`;
relative approximation error survives this scaling unchanged. The reference
normalizer stores the equivalent unary shift `s = 10-k`.

On interval `m=a+t`, `0 <= t <= h=2^-r`, discover a local relative minimax
polynomial at 50-digit precision. Round its coefficients to signed integers
divided by `2^w`. Numerical minimax output is used only to propose a rational
bound eta. For every interval, exact Resolve proves

    p(t)>0,
    (1-eta)^2 <= (a+t) p(t)^2 <= (1+eta)^2,
    -1 <= p'(t) <= 0.

Since eta is between zero and one and p is positive, these are relative-error
certificates without numerical square-root comparison. The positive Horner
intermediates are also bounded away from zero and one, including an arithmetic
allowance. No finite WorkingPrecision is supplied to Resolve. A timeout never
counts as a certificate.

Each controlled partial-product row drops less than one output ulp. With
`R=min(d,20)-r` rows and degree L, an absolute Horner bound is

    H = R 2^-w sum(h^j, j=0..L-1).

Using only d<20 input fraction bits adds at most `2^-d` absolute polynomial
error, by the derivative certificate. Using all 20 adds zero. Multiplication by
sqrt(m)<=2 converts these to a conservative relative bound. Coefficient rounding
is already included in the polynomial certified by Resolve.

These bounds concern the specified arithmetic model. Finite-register capacities
are additionally checked exhaustively over the actual benchmark inputs, including
the Newton registers. They do not establish the entire generated measurement
circuit's action on all branches.

## A: better seed followed by one Newton step

For relative seed error e, symbolic expansion gives

    e_new = -3 e^2/2 - e^3/2.

The code bounds omitted square bit products with an exact rational sum, and the
subsequent multiplication row errors separately. It does not replace row
truncation with rounding a whole exact product. It considers two architectures:

| Seed | Local degree | Intervals | Coefficient bits | Input fraction bits |
|---|---:|---:|---:|---:|
| Quadratic | 2 | 12 | 24 | 16 |
| Linear | 1 | 96 | 24 | 16 |

The continuum bound chooses uniform Newton widths 36 and 32 respectively to
reserve roughly half the tolerance. These are conservative sufficient widths,
not proven minimum widths. Both are more costly than the incumbent. Additional
28-bit probes avoid unfairly rejecting this family solely because its continuum
bound is conservative: the quadratic still costs more; the linear fails at
exactly x=1,048,575. No impossibility claim is made for other partitions, sparse
coefficients, coefficient widths or mixed Newton widths.

## B: approximate the phase directly

The same certificates allow direct evaluation to the final phase tolerance,
eliminating all three Newton arithmetic stages and their registers. The tested
architectures are a cubic on 48 intervals and a quadratic on 192 intervals.
Both use 36 fractional coefficient bits and all 20 input fraction bits. The
192-interval choice trades a larger lookup for one fewer Horner multiplication.

All coefficients are loaded together with init_lookup into the coefficient
registers, then removed with del_lookup after reversing Horner. Lookup addresses
are the normalized input's top r+2 bits; the unused lower quarter of addresses is
padded, and that padding is charged. Both interval endpoints are covered by the
polynomial certificates, although the circuit assigns boundaries to one interval.

The certified total relative bounds are 2.9363036447e-8 (cubic) and
3.6824530980e-8 (quadratic), below 2^-24 = 5.9604644775e-8. At least 50.7% and
38.2% of tolerance respectively remains under these conservative bounds. The
comparison includes arithmetic scratch and carry ancillas allocated internally
by kickmix, not just the explicitly requested registers.

## C: propagate phases through the final multiplication

For bits a,b, the rational identity `ab=(a+b-(a XOR b))/2` converts a controlled
phase into three parity phases. For three bits,

    abc = (a+b+c-(a XOR b)-(a XOR c)-(b XOR c)+(a XOR b XOR c))/4.

The code expands each retained final-multiplication partial product into these
parities, combines equal supports, and reduces their coefficients modulo one
turn. This preserves the *truncated* product actually computed, rather than
silently restoring dropped terms. Boolean arithmetic and arbitrary-angle phase
arithmetic use distinct coefficient domains. The zero-state phase is preserved,
not discarded as a global phase.

With normalization exponent fixed to s=0, this construction needs 403 fine
rotations. That is a diagnostic, not a circuit for the full task and not a lower
bound on every possible phase synthesis. Handling all eleven coherent exponent
values with equality indicators gives 8,875 fine rotations. An inactive quantum
indicator does not make an emitted rotation free: the operations are applied
coherently, not skipped by a classical branch. At half the smallest tolerance
reserved for synthesis, this version costs about 16.40 million, compared with
1.94 million for the incumbent. This rejects this expansion as a first
implementation; it does not rule out other phase synthesis methods.

## Comparison and evidence levels

Objectives below use exhaustive-model minimum absolute margins and the verifier's
rotation price for A/B. C uses an explicitly hypothetical half-tolerance margin.
Counts are compiled circuit counts. Costs remain predictions until the official
circuit verification succeeds.

| Probe | Q | Toffolis | Fine rotations | Predicted Q*C | Model result |
|---|---:|---:|---:|---:|---|
| Incumbent | 256 | 7,327 | 36 | 1,941,281 | Full local verifier passed previously |
| A quadratic, Newton 36 | 300 | 10,497 | 44 | 3,237,997 | All inputs pass |
| A quadratic, Newton 28 | 260 | 8,361 | 36 | 2,244,985 | All inputs pass; no continuum 28-bit guarantee |
| A linear, Newton 32 | 256 | 8,256 | 40 | 2,182,382 | All inputs pass |
| A linear, Newton 28 | 236 | 7,180 | 36 | Ineligible | One input fails |
| B cubic, 48 intervals | 256 | 6,638 | 43 | 1,773,730 | All inputs pass |
| B quadratic, 192 intervals | 220 | 4,449 | 43 | 1,043,220 | All inputs pass |
| C, all exponents | 217 | 5,099 | 8,875 | 16,399,036 | Exact phase identities; full circuit unverified |

Each of six A/B variants checks all 3,133,452 inputs in the integer model.
None has a register-capacity violation. Small circuit/model cross-checks cover
81 distinct inputs with two seeded shots for the quadratic seed and 93 with two
shots for each other design. They check exact phase, restored x and clean scratch.
These are sampled compiler checks, not an all-input circuit verification.

The direct quadratic's predicted improvement is about 1.99 times the supplied
baseline and 1.86 times the incumbent. With the other two tasks at baseline,
its predicted aggregate improvement would be only about 1.26. This does not
match an aggregate leaderboard score of 8.43991.

## Required regimes

| Regime | Evidence |
|---|---|
| General symbolic | arbitrary relative error e, normalization scale, arbitrary real phase angle; parity identity tested for sizes 1 through 6 |
| Exact special | zero seed error; exact bit-product identities and symbolic diagonal unitarity |
| Limit | e_new/e^2 tends to -3/2 as e tends to zero |
| Numeric reference | high-precision polynomial samples; exhaustive bit-exact domain model; sampled real compiler circuits |
| Failure/edge | inaccurate polynomial rejected; inexact/empty certificate inputs rejected; row-versus-whole-product witness; one failing linear-seed input; modulo-one and zero-state phases |

## Reproduction and selection gate

Run from this directory:

```bash
WolframKernel -noinit -noprompt -script run-feasibility.wls
WolframKernel -noinit -noprompt -script verify-feasibility.wls
python3 resource_probe.py --data "$QF_INVSQRT_DATA"
```

DOCUMENTATION.md records the documentation read and tested signatures. WL
certificates and resource/model evidence are retained as JSON with source hashes.
No runtime/performance speedup is claimed.

## Reviewed decision

Choose **B: direct quadratic phasing on 192 intervals** for the next production
implementation. It has the lowest predicted objective among the tested designs,
exact continuum bounds, a passing exhaustive integer model, no register-capacity
violations, and matching sampled compiled circuits. It removes the Newton square
and both multiplications while paying the joint coefficient lookup and full
uncomputation costs. The conservative symbolic error bound also leaves positive
synthesis margin, independently of the exhaustive-model result.

The tested A designs provide no improvement over the incumbent when eligible.
The tested C expansion is dominated by its fine-rotation cost. These decisions
apply to these bounded probes; they do not establish optimality or exclude other
designs in either family.

Fresh [wl-verify round 1](reviews/correctness-1.md) and
[wl-quality round 1](reviews/quality-1.md) both finished with **OPEN ISSUES: 0**.
The correctness reviewer independently reproduced all 1,044 exact interval
certificates, the arithmetic bounds, complete resource/model results and phase
identities. The quality reviewer passed all seven axes, including symbolic-first
derivation, nontrivial examples, invariants, documented APIs and all five regimes.
Neither review required an executable-code change. Independent C samples also
matched the separate official reference model across all eleven exponents.

The next implementation milestone is to preserve these coefficients and error
certificates, build a standalone reproducible circuit generator, and run the
unmodified 200,540,928-shot official verifier with its actual resource accounting.
Until that succeeds, the 1,043,220 objective is a prediction and the existing
verified incumbent remains the submission candidate. No Kaggle score is claimed.
