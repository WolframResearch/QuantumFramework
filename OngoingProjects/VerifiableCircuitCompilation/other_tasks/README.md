# DQI arithmetic and DFTHC mixing feasibility

This local research compares arithmetic fragments and component mixing. It
does not implement a new DQI decoder or replace the inverse-square-root
incumbent. No global-optimality or leaderboard score is claimed.

## Full unchanged baseline profiles

| Task | Official objective, lower is better | Full default verification |
|---|---:|---|
| DQI | 24,439,458,623,294.45 | PASS, 102,400 shots, zero errors |
| DFTHC-SOS | 3,355,599,588.3155365 | PASS, 3,976 addresses, 64 shots each |

DQI has 1,893 qubits, 5,725,107.220107422 mean executed Toffolis and
5,725,461 static Toffolis. Even zero observed errors gives a nonzero
confidence bound: epsilon-plus 0.0008660587301657926 and repetition factor
1.1912604648115832, at the unchanged confidence setting delta=2^-128.

The native flame chart attributes about 36.84% of static Toffolis to forward
EEA, 36.84% to reverse EEA and 25.23% to Chien/Forney. SVG widths are rounded
display values, so these stage percentages are approximate. Exact displayed
group hit counts identify 58,211 degree-12 multiplications, each 63 Toffolis:
64.0524% of all static Toffolis. Static counts include classical conditions;
they are not interchangeable with the measured mean.

DFTHC uses 93,237 Toffolis and 730 total qubits. Angle loading contributes
79,480 Toffolis; the other 13,757 come from fixed-model preparation and
Givens rotations. Shared temporary space is dominated by the preparation's
544 qubits, rather than the angle circuits' 11 temporary qubits.

## Exact DQI arithmetic comparison

The native field is GF(4096) in the polynomial basis defined by
`x^12+x^6+x^4+x+1`, with primitive element represented by integer 2.
`FieldArithmetic.wl` uses that exact field, finds compatible subfields and
constructs all evaluation and interpolation maps. A multiplication identity
is proved by equality of all 1,728 bilinear coefficient-tensor entries.
This is an algebraic certificate for arbitrary inputs, not a random sample.
Outputs are converted back to the native bit basis, so classical MBUC masks
keep their original meaning.

| Construction | Multiply Q / Toffoli | Controlled multiply Q / Toffoli | Inverse phase Q / Toffoli |
|---|---:|---:|---:|
| Native Karatsuba | 36 / 63 | 49 / 75 | 72 / 315 |
| Quadratic over GF(64) | 54 / 63 | 67 / 75 | 90 / 315 |
| Cubic interpolation over GF(16) | 48 / 45 | 61 / 57 | 84 / 225 |
| Quartic interpolation over GF(8) | 46 / 42 | 59 / 54 | 82 / 210 |

The GF(8) route uses seven small products, each implemented by a six-product
identity proved symbolically modulo `x^3+x+1`. Linear maps use CNOT/SWAP;
all nonlinear gates, scratch qubits, measurements and cleanup are included
in these compiled fragment counts. Constant-field and Frobenius transforms
are linear. The inversion keeps the native five-multiplication Itoh-Tsujii
structure and zero-Toffoli measurement cleanup. Zero maps to zero by convention;
this is an extension of inversion, not a multiplicative inverse of zero.

Each of the twelve compiled fragments was simulated on 8,192 basis shots,
covering every field element, independent partner words, dirty XOR targets,
both control values, input restoration, clean scratch and correct phase.
Additional tests cover zero and multibit phase masks. These tests sample
measurement branches; they do not exhaust every possible branch.

A conditional whole-decoder substitution model gives factors 0.9812,
1.2086 and 1.2581 for the GF(64), GF(16) and GF(8) routes respectively.
It assumes unchanged repetition factor and an added global scratch peak of
18, 12 or 10 qubits. It uses static multiplication counts. No modified
decoder has been compiled or passed the full verifier: the actual global
peak, executed counts and accuracy remain integration questions.

## DFTHC component mixing

For each outer factor, the component index can be rotated by a real
orthogonal matrix. Lowercase weights and native identity offsets must rotate
together. `ComponentMixing.wl` proves invariance of the bilinear two-body
tensor, the offset contraction in the one-body operator, and the squared
offset norm entering the SOS energy. The decoded vectors, fixed symmetry
shift and factorization dimensions remain unchanged. These identities hold
for arbitrary tensor entries; pair rotations generate the tested orthogonal
transformations.

The LCU normalization depends on absolute coefficient sums and is not
invariant. On every fixed sign region, a pair objective is

```
constant + A*cos(2*theta) + B*sin(2*theta).
```

WL derives this form. `component_mixing.py` evaluates every sign boundary and
stationary point in a pi/2 period for each attempted pair. Floating evaluation
is a numerical instantiation, not an exact floating-point optimality proof.
The bounded search performs two rounds of 60 derivative-prioritized pairs
per outer factor: 7,800 attempted rotations, 5,944 accepted. It is not a
global optimizer. A separate spectral rotation serves as a control.

| Mixing | Raw Lambda | Full-resource cost | Improvement |
|---|---:|---:|---:|
| Identity | 308.4559397331 | 3,355,599,588.32 | 1.000000 |
| Spectral control | 334.7624281927 | 3,496,650,989.27 | 0.959661 |
| Pair coordinate | 300.7846945186 | 3,313,336,907.31 | 1.012755 |

The best tested mixing passes the full official verifier with unchanged angle
circuits and 254,464 address shots. Its effective normalization is
48.68043460482008. One- and two-body errors remain approximately 0.0031852
and 0.0186396 Ha. Relative augmented-Gram drift is 3.22e-16, and all official
thresholds, fixed shifts and float32 reference conventions are retained.
The normalization floor applies to raw Lambda, not lambda-eff. Candidate
factorizations are stored only in ignored `local/`; no competition source
or integrals are vendored.

## Required regimes

| Regime | DQI | DFTHC |
|---|---|---|
| General symbolic | Complete bilinear tensors; exact field embeddings | Arbitrary component tensor/offset invariants |
| Exact special | GF(8) six-product identity and all field inverses | A two-component example gives exact 2x LCU-H2 gain |
| Limit/decoupling | Zero control, zero classical phase mask | Identity rotation and limiting sign flip |
| Independent numeric | Compiled native-field simulation with dirty targets | Dense pair-search references, Gram tests and official metrics |
| Failure/edge | Unsupported/inexact degrees and zero inverse convention | Complex/malformed rotation rejection; spectral control worsens cost |

The relevant quantum invariant is a reversible XOR map or a diagonal phase,
with restored inputs and clean scratch. Hamiltonian Hermiticity is relevant
to DFTHC: real linear combinations preserve symmetry, and the reconstructed
two-body tensor and one-body operator are preserved by the proved identities.
No new quantum Hamiltonian evolution is performed by the arithmetic probes.

## Reproduction and provenance

Set `QF_INVSQRT_DATA` to the external competition folder. From the parent
project directory, in the pinned kickmix environment:

```sh
python other_tasks/profile_baselines.py dqi --data "$QF_INVSQRT_DATA" --workers 1
python other_tasks/profile_baselines.py dfthcsos --data "$QF_INVSQRT_DATA"
WolframKernel -noinit -noprompt -script other_tasks/run-symbolic.wls
WolframKernel -noinit -noprompt -script other_tasks/verify-symbolic.wls
python other_tasks/field_probes.py
python other_tasks/profile_stages.py --data "$QF_INVSQRT_DATA"
python other_tasks/component_mixing.py --data "$QF_INVSQRT_DATA"
python -m unittest other_tasks.test_feasibility -v
```

Reports record input/source SHA-256 hashes. `dqi-baseline.json` and
`dfthcsos-baseline.json` contain actual full runs. `field-results.json`
contains actual fragment resources, `dqi-profile.json` separates estimates,
and `mixing-results.json` contains official metrics and a full best-mixing run.
Run durations are observations, not general speed benchmarks.

Installed WL 15 documentation was read for `FiniteField`, `FiniteFieldElement`
and `PolynomialExtendedGCD`, including custom irreducible polynomials,
coefficient/index enumeration and finite-field linear algebra. Live probes
confirmed `Inverse[...,Modulus->2]` and exact finite-field matrix inversion.
The pinned kickmix implementation at commit
`638777fac4379bc62e2a1e4a52cf1a1753895d48` and its live method documentation
were inspected for multiplication, inversion scaffolds, phase conventions,
Frobenius maps, cleanup and native flame-chart profiling. No new QF APIs are
used. Bare wolframscript is not configured on the exercised machine; the
documented direct WolframKernel invocation was used instead.
