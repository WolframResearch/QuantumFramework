# GF(8) DQI decoder integration

This generator compiles the external, unchanged Proos-Zalka reference decoder
using the exact native-basis GF(8) multiplication plan from `other_tasks`.
The EEA iteration count, syndrome and classical measurement registers, native
field, Forney phases, decoding algorithm and verifier distribution are fixed.
The inverse-square-root and DFTHC artifacts are unchanged. No push is part of
this stage.

## Full official result

The default full official verifier passes all 102,400 shots with zero value
errors and zero phase errors. These are actual measurements, replacing the
earlier conditional substitution model.

| Quantity | Unchanged baseline | GF(8) decoder |
|---|---:|---:|
| Official cost | 24,439,458,623,294.45 | 19,425,188,366,155.996 |
| Qubits | 1,893 | 1,903 |
| Mean executed Toffolis | 5,725,107.220107422 | 4,502,782.126171875 |
| Static Toffolis | 5,725,461 | 4,503,030 |
| Total operations | 86,051,014 | 117,029,650 |
| Reaction depth | 170,757 | 174,852 |
| Repetition penalty | 1.1912604648115832 | 1.1912604648115832 |

The measured improvement is **1.2581323878369532x**, a **20.517108559676135%**
cost reduction. More Clifford operations and ten additional qubits are the
tradeoff for fewer Toffolis. Wall-clock speed is not the scoring objective,
and the recorded run duration is not a portable speed claim.

All 58,211 static multiplications are replaced: 32,408 controlled products,
888 division scaffolds with six products each, and 4,095 inverse scaffolds
with five products each. Native zero-Toffoli cleanup and phase masks remain.
The actual peak and static Toffolis match the feasibility model; the actual
executed mean comes from the candidate's own hash-derived official samples.

The complete terminal report is `verified.json`; `build.json` describes the
earlier build stage and deliberately does not certify verification. This is
not a new submitted Kaggle score.

## Construction and contract

`GF8Builder` is a narrow adapter for the three computation entry points used
by the reference: XOR multiplication, scaffolded inversion and scaffolded
division. It delegates other operations to an unpatched kickmix builder.
The degree-12 multiplication is quartic interpolation over GF(8): seven
subproducts, each using the exact six-product identity. Binary input/output
maps preserve the native polynomial basis `x^12+x^6+x^4+x+1`.

The inverse uses the native five-multiplication Itoh-Tsujii schedule. Its four
ordered scaffold words hold powers `{3,15,1020,2046}` of the original input;
the output holds power 4094. Native `del_gf2_inverse_with_scaffold` and
`del_gf2_div_with_scaffold` can therefore clear these same values without
rebuilding the chain. The native multiplication measurement cleanup also
remains in use. No alternative field basis is exposed to classical phase
masks. Zero inversion is extended to zero, matching native behavior.

`run-certificates.wls` independently regenerates the selected multiplication
plan and the exact Frobenius maps consumed by the Python compiler. It proves
the symbolic chain-composition identity, scaffold exponents, inverse exponent,
all nonzero inverse values, map invertibility, inverse maps and field period.
The compiler checks certificate/source bindings before construction.

The adapter requires disjoint native-width registers in the supported
decoder use. It is not a complete replacement for every kickmix arithmetic
method. The reference source is imported from the external data bundle and
must match the previously profiled SHA-256. No reference or verifier source
is copied into this repository.

## Reproduce

Export `QF_INVSQRT_DATA` with the external competition bundle path. From the parent
project directory, with the pinned kickmix environment and a direct WL kernel:

```sh
WolframKernel -noinit -noprompt -script dqi_gf8/run-certificates.wls
python -m unittest dqi_gf8.test_integration -v
python dqi_gf8/generator.py --data "$QF_INVSQRT_DATA"
python dqi_gf8/verify.py --data "$QF_INVSQRT_DATA" --workers 1
```

The full runner fixes the official defaults: m=4095, n=70, degree=12, r=2016,
102,400 shots, batches of 256 and log2(delta)=-128. The fixed salt matches
the baseline; each circuit's official hash still determines its own samples.
The actual mean executed Toffoli count, global qubit peak and confidence
penalty determine the cost, rather than static fragment substitution.
Reports pin the circuit, sources and official inputs. Durations are observed
run times, not speed benchmarks. The generated binary is in ignored `local/`.

The runtime uses WL 15.0 and unpatched kickmix 0.1 at source commit
`638777fac4379bc62e2a1e4a52cf1a1753895d48`. Installed `FiniteField` and
`FiniteFieldElement` documentation and the native inversion/division source
were consulted for polynomial representation, matrix inversion and the
ordered scaffold contract. No QuantumFramework API is used in this adapter.

## Regimes and evidence

| Regime | Evidence |
|---|---|
| General symbolic | Complete multiplication coefficient tensor and arbitrary positive-integer chain composition |
| Exactly solvable special | GF(8) six-product identity; every native field element through inversion; exact scaffold exponents |
| Limit/decoupling | False quantum-control branch and disabled classical condition; zero inverse extension |
| Independent numeric | Native GF arithmetic versus compiled XOR, Frobenius, division and phase results; small and full official decoder protocols |
| Failure/edge | Wrong field/width, invalid target, singular map and invalid decoder dimensions rejected |

The six integration tests cover 8,192 multiplication shots for each control
mode, all 4,096 inverse inputs with independent classical masks, every
divisor with both enabled/disabled branches, all six Frobenius maps on every
field word, and a 2,048-shot reduced decoder. They assert restored inputs,
correct phase and clean scratch. Measurement branches are sampled rather
than exhaustively enumerated. Relevant quantum invariants are reversible
XOR behavior and the intended diagonal phase, with restored syndrome and
scratch registers; no Hamiltonian evolution or density operator is built.

The review graph is certificates -> integration checks -> full official run
-> fresh wl-verify -> fresh wl-quality. Each review permits at most three
fresh rounds. Any executable quality fix returns to correctness and, when it
changes the generated circuit, full verification. This is a code deliverable;
there is no notebook build.

## Completed review gates

- [wl-verify, round 1](reviews/correctness-1.md): zero open issues. The full
  candidate was independently rebuilt byte-for-byte and its complete
  102,400-shot official run reproduced the exact cost and zero errors.
  Additional polynomial, inverse-scaffold, mask and decoder-boundary checks
  passed. The baseline comparison uses the retained historical full report.
- [wl-quality, round 1](reviews/quality-1.md): all seven axes passed, zero open
  issues. Documented APIs, functional WL, exact derivations, nontrivial
  quantum invariants and all five regimes were checked with a fresh kernel.
  Full-run records were audited in this quality review rather than rerun.

[evidence-manifest.json](evidence-manifest.json) binds the final local source,
results and redacted reviews. Private-path scanning is clean; the circuit is
ignored and external competition inputs remain external. Nothing was
committed or pushed.
