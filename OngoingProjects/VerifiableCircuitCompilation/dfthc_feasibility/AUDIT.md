# Released designs: applicability audit

Audit date: 2026-10-24. Published counts below are attributed source claims,
not local competition-verification results. No external source is vendored.

## Original released circuits

Pinned [andigu/uncompute](https://github.com/andigu/uncompute/tree/64b5a13aa249f1360846657392e6c66817b26d0f)
at commit `64b5a13aa249f1360846657392e6c66817b26d0f` (Apache-2.0).
The task READMEs, discovered-design READMEs, and DQI circuit/constants source
were inspected for contracts and construction; this is not a reproduction of
their large sampling experiments. Exact fetched blob IDs are in
`released-sources.json`.

| Design | Released contract | Competition implication |
|---|---|---|
| Coulomb interpolator | Six 8-bit coordinates, squared-distance computation, a 25-bit output accumulator, and error allowance 32*2^-24; source reports 1140 Toffolis, Q142. | Our task phases a 22-bit integer over the official domain at a different tolerance. Reuse normalization, fitting and erasure ideas; the published coefficient table and score cannot be substituted directly. |
| DQI truncated BM | Modulus 4331, fixed SHA-derived classical mask, quantized error-weight law, 100000 tests with 100 failures permitted. Source reports 442397 Toffolis, Q1837. | Our modulus is 4179, the mask is a classical input, and the verifier uses a different sampling/confidence penalty. The old reported score is not transferable. |
| Chemistry seven-frame walk | A 20-orbital Fe2S2 Pauli SELECT/alias table; source reports 14978 Toffolis, Q126. | Our task is a 76-orbital FeMoCo DFTHC angle loader, with fixed shifts and a full block-encoding cost. Frame compression does not directly replace its Givens-angle table. |

The original DQI source has a 42-product field multiplier, already the count of
our GF(8) backend. Its two-level norm-tower inverter costs 66 nonlinear products
in its basis, versus our integrated five-multiply inverse's 210. That is an
arithmetic lead, not a certified replacement in our basis or workspace. The
BM/reciprocal-derivative polynomial design offers much larger potential changes,
but its approximate cutoffs and one-entry history require our verifier's tests.

## Official starter kit

[Starter kit](https://www.kaggle.com/datasets/sohier/verifiable-circuit-compilation-starter-kit),
version 4, release note revision `f02280140eb100d604a077f4c774f85d68e1269f`.
The public manifest has eight research baselines for different systems, with
only one FeMoCo76o113e file. The seven other molecules/orbital spaces are not
candidate submissions and were not downloaded or numerically tested. The
FeMoCo `(65,60,60)` research NPZ has exactly the same five parameter arrays as
the provided competition sample; its stored training diagnostics are not a
separate better competition baseline.

Read the factorization starter, solution IO, training objective and loop,
system registry, and Givens-tree/angle encoder contracts. The competition
training starter freezes `Hsym_lower` and `H1_I`, uses the two strict floors,
and recommends a 10 mHa training target because the generic error penalty can
drift beyond the 20 mHa gate. Its shape adaptation copies leading slices into
randomly initialized arrays: this does not certify a compressed fit. The
[research paper's](https://arxiv.org/html/2502.15882v1) FeMoCo `(15,57,19)` used a different accuracy/optimization
contract; under our fixed-shift Frobenius criterion its rank cap 285 is
numerically ruled out by the present spectrum.

## Directly relevant public competition release

[VQCC FeMoCo Angle Words](https://www.kaggle.com/datasets/dmitriigluzdov/vqcc-femoco-public-factorization),
version 2, Apache-2.0. NPZ SHA256:
`9afd0ca5e90dd1740fe713e4ac1c26736184ea9bc4d4b57a6eeefafd608ee193`.
It supplies `(65,60,30)` weights/offsets and 3976*75 ten-bit angle words.
We independently reproduced its claimed metric acceptance with the unmodified
official evaluator: error_1=0.01391173198501434 Ha,
error_2=0.019980679640672218 Ha, Lambda=290.7542054993138,
E_gap=3.9159356215095613, lambda_eff=47.55854233515931.
The 20 mHa two-body margin is only 0.0193203593 mHa. Its gap is only
0.0010376215 Ha above the strict floor. Preserve these words/factors until
any proposed change is evaluated again.

Credits stated by the release: k256.dev baseline refinement; kozykappa ensemble
compression and offsets; Artem1410/anazemcev normalization tuning. This is
external public work, not an improvement discovered by this project.

[Quantum Circuit Compilation from Source](https://www.kaggle.com/code/dmitriigluzdov/quantum-circuit-compilation-from-source),
public version 8, kernel ID 140416287. Retrieved notebook JSON source SHA256:
`e4793278e7fd7f6828cb38a9d0983d2061c91a101e68a9fbb19c3ec44ae962a4`.
The audit read the pipeline, DFTHC loader/gauge/erasure dependency closure,
DQI field/mask adapter entry points, and results/provenance description.
It did not independently reproduce the entire notebook or all three tasks.

This newer notebook explicitly ports BM to the competition modulus and
arbitrary classical records. Thus the original paper's incompatibilities
have already been addressed in a public candidate; it deserves a future
standalone audit and full verification before we redo that port ourselves.
Its claimed full-run results are not imported as our verified score.

For DFTHC it combines signed Givens words, three nine-bit difference banks,
four-address Boolean interpolation, truncated ANF loading, five-wire pebbling,
and joint measurement phase correction. We ran the pinned released modules
in memory with fixed public factors, without column-gauge optimization.
The signed gauge changes vector components by at most 7.78e-16 up to row sign;
the official metric gates still pass. All 25 circuits compile within Q128.
Static resource projection: angle Toffolis 56450 (both passes), total Toffolis
65244, total block-encoding Q661, conditional cost 2051023203.3721037,
or 1.615455594 times cheaper than our incumbent DFTHC result.
This is an audited resource projection, **not full circuit verification**.
Measurement phase repair, arbitrary dirty target words, chained batches and
the official forward/reverse shot coverage remain required before adoption.

The public credits assign the signed Givens/four-address/five-wire design to
kozykappa, banked lookup to Sumit Chongder and kozykappa, and ANF/gauge work
to Artem1410/anazemcev. These credits and Apache-2.0 notices must accompany
any later port. No source is copied into the project during this audit.
