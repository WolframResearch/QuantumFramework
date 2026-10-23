# Inverse-square-root phasing: QF / kickmix prototype

The subsequent [WL-backed architectural comparison](feasibility/COMPARISON.md)
compares better Newton seeds, direct polynomial phasing, and phase propagation.
Its exact certificates, resource counts, full-domain integer models, and
independent review records are separate from the verified incumbent below.
The reviewed choice is direct quadratic phasing on 192 intervals, with a predicted
inverse-square-root improvement of 1.99 times the supplied baseline. Both fresh
WL reviews passed; full repeated-shot verification of that new design remains.

This experiment searches the precision parameters of the supplied Newton–Raphson
circuit generator. The competition verifier and generator are imported unchanged
from the downloaded bundle. QF independently checks small arithmetic and
measurement/uncomputation fragments. The precision search itself uses Python and
the host's integer model; a measured improvement is not evidence that QF caused it.

## The contract

For every integer `1 <= x <= 3,133,452`, the submitted circuit must implement
`|x>|0> -> exp(2 pi i / sqrt(x)) |x>|0>`. The only named register is `x`,
containing 22 distinct little-endian qubits. Input restoration, scratch cleanup,
and phase are checked on **64 shots per x: 200,540,928 shots**. Allowed phase
error is `2^-24 / sqrt(x)` turns, modulo one turn.

Cost is `Q * (Toffolis + T/4)`, including the verifier's rotation-synthesis price.
The price depends on the smallest remaining phase-error margin. Depth is
reported but is not the objective. These statements come from the downloaded
`verifiers/invsqrt/verify_submission_invsqrt.py`, not from an inferred interpretation
of the competition title.

## What runs

- `experiment.py baseline` verifies the supplied binary baseline.
- `experiment.py search` samples 1,280 precision settings, orders survivors by
  a resource proxy, then checks the first 24 over the entire input domain using
  the verifier's exact uint64 phase/tolerance arithmetic. This is a bounded search,
  not a global optimality claim. It writes `candidate.kmb` and screening reports.
- `experiment.py candidate` runs the unmodified verifier at full benchmark size.
  A `MODEL_PASS` alone is never presented as a circuit `PASS`.
- `fragments.py` exports real kickmix AND/phase/uncompute and adder fragments.
- `KickmixBridge.wl` maps a restricted instruction subset to QF branch operators.
- `verify-qf.wls` checks phases, coherence, Kraus normalization, scratch cleanup,
  little-endian arithmetic, symbolic Newton error, and negative controls.
- `test_experiment.py` checks contract rejection and model/simulator agreement.

The bridge supports small fragments (at most eight qubits), not full benchmark
simulation or general kickmix import. `HMR` takes explicit outcomes and retains
the `1/sqrt(2)` Kraus factors and classical corrections. `R` is represented by
the clean-wire projector `|0><0|`; it asserts a clean ancilla, rather than
implementing an arbitrary physical reset. Supported gate phases use kickmix's
half-turn convention. Other instruction families are rejected.

## Environment

The exercised runtime is Python 3.13.1, NumPy 2.5.4, kickmix 0.1,
Wolfram Language 15.0, and this checkout's QF 2.1.3. Google Qualtran/kickmix source
is pinned to commit `638777fac4379bc62e2a1e4a52cf1a1753895d48`.
The macOS ARM build used Homebrew GCC 15.2 and the build dependencies in
`requirements-build.txt`; Apple's installed Clang 15 failed to compile the
upstream C++20 aggregate initialization. The downloaded Linux `kickmix.so` is
not used on macOS. Source code was not patched.

Install the build requirements into an isolated environment, then install the
pinned official repository's `kickmix/` directory with
`pip install --no-build-isolation /path/to/Qualtran/kickmix`.
Activate that environment before running the commands below.
Competition data remains external; JSON reports record SHA-256 hashes of the
baseline, generator, and verifier. No verifier or competition data is vendored here.

From this directory:

```bash
QF_INVSQRT_DATA=/path/to/extracted/competition
QF_INVSQRT_PYTHON=python3
"$QF_INVSQRT_PYTHON" experiment.py baseline --data "$QF_INVSQRT_DATA" --workers 4
"$QF_INVSQRT_PYTHON" experiment.py search --data "$QF_INVSQRT_DATA"
"$QF_INVSQRT_PYTHON" experiment.py candidate --data "$QF_INVSQRT_DATA" --workers 4
"$QF_INVSQRT_PYTHON" fragments.py
export QF_INVSQRT_DATA
"$QF_INVSQRT_PYTHON" -m unittest -v test_experiment.py
WolframKernel -noinit -noprompt -script verify-qf.wls
```

The WL runner registers existing Arrays, TensorNetworks, IBMQuantumPlatform,
and ImageMetadataTools installations for that kernel. It does not install or
silently replace missing dependencies. Fresh-kernel failures remain failures.

## Evidence and limits

Both full-size verifier runs passed with zero failures:

| Quantity | Supplied baseline | Candidate |
|---|---:|---:|
| Qubits | 260 | 256 |
| Toffolis | 7,735 | 7,327 |
| Fine rotations | 37 | 36 |
| Toffoli-equivalent cost | 7,987.593207 | 7,583.127187 |
| Objective `Q*C` | 2,076,774.233715 | 1,941,280.559862 |
| Verified shots | 200,540,928 | 200,540,928 |

This is a **6.524237% reduction in inverse-square-root cost**, or a task
improvement factor of **1.069796**. Holding the other two tasks at baseline
would give an aggregate improvement factor of **1.022744**.
The supplied generator also reproduces the baseline binary byte for byte in
this runtime, confirming the baseline calibration.

The only precision changes are `delta_bits: 14 -> 12`, `product_bits: 28 -> 27`,
and `result_bits: 29 -> 28`. Coefficient and square widths remain 17 and 28.
The candidate's worst observed phase error uses 92.51% of its allowed tolerance;
rotation synthesis is priced using the remaining margin. This is an initial
improvement to the reference algorithm, not a new inverse-square-root algorithm.

`results/baseline-verified.json` and `results/candidate-verified.json` hold full
verifier results. The salt is a fixed,
declared local validation salt. Kaggle uses its own secret salts: this local
run is not a leaderboard submission or a universal proof over measurement
outcomes. Candidate and baseline use the same verifier settings and shot count.

`results/search.json` retains rejected candidates, including failures that a
sample-only screen missed. `results/qf-verification.json` records the independent
matrix checks. Generated binary circuits are ignored by Git and can be regenerated
from the precision configuration in `results/candidate-model.json`.

Independent WL correctness review ran for its three-round limit and reproduced
the quantum/mathematical claims. Its last finding was that wrong argument shapes
stayed unevaluated; an explicit `Failure` fallback and regression checks fix it.
That final defensive fix was locally tested after the fresh-review limit, rather
than receiving a fourth independent round. The quality review passed all seven
axes. No change from these bridge reviews altered the verified candidate circuit.

Primary sources: [competition](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation),
[companion paper](https://arxiv.org/html/2610.23832v2#S3.SS1),
[official kickmix source](https://github.com/quantumlib/Qualtran/tree/638777fac4379bc62e2a1e4a52cf1a1753895d48/kickmix).
