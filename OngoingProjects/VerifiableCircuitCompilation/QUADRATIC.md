# Fixed quadratic phase circuit: full verification

The standalone direct quadratic design **passed the unmodified official local
verifier on all 200,540,928 shots**, covering every integer x from 1 through
3,133,452 on 64 shots each. Input restoration, clean scratch, and phase accuracy
all passed. The verifier's default 22-bit width, beta=1, relative tolerance
2^-24, and batch size 256 were retained. The run used four workers and the same
declared local salt as the earlier baseline and Newton incumbent.

| Quantity | Supplied baseline | Initial Newton incumbent | Fixed quadratic |
|---|---:|---:|---:|
| Qubits | 260 | 256 | 220 |
| Raw Toffolis | 7,735 | 7,327 | 4,449 |
| Fine rotations | 37 | 36 | 43 |
| Toffoli-equivalent cost C | 7,987.593207 | 7,583.127187 | 4,741.910048 |
| Actual objective Q*C | 2,076,774.233715 | 1,941,280.559862 | **1,043,220.210631** |
| Verified shots | 200,540,928 | 200,540,928 | 200,540,928 |

This is a **1.990734-fold improvement over the supplied baseline** (49.7673%
cost reduction), and a **1.860854-fold improvement over the initial Newton
incumbent** (46.2612% reduction). Rotation-synthesis costs use the official
verifier's measured remaining phase-error margin. The actual objective equals
the earlier feasibility estimate exactly.

Worst measured error used 61.2532481% of the allowed tolerance, at x=1,052,655.
The full verifier's worst ratio and input agree with the exhaustive integer
model. The completed run took 205.63 seconds on the exercised local environment;
this is an observed run duration, not a general performance guarantee.

## Implementation and coefficient integrity

`quadratic.py` is a standalone generator using kickmix and the checked-in
`feasibility/B-quadratic-192.json` certificate. Building imports neither the
competition's generator/verifier nor the feasibility resource-probe module.
There are no coefficient-fitting or tuning options. The fixed design uses
192 intervals, degree two, 36-bit coefficients and all 20 fractional input bits.
Both the certificate document and the integer coefficient table are hash-pinned.

After normalization, M=x*2^(2s)/2^20 lies in [1,4). The circuit evaluates the
local quadratic with row-truncated Horner arithmetic. Its 36-bit result, extended
by ten clean bits and shifted by s, supplies a 46-bit phase word. Applying each
bit's phase gives approximately 2^(s-10)/sqrt(M)=1/sqrt(x) turns. Arithmetic,
coefficient lookup and normalization are then uncomputed.

The standalone circuit is byte-for-byte identical to the reviewed feasibility
circuit. No certified coefficient, arithmetic width, lookup strategy or phase
angle changed. The canonical binary's SHA-256 is

```text
76312120708d23ca141a9acc55037830363fcacebbd2492bd2145358f31f0bdd
```

The original certificate document's SHA-256 remains

```text
2b768e13c7054d0e58bc4edf7afc0ff64559685aa309ad250b275fd9fd6a4c8d
```

Seven generator tests passed: frozen coefficient integrity with mutation
rejection, byte identity with the reviewed compiler, resource/register counts,
normalization/interval-boundary model agreement, a 65,536-shot small-domain
official verifier run, binary round-trip integrity, and building in an isolated
directory containing only the generator and its certificate JSON.

## Reproduce and inspect

With the kickmix environment activated, run from this project directory:

```bash
python3 quadratic.py
QF_INVSQRT_DATA=/path/to/extracted/competition python3 -m unittest -v test_quadratic.py
python3 verify_quadratic.py --data /path/to/extracted/competition --workers 4
```

The runner hash-checks the external official verifier before use, preserves its
full problem specification, verifies the serialized circuit, and checks that the
coefficients, verifier, binary and executable sources did not change during the
run. The build manifest records the generator, coefficient and circuit hashes.
No competition source or downloaded data is copied into this project.

- [Build manifest](results/quadratic.build.json)
- [Actual full-verifier result and comparisons](results/quadratic-verified.json)
- Binary artifact: `results/quadratic.kmb` (generated locally and ignored by Git).
- [Prior symbolic certificates and independent reviews](feasibility/COMPARISON.md)

This is a local full-verifier PASS using the declared salt. It is not a Kaggle
submission or a proof over every possible measurement outcome. With the other
two competition tracks at baseline, the measured inverse-square-root improvement
would give an aggregate factor of 1.257972. The historical feasibility decision
and old verifier reports remain unchanged; the JSON linked here records the new
completed run.
