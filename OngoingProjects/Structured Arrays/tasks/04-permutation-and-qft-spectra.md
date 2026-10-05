# Brief 4: spectra of permutation operators and of the QFT (investigate first)

Read `README.md` in this folder first, and brief 2: this brief depends on the eigenvalue order decided there.

## Goal

A permutation operator P|x⟩ = |π(x)⟩ has as eigenvalues, for each cycle of π of length ℓ, the ℓ-th roots of unity, with eigenvectors that are discrete Fourier sums over the cycle. The quantum Fourier transform F on n qubits satisfies F⁴ = 1, so its eigenvalues lie in {1, −1, i, −i}, with multiplicities known in closed form. QF finds both spectra with a general exact eigensolver. The task: decide whether QF should recognize these operators, and if so, give their spectra from the structure.

## What QF does today

Measured on main at `f6299b43` by `baseline.wls permutation` (raw output in `baseline-out/permutation.txt`). `"Eigenvalues"` goes through `eigensystem` on the exact matrix, as described in brief 2.

| Operator | QF `"Eigenvalues"` | structured head | |
|---|---|---|---|
| `QuantumOperator["X"[16]]` (cyclic shift) | 0.10 s | 0.04 s | `PermutationMatrix` |
| `"X"[64]` | 2.3 s | 0.24 s | |
| `"X"[128]` | 33 s, prints `Eigensystem::arhm` | 1.1 s | |
| `QuantumCircuitOperator["Fourier"[3]]["QuantumOperator"]` | 0.72 s | 0.0004 s | `FourierMatrix` |
| the same for 4 qubits | 9.9 s | 0.0004 s | |

QF's `"X"[d]` matrix is `PermutationMatrix[RotateRight[Range[d]]]` (checked at d = 16, 64, 128), built by `pauliMatrix[1, d]` (`Utilities.m:347`). The survey checked that QF's QFT operator equals `FourierMatrix[2^n]`. Even the structured `PermutationMatrix` takes 1 s at d = 128; reading the cycles of π with `PermutationCycles` and writing down roots of unity is linear in d.

## Investigation (do this first, report before changing code)

1. **Is it worth it?** Find who takes spectra of permutations or of the QFT in QF: doc pages, tests, tutorials, `QuantumMeasurementOperator` on such an operator, `"Diagonalize"`. Reversible circuits (X, CNOT, Toffoli, SWAP networks, `"X"[d]` on qudits) are all permutations. Report the uses found; if there are none, the recommendation can be to stop here.
2. **Recognizing a permutation** must be cheap and exact: a square exact `SparseArray` whose stored values are all 1, with exactly one in each row and each column. Measure the test on a 2¹⁶-dimensional permutation and on a non-permutation of the same size.
3. **Recognizing the QFT** from its matrix needs an exact comparison with `FourierMatrix`, which costs a `Simplify` of d² entries. Find whether anything cheaper identifies it: the label QF gives the operator, or the circuit it came from. Labels can be changed by the user, so a label alone cannot decide the route.
4. **Eigenvectors.** Eigenvalues of both are closed forms. Eigenvectors of a permutation are closed forms too; exact eigenvectors of the Fourier matrix are not simple (the comment at `Utilities.m:166-187` describes the 7 × 7 case). Decide whether the fast route covers `"Eigenvalues"` only, or `"Eigensystem"` as well.

## Likely change, to confirm with the investigation

A pattern-matched case in `eigensystem` (`Utilities.m:212`), next to brief 2's diagonal case: a permutation matrix gets eigenvalues from its cycles and eigenvectors e^(−2πi jk/ℓ) on each cycle, normalized, in the order brief 2 fixes. The QFT, only if step 3 finds a cheap and safe test.

## What must not change

The spectrum of every operator that is not recognized, identical to the last bit, and the roundoff rules of `chopEigensystem`.

## Tests to write

- `"Eigenvalues"` of `"X"[128]` inside a time limit, equal as a set to the 128th roots of unity.
- A permutation with several cycles (a SWAP network, a Toffoli): eigenvalues as a set from its cycle lengths; every eigenvector satisfies P v = λ v; the eigenvectors are orthonormal.
- Regression on both copies: a matrix that is almost a permutation (one entry 2, or one row with two entries) keeps today's result.

## Traps

- Exact roots of unity have many forms (`(-1)^(2/3)`, `E^(2 I Pi/3)`); compare values.
- The order of eigenvalues is brief 2's decision; do not invent a second one here.
