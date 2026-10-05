# Brief 1: change of basis with a unitary basis matrix

Read `README.md` in this folder first: it has the workflow, the house rules and the traps every brief shares.

## Goal

Rewriting a state in a basis with matrix B solves |ψ'⟩ = B⁻¹|ψ⟩. When B is unitary, B⁻¹ = B†, which costs nothing. QF instead computes the exact inverse of B, and for an exact basis matrix of roots of unity that inverse does not finish in two minutes at dimension 16. The task: when a basis matrix is exactly unitary, use its conjugate transpose.

## What QF does today

Measured on main at `f6299b43` by `baseline.wls basis` and `supplement.wls` (raw output in `baseline-out/basis.txt` and `baseline-out/supplement.txt`).

All inverses of basis matrices go through one function, `MatrixInverse` (`QuantumFramework/Kernel/Utilities.m:421-425`):

```wl
MatrixInverse[matrix_] := If[
    SquareMatrixQ[matrix],
    Quiet[Check[Inverse[matrix], PseudoInverse[matrix], Inverse::sing], Inverse::sing],
    PseudoInverse[matrix]
]
```

Its callers:

| Where | What it inverts |
|---|---|
| `QuantumState/QuantumState.m:187-188` | the new and the old basis matrices when a pure state changes basis, `QuantumState[qs, newBasis]` |
| `QuantumState/QuantumState.m:171` (`doubledReducedMatrix`), used at `:201-202` | the same for a mixed state |
| `QuditBasis/Properties.m:225` | the `"Inverse"` property of a `QuditBasis` |
| `QuantumCircuitOperator/TensorNetwork.m:33` and `:103` | the input basis of an operator whose input is not computational, when a circuit is contracted |

QF's Fourier basis matrix is exact, with entries such as e^(iπ/8)/4, and equals `FourierMatrix[d]` (checked for d = 16, 32, 64).

| d | `QuantumState[qs, QuantumBasis["Fourier"[d]]]` | `Inverse[m]` alone | `ConjugateTranspose[m]` | exact `UnitaryMatrixQ[m]` | `Inverse[N[m]]` |
|---|---|---|---|---|---|
| 8 | 0.56 s | 0.55 s | 0.0002 s | 0.003 s | |
| 16 | did not finish in 120 s | did not finish in 120 s | 0.0009 s | 0.015 to 0.027 s | 0.0001 s |
| 32 | | | | 0.11 s | 0.0008 s |
| 64 | | | | 0.97 s | 0.007 s |

`QuditBasis["Fourier"[16]]["Inverse"]` also did not finish in 120 s. At d = 8, QF's result and `ConjugateTranspose[m] . v` agree to 1.1 × 10⁻¹⁶, so the conjugate transpose is the same answer. The whole cost is the exact inverse: the machine inverse of the same matrix is instant.

Which named bases, at their default arguments, have a unitary matrix (numerically, to 10⁻¹⁰):

- Unitary: Computational, PauliX, PauliY, PauliZ, Identity, I, X, Y, Z, JX, JY, JZ, J, JI, Bell, Fourier, Ivanovic.
- Square but not unitary (frames and operator bases): Schwinger, Dirac, Wigner, WignerMIC, Pauli, GellMann, GellMannMIC, Bloch, BlochSphere, GellMannBloch, GellMannBlochMIC, Wootters, Feynman, Tetrahedron, RandomMIC, RandomHaarMIC, RandomBlochMIC, QBismSIC, HesseSIC, HoggarSIC.

A Fourier basis built over another orthonormal basis (`QuantumBasis["Fourier"[QuditBasis["PauliX"[3]]]]`) is unitary too, so testing the matrix covers cases a test of the basis name would miss. Exact `UnitaryMatrixQ` on the 64 × 64 HoggarSIC matrix returns `False` in 0.05 s, and on a symbolic matrix it returns `False` at once.

## Proposed change

In `MatrixInverse`, before `Inverse`: a square matrix whose entries are exact numbers and that `UnitaryMatrixQ` proves unitary returns `ConjugateTranspose[matrix]`. Every caller then gains at once. Use pattern-matched definitions, not another `If`:

```wl
MatrixInverse[matrix_ ? exactUnitaryMatrixQ] := ConjugateTranspose[matrix]
MatrixInverse[matrix_] := (the present definition)
```

with `exactUnitaryMatrixQ` true for a square matrix of exact numeric entries (`MatrixQ[m, NumericQ]` and `Precision[m] === Infinity`) on which `UnitaryMatrixQ` returns `True`. Keep a `SparseArray` sparse (`ConjugateTranspose` does).

## What must not change

- The result for every basis that is not exactly unitary: frames, operator bases, machine and symbolic matrices. They keep `Inverse` and `PseudoInverse` as today.
- The value of every exact result. The conjugate transpose of an exact unitary matrix is its exact inverse, but the two can be written differently (`E^(-I Pi/8)/4` against a form `Inverse` simplified). Tests compare values with `Simplify[a - b] === 0`-style checks, not with `===`.
- The existing `Quiet` in `MatrixInverse`. It breaks the house rule, but changing it is a separate decision; do not widen it.

## Decisions to make first

1. **Exact matrices only, or machine matrices too?** Recommendation: exact only. The machine inverse is already instant (0.007 s at d = 64), and the conjugate transpose of a machine matrix that is unitary only to roundoff differs from its inverse by that roundoff, so the gain is small and the results would move.
2. **Is the test cost acceptable for non-unitary exact bases?** The exact test grows about eightfold per doubling of d (1 s at d = 64). For a frame basis it is spent and then `Inverse` runs anyway. Recommendation: accept it, after measuring it on the largest exact frame basis QF users build, and reporting the number.

## Tests to write

A new file, `Tests/QuantumBasisUnitaryInverse.wlt` (or tests added to `Tests/QuantumState.wlt`). Each must fail on today's main and pass after the change, except the last group, which must pass on both:

- `QuantumState[qs, QuantumBasis["Fourier"[16]]]` for an exact symbolic state finishes inside `TimeConstraint -> 10` and equals `ConjugateTranspose[FourierMatrix[16]] . v` (compare with `Simplify` of the difference).
- The round trip: changing to the 16-dimensional Fourier basis and back to the computational one returns the state.
- A mixed state in the 16-dimensional Fourier basis (the `doubledReducedMatrix` path) equals F† ρ F.
- `QuditBasis["Fourier"[16]]["Inverse"]` finishes inside a time limit and its matrix times the basis matrix is the identity.
- A circuit containing an operator with a Fourier input basis of dimension 16 (the `TensorNetwork.m` path) gives the same state as the explicit matrix product.
- Regression, must pass before and after: a change of basis into a frame basis (QBismSIC or Tetrahedron) and into a machine Fourier basis gives the same result as today, compared to the last bit with `Hash[BinarySerialize[...]]` on both copies.

Then run `Tests/QuantumState.wlt`, `Tests/QuantumBasis.wlt`, `Tests/QuditBasis.wlt` and `Tests/QuantumCircuitOperator.wlt`, which already use Fourier bases, and the full suite.

## How to verify and land

Follow the README: two copies of current main, the new tests on both, the full suite on the changed copy, one commit. Report the before and after times of the table above at d = 8, 16, 32, 64.

## Traps

- `UnitaryMatrixQ` on a symbolic matrix returns `False` without assumptions, even for a matrix that is unitary for every value of its symbols. That is the safe direction (it keeps `Inverse`); do not add assumptions.
- `Conjugate` of a symbolic entry stays `Conjugate[...]`. Only exact numeric matrices take the new route, so this does not arise.
