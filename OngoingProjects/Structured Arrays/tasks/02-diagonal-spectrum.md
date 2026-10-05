# Brief 2: the spectrum of a diagonal operator

Read `README.md` in this folder first: it has the workflow, the house rules and the traps every brief shares.

## Goal

The eigenvalues of a diagonal matrix are its diagonal entries, and its eigenvectors are the unit vectors. QF runs a full eigensolver on it instead: at 12 qubits `qo["Eigenvalues"]` of a diagonal Ising operator takes 46 s, where reading the diagonal takes milliseconds. The task: give a diagonal operator its spectrum directly, in an order decided first (decision 1).

A smaller second part: the exact non-integer power of an exact diagonal operator.

## What QF does today

Measured on main at `f6299b43` by `baseline.wls spectrum`, `baseline.wls power` and `supplement2.wls` (raw output in `baseline-out/`).

`"Eigenvalues"`, `"Eigenvectors"` and `"Eigensystem"` (`QuantumFramework/Kernel/QuantumOperator/Properties.m:520-524`) call `eigenvalues`, `eigenvectors` and `eigensystem` (`Utilities.m:210-255` and `:329-335`) on `qo["MatrixRepresentation"]`, with `"Sort" -> False` and `"Normalize" -> True`. With the default `Chop -> False`, `eigensystem` calls `Eigensystem` on the stored `SparseArray`. WL converts it to a dense matrix, printing `Eigensystem::arh`, and solves the full dense problem; `Simplify` and the roundoff rules (`chopEigensystem`, `:282-286`) follow. The same `eigensystem` also serves `"Projectors"`, `"Diagonalize"` and the eigenbasis of a measurement operator (`orthonormalEigensystem`, `Properties.m:635`).

| Diagonal Ising operator | `"Eigenvalues"` | `"Eigensystem"` | bare `Eigensystem` of the sparse matrix | `Simplify` of its result |
|---|---|---|---|---|
| 8 qubits | 0.24 s | 0.16 s | 0.23 s | 0.004 s |
| 10 qubits | 3.1 s | 2.5 s | 3.4 s | 0.08 s |
| 12 qubits | 46 s | 39 s | 37 s | 0.9 s |

The eigensolver is the whole cost. The survey measured `Eigenvalues` of a structured `DiagonalMatrix` at 12 qubits in 0.002 s.

The order QF returns today follows WL's eigensolver, not the diagonal, and differs between exact and machine matrices:

| Matrix | Eigenvalues | Eigenvectors, as unit vectors |
|---|---|---|
| exact diag(2, 1, 2) | 2, 2, 1 | e₃, e₁, e₂ |
| machine diag(2., 1., 2.) | 2., 2., 1. | e₁, e₃, e₂ |
| exact diag(1, -1, i, -i, 1) | -1, i, -i, 1, 1 | e₂, e₃, e₄, e₅, e₁ |
| `"Z"[3]`, diagonal (1, e^(-2πi/3), e^(2πi/3)) | -(-1)^(1/3), (-1)^(2/3), 1 | |

The power: a machine diagonal raised to a non-integer power is already computed entry by entry (`qo^0.5` of a 12-qubit diagonal in 0.03 s, exact to the last bit), because `inexactMatrixPower` (`Utilities.m:449-490`) sends a diagonal to the diagonal branch of `scalarMatrixFunction` (`:513-524`). An exact diagonal takes `MatrixPower` (`:452`): `qo^(1/2)` of `DiagonalMatrix[Range[256]]` takes 2.2 s, while `Sqrt[qo]` of the same operator takes 0.004 s.

## Proposed change

1. In `eigensystem` (`Utilities.m:212`), a matrix that passes `diagonalMatrixQ` (strict, no tolerance; defined in `Utilities.m` by commit `191f1ac3`) skips `Eigensystem`: its values are its diagonal entries and its vectors are the unit vectors, as a `SparseArray`, in the order decision 1 fixes. Everything after the solve stays as it is: `chopEigensystem`, the `"Sort"`, `"Normalize"` and `"Orthogonalize"` options. A pattern-matched helper, not another branch in the `Which`.
2. In `matrixFunction[Power, ...]` (`Utilities.m:449-452`), a diagonal matrix with a non-integer exponent goes to `scalarMatrixFunction`'s diagonal branch whatever its precision, as `Sqrt` already does.

## What must not change

- The roundoff rules for inexact eigenvalues. `chopEigensystem` writes 0 for an eigenvalue below the roundoff of the largest (`diag(10^10, 10^-5)` at machine precision has the eigenvalues 10¹⁰ and 0, by design). The tests `InexactEigensystem-*` in `Tests/QuantumOperator.wlt` (around lines 185-260) use diagonal matrices with tiny entries and must keep passing unchanged.
- Spectra of non-diagonal operators: identical to the last bit.
- Symbolic diagonals keep exact closed forms.
- `"MatrixRepresentation"` is the matrix in the computational basis. An operator stored diagonal in its own basis, such as `"JX"`, is not diagonal there and keeps the eigensolver.

## Decisions to make first

1. **The order of eigenvalues and eigenvectors.** WL's order is not the diagonal's, and it differs between exact and machine matrices (table above), so reproducing it exactly would mean reimplementing two solvers' tie rules.
   - Option (a): the diagonal's own order, e₁, e₂, e₃, ..., for exact and machine matrices alike, documented on the reference page. Simple and deterministic, but it changes the order wherever WL's differed: ties, and phases of equal modulus such as the clock matrix `"Z"[3]`.
   - Option (b): sort the diagonal by a stated rule that matches WL on the cases that matter downstream.

   Recommendation: (a), but only after the audit in step 2 shows which callers depend on the order. Mads decides.
2. **Audit the consumers before changing the order.** List every caller of `eigensystem` and of the eigen properties: `"Projectors"`, `"Diagonalize"`, `orthonormalEigensystem`, `QuantumMeasurementOperator`'s outcomes (`QuantumMeasurementOperator[QuantumOperator["Z"]]["Eigenvalues"]` gives the outcomes 1 then -1 today), `QuantumEvolve`, and any doc page or test that prints a spectrum. For each, state whether a new order changes what a user sees. That list is the input to decision 1.

## Tests to write

A new file, `Tests/QuantumOperatorDiagonalSpectrum.wlt`. Each must fail on today's main and pass after the change, except the regression group:

- `"Eigenvalues"` of a 12-qubit diagonal Ising operator finishes inside `TimeConstraint -> 10` with no message, and its eigenvalues are the diagonal entries as a set (compare sorted lists against a bound, not `===`).
- `"Eigensystem"` of the same: every eigenvector v with eigenvalue λ satisfies M v = λ v, and the eigenvectors are orthonormal.
- The order fixed by decision 1, on exact, machine and complex diagonals with ties (the cases in the table above).
- `qo^(1/2)` of an exact 256-entry diagonal finishes inside a time limit and equals `Sqrt[qo]`.
- Regression, must pass before and after: the `InexactEigensystem-*` tests; a symbolic diag(a, b, c); a non-diagonal Hermitian matrix whose spectrum must be identical to the last bit on both copies.

## How to verify and land

Follow the README. Two commits, so each can be reverted alone: the spectrum, and the power. Report the times of the table above before and after.

## Traps

- Do not use a structured `DiagonalMatrix` head: in WL 15.0, `Eigensystem` of `DiagonalMatrix[list, TargetStructure -> "Structured"]` returns malformed eigenvectors. A bare `DiagonalMatrix[list]` with more than 1000 entries is structured by default.
- `Eigensystem::arh` is printed from inside QF's eigensolver call, and it was not seen in `$MessageList` around the property call in `baseline.wls`. Check how a test can assert "no message" before relying on it.
- Exact eigenvalues such as `(-1)^(2/3)` and `E^(2 I Pi/3)` are the same number written differently; compare values, not forms.
