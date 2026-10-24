# Mathematical contract and documentation

The physical object is the shifted two-electron coefficient tensor on 76 real
spatial orbitals. This is an audit of its factorization and block encoding
resources, rather than simulation of the 113-electron many-body state.

For real symmetric matrices use the orthonormal coordinate map svec: diagonal
entries have weight 1, off-diagonal entries sqrt(2). Then
`svec(X).svec(Y) = Tr[X.Y]` and tensor Frobenius norms are preserved. The
pair dimension is 76*77/2 = 2926. The official shifted target is
`tei + (S tensor I + I tensor S)/2`; S and beta are fixed by the verifier.

Write its pair matrix as H. Each DFTHC factor gives a vector x_rc=svec(X_rc),
so its reconstruction is `Sum[x_rc x_rc^T]`: positive semidefinite (PSD).
Within rho, all x_rc lie in the span of B rank-one projectors u_b u_b^T.
Consequently its rank is at most `Min[2926,R Min[B,C]]`.

The best unrestricted PSD rank-k approximation keeps the k largest positive
eigenvalues of H, sets the others to zero, and drops every negative eigenvalue.
To see the minimum, diagonalize H orthogonally. Von Neumann's trace inequality
bounds Tr[H.Y] by the similarly ordered eigenvalue inner product. Thus
`||H-Y||_F^2` is bounded below by the scalar squared eigenvalue differences.
Minimizing over at most k nonzero, nonnegative eigenvalues keeps those largest
positive eigenvalues. The residual squared is the sum of discarded positive
eigenvalue squares plus *all* negative eigenvalue squares. This minimum is
achieved by that spectral truncation. It is a necessary bound for DFTHC;
the shared-projector structure can make the actual best fit much worse.

Machine eigenvalues are numerically crosschecked in a separate SciPy driver,
and WL checks eigenpair residuals and orthogonality. The 1e-8 numerical guard
is deliberately explicit, and is not a rigorous interval rounding certificate.
A shape is only numerically rejected if its bound minus the guard exceeds
0.020 Ha. A result within the guard returns REFINE, never feasibility success.

The cost sweep holds lambda_eff at the incumbent value, 48.68043460482008,
solely for comparison. Lower rank and lower precision do not establish this
normalization or either residual. Ordinary unary loading is a reference
construction: data-dependent pruning and other loaders change its cost.
There is no claim that fixing lambda gives a global lower bound on future cost.

The released four-address design uses a Boolean identity independently checked
in WL for all four low-address cases. Row signs leave u u^T unchanged. Signed
10-bit words for [-pi/2,pi/2) have bit9=bit8, checked on all 512 words. These
facts alone do not verify an emitted circuit or its measurement phase repair.
The public module checks its table and row-sign identities; full official
circuit verification is still required before using it as a submission.

## Authoritative APIs read

Installed Wolfram 15.0 documentation notebooks were converted with nb-reader
to temporary Markdown and read, including Details, Options and relevant examples.

| API | Documented behavior used |
|---|---|
| [Eigenvalues](https://reference.wolfram.com/language/ref/Eigenvalues.html) | Numerical eigenvalues sort by absolute magnitude; the tail explicitly sorts positive values algebraically. |
| [Eigensystem](https://reference.wolfram.com/language/ref/Eigensystem.html) | Eigenvectors are rows; both the eigen-equation and orthogonality are checked. |
| [Norm](https://reference.wolfram.com/language/ref/Norm.html) | Matrix default is operator norm; all matrix residuals explicitly use `"Frobenius"`. |
| [BinaryReadList](https://reference.wolfram.com/language/ref/BinaryReadList.html) | `BinaryReadList[file,"Real64",ByteOrdering->-1]` reads the little-endian array and closes a file opened by the call. |
| [PolynomialMod](https://reference.wolfram.com/language/ref/PolynomialMod.html) | Reduces the four-address Boolean interpolation coefficients modulo 2. |

Core WL symbols only; no paclet or QuantumFramework API is used in this
coefficient-tensor audit. The official Python cost implementation is read and
independently reproduced in WL, including its special power-of-two B handling
and the distinct inner scratch exponent used for qubit accounting.

## Required regimes

| Regime | Evidence |
|---|---|
| General symbolic | Pair inner product, tensor norm, projector Gram, row-sign and four-address identities; spectral argument above for arbitrary dimension. |
| Exact special | Diagonal spectrum {a,b,-c}, a>=b>=0, c>=0 gives rank-one error sqrt(b^2+c^2); mixed spectrum {5,3,-2}. |
| Limiting | A discarded positive mode tends to zero, leaving the negative-spectrum floor; large allowed rank still retains this floor. |
| Numerical reference | Official full N^4 baseline residual versus packed residual; full 2926 eigensystem versus independent SciPy; complete official cost sweep crosscheck. |
| Failure/edge | Rank zero, empty/degenerate spectra, negative eigenvalues, invalid spectrum/rank/shape, power-of-two B, and explicit Q128 acceptance/Q133 rejection. Every retained cost row also asserts `within_cap`. |
