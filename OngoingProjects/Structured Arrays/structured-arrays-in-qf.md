---
Template: Default
---

# Structured arrays in QuantumFramework: where they pay

QuantumFramework (QF) stores every operator as a `SparseArray`, even when the operator was built as a diagonal, a permutation, a Fourier matrix or a block-diagonal matrix. Wolfram Language structured arrays (`DiagonalMatrix`, `PermutationMatrix`, `FourierMatrix`, `BlockDiagonalMatrix` with `TargetStructure -> "Structured"`) keep that structure, and WL's linear algebra uses it. Each section below takes one operation QF performs, writes it as a formula, runs QF's current code, runs the same operation on a structured array, and times both.

The sizes are chosen where QF's current code takes at least about a second. Below that, no difference is worth a change to QF. QF's own property cache stays on, and each QF timing is the first call on a freshly built operator. Each structured timing is the fastest of three runs, each after clearing the system cache.

## Setup

Load QF from the working tree:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
{$Version, $ProcessorCount}
```

`timeQF` times one call; `timeNew` takes the fastest of three calls; `compare` sets the two side by side, with their ratio:

```wl
SetAttributes[{timeQF, timeNew}, HoldFirst];
timeQF[e_] := (ClearSystemCache[]; First[AbsoluteTiming[e;]]);
timeNew[e_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[e;]], 3]];
compare[q_, s_] := Grid[{{"QF", "structured", "QF / structured"}, {Quantity[q, "Seconds"], Quantity[s, "Seconds"], NumberForm[q/s, 2]}}, Frame -> All];
```

Several sections use the Ising energy $h_x = \sum_{i<j} J_{ij} z_i z_j$ of every computational basis state $x$, with $z_i = \pm 1$ the spins of $x$. This is the diagonal of $H = \sum_{i<j} J_{ij} Z_i Z_j$:

```wl
zz[n_] := With[{J = RandomReal[{-1, 1}, {n, n}], z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
SeedRandom[7];
Length[zz[3]]
```

## 1. Evolution under a diagonal Hamiltonian

The time-evolution operator of an Ising Hamiltonian, the cost unitary of QAOA:

$$ U(\theta) = e^{-i \theta H}, \qquad H = \sum_{i<j} J_{ij} Z_i Z_j $$

$H$ is diagonal, so $U$ is diagonal with entries $e^{-i\theta h_x}$. QF computes `Exp[-I θ qo]` as `MatrixExp[Log[E] op["Matrix"]]` on the stored `SparseArray` (`QuantumOperator.m:577`). The structured array exponentiates each diagonal entry.

```wl
n = 14; θ = 0.3; h = zz[n];
qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]];
```

```wl
compare[
    timeQF[Exp[-I θ qo]],
    timeNew[MatrixExp[DiagonalMatrix[-I θ h, TargetStructure -> "Structured"]]]
]
```

Both give the same operator; the difference of the two results applied to a random state is:

```wl
Norm[Exp[-I θ qo]["Matrix"] . v - MatrixExp[DiagonalMatrix[-I θ h, TargetStructure -> "Structured"]] . v]
```

The work QF does grows with the matrix, while the structured array does one exponential per basis state. Each qubit added doubles the structured cost and multiplies QF's by far more, so the gap widens with every qubit.

## 2. A function of a diagonal operator

Any function of a diagonal operator acts entry by entry:

$$ \cos H = \operatorname{diag}\big(\cos h_x\big) $$

QF sends `Cos[qo]` through `matrixFunction` (`Utilities.m:313-317`). When the stored matrix passes `DiagonalMatrixQ` (`Utilities.m:330`), QF applies the function to each diagonal entry, which is the same computation a structured `DiagonalMatrix` performs, so here QF already uses the structure without storing it:

```wl
n = 12; h = zz[n];
qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
compare[
    timeQF[Cos[qo]],
    timeNew[DiagonalMatrix[Cos[h], TargetStructure -> "Structured"]]
]
```

For a diagonal operator, $\operatorname{diag}(\cos h_x)$ is the exact answer, so it is also the check; the largest entry-wise error of QF's result:

```wl
Max[Abs[Diagonal[Cos[qo]["Matrix"]] - Cos[h]]]
```

`DiagonalMatrixQ` treats an off-diagonal entry at the level of rounding error, relative to the size of the matrix, as zero. Dropping such an entry changes $f(M)$ by about as much as floating-point evaluation of $f(M)$ already does, while any entry large enough to be a real coupling keeps the matrix off this route. The test reads only the stored nonzero entries of a `SparseArray`, so it costs a fraction of a millisecond even for twelve qubits:

```wl
timeNew[DiagonalMatrixQ[qo["Matrix"]]]
```

`MatrixFunction` applied to a structured `DiagonalMatrix` computes the right values but returns an ordinary dense matrix, so the structure is lost; writing the diagonal of $f(H)$ directly keeps it:

```wl
Head /@ {MatrixFunction[Cos, DiagonalMatrix[h, TargetStructure -> "Structured"]], DiagonalMatrix[Cos[h], TargetStructure -> "Structured"]}
```

## 3. Spectrum of a permutation

The cyclic shift $X_d \lvert k\rangle = \lvert k+1 \bmod d\rangle$ is a permutation with a single cycle of length $d$, so its eigenvalues are the $d$-th roots of unity:

$$ \operatorname{spec}(X_d) = \{ e^{2\pi i k/d} \}_{k=0}^{d-1} $$

QF's `"Eigenvalues"` (`QuantumOperator/Properties.m:506`) computes the full `Eigensystem` of the exact matrix and simplifies it (`Utilities.m:162`). `PermutationMatrix` reads the spectrum from the cycle structure. QF's `"X"[d]` is `PermutationMatrix[RotateRight[Range[d]]]`:

```wl
d = 64;
qo = QuantumOperator["X"[d]];
P = PermutationMatrix[RotateRight[Range[d]], TargetStructure -> "Structured"];
Normal[qo["Matrix"]] === Normal[P]
```

```wl
compare[timeQF[qo["Eigenvalues"]], timeNew[Eigenvalues[P]]]
```

The two spectra are the same set:

```wl
Sort[N[qo["Eigenvalues"]]] == Sort[N[Eigenvalues[P]]]
```

## 4. Spectrum of the quantum Fourier transform

The QFT on $n$ qubits is the $2^n$-point Fourier matrix $F_{jk} = \omega^{jk}/\sqrt{2^n}$, $\omega = e^{2\pi i/2^n}$. Since $F^4 = I$, its eigenvalues lie in $\{1, -1, i, -i\}$, with known multiplicities. QF builds the operator of its `"Fourier"` circuit and takes the same `"Eigenvalues"` route as above; `FourierMatrix` knows the answer:

```wl
n = 3;
qo = QuantumCircuitOperator["Fourier"[n]]["QuantumOperator"];
F = FourierMatrix[2^n, TargetStructure -> "Structured"];
Simplify[Normal[qo["Matrix"]] - Normal[F]] === ConstantArray[0, {2^n, 2^n}]
```

```wl
compare[timeQF[qo["Eigenvalues"]], timeNew[Eigenvalues[F]]]
```

```wl
Sort[N[qo["Eigenvalues"]]] == Sort[N[Eigenvalues[F]]]
```

## 5. Spectrum of a diagonal operator

The eigenvalues of a diagonal operator are its diagonal entries, $\operatorname{spec}(H) = \{h_x\}$. QF runs the full eigensolver on the stored matrix, converting it to dense first:

```wl
n = 10; h = zz[n];
qo = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
compare[
    timeQF[qo["Eigenvalues"]],
    timeNew[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]]
]
```

```wl
Max[Abs[Sort[qo["Eigenvalues"]] - Sort[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]]]]
```

A second `"Eigenvalues"` call on the same operator is answered from QF's property cache, so this gain is on the first call.

## 6. Change of basis into the Fourier basis

Rewriting a state in a new basis with matrix $B$ solves

$$ \lvert\psi'\rangle = B^{-1} \lvert\psi\rangle $$

For the Fourier basis, $B = F_d$ and $F_d^{-1} = F_d^\dagger$. QF inverts the exact basis matrix with `MatrixInverse` (`QuantumState.m:187-188`, `Utilities.m:295-299`), an exact inverse of a matrix of roots of unity. QF's Fourier basis matrix is `FourierMatrix[d]`:

```wl
d = 8;
qb = QuantumBasis["Fourier"[d]];
Simplify[Normal[qb["Output"]["ReducedMatrix"]] - Normal[FourierMatrix[d]]] === ConstantArray[0, {d, d}]
```

```wl
v = Normalize[RandomComplex[{-1 - I, 1 + I}, d]];
qs = QuantumState[v];
F = FourierMatrix[d, TargetStructure -> "Structured", WorkingPrecision -> MachinePrecision];
compare[timeQF[QuantumState[qs, qb]], timeNew[Inverse[F] . v]]
```

```wl
Norm[Normal[QuantumState[qs, qb]["StateVector"]] - Inverse[F] . v]
```

The exact inverse grows so fast with $d$ that QF's change of basis into a 16-dimensional Fourier basis does not finish in minutes, while the structured array still answers in milliseconds.

## 7. Quantum Fourier transform of a state

Applying the QFT to a state vector is a discrete Fourier transform of its amplitudes, $\lvert\psi\rangle \mapsto F \lvert\psi\rangle$. QF contracts the gates of its `"Fourier"` circuit as a tensor network. `FourierMatrix` applied to a vector runs a fast Fourier transform:

```wl
n = 16;
qc = QuantumCircuitOperator["Fourier"[n]];
v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]];
qs = QuantumState[v];
F = FourierMatrix[2^n, TargetStructure -> "Structured", WorkingPrecision -> MachinePrecision];
compare[timeQF[qc[qs]], timeNew[F . v]]
```

```wl
Norm[Normal[qc[qs]["StateVector"]] - F . v]
```

Both costs roughly double with each added qubit, so the gap stays about a constant factor as $n$ grows. The factor itself is the difference between contracting the circuit gate by gate and doing the whole transform in one fast Fourier pass over the amplitudes.

## 8. Controlled evolution

A Hamiltonian that acts only when a control qubit is in $\lvert 1\rangle$ generates a controlled evolution:

$$ e^{-i t\, \lvert 1\rangle\langle 1\rvert \otimes H} = \lvert 0\rangle\langle 0\rvert \otimes I + \lvert 1\rangle\langle 1\rvert \otimes e^{-i t H} $$

The generator is block diagonal, with a zero block and $H$. QF builds it with `QuantumTensorProduct` and exponentiates the stored `SparseArray`. A structured `BlockDiagonalMatrix` exponentiates only the $H$ block:

```wl
SeedRandom[5]; n = 7; t = 0.7; k = 2^(n - 1);
Hm = (# + ConjugateTranspose[#])/2 &[RandomComplex[{-1 - I, 1 + I}, {k, k}]];
qo = QuantumTensorProduct[QuantumOperator[{{0, 0}, {0, 1}}, {1}], QuantumOperator[Hm, Range[2, n]]];
B = BlockDiagonalMatrix[{ConstantArray[0., {k, k}], -I t Hm}, TargetStructure -> "Structured"];
v = Normalize[RandomComplex[{-1 - I, 1 + I}, 2 k]];
compare[timeQF[Exp[-I t qo]], timeNew[MatrixExp[B]]]
```

```wl
Norm[Exp[-I t qo]["Matrix"] . v - MatrixExp[B] . v]
```

The fair baseline here is the dense matrix, which WL exponentiates quickly. The structured array beats it only once $H$ is large, because it skips the zero block:

```wl
Dn = Developer`ToPackedArray[Normal[-I t qo["Matrix"]]];
compare[timeNew[MatrixExp[Dn]], timeNew[MatrixExp[B]]]
```

At seven qubits QF, the dense matrix and the structured array all finish in a fraction of a second, so nothing here is worth changing. The structured array pulls ahead of the dense exponential only from about eleven qubits, where skipping the zero block saves real time; this is the one row where the gain is a constant factor rather than a change in how the cost grows.

## 9. Dephasing of a qubit register

Each qubit of a register loses phase coherence at rate $\gamma$ while evolving under a diagonal $H$:

$$ \dot\rho = \mathcal L \rho = -i[H, \rho] + \gamma \sum_k \left( Z_k \rho Z_k - \rho \right) $$

For diagonal $H$, $\mathcal L$ is diagonal on the $4^n$ entries of $\rho$: each coherence $\rho_{xy}$ picks up a phase $e^{-i(h_x - h_y)t}$ and decays with the number of qubits where $x$ and $y$ differ. QF builds $\mathcal L$ with its `"Liouvillian"` named operator (`NamedOperators.m:904-925`) and exponentiates it as a dense matrix:

```wl
n = 5; t = 0.7; h = zz[n];
H = QuantumOperator[SparseArray[Band[{1, 1}] -> h], Range[n]];
zs = Table[QuantumOperator[SparseArray[KroneckerProduct @@ ReplacePart[ConstantArray[IdentityMatrix[2], n], k -> {{1, 0}, {0, -1}}]], Range[n]], {k, n}];
L = QuantumOperator["Liouvillian"[H, zs, ConstantArray[0.2, n]]];
m = Normal[L["Matrix"]];
Max[Abs[m - DiagonalMatrix[Diagonal[m]]]]
```

The largest off-diagonal entry is zero, so $\mathcal L$ is its diagonal:

```wl
l = Diagonal[m];
v = Normalize[RandomComplex[{-1 - I, 1 + I}, Length[l]]];
compare[timeQF[Exp[t L]], timeNew[MatrixExp[t DiagonalMatrix[l, TargetStructure -> "Structured"]]]]
```

```wl
Norm[Exp[t L]["Matrix"] . v - MatrixExp[t DiagonalMatrix[l, TargetStructure -> "Structured"]] . v]
```

Open-system evolution is where the gap opens earliest: the superoperator lives on $4^n$ entries, so QF's dense exponential becomes slow already at a handful of qubits.

## 10. Whole jobs: phase estimation and QAOA

A single operation can be fast or slow in isolation; what a user waits for is a whole job. Two common jobs are timed here, each against every method QF offers for applying a circuit: its default tensor-network contraction, the same contraction along a greedy contraction path, and the gate-by-gate `"Schrodinger"` method. `table` lists the times:

```wl
table[rows_] := Grid[Prepend[{#1, Quantity[#2, "Seconds"]} & @@@ rows, {"method", "time"}], Frame -> All, Alignment -> Left];
greedy = {"TensorNetwork", "Path" -> "Greedy"};
```

### Phase estimation

QF's `"PhaseEstimation"` circuit prepares the counting register, applies the controlled powers of $U$, and ends with the inverse QFT and a measurement of each counting qubit. Dropping the measurements leaves the state to compare; dropping the inverse QFT too leaves the part that stays with QF, after which the inverse QFT is applied as an inverse fast Fourier transform along the counting register:

```wl
invQFT[v_, n_] := With[{m = ArrayReshape[v, {2^n, Length[v]/2^n}]}, Flatten[Transpose[InverseFourier /@ Transpose[m]]]];
SeedRandom[3]; u = QuantumOperator["RandomUnitary"[{2, 2}]];
n = 12;
els = QuantumCircuitOperator["PhaseEstimation"[u, n]]["Elements"];
body = QuantumCircuitOperator[Drop[els, -n]];
front = QuantumCircuitOperator[Drop[els, -(n + 1)]];
Last[body["Elements"]]["Label"]
```

The block dropped from `front` is the inverse QFT. The whole job, each way:

```wl
table[{
    {"QF, default contraction", timeQF[body[]]},
    {"QF, greedy path", timeQF[body[Method -> greedy]]},
    {"QF, Schrodinger", timeQF[body[Method -> "Schrodinger"]]},
    {"QF without inverse QFT, then FFT", timeQF[invQFT[Normal[front[]["StateVector"]], n]]},
    {"QF greedy without inverse QFT, then FFT", timeQF[invQFT[Normal[front[Method -> greedy]["StateVector"]], n]]}}]
```

```wl
Norm[invQFT[Normal[front[Method -> greedy]["StateVector"]], n] - Normal[body[]["StateVector"]]]
```

The inverse QFT is a large share of the default contraction, and there are two ways to remove that cost: the fast Fourier transform, or a greedy contraction path, which orders the same gates so that the intermediate tensors stay small. Each recovers a similar factor on its own, and together they give the fastest run. With a two-qubit $U$ the whole job stays short at every size tried, so for phase estimation the saving is real but small in absolute terms; it matters when the job is repeated many times. A larger $U$ makes its controlled powers more expensive, which should shrink the QFT's share further; that case is not timed here.

### QAOA

A depth-one QAOA scan evaluates the energy $\langle C\rangle$ of the Ising cost function on a grid of angles $(\gamma, \beta)$. The circuit is a layer of Hadamards, one $e^{-i\gamma J_{ij} Z_i Z_j}$ gate for each of the $n(n-1)/2$ pairs, and one $R_X(2\beta)$ on each qubit. The whole cost layer is the diagonal $e^{-i\gamma h}$, so the structured route replaces those $n(n-1)/2$ gates by one diagonal multiply and hands only the mixer to QF:

```wl
zzJ[J_, n_] := With[{z = 1 - 2 Tuples[{0, 1}, n]},
    Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
n = 12; SeedRandom[5]; J = RandomReal[{-1, 1}, {n, n}]; h = zzJ[J, n];
qc = QuantumCircuitOperator[Join[
        Table[QuantumOperator["H", {i}], {i, n}],
        Flatten @ Table[QuantumOperator["R"[2 γ J[[i, j]], "ZZ"], {i, j}], {i, n}, {j, i + 1, n}],
        Table[QuantumOperator["RX"[2 β], {i}], {i, n}]],
    "Parameters" -> {γ, β}];
mix = QuantumCircuitOperator[Table[QuantumOperator["RX"[2 β], {i}], {i, n}], "Parameters" -> {β}];
grid = Tuples[{Subdivide[0.1, 1.0, 4], Subdivide[0.1, 1.0, 4]}];
v0 = ConstantArray[1/Sqrt[2.^n], 2^n];
energy[v_] := Re[Conjugate[v] . (h v)];
qfScan[m_] := energy[Normal[qc[#[[1]], #[[2]]][Method -> m]["StateVector"]]] & /@ grid;
structuredScan[m_] := energy[Normal[mix[#[[2]]][QuantumState[DiagonalMatrix[Exp[-I #[[1]] h], TargetStructure -> "Structured"] . v0], Method -> m]["StateVector"]]] & /@ grid;
Length[grid]
```

The QF `"R"[θ, "ZZ"]` gate is $e^{-i\theta Z Z/2}$, so each gate above is $e^{-i\gamma J_{ij} Z_i Z_j}$ and their product is exactly $e^{-i\gamma h}$. The scan, each way (the `"Schrodinger"` method is far slower than both contractions for this circuit and is left out to keep the page quick to run):

```wl
table[{
    {"QF, default contraction", timeQF[qfScan[Automatic]]},
    {"QF, greedy path", timeQF[qfScan[greedy]]},
    {"structured cost layer, QF default mixer", timeQF[structuredScan[Automatic]]},
    {"structured cost layer, QF greedy mixer", timeQF[structuredScan[greedy]]}}]
```

```wl
Max[Abs[qfScan[greedy] - structuredScan[greedy]]]
```

Here a contraction path helps much less than in phase estimation: the cost of this circuit is the number of two-qubit gates, and reordering them does not remove any. The structured cost layer removes them all, so it wins whichever method applies the mixer, and the time that remains is the mixer's single-qubit layer run through QF. Over the sizes tried the factor between QF's default and the structured route stays roughly constant, so the time saved grows in step with the job itself.

## Where structured arrays do not help

`HadamardMatrix` has no structured form that applies $H^{\otimes n}$ quickly: it accepts only `"Dense"`, `"Orthogonal"` and `"Unitary"` as target structures, and the last two wrap the dense matrix.

```wl
Head /@ {HadamardMatrix[8, Method -> "BitComplement", TargetStructure -> "Orthogonal"], HadamardMatrix[8, Method -> "BitComplement", TargetStructure -> "Unitary"]}
```

`Eigensystem` of a `DiagonalMatrix` head returns malformed eigenvectors in WL 15.0, while `Eigenvalues` is correct. A bare `DiagonalMatrix[list]` with more than 1000 entries is structured by default, so this also happens with no option given:

```wl
Eigensystem[DiagonalMatrix[{1., 3., 4.}, TargetStructure -> "Structured"]]
```

The structured array pays whenever the operator's structure already fixes the answer and QF recomputes it from scratch: the exponential and spectrum of a diagonal are its entries, a permutation's spectrum comes from its cycles, $F^{-1} = F^\dagger$, and $F$ on a vector is a fast Fourier transform. It stops paying where the work is genuinely dense, as for the $H$ block of a controlled evolution, and wherever WL has no structured form, as for $H^{\otimes n}$. In whole jobs the gain follows the share of the work the structure covers: large for QAOA, whose cost layer is entirely diagonal, and modest for phase estimation, where the QFT is one block among the controlled powers and a greedy contraction path already recovers much of the same time.
