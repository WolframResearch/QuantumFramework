## Part 4. Composite systems: tensor product and partial trace

Two quantum systems join by the tensor product, and the reverse question, what one subsystem looks like
once the other is ignored, is answered by the partial trace. This Part builds each direction from
primitives and then confirms that the framework object agrees: the tensor product is a Kronecker product,
and the partial trace is a reshape-and-contract.

### 4.1 [BSc] How do I form the tensor product of states and of operators?

Two subsystems combine by the tensor product: states as $|\psi\rangle\otimes|\phi\rangle$, operators as
$A\otimes B$, and the two are compatible, $(A\otimes B)(|\psi\rangle\otimes|\phi\rangle) =
(A|\psi\rangle)\otimes(B|\phi\rangle)$. The joint system carries one amplitude for every pair of subsystem
labels, so the composite lives in the $d_1 d_2$-dimensional product space, and in the computational basis
the tensor product is exactly the Kronecker product.

**WL** : the tensor product of two single-qubit amplitude vectors $\{a,b\}$ and $\{c,d\}$, flattened into
one length-4 vector, one amplitude per computational label $00, 01, 10, 11$.

```wl
Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

Operators combine the same way; here $X\otimes Z$ as the Kronecker product of the two matrices.

```wl
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]
```

The compatibility rule is what makes this a tensor product: acting with $X\otimes Z$ on the product
state is the same as acting with $X$ and $Z$ on the two factors separately, for every $a, b, c, d$.

```wl
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]] . Flatten[KroneckerProduct[{a, b}, {c, d}]] == Flatten[KroneckerProduct[PauliMatrix[1] . {a, b}, PauliMatrix[3] . {c, d}]]
```

**QF** : `QuantumTensorProduct` performs the same join but keeps the factor structure. The result is a
two-qubit `QuantumState`, shown as its summary box.

```wl
tp = QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]
```

Its `"StateVector"` reads the amplitudes back out.

```wl
tp["StateVector"] // Normal
```

The object route and the hand computation agree.

```wl
Normal[tp["StateVector"]] == Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

It joins operators the same way; here $X\otimes Z$ kept as a `QuantumOperator`, shown as its summary box.

```wl
xz = QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]
```

Its matrix is the Kronecker product built by hand.

```wl
Normal[xz["Matrix"]] == KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]
```

The compatibility rule holds at the object level too: the joint operator applied to the joint state is
the tensor product of the two factor images, equal as states.

```wl
xz[tp] == QuantumTensorProduct[QuantumOperator["X"][QuantumState[{a, b}]], QuantumOperator["Z"][QuantumState[{c, d}]]]
```

Both give the amplitude vector $\{ac, ad, bc, bd\}$, the pairwise products of the two factors' amplitudes,
and the operator $X\otimes Z$, and both confirm the compatibility rule for arbitrary amplitudes: the
composite factorizes into its parts, and the framework object additionally knows it is a two-qubit system
(its factors, order, and basis), a matrix being one representation read from it.

### 4.2 [BSc] How do I take a partial trace of a two-party state and obtain a reduced density matrix?

The state of subsystem $A$ alone is the reduced density matrix $\rho_A = \mathrm{Tr}_B[\rho_{AB}]$, formed
by tracing out $B$: an average over everything in $B$ that is not observed, keeping only what $A$ can
still predict on its own. For a product state $\rho_A$ stays pure; for an *entangled* state it comes out
mixed, which is the operational signature of entanglement. Take the one-parameter family
$|\psi(\lambda)\rangle = \cos\lambda\,|00\rangle + \sin\lambda\,|11\rangle$, which runs from the product
state $|00\rangle$ at $\lambda = 0$ to the Bell state $|\Phi^+\rangle$ at $\lambda = \pi/4$; every two-qubit
pure state is one of these up to local unitaries, its Schmidt form (Part 11, 11.3).

**WL** : the amplitudes of $|\psi(\lambda)\rangle$.

```wl
psi = {Cos[\[Lambda]], 0, 0, Sin[\[Lambda]]}
```

Its density matrix is the outer product $|\psi\rangle\langle\psi|$.

```wl
rho = KroneckerProduct[psi, Conjugate[psi]];
```

Reshape the $4\times4$ matrix into a $2\times2\times2\times2$ tensor, whose four legs are
$\mathrm{ket}_A, \mathrm{ket}_B, \mathrm{bra}_A, \mathrm{bra}_B$, and trace out $B$ by contracting its two
legs, indices $\{2,4\}$; $\lambda$ is a real angle.

```wl
rhoA = Simplify[TensorContract[ArrayReshape[rho, {2, 2, 2, 2}], {{2, 4}}], \[Lambda] \[Element] Reals]
```

Its purity $\mathrm{Tr}[\rho_A^2]$ measures how much of the pair's structure the reduction kept.

```wl
purityA = Simplify[Tr[rhoA . rhoA], \[Lambda] \[Element] Reals]
```

At the two ends of the family, the product state and the Bell state.

```wl
purityA /. {{\[Lambda] -> 0}, {\[Lambda] -> Pi/4}}
```

**QF** : `QuantumPartialTrace` traces out the listed qubit in one call and returns the reduced state as an
object, shown as its summary box.

```wl
redA = QuantumPartialTrace[QuantumState[psi], {2}]
```

Its `"Purity"` is the same function of $\lambda$.

```wl
Simplify[redA["Purity"], \[Lambda] \[Element] Reals]
```

The hand computation and the framework agree.

```wl
Simplify[Normal[redA["DensityMatrix"]] == rhoA, \[Lambda] \[Element] Reals]
```

The two ends of the family are named states: the Bell pair `"PhiPlus"` and a product ket such as `"01"`.

```wl
QuantumPartialTrace[QuantumState[#], {2}]["Purity"] & /@ {"PhiPlus", "01"}
```

Both give $\rho_A = \mathrm{diag}(\cos^2\lambda, \sin^2\lambda)$ with purity $\cos^4\lambda + \sin^4\lambda$:
the reduced state carries the Schmidt weights on its diagonal, pure at $\lambda = 0$ where the pair is a
product and the maximally mixed $I/2$ at $\lambda = \pi/4$ where the pair is maximally entangled, so the
mixedness of the reduced state is exactly the entanglement discarded along with $B$.

**A heterogeneous register.** The reshape-and-contract recipe does not care that the qudits differ in
size. In the two-qubit case each qubit owned a ket leg and a bra leg; the same holds for qudits of any
dimensions. Take three qudits of dimensions $(d_1,d_2,d_3) = (3,2,5)$ in the same family,
$|\psi(\lambda)\rangle = \cos\lambda\,|0,0,0\rangle + \sin\lambda\,|1,1,1\rangle$, living in
$\mathbb C^{3}\otimes\mathbb C^{2}\otimes\mathbb C^{5} = \mathbb C^{30}$. The joint density matrix is
$30\times30$; reshaped to a rank-6 array of shape $(3,2,5,3,2,5)$, its first three axes are the ket (row)
legs of particles $1,2,3$ and the last three the bra (column) legs, so particle $p$ owns ket-leg $p$ and
bra-leg $p+3$, each of size $d_p$. Tracing out a set of particles contracts their leg pairs, and the legs
that survive, grouped into rows and columns, are the reduced density matrix on the remaining qudits, of
size $\prod_{\text{kept}} d_q$.

**WL** : build the two basis kets from the dimension list, one level $k$ on every qudit.

```wl
kets = Table[Flatten[KroneckerProduct @@ (UnitVector[#, k] & /@ {3, 2, 5})], {k, 2}];
```

Combine them with the Schmidt weights.

```wl
psi = Cos[\[Lambda]] kets[[1]] + Sin[\[Lambda]] kets[[2]];
```

Reshape its density matrix to the $(3,2,5,3,2,5)$ tensor.

```wl
rhoT = ArrayReshape[KroneckerProduct[psi, Conjugate[psi]], {3, 2, 5, 3, 2, 5}];
```

Remove the dimension-$2$ qudit (particle $2$) by contracting legs $\{2,5\}$, then group the surviving ket
legs into rows and bra legs into columns.

```wl
rho13 = Flatten[TensorContract[rhoT, {{2, 5}}], {{1, 2}, {3, 4}}];
```

Its size is the product of the kept dimensions.

```wl
Dimensions[rho13]
```

It is a state: the input is normalized, and the reduction has unit trace, is Hermitian, and has a
nonnegative spectrum, the two Schmidt weights and zeros.

```wl
Simplify[{Norm[psi], Tr[rho13], HermitianMatrixQ[rho13], Eigenvalues[rho13]}, \[Lambda] \[Element] Reals]
```

Its purity is the same function of $\lambda$ as the qubit pair's.

```wl
Simplify[Tr[rho13 . rho13], \[Lambda] \[Element] Reals]
```

Remove two qudits at once: particle $p$ owns legs $p$ and $p+3$, so each pair of particles contracts two
leg pairs, and the pairs $\{1,2\}, \{1,3\}, \{2,3\}$ leave the dimension-$5$, dimension-$2$, and
dimension-$3$ qudit.

```wl
pairs = Simplify[TensorContract[rhoT, Transpose[{#, # + 3}]] & /@ Subsets[Range[3], {2}], \[Lambda] \[Element] Reals]
```

**QF** : the framework needs only the amplitudes and the dimension list to build the same register.

```wl
qs = QuantumState[psi, {3, 2, 5}]
```

Trace out the dimension-$2$ qudit once, as an object that keeps the $(3,5)$ factorization.

```wl
red2 = QuantumPartialTrace[qs, {2}]
```

Its dimensions are reported as $(3,5)$ rather than the flattened $15$.

```wl
red2["Dimensions"]
```

Its purity is the hand value.

```wl
Simplify[red2["Purity"], \[Lambda] \[Element] Reals]
```

The hand matrix and the framework's agree.

```wl
Simplify[Normal[red2["DensityMatrix"]] == rho13, \[Lambda] \[Element] Reals]
```

The dimensions left by each single removal.

```wl
QuantumPartialTrace[qs, #]["Dimensions"] & /@ Subsets[Range[3], {1}]
```

The three pair removals equal the hand contractions.

```wl
Simplify[(Normal[QuantumPartialTrace[qs, #]["DensityMatrix"]] & /@ Subsets[Range[3], {2}]) == pairs, \[Lambda] \[Element] Reals]
```

The purity of all six reductions, the three single removals followed by the three pairs.

```wl
Simplify[QuantumPartialTrace[qs, #]["Purity"] & /@ Subsets[Range[3], {1, 2}], \[Lambda] \[Element] Reals]
```

Both routes give reductions of dimension $\prod_{\text{kept}} d_q$, $(2,5), (3,5), (3,2)$ for the single
removals and $5, 2, 3$ for the pairs, and every one of the six carries the same spectrum, the two Schmidt
weights $\cos^2\lambda, \sin^2\lambda$ padded with zeros, hence the same purity $\cos^4\lambda + \sin^4\lambda$:
the state is entangled across every cut with the same Schmidt weights, and which qudit sizes sit on
either side of the cut does not enter. The contraction is the same construction at any sizes; only the
leg dimensions change, and where the hand computation sees a flattened square matrix, the framework keeps
the qudit factorization and reports its dimensions accordingly.
