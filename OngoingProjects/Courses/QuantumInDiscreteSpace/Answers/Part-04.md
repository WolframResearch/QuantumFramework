## Part 4. Composite systems: tensor product and partial trace

Two quantum systems join by the tensor product, and the reverse question, what one subsystem looks like once
we ignore the other, is answered by the partial trace. In other words, this Part is about going up (stacking
two systems into one) and coming back down (reading a piece of a joint state on its own). We build each
direction from primitives and then confirm the framework agrees: the tensor product is a Kronecker product,
and the partial trace is a reshape-and-contract.

### 4.1 [BSc] How do I form the tensor product of states and of operators?

Two subsystems combine by the tensor product: states as $|\psi\rangle\otimes|\phi\rangle$, operators as
$A\otimes B$, and the two are compatible, $(A\otimes B)(|\psi\rangle\otimes|\phi\rangle) =
(A|\psi\rangle)\otimes(B|\phi\rangle)$. In other words, the joint system carries one amplitude for every pair
of subsystem labels, so the composite lives in the $d_1 d_2$-dimensional product space. In the computational
basis this tensor product is exactly the Kronecker product.

Form the tensor product of two single-qubit amplitude vectors $\{a,b\}$ and $\{c,d\}$, flattened into one
length-4 vector:

```wl
Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

As one can see, the four amplitudes are the pairwise products $\{ac, ad, bc, bd\}$, one for each
computational label $00, 01, 10, 11$.

Operators combine the same way. Form $X\otimes Z$ as the Kronecker product of their matrices:

```wl
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]
```

The framework performs the same join but keeps the factor structure. Combine the two states with
`QuantumTensorProduct`:

```wl
tp = QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]
```

The result is a two-qubit `QuantumState`; its summary box reports a $2\times2$-dimensional system, not just a
bare vector. Read its amplitude vector back out:

```wl
tp["StateVector"] // Normal
```

Confirm the object route and the hand computation agree:

```wl
tp["StateVector"] // Normal == Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

It joins operators the same way; here $X\otimes Z$ kept as a `QuantumOperator`:

```wl
QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]
```

Therefore both routes give the amplitude vector $\{ac, ad, bc, bd\}$ and the operator $X\otimes Z$: the
composite factorizes into its parts, and the framework object additionally knows it is a two-qubit system
(its factors, order, and basis), a matrix being just one representation read from it.

### 4.2 [BSc] How do I take a partial trace of a two-party state and obtain a reduced density matrix?

The state of subsystem $A$ alone is the reduced density matrix $\rho_A = \mathrm{Tr}_B[\rho_{AB}]$, formed by
tracing out $B$. In other words, we average over everything we cannot see in $B$, keeping only what $A$ can
still predict on its own. For a product state $\rho_A$ stays pure; for an *entangled* state it comes out
mixed, which is the operational signature of entanglement.

Take the Bell state $|\Phi^+\rangle = (|00\rangle + |11\rangle)/\sqrt2$ and trace out qubit $B$ by hand.
First define its amplitudes:

```wl
bell = {1, 0, 0, 1}/Sqrt[2]
```

Its density matrix is the outer product $|\Phi^+\rangle\langle\Phi^+|$:

```wl
rho = KroneckerProduct[bell, Conjugate[bell]]
```

Reshape the $4\times4$ matrix into a $2\times2\times2\times2$ tensor, whose four legs are
$\mathrm{ket}_A, \mathrm{ket}_B, \mathrm{bra}_A, \mathrm{bra}_B$, and trace out $B$ by contracting its two
legs, indices $\{2,4\}$:

```wl
rhoA = TensorContract[ArrayReshape[rho, {2, 2, 2, 2}], {{2, 4}}]
```

As one can see, the reduced state is $\tfrac12 I$: maximally mixed, retaining none of the pair's structure.

The framework traces out the listed qubit in one call. Reduce the same Bell state, keeping the object form:

```wl
redA = QuantumPartialTrace[QuantumState["PhiPlus"], {2}]
```

Its purity is below one, the quantitative signature of the entanglement discarded with $B$:

```wl
redA["Purity"]
```

Confirm the hand computation and the framework agree:

```wl
redA["DensityMatrix"] // Normal == rhoA
```

A *product* input, by contrast, reduces to a pure state. The product ket $|01\rangle$ is written directly as
`QuantumState["01"]`; trace out qubit 2 and read the purity:

```wl
QuantumPartialTrace[QuantumState["01"], {2}]["Purity"]
```

As expected, purity $1$: a product state carries no entanglement, so nothing is lost in the reduction. The
mixedness of the reduced state is exactly the entanglement discarded along with $B$.

**A heterogeneous register.** The reshape-and-contract recipe does not care that the qudits differ in size.
Recall that the two-qubit case reshaped the density matrix so each qubit owned a ket leg and a bra leg; the
same holds for qudits of any dimensions. Take three qudits of dimensions $(d_1,d_2,d_3) = (3,2,5)$ in the
entangled state $|\psi\rangle = \tfrac1{\sqrt2}(|0,0,0\rangle + |1,1,1\rangle)$, living in
$\mathbb C^{3}\otimes\mathbb C^{2}\otimes\mathbb C^{5} = \mathbb C^{30}$.

All the bookkeeping is in the indices, so let us be explicit about them. The joint density matrix is
$30\times30$; reshape it to a rank-6 array of shape $(3,2,5,3,2,5)$. The first three axes are the ket (row)
legs of particles $1,2,3$; the last three are the bra (column) legs. In short, particle $p$ owns two legs of
size $d_p$: ket-leg $p$ and bra-leg $p+3$. Tracing out particle $p$ then means contracting those two legs
together, and the legs that survive, reshaped into a square matrix of size $\prod_{\text{kept}} d_q$, are the
reduced density matrix on the remaining qudits.

Encode $|\psi\rangle$ as a length-30 vector, each term a Kronecker product of the three qudit levels:

```wl
psi = Normalize[Flatten[KroneckerProduct[{1, 0, 0}, {1, 0}, {1, 0, 0, 0, 0}]] + Flatten[KroneckerProduct[{0, 1, 0}, {0, 1}, {0, 1, 0, 0, 0}]]];
```

Reshape its density matrix to the $(3,2,5,3,2,5)$ tensor:

```wl
rhoT = ArrayReshape[KroneckerProduct[psi, Conjugate[psi]], {3, 2, 5, 3, 2, 5}];
```

The framework needs only the amplitudes and the dimension list to build the same register:

```wl
qs = QuantumState[psi, {3, 2, 5}]
```

Now trace out one qudit by contracting its single ket-bra leg pair. Remove the dim-$2$ qudit (particle $2$)
by contracting legs $\{2,5\}$:

```wl
rho13 = ArrayReshape[TensorContract[rhoT, {{2, 5}}], {15, 15}];
```

The kept legs $(3,5)$ flatten to a $15\times15$ matrix:

```wl
Dimensions[rho13]
```

and its purity is below one:

```wl
Tr[rho13 . rho13]
```

The framework traces the same qudit from its indices, and because it keeps the $(3,5)$ factorization it
reports the reduced dimensions as $(3,5)$ rather than the flattened $15$:

```wl
QuantumPartialTrace[qs, {2}]["Dimensions"]
```

with the same purity:

```wl
QuantumPartialTrace[qs, {2}]["Purity"]
```

The leg pattern is uniform: particle $1$ (dim $3$) is legs $\{1,4\}$, particle $3$ (dim $5$) is $\{3,6\}$.
Trace out each single qudit in turn and read the dimensions each reduction leaves behind:

```wl
Table[QuantumPartialTrace[qs, {p}]["Dimensions"], {p, 3}]
```

As one can see, removing one qudit leaves the other two: $(2,5), (3,5), (3,2)$ as we drop particles
$1, 2, 3$.

Next trace out two qudits by contracting both leg pairs at once. Remove the dim-$3$ and dim-$2$ qudits
(particles $\{1,2\}$, legs $\{1,4\}$ and $\{2,5\}$), leaving only the dim-$5$ qudit:

```wl
TensorContract[rhoT, {{1, 4}, {2, 5}}]
```

The framework agrees:

```wl
QuantumPartialTrace[qs, {1, 2}]["DensityMatrix"] // Normal
```

Remove instead the dim-$3$ and dim-$5$ qudits (particles $\{1,3\}$, legs $\{1,4\}$ and $\{3,6\}$), leaving
the dim-$2$ qudit maximally mixed:

```wl
TensorContract[rhoT, {{1, 4}, {3, 6}}]
```

and again the framework matches:

```wl
QuantumPartialTrace[qs, {1, 3}]["DensityMatrix"] // Normal
```

Therefore every reduction has dimension $\prod_{\text{kept}} d_q$: $15, 10, 6$ when one qudit is removed,
$5, 2, 3$ when two are. Each comes out mixed at purity $\tfrac12$, because
$(|0,0,0\rangle + |1,1,1\rangle)/\sqrt2$ is entangled across every cut. The contraction is the same
construction at any sizes; only the leg dimensions change, and where the hand computation sees a flattened
square matrix, the framework keeps the qudit factorization and reports its dimensions accordingly.
