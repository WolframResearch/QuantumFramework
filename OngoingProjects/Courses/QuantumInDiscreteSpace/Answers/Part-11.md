## Part 11. Composite systems and entanglement

A joint state that is not a product carries correlations no local description can reproduce. This Part builds
the standard entangled states, reads their correlations off a measurement, quantifies how entangled a *pure*
state is, then moves to the harder mixed-state question, where separability is no longer visible in a
spectrum and needs a criterion (positive partial transpose), a measure (concurrence, negativity), or an
observable that certifies it (a witness). The last questions push on the limits: which conversions local
operations allow, how entanglement is shared among three parties, states that are entangled yet
undistillable, and the exchange symmetry of identical particles. Each answer builds every state and operator
it uses from scratch, out of the single-qubit kets $|0\rangle$ and $|1\rangle$.

### 11.1 [BSc] How do I build the Bell, GHZ, and W states?

The Bell states are the four maximally entangled two-qubit states and form an orthonormal basis of the
four-dimensional joint space. For three qubits there are two inequivalent kinds: the GHZ state, an equal
superposition of the all-zeros and all-ones configurations, and the W state, a single excitation delocalized
over the sites. Both generalize to any number of qubits.

**WL** : write the kets as the two basis vectors and combine them with `KroneckerProduct`. The seed is
$|\Phi^+\rangle = (|00\rangle+|11\rangle)/\sqrt2$.

```wl
phiPlus = With[{zero = {1, 0}, one = {0, 1}}, Normalize @ Flatten[KroneckerProduct[zero, zero] + KroneckerProduct[one, one]]]
```

A Pauli acting on the second qubit alone carries that seed to each of the other three, so the whole basis
comes from one construction (with `PauliMatrix[0]` the identity).

```wl
bell = Table[KroneckerProduct[IdentityMatrix[2], PauliMatrix[k]] . phiPlus, {k, 0, 3}]
```

The Paulis are unitary, so the four images of an orthonormal seed are again orthonormal; their Gram matrix
makes that explicit.

```wl
Outer[Conjugate[#1] . #2 &, bell, bell, 1] == IdentityMatrix[4]
```

The GHZ state of $n$ qubits is the all-zeros ket plus the all-ones ket, each a Kronecker product of $n$
copies.

```wl
ghz = With[{n = 3, zero = {1, 0}, one = {0, 1}}, Normalize @ Flatten[KroneckerProduct @@ ConstantArray[zero, n] + KroneckerProduct @@ ConstantArray[one, n]]]
```

The W state puts one excitation on each site in turn, which is exactly the set of distinct arrangements of
the multiset of $n-1$ zeros and one excitation, so `Permutations` generates the terms.

```wl
w = With[{n = 3, zero = {1, 0}, one = {0, 1}}, Normalize @ Total[Flatten /@ (KroneckerProduct @@@ Permutations[Append[ConstantArray[zero, n - 1], one]])]]
```

**QF** : all three are named states. The Bell family:

```wl
bellQF = QuantumState /@ {"PhiPlus", "PsiPlus", "PsiMinus", "PhiMinus"}
```

Compare them with the Pauli-generated list. State equality in the framework is equality of rays, up to a
global phase, which matters here: the $\sigma_y$ image carries a factor $i$, and a global phase is
unobservable (Part 1, 1.4), so it is the same state.

```wl
MapThread[QuantumState[#1] == #2 &, {bell, bellQF}]
```

The GHZ state carries its qubit count as an argument:

```wl
ghzQF = QuantumState["GHZ"[3]]
```

```wl
Normal[ghzQF["StateVector"]] == ghz
```

and so does the W state:

```wl
wQF = QuantumState["W"[3]]
```

```wl
Normal[wQF["StateVector"]] == w
```

The object also renders itself in Dirac notation, which shows the single-excitation structure directly:

```wl
wQF["Formula"]
```

Both routes give the same three families. The framework objects additionally know their factorization into
qubits, so they can be partially traced, measured on one wire, or fed to a circuit, whereas a bare amplitude
vector has to be reshaped by hand each time.

### 11.2 [BSc] How do I measure a Bell pair jointly and exhibit the perfect correlations?

Measuring both qubits of a Bell pair in the computational basis gives $00$ or $11$ with probability one half
each and never $01$ or $10$: the outcomes are individually random but perfectly correlated. The correlation
survives a change of basis, which is what the three parity correlators
$\langle\sigma_k\otimes\sigma_k\rangle$ record.

**WL** : build the pair from the kets.

```wl
phiPlus = With[{zero = {1, 0}, one = {0, 1}}, Normalize @ Flatten[KroneckerProduct[zero, zero] + KroneckerProduct[one, one]]];
```

The joint Born probabilities of the four computational outcomes are the squared moduli of the four
amplitudes.

```wl
Abs[phiPlus]^2
```

Each qubit on its own is unbiased: summing the joint distribution over the second label leaves the marginal
of the first.

```wl
Total[Partition[Abs[phiPlus]^2, 2]]
```

The correlators are expectation values of the three parity observables $\sigma_k\otimes\sigma_k$.

```wl
Table[Conjugate[phiPlus] . KroneckerProduct[PauliMatrix[k], PauliMatrix[k]] . phiPlus, {k, 3}]
```

Measuring both qubits in the $X$ basis instead means rotating by a Hadamard on each side first and reading
the same four probabilities.

```wl
Abs[KroneckerProduct[HadamardMatrix[2], HadamardMatrix[2]] . phiPlus]^2
```

**QF** : measuring the wire list is the joint computational measurement; it returns a `QuantumMeasurement`
object.

```wl
qm = QuantumMeasurementOperator[{1, 2}][QuantumState["PhiPlus"]]
```

Its `"Probabilities"` are keyed by the outcome pair, so the labels travel with the numbers:

```wl
qm["Probabilities"]
```

The correlators come from measuring the two-qubit parity observables, whose names are the two letters:

```wl
Table[QuantumMeasurementOperator[QuantumOperator[StringRepeat[s, 2]]][QuantumState["PhiPlus"]]["Mean"], {s, {"X", "Y", "Z"}}]
```

```wl
Values[qm["Probabilities"]] == Abs[phiPlus]^2
```

Both give one half on $00$ and on $11$ and zero on the mixed outcomes, marginals of one half on each side,
and correlators $+1$ for $XX$ and $ZZ$ and $-1$ for $YY$. The $X$-basis measurement is just as sharply
correlated as the $Z$-basis one, which is what no product state can imitate: a product state cannot be
perfectly correlated in two complementary bases at once.

### 11.3 [BSc] How do I compute the Schmidt decomposition and Schmidt rank of a bipartite pure state?

Reshaping the amplitudes of a bipartite state into the coefficient matrix $C$, with one index per party,
turns the Schmidt decomposition into a singular value decomposition: the singular values of $C$ are the
Schmidt coefficients, their squares sum to one, and the left and right singular vectors are the two local
bases in which the state is a single diagonal sum. The Schmidt rank is the number of nonzero coefficients,
and it is one exactly for a product state.

**WL** : take a fully general two-qubit state and reshape its four amplitudes into the coefficient matrix.

```wl
c = ArrayReshape[{\[FormalA], \[FormalB], \[FormalC], \[FormalD]}, {2, 2}]
```

The Schmidt coefficients are the singular values, that is the square roots of the eigenvalues of
$C\,C^\dagger$.

```wl
Simplify @ Sqrt @ Eigenvalues[c . ConjugateTranspose[c]]
```

The determinant sits under the inner square root, so the smaller coefficient vanishes exactly when the
determinant does. The Schmidt rank is therefore the rank of the coefficient matrix.

```wl
{MatrixRank[c], Det[c]}
```

Imposing a vanishing determinant collapses the rank to one, the product case.

```wl
MatrixRank[c /. \[FormalD] -> \[FormalB] \[FormalC]/\[FormalA]]
```

For a worked instance take the state whose three terms are $|00\rangle$, $|01\rangle$ and $|10\rangle$,
whose two Schmidt coefficients come out unequal.

```wl
psi = With[{zero = {1, 0}, one = {0, 1}}, Normalize @ Total[Flatten /@ {KroneckerProduct[zero, zero], KroneckerProduct[zero, one], KroneckerProduct[one, zero]}]];
```

```wl
SingularValueList[ArrayReshape[psi, {2, 2}]]
```

```wl
Simplify @ Total[SingularValueList[ArrayReshape[psi, {2, 2}]]^2]
```

**QF** : `"SchmidtBasis"` returns the same state re-expressed in its own Schmidt bases, as a state object
carrying those bases.

```wl
sb = QuantumState[psi]["SchmidtBasis"]
```

In that basis the amplitude vector is diagonal: the coefficients sit on the matched pairs and everything else
is zero.

```wl
Normal[sb["StateVector"]]
```

The Schmidt rank is the number of surviving amplitudes.

```wl
Count[Normal[sb["StateVector"]], Except[0]]
```

`"SchmidtDecompose"` renders the decomposition as the sum of Kronecker products itself:

```wl
QuantumState[psi]["SchmidtDecompose"]
```

A product input collapses to a single term, the rank-one case:

```wl
Count[Normal[QuantumState[With[{zero = {1, 0}, one = {0, 1}}, Flatten @ KroneckerProduct[Normalize[zero + one], zero]]]["SchmidtBasis"]["StateVector"]], Except[0]]
```

```wl
DeleteCases[Normal[sb["StateVector"]], 0] == SingularValueList[ArrayReshape[psi, {2, 2}]]
```

Both give the coefficients $\sqrt{(3\pm\sqrt5)/6}$ and rank $2$, dropping to rank $1$ on a product state.
Reshape-plus-SVD returns bare numbers, while the framework returns the state itself in its Schmidt bases, so
the local bases that realize the decomposition come along with the coefficients.

### 11.4 [BSc] How do I quantify pure-state entanglement by the entropy of entanglement?

For a pure bipartite state the entanglement is measured by the von Neumann entropy of either reduced state,
which in terms of the Schmidt coefficients is the Shannon entropy of their squares. It is zero for a product
state and $\log_2 d$ for a maximally entangled one. Every two-qubit pure state is, up to local unitaries,
$\cos\theta\,|00\rangle + \sin\theta\,|11\rangle$, so that one-parameter family carries the whole answer.

**WL** : build the family from the kets.

```wl
ket = With[{zero = {1, 0}, one = {0, 1}}, Cos[\[Theta]] Flatten[KroneckerProduct[zero, zero]] + Sin[\[Theta]] Flatten[KroneckerProduct[one, one]]]
```

Reduce it by tracing out qubit $B$, contracting legs $2$ and $4$ of the reshaped density matrix as in Part 4,
4.2.

```wl
rhoA = Simplify[TensorContract[ArrayReshape[KroneckerProduct[ket, Conjugate[ket]], {2, 2, 2, 2}], {{2, 4}}], \[Theta] \[Element] Reals]
```

Its two eigenvalues are the Schmidt weights; feed them to the entropy sum.

```wl
s = Simplify[-Total[# Log2[#] & /@ Eigenvalues[rhoA]], 0 < \[Theta] < Pi/2]
```

Maximizing over the family locates the maximally entangled point.

```wl
Simplify @ Maximize[{s, 0 < \[Theta] < Pi/2}, \[Theta]]
```

**QF** : trace out qubit $2$ and ask the reduced state for its `"VonNeumannEntropy"`, which comes back as a
`Quantity` in bits.

```wl
redA = QuantumPartialTrace[QuantumState[ket], {2}]
```

```wl
Simplify[redA["VonNeumannEntropy"], 0 < \[Theta] < Pi/2]
```

```wl
Simplify[redA["VonNeumannEntropy"] == Quantity[s, "Bits"], 0 < \[Theta] < Pi/2]
```

The packaged monotone reads the same number off the state directly, once the Schmidt weights are explicit
numbers. Unlike `"Concurrence"`, its pure-state branch still selects the weights with a bare `# > 0`, which
no symbolic weight satisfies, so a symbolic angle silently returns zero bits instead of refusing. Give it a
definite angle:

```wl
QuantumEntanglementMonotone[QuantumState[ket /. \[Theta] -> Pi/4], "EntanglementEntropy"]
```

Both give the entropy of the two weights $\cos^2\theta$ and $\sin^2\theta$, rising from zero at
$\theta = 0$, where the state is the product $|00\rangle$, to exactly one bit at $\theta = \pi/4$, one Bell
pair's worth of entanglement. The reduced state is where the entanglement is read: a pure whole with a mixed
part is the signature, and the amount of mixedness is the amount of entanglement.

### 11.5 [BSc] How do I verify no-signalling through the invariance of reduced states?

Whatever Bob does to his qubit alone, Alice's reduced state is unchanged: for a unitary because the unitary
and its adjoint meet under the trace over Bob's factor, and for a measurement Bob does not report because the
non-selective Lüders sum (Part 5, 5.4) contracts to the same thing. Since Alice's statistics are fixed by her
reduced state, no choice Bob makes can carry information to her. The claim is linear in the joint state, so a
generic $4\times4$ array with nothing assumed about it already settles it for every state.

**WL** : take a generic two-qubit operator and its reduction over $B$.

```wl
r = Array[Subscript[\[FormalR], ##] &, {4, 4}];
```

```wl
rhoA = TensorContract[ArrayReshape[r, {2, 2, 2, 2}], {{2, 4}}]
```

A general single-qubit unitary is the $Z$-$Y$-$Z$ Euler product of Part 9, 9.4, with three free angles.

```wl
u = MatrixExp[-I \[Phi] PauliMatrix[3]/2] . MatrixExp[-I \[Theta] PauliMatrix[2]/2] . MatrixExp[-I \[Lambda] PauliMatrix[3]/2];
```

Act with it on qubit $B$ only and reduce again.

```wl
FullSimplify[TensorContract[ArrayReshape[KroneckerProduct[IdentityMatrix[2], u] . r . ConjugateTranspose[KroneckerProduct[IdentityMatrix[2], u]], {2, 2, 2, 2}], {{2, 4}}] == rhoA, {\[Theta], \[Phi], \[Lambda]} \[Element] Reals]
```

Now let Bob measure along an arbitrary axis without reporting the outcome, using the two projectors
$(I \pm \hat n\cdot\vec\sigma)/2$, and reduce the resulting mixture.

```wl
lueders = Sum[KroneckerProduct[IdentityMatrix[2], proj] . r . KroneckerProduct[IdentityMatrix[2], proj], {proj, (IdentityMatrix[2] + # FromSphericalCoordinates[{1, \[Theta]1, \[Phi]1}] . PauliMatrix[{1, 2, 3}])/2 & /@ {1, -1}}];
```

```wl
FullSimplify[TensorContract[ArrayReshape[lueders, {2, 2, 2, 2}], {{2, 4}}] == rhoA, {\[Theta]1, \[Phi]1} \[Element] Reals]
```

**QF** : the same generic array as a two-qubit state object, with the framework's own general single-qubit
unitary `"U3"` placed on wire $2$.

```wl
qr = QuantumState[r, {2, 2}]
```

```wl
FullSimplify[Normal[QuantumPartialTrace[QuantumOperator["U3"[\[Theta], \[Phi], \[Lambda]], {2}][qr], {2}]["DensityMatrix"]] == Normal[QuantumPartialTrace[qr, {2}]["DensityMatrix"]], {\[Theta], \[Phi], \[Lambda]} \[Element] Reals]
```

For the unreported measurement, build the two projectors as operators on wire $2$ and sum their action on
the state.

```wl
lu = Total[#[qr] & /@ ((QuantumOperator["I", {2}] + # QuantumOperator[FromSphericalCoordinates[{1, \[Theta]1, \[Phi]1}] . PauliMatrix[{1, 2, 3}], {2}])/2 & /@ {1, -1})]
```

```wl
FullSimplify[Normal[QuantumPartialTrace[lu, {2}]["DensityMatrix"]] == Normal[QuantumPartialTrace[qr, {2}]["DensityMatrix"]], {\[Theta]1, \[Phi]1} \[Element] Reals]
```

Every check returns `True`, for an arbitrary joint state, an arbitrary local unitary, and an arbitrary
measurement axis. Entanglement therefore correlates outcomes without transmitting anything: Alice sees
exactly the same statistics whether Bob acts, measures, or does nothing, and the correlations only appear
when the two records are later compared.

### 11.6 [MSc] How do I test a mixed state for entanglement by the positive-partial-transpose (Peres-Horodecki) criterion?

Transposing only Bob's factor of a separable state sends each of its product terms to another product term,
so the result is still a legitimate state and its spectrum stays nonnegative. A negative eigenvalue of the
partial transpose therefore certifies entanglement. Take the family
$\rho(q) = q\,|\Phi^+\rangle\langle\Phi^+| + (1-q)\,|01\rangle\langle01|$, separable at $q = 0$ and a Bell
state at $q = 1$.

**WL** : build the family from the kets.

```wl
rho = With[{zero = {1, 0}, one = {0, 1}}, With[{bell = Normalize @ Flatten[KroneckerProduct[zero, zero] + KroneckerProduct[one, one]], ket01 = Flatten @ KroneckerProduct[zero, one]}, q KroneckerProduct[bell, bell] + (1 - q) KroneckerProduct[ket01, ket01]]]
```

The partial transpose on $B$ reshapes the matrix into its four legs and exchanges the two belonging to $B$,
positions $2$ and $4$.

```wl
pt = ArrayReshape[Transpose[ArrayReshape[rho, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]
```

```wl
eigs = Simplify @ Eigenvalues[pt]
```

Demanding that all four be nonnegative locates the entire positive-partial-transpose region of the family.

```wl
Reduce[And @@ Thread[eigs >= 0] && 0 <= q <= 1, q, Reals]
```

**QF** : `QuantumPartialTranspose` transposes the listed qudits and hands back a state object, whose spectrum
reads off directly.

```wl
ptQF = QuantumPartialTranspose[QuantumState[rho, {2, 2}], {2}]
```

```wl
ptQF["Eigenvalues"] === eigs
```

The packaged test wraps the same criterion into a boolean; its third argument selects the measure that
decides, the second being the bipartition.

```wl
Table[QuantumEntangledQ[QuantumState[rho /. q -> v, {2, 2}], Automatic, "Negativity"], {v, {0, 1/4, 1/2, 1}}]
```

Both give the spectrum $\{q/2,\ q/2,\ (1-q\pm\sqrt{1-2q+2q^2})/2\}$, whose third entry is negative for every
$q > 0$: the only member of the family with a positive partial transpose is the separable endpoint $q = 0$.
So an admixture of a Bell state, however small, is already detected. Note the logic is one-directional in
general: a negative eigenvalue proves entanglement, but a positive partial transpose proves separability only
in small dimensions (11.12a and 11.12b).

### 11.7 [MSc] How do I compute the concurrence and the negativity of a mixed two-qubit state?

Both turn the criterion of 11.6 into a number. The negativity is half the amount by which the trace norm of
the partial transpose exceeds one, that is the total weight of the negative part of its spectrum. The
concurrence is Wootters' formula: with $\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge \lambda_4$ the square roots
of the eigenvalues of $\rho$ times its spin flip
$(\sigma_y\otimes\sigma_y)\,\rho^*\,(\sigma_y\otimes\sigma_y)$, the concurrence is
$\max\{0,\ \lambda_1-\lambda_2-\lambda_3-\lambda_4\}$. Use the same family as in 11.6.

**WL** : the state, and its spin flip.

```wl
rho = With[{zero = {1, 0}, one = {0, 1}}, With[{bell = Normalize @ Flatten[KroneckerProduct[zero, zero] + KroneckerProduct[one, one]], ket01 = Flatten @ KroneckerProduct[zero, one]}, q KroneckerProduct[bell, bell] + (1 - q) KroneckerProduct[ket01, ket01]]];
```

```wl
flipped = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]] . Conjugate[rho] . KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]
```

The four $\lambda_i$ are the square roots of the spectrum of the product.

```wl
lambdas = Simplify[Sqrt @ Eigenvalues[rho . flipped], 0 < q < 1]
```

```wl
concurrence = Simplify[First[#] - Total[Rest[#]] &@ ReverseSort[lambdas], 0 < q < 1]
```

The negativity is half the excess of the partially transposed trace norm over one, so it is half the sum of
the absolute eigenvalues minus one.

```wl
negativity = Simplify[(Total[Abs[Eigenvalues @ ArrayReshape[Transpose[ArrayReshape[rho, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]]] - 1)/2, 0 < q < 1]
```

**QF** : `QuantumEntanglementMonotone` names both measures.

```wl
negQF = QuantumEntanglementMonotone[QuantumState[rho, {2, 2}], "Negativity"]
```

The framework reaches the value through the trace norm, so its expression carries nested radicals; the two
agree as numbers, which is what `==` decides.

```wl
Simplify[negQF == negativity, 0 < q < 1]
```

The logarithmic version, the log of that trace norm, is the additive relative of the same quantity:

```wl
Simplify[QuantumEntanglementMonotone[QuantumState[rho, {2, 2}], "LogNegativity"], 0 < q < 1]
```

The mixed-state concurrence resolves symbolically, so it matches the closed form found above for the whole
family at once rather than value by value.

```wl
Simplify[QuantumEntanglementMonotone[QuantumState[rho, {2, 2}], "Concurrence"], 0 < q < 1]
```

```wl
Simplify[QuantumEntanglementMonotone[QuantumState[rho, {2, 2}], "Concurrence"] == concurrence, 0 < q < 1]
```

Both give concurrence exactly $q$, linear in the Bell admixture and reaching $1$ on the Bell state, and
negativity $(\sqrt{1-2q+2q^2} - (1-q))/2$, which is likewise positive for every $q > 0$ but smaller, reaching
one half at $q = 1$. The two measures vanish together and order the family the same way, yet they are
different functions: neither is a rescaling of the other, and only the concurrence reaches $1$ on a maximally
entangled state.

### 11.8 [MSc] How do I construct an entanglement witness for a given entangled state?

A witness is a Hermitian operator $W$ whose expectation value is nonnegative on every separable state and
negative on the target, so one measured number certifies entanglement. When the partial transpose of the
target has a negative eigenvalue with eigenvector $|e\rangle$, the choice
$W = (|e\rangle\langle e|)^{T_B}$ works: transposing back inside the trace turns
$\mathrm{Tr}[W\rho]$ into $\langle e|\rho^{T_B}|e\rangle$, which is that negative eigenvalue, while on a
product state $W$ cannot go negative. Build it for the Bell state $|\Phi^+\rangle$.

**WL** : the target state, from the kets.

```wl
rho = With[{zero = {1, 0}, one = {0, 1}}, KroneckerProduct[#, #] &@ Normalize @ Flatten[KroneckerProduct[zero, zero] + KroneckerProduct[one, one]]];
```

```wl
es = Eigensystem[ArrayReshape[Transpose[ArrayReshape[rho, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]]
```

Pick the eigenvector belonging to the negative eigenvalue.

```wl
e = Normalize @ First @ Pick[es[[2]], Negative /@ es[[1]]]
```

Transposing its projector back on the $B$ factor gives the witness.

```wl
wit = ArrayReshape[Transpose[ArrayReshape[KroneckerProduct[e, Conjugate[e]], {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]
```

On the target it returns the negative eigenvalue itself.

```wl
Tr[wit . rho]
```

On a general product state, each factor carried by its own pair of Bloch angles, the expectation value has a
closed form.

```wl
expW = FullSimplify @ ComplexExpand @ Re[Conjugate[#] . wit . # &@ Flatten @ KroneckerProduct[{Cos[\[Theta]a/2], Exp[I \[Phi]a] Sin[\[Theta]a/2]}, {Cos[\[Theta]b/2], Exp[I \[Phi]b] Sin[\[Theta]b/2]}]]
```

Searching that expression for a negative value over the whole product-state manifold comes up empty, so no
product state, and by convexity no separable state, is flagged.

```wl
Reduce[expW < 0 && 0 <= \[Theta]a <= Pi && 0 <= \[Theta]b <= Pi, {\[Theta]a, \[Phi]a, \[Theta]b, \[Phi]b}, Reals]
```

**QF** : the negative eigenvector here is the singlet, so the witness is the partial transpose of the
singlet, read as an operator off the transposed state object.

```wl
witQF = QuantumState["PsiMinus"]["Transpose", {2}]["Operator"]
```

```wl
Normal[witQF["Matrix"]] == wit
```

Measuring that observable on the target gives its expectation value:

```wl
QuantumMeasurementOperator[witQF][QuantumState["PhiPlus"]]["Mean"]
```

and on a general product state, assembled with `QuantumTensorProduct`, the same closed form:

```wl
FullSimplify[QuantumMeasurementOperator[witQF][QuantumTensorProduct[QuantumState[{Cos[\[Theta]a/2], Exp[I \[Phi]a] Sin[\[Theta]a/2]}], QuantumState[{Cos[\[Theta]b/2], Exp[I \[Phi]b] Sin[\[Theta]b/2]}]]]["Mean"] == expW, {\[Theta]a, \[Phi]a, \[Theta]b, \[Phi]b} \[Element] Reals]
```

Both give the same witness, with expectation value $-1/2$ on the Bell state, and on product states the
closed form $(1 - \hat n_a\cdot\hat n_b)/4$, which is nonnegative because a dot product of unit vectors never
exceeds one, and vanishes exactly on the aligned product states. So a single Hermitian observable, measurable
without any tomography, separates this entangled state from every separable one; the price is that the
witness is tailored to the target, and a different entangled state generally needs a different one.

### 11.9 [MSc] How do I decide whether a pure state can be converted to another by LOCC (Nielsen's theorem)?

Nielsen's theorem: one pure state can be turned into another with certainty by local operations and
classical communication if and only if the first state's vector of Schmidt weights, the squared Schmidt
coefficients sorted decreasing, is majorized by the second's. Majorized means every partial sum of the first
is at most the corresponding partial sum of the second. Flatter weights mean more entanglement, and LOCC can
only steepen them. Take two states of two qutrits with weights $\{1/2, 3/10, 1/5\}$ and $\{3/5, 1/5, 1/5\}$.

**WL** : a state with prescribed Schmidt weights is the sum of the matched basis kets $|kk\rangle$ weighted
by the square roots.

```wl
psi = Normalize @ Total[Sqrt[{1/2, 3/10, 1/5}] (Flatten[KroneckerProduct[#, #]] & /@ IdentityMatrix[3])];
```

```wl
phi = Normalize @ Total[Sqrt[{3/5, 1/5, 1/5}] (Flatten[KroneckerProduct[#, #]] & /@ IdentityMatrix[3])];
```

Recover the weights from the states themselves, as the squared singular values of the reshaped amplitudes
(11.3).

```wl
{lp, lq} = Sort[SingularValueList[ArrayReshape[#, {3, 3}]]^2, Greater] & /@ {psi, phi}
```

Majorization is the comparison of the two partial-sum sequences, checked in both directions.

```wl
{And @@ Thread[Accumulate[lp] <= Accumulate[lq]], And @@ Thread[Accumulate[lq] <= Accumulate[lp]]}
```

The entropies order consistently, since the allowed direction can only lose entanglement.

```wl
N[-Total[# Log2[#] & /@ #] & /@ {lp, lq}]
```

**QF** : the Schmidt weights are the probabilities the state reports in its own Schmidt basis.

```wl
qpsi = QuantumState[psi, {3, 3}]
```

```wl
{wp, wq} = Sort[Values[#["SchmidtBasis"]["Probability"]], Greater] & /@ {qpsi, QuantumState[phi, {3, 3}]}
```

```wl
{And @@ Thread[Accumulate[wp] <= Accumulate[wq]], And @@ Thread[Accumulate[wq] <= Accumulate[wp]]}
```

```wl
{wp, wq} == {lp, lq}
```

Both directions come back true then false: the partial sums $1/2, 4/5, 1$ of the first state never exceed
$3/5, 4/5, 1$ of the second, so the first can be converted into the second with certainty, while the reverse
fails at the very first sum, since $3/5$ exceeds $1/2$. The conversion is one-way. The entropies, about
$1.485$ and $1.371$ bits, order the same way, as they must, but that ordering is only necessary:
majorization is a whole chain of inequalities, and two states can be ordered by entropy while neither
majorizes the other, in which case no certain conversion exists in either direction.

### 11.10 [MSc] How do I verify monogamy of entanglement (the CKW inequality) using the three-qubit tangle?

Entanglement cannot be freely shared: if $A$ is strongly entangled with $B$ it has little left for $C$. The
Coffman-Kundu-Wootters inequality makes this quantitative. The tangle across the $A$ versus $BC$ cut is
$\tau_{A|BC} = 4\det\rho_A$, which for a qubit is the same as $2(1-\mathrm{Tr}\rho_A^2)$, and it is at least
the sum of the squared Wootters concurrences of the two pairs, $C_{AB}^2 + C_{AC}^2$. The W and GHZ states
sit at the two extremes.

**WL** : both states at once, so the contrast is carried through every step.

```wl
states = With[{zero = {1, 0}, one = {0, 1}}, <|"W" -> Normalize @ Total[Flatten /@ (KroneckerProduct @@@ Permutations[{zero, zero, one}])], "GHZ" -> Normalize @ Flatten[KroneckerProduct[zero, zero, zero] + KroneckerProduct[one, one, one]]|>]
```

The tangle across the cut needs the single-qubit reduction, obtained by contracting the ket-bra leg pairs of
qubits $2$ and $3$.

```wl
tangle = 4 Det @ ArrayReshape[TensorContract[ArrayReshape[KroneckerProduct[#, Conjugate[#]], ConstantArray[2, 6]], {{2, 5}, {3, 6}}], {2, 2}] & /@ states
```

The pairwise concurrences come from the two-qubit reductions: dropping qubit $3$ leaves the pair $AB$,
dropping qubit $2$ leaves $AC$, and each is fed to Wootters' formula.

```wl
concurrences = Table[With[{r = ArrayReshape[TensorContract[ArrayReshape[KroneckerProduct[#, Conjugate[#]], ConstantArray[2, 6]], {{k, k + 3}}], {4, 4}], y = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]}, With[{l = ReverseSort @ Sqrt @ Chop @ Re @ Eigenvalues[r . y . Conjugate[r] . y]}, Max[0, First[l] - Total[Rest[l]]]]], {k, {3, 2}}] & /@ states
```

The slack in the inequality is the tangle minus what the two pairs account for.

```wl
residual = tangle - Total /@ (concurrences^2)
```

```wl
NonNegative /@ residual
```

**QF** : the tangle across a cut is the square of the concurrence for that bipartition, which the framework
computes when the bipartition is given explicitly.

```wl
qstates = <|"W" -> QuantumState["W"[3]], "GHZ" -> QuantumState["GHZ"[3]]|>
```

```wl
tangleQF = QuantumEntanglementMonotone[#, {{1}, {2, 3}}, "Concurrence"]^2 & /@ qstates
```

The pairwise concurrences are the same monotone on the traced-down two-qubit states.

```wl
concQF = Table[QuantumEntanglementMonotone[QuantumPartialTrace[#, {k}], "Concurrence"], {k, {3, 2}}] & /@ qstates
```

```wl
{tangleQF == tangle, concQF == concurrences}
```

For the W state the tangle is $8/9$ and both pairwise concurrences are $2/3$, so the two squares add to
exactly $8/9$: the inequality is saturated, and W distributes all of qubit $A$'s entanglement into pairwise
sharing. For GHZ the tangle is $1$ while both pairwise concurrences vanish, so the slack is a full $1$: every
pair of GHZ qubits is separable once the third is traced out, and all the entanglement is irreducibly
three-way. That leftover slack is a quantity in its own right, the three-tangle of 11.11.

### 11.11 [MSc] How do I distinguish the GHZ and W classes of tripartite entanglement under SLOCC?

Three qubits have two inequivalent kinds of genuine entanglement, and the residual of the monogamy
inequality separates them: the three-tangle $\tau_3 = \tau_{A|BC} - C_{AB}^2 - C_{AC}^2$ is positive on the
GHZ class and zero on the W class. It equals four times the modulus of Cayley's hyperdeterminant of the
amplitude tensor, and for a $2\times2\times2$ tensor that hyperdeterminant has a short form: writing the
tensor as two square slices $M_0$ and $M_1$, it is the discriminant, in $x$, of the quadratic
$\det(M_0 + x\,M_1)$. Being a degree-four invariant, it can only be rescaled by an invertible local map, so
its vanishing is a property of the whole SLOCC class.

**WL** : both states, as amplitude vectors to be reshaped into $2\times2\times2$ tensors.

```wl
states = With[{zero = {1, 0}, one = {0, 1}}, <|"W" -> Normalize @ Total[Flatten /@ (KroneckerProduct @@@ Permutations[{zero, zero, one}])], "GHZ" -> Normalize @ Flatten[KroneckerProduct[zero, zero, zero] + KroneckerProduct[one, one, one]]|>]
```

The determinant of the slice combination is a quadratic in $x$; its discriminant is the hyperdeterminant.

```wl
tau3 = 4 Abs[With[{c = Coefficient[Det[#[[1]] + \[FormalX] #[[2]]], \[FormalX], {0, 1, 2}]}, c[[2]]^2 - 4 c[[1]] c[[3]]] &@ ArrayReshape[#, {2, 2, 2}]] & /@ states
```

Now apply an invertible but non-unitary local map, one factor per qubit, and renormalize. This moves inside
the SLOCC class while leaving the unitary orbit.

```wl
slocc = Normalize[KroneckerProduct[{{2, 1}, {0, 1}}, {{1, 0}, {1, 3}}, {{1, 2}, {1, 1}}] . #] & /@ states;
```

```wl
4 Abs[With[{c = Coefficient[Det[#[[1]] + \[FormalX] #[[2]]], \[FormalX], {0, 1, 2}]}, c[[2]]^2 - 4 c[[1]] c[[3]]] &@ ArrayReshape[#, {2, 2, 2}]] & /@ slocc
```

**QF** : assemble the same invariant as the monogamy residual, every piece a framework monotone.

```wl
qstates = <|"W" -> QuantumState["W"[3]], "GHZ" -> QuantumState["GHZ"[3]]|>
```

```wl
tau3QF = QuantumEntanglementMonotone[#, {{1}, {2, 3}}, "Concurrence"]^2 - Total[Table[QuantumEntanglementMonotone[QuantumPartialTrace[#, {k}], "Concurrence"], {k, {3, 2}}]^2] & /@ qstates
```

```wl
tau3QF == tau3
```

Both give a three-tangle of $1$ for GHZ and exactly $0$ for W. Under the invertible local map the GHZ value
drops to $36/5041$, rescaled but still nonzero, while the W value stays exactly $0$: no invertible local map
takes one to the other, in either direction, so the two states head genuinely different classes. The
hyperdeterminant and the monogamy residual are two readings of the same invariant, one from the amplitudes
and one from the pairwise entanglement they leave behind.

### 11.12a [MSc] How do I construct a $3\times3$ PPT-entangled (bound) state, entangled yet undistillable because every LOCC operation preserves a positive partial transpose while a Bell state has a negative one?

Beyond the smallest dimensions a positive partial transpose no longer implies separability: there are
entangled states whose partial transpose stays positive, and because distillation is local and cannot turn a
positive partial transpose negative, no amount of processing extracts a single Bell pair from them. They are
*bound* entangled. The standard construction uses an unextendible product basis, a set of orthonormal
product vectors whose orthogonal complement contains no product vector at all. The Tiles basis in $3\times3$
has five members, and the normalized projector onto its complement is bound entangled.

**WL** : the five Tiles vectors, each a Kronecker product of two qutrit combinations built from the qutrit
kets.

```wl
tiles = Normalize /@ (Flatten[KroneckerProduct @@ #] & /@ With[{u = UnitVector[3, # + 1] &}, {{u[0], u[0] - u[1]}, {u[2], u[1] - u[2]}, {u[0] - u[1], u[2]}, {u[1] - u[2], u[0]}, {u[0] + u[1] + u[2], u[0] + u[1] + u[2]}}])
```

They are orthonormal, so their projector has rank five and the complement has rank four.

```wl
Outer[#1 . #2 &, tiles, tiles, 1] == IdentityMatrix[5]
```

Normalizing the projector onto that four-dimensional complement gives the candidate state.

```wl
rho = (IdentityMatrix[9] - Total[KroneckerProduct[#, #] & /@ tiles])/4
```

```wl
{Tr[rho], HermitianMatrixQ[rho], PositiveSemidefiniteMatrixQ[rho]}
```

Its partial transpose on the second qutrit has no negative eigenvalue, so the criterion of 11.6 says
nothing.

```wl
Eigenvalues @ ArrayReshape[Transpose[ArrayReshape[rho, {3, 3, 3, 3}], {1, 4, 3, 2}], {9, 9}]
```

It is nevertheless entangled, and the reason is unextendibility: solve for a normalized product vector
orthogonal to all five tiles.

```wl
Solve[Join[Thread[tiles . Flatten[KroneckerProduct[Array[Subscript[\[FormalA], #] &, 3], Array[Subscript[\[FormalB], #] &, 3]]] == 0], {# . # == 1 & @ Array[Subscript[\[FormalA], #] &, 3], # . # == 1 & @ Array[Subscript[\[FormalB], #] &, 3]}], Join[Array[Subscript[\[FormalA], #] &, 3], Array[Subscript[\[FormalB], #] &, 3]]]
```

An independent numerical certificate is the realignment (computable cross-norm) criterion: reshape the
density matrix by exchanging its two middle legs, then the singular values sum to more than one only for
entangled states.

```wl
N[Total @ SingularValueList @ ArrayReshape[Transpose[ArrayReshape[rho, {3, 3, 3, 3}], {1, 3, 2, 4}], {9, 9}] - 1]
```

**QF** : the same matrix as a state of two qutrits.

```wl
qs = QuantumState[rho, {3, 3}]
```

```wl
QuantumPartialTranspose[qs, {2}]["Eigenvalues"]
```

The negativity, which measures exactly the negative part of that spectrum, is therefore zero:

```wl
QuantumEntanglementMonotone[qs, "Negativity"]
```

yet the entanglement test, which uses the realignment criterion by default, returns `True`:

```wl
{QuantumEntangledQ[qs], N @ QuantumEntanglementMonotone[qs, "Realignment"]}
```

Both routes agree: the partial transpose has four eigenvalues of $1/4$ and five zeros, nothing negative; no
product vector is orthogonal to all five tiles, so the solution set is empty; and the realignment sum exceeds
one by about $0.0874$. So this state is entangled with vanishing negativity, an entanglement that exists but
cannot be concentrated.

### 11.12b [MSc] Why can no bound-entangled state exist for two qubits?

Bound entangled means entangled and undistillable. A positive partial transpose is what makes a state
undistillable, since no local operation with classical communication can turn a positive partial transpose
negative and a Bell state's is negative, so the candidates are exactly the states that survive the partial
transpose. The concurrence of 11.7 vanishes exactly on the separable states. A two-qubit bound-entangled state
would therefore have to have a positive partial transpose *and* a nonzero concurrence, and both of those are
things a machine computes, so the question is a search rather than an argument: generate states, keep the ones
whose partial transpose is still a state, and read off their concurrence. It is zero on every one of them. The
states that pass the filter are not entangled, so there is nothing left to be bound. (One dimension up the
search finds something, which is 11.12a.)

**WL** : generate two hundred random mixed two-qubit states. A Ginibre matrix does it: $g g^\dagger$ normalized
is positive, Hermitian and of unit trace by construction, and generic, so the sample sweeps the interior of the
state space rather than a family.

```wl
SeedRandom[11]; randomStates = Table[#/Tr[#] &@With[{g = RandomComplex[{-1 - I, 1 + I}, {4, 4}]}, g . ConjugateTranspose[g]], 200];
```

Keep the ones whose partial transpose is still a physical state. The transpose acts on the second factor's
indices, so reshaping to $2\times2\times2\times2$ and swapping the pair puts it in place, and being a state is
positive semidefiniteness, the trace being one already.

```wl
ptB[r_] := ArrayReshape[Transpose[ArrayReshape[r, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]; filtered = Select[randomStates, PositiveSemidefiniteMatrixQ[ptB[#]] &];
```

To be bound entangled those states would need a nonzero concurrence, so compute it: the square roots of the
eigenvalues of $\rho$ times its spin flip, largest minus the other three, clipped at zero.

```wl
conc[r_] := With[{y = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]}, Max[0, 2 Max[#] - Total[#]] &@Sqrt@Chop@Re@Eigenvalues[r . y . Conjugate[r] . y]]; Chop[conc /@ filtered]
```

Fifty-one of the two hundred pass the filter and the concurrence is $0$ on every one. The other $149$ all have
positive concurrence, between $0.0025$ and $0.5225$, so the two conditions never come apart in either
direction: positive partial transpose means separable here, and entangled means the partial transpose has gone
negative.

```wl
MinMax@Chop[conc /@ Complement[randomStates, filtered]]
```

Two hundred states is a search, not a proof, and what makes the empty result more than luck is that the two
quantities are tied together by an identity. It is first a statement about the four numbers
$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge \lambda_4$ that the concurrence is built from: expanding
$\prod_k(2\lambda_k - \sum_j \lambda_j)$ and re-collecting it in the power sums kills every odd symmetric
function, which is what lets the square roots disappear.

```wl
With[{l = Array[Subscript[\[FormalL], #] &, 4]}, Expand[Times @@ (2 l - Total[l]) - (Total[l^2]^2 - 2 Total[l^4] + 8 Times @@ l)] === 0]
```

The sign argument is decidable too, so it need not be talked through: over every ordered nonnegative quadruple
whose product of factors is nonnegative, the leading factor is nonpositive, which is the concurrence clipping
to zero.

```wl
With[{l = Array[Subscript[\[FormalL], #] &, 4]}, Resolve[ForAll[l, GreaterEqual @@ Append[l, 0] && Times @@ (2 l - Total[l]) >= 0, First[l] - Total[Rest[l]] <= 0], Reals]]
```

That product is $16\det\rho^{T_B}$. The $\lambda_k^2$ are the eigenvalues of $\rho\tilde\rho$, so the power sums
are $\mathrm{Tr}[\rho\tilde\rho]$ and $\mathrm{Tr}[(\rho\tilde\rho)^2]$, and
$\prod_k \lambda_k = \det\rho$ because the spin flip preserves the determinant. No eigenvalues are left, so what
remains is a polynomial identity in the sixteen real parameters of a Hermitian $4\times4$ matrix, and the Pauli
basis sweeps those exactly once. Put `Conjugate` on the numeric Paulis rather than on the symbolic real
coefficients, where it would stay unevaluated and silently break the check.

```wl
paulis = Flatten[Table[KroneckerProduct[PauliMatrix[i], PauliMatrix[j]], {i, 0, 3}, {j, 0, 3}], 1]; {rho, rhoStar} = With[{c = Array[Subscript[\[FormalC], #] &, 16]}, {c . paulis, c . Conjugate[paulis]}];
```

```wl
rt = With[{y = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]}, rho . y . rhoStar . y]; Expand[16 Det[ptB[rho]] - (Tr[rt]^2 - 2 Tr[rt . rt] + 8 Det[rho])] === 0
```

So the search result holds for every two-qubit state and not just the sample: a positive partial transpose makes
the determinant nonnegative, the identity and the quantifier elimination then force the concurrence to zero, and
a vanishing concurrence is separability (Wootters, the one result borrowed here). The filter can never catch an
entangled state.

**QF** : the framework generates the states and names both quantities, so the same search is three lines.

```wl
SeedRandom[11]; qfStates = Table[QuantumState["RandomMixed"[2]], 200];
```

`"PhysicalQ"` asks whether the partial transpose is still a state, which is the filter.

```wl
qfFiltered = Select[qfStates, QuantumPartialTranspose[#, {2}]["PhysicalQ"] &];
```

```wl
QuantumEntanglementMonotone[#, "Concurrence"] & /@ qfFiltered
```

The same fifty-one states, all with concurrence $0$. Not merely the same count: `"RandomMixed"` is the Ginibre
construction above, so under one seed the two routes generate the identical sample and the deviation is exactly
zero.

```wl
SeedRandom[11]; Max@Abs@Flatten[Normal@QuantumState["RandomMixed"[2]]["DensityMatrix"] - First[randomStates]]
```

The framework's own entanglement test agrees whichever of the two criteria decides it, which is the same
statement in its own terms.

```wl
Tally@Table[QuantumEntangledQ[q, Automatic, #] & /@ {"Negativity", "Concurrence"}, {q, Take[qfStates, 30]}]
```

At $3\times3$ the search stops coming back empty, which is the point of 11.12a. Rebuild the
unextendible-product-basis state in the framework's own terms. Two-qutrit computational kets are named by their
digit string, so each Tiles vector is a difference of two of them and the fifth is the uniform superposition
over all nine.

```wl
qt[s_] := QuantumState[s, {3, 3}]; tilesQF = #["Normalized"] & /@ {qt["00"] - qt["01"], qt["21"] - qt["22"], qt["02"] - qt["12"], qt["10"] - qt["20"], Total[qt /@ StringJoin @@@ Tuples[{"0", "1", "2"}, 2]]};
```

They are orthonormal, which the framework checks with bras applied to kets rather than with dot products.

```wl
Chop@Outer[#1["Dagger"][#2]["Number"] &, tilesQF, tilesQF, 1] === IdentityMatrix[5]
```

Each state carries its own `"Projector"`, and `"I"[{3, 3}]` is the two-qutrit identity, so the complement
projector is operator arithmetic and its trace is $9 - 5 = 4$.

```wl
boundQF = QuantumState[((QuantumOperator["I"[{3, 3}]] - Total[QuantumOperator[#["Projector"]] & /@ tilesQF])/4)["Matrix"], {3, 3}];
```

Now the same two tests, with the realignment criterion standing in for the concurrence, put it exactly where no
two-qubit state goes: partial transpose still a state, entanglement detected anyway.

```wl
{boundQF["PhysicalQ"], QuantumPartialTranspose[boundQF, {2}]["PhysicalQ"], QuantumEntangledQ[boundQF], N@QuantumEntanglementMonotone[boundQF, "Realignment"], QuantumEntanglementMonotone[boundQF, "Negativity"]}
```

The search settles it. Of two hundred random two-qubit states, $51$ have a partial transpose that is still a
state and all $51$ return concurrence exactly $0$; the remaining $149$ all return positive concurrence, from
$0.0025$ to $0.5225$. Pure Wolfram Language and the framework agree state for state, since `"RandomMixed"` is
the normalized Ginibre construction and the deviation between the two samples is $0$, and `QuantumEntangledQ`
returns the same verdict under `"Negativity"` and `"Concurrence"` on all $30$ states tested. The states with a
positive partial transpose are simply not entangled, so the conjunction that defines bound entanglement is never
satisfied. That the sample was not lucky follows from one identity checked symbolically rather than sampled:
`Expand` confirms $16\det\rho^{T_B} = \prod_k(2\lambda_k - \sum_j\lambda_j)$ for all sixteen real parameters of a
Hermitian $4\times4$ matrix at once, and `Resolve` confirms $\lambda_1 \le \lambda_2+\lambda_3+\lambda_4$ for
every ordered nonnegative quadruple keeping that product nonnegative, so a positive partial transpose forces
zero concurrence for every two-qubit state, with Wootters' theorem turning that into separability. One dimension
up the two tests separate: the $3\times3$ unextendible-product-basis state, assembled from `QuantumState` qutrit
kets and their `"Projector"`s, passes `"PhysicalQ"` after the partial transpose and has zero negativity, yet
`QuantumEntangledQ` returns `True` and the realignment monotone $0.0874$. That is a bound-entangled state, and
it is what 11.12a constructs.

### 11.13 [MSc] How do I build the symmetric and antisymmetric subspaces of identical particles and project onto them?

Identical particles occupy definite-symmetry states: bosons live in the symmetric subspace, fermions in the
antisymmetric one. For two particles with $d$ modes each, the exchange operator (the swap) squares to the
identity, so its eigenvalues are $\pm1$ and the projectors onto its two eigenspaces are $(I \pm P)/2$, of
dimensions $d(d+1)/2$ and $d(d-1)/2$. Antisymmetrizing a doubly occupied mode gives zero, which is the Pauli
exclusion principle. Take three modes.

**WL** : the swap is the permutation matrix that exchanges the two labels of every basis pair, generated from
the label list rather than written out.

```wl
swap = With[{d = 3}, PermutationMatrix[FindPermutation[Tuples[Range[d], 2], Reverse /@ Tuples[Range[d], 2]], d^2]]
```

That returns a structured array rather than an explicit list of rows, so display it to see the entries.

```wl
% // MatrixForm
```

The two projectors:

```wl
{sym, anti} = (IdentityMatrix[9] + # swap)/2 & /@ {1, -1};
```

They are idempotent, mutually orthogonal, and complete, which is what makes them a splitting of the space.

```wl
{sym . sym == sym, anti . anti == anti, sym . anti == 0 anti, sym + anti == IdentityMatrix[9]}
```

Their traces are the two subspace dimensions, matching the formulas above.

```wl
With[{d = 3}, {{Tr[sym], Tr[anti]}, {d (d + 1)/2, d (d - 1)/2}}]
```

Project a product state of two *different* modes into each sector.

```wl
{sym . #, anti . #} &@ Flatten @ KroneckerProduct[UnitVector[3, 1], UnitVector[3, 2]]
```

Project a *doubly occupied* mode into the antisymmetric sector.

```wl
anti . Flatten @ KroneckerProduct[UnitVector[3, 1], UnitVector[3, 1]]
```

**QF** : the swap on two qutrits is a named operator, and the projectors are built with operator arithmetic,
so each stays a `QuantumOperator` that knows it acts on two wires of dimension three.

```wl
symQF = (QuantumOperator["I"[{3, 3}]] + QuantumOperator["SWAP"[3]])/2
```

```wl
antiQF = (QuantumOperator["I"[{3, 3}]] - QuantumOperator["SWAP"[3]])/2
```

```wl
{Tr @ Normal @ symQF["Matrix"], Tr @ Normal @ antiQF["Matrix"]}
```

Applying an operator to a state returns a state, so the projections are state objects:

```wl
Normal /@ {symQF[QuantumState[Flatten @ KroneckerProduct[UnitVector[3, 1], UnitVector[3, 2]], {3, 3}]]["StateVector"], antiQF[QuantumState[Flatten @ KroneckerProduct[UnitVector[3, 1], UnitVector[3, 2]], {3, 3}]]["StateVector"]}
```

```wl
Normal[antiQF[QuantumState[Flatten @ KroneckerProduct[UnitVector[3, 1], UnitVector[3, 1]], {3, 3}]]["StateVector"]]
```

```wl
Normal[symQF["Matrix"]] == sym
```

Both give dimensions $6$ and $3$ for the symmetric and antisymmetric sectors of three modes, the two
projections of a distinct-mode product state onto the symmetric and antisymmetric combinations of
$|0,1\rangle$ and $|1,0\rangle$, and exactly zero for the antisymmetrized doubly occupied mode. So exchange
symmetry is a linear constraint, not an extra postulate about the dynamics: the state space of two identical
particles is one of the two eigenspaces of the swap, and the fermionic one simply has no vector with two
particles in the same mode.
