---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 2)

A companion answer key to `Question-List.md`. For each question it gives **two** worked answers:

1. **WL** : native Wolfram Language only (plain vectors and matrices, `PauliMatrix`, `Conjugate`,
   `Normalize`, `Eigenvectors`, `MatrixExp`, and friends). Nothing beyond the built-in language.
2. **QF** : the same task done through the **QuantumFramework** paclet, leaning on its objects and
   property downvalues (`QuantumState[...]["..."]`) rather than rebuilding matrices by hand. A QF
   result is a rich object, a state or an operator, not a bare array: it can act on a target, compose
   with others, report its spectrum, change basis, and more. The answers keep results as objects and
   treat a matrix or amplitude vector as just one of the object's many representations, not the thing
   itself.

Three habits run through the answers. They stay **symbolic** wherever possible (general amplitudes
$\alpha,\beta$, general angles $\theta,\phi$), so each operation is exercised in full generality
rather than on a lucky special case. They write the computation **explicitly** rather than hiding it
in a helper function, because the code is itself part of the explanation. And they prefer **built-in
features** (a `Normalize`, a named state) over hand-rolled constructions, so nothing is hardcoded:
every cell computes its result. Each cell does **one** thing and is preceded by a sentence saying
what, so the notebook reads as a sequence of single, captioned computations. Run the WL answers in a
bare kernel; the QF answers need the paclet loaded once.

## Setup

The QF answers load the paclet once:

```wl
Needs["Wolfram`QuantumFramework`"]
```

---

## Part 2. Observables, spectra, and approximation

### 2.1 [BSc] How do I compute a commutator and an anticommutator of two operators?

The commutator $[A,B]=AB-BA$ measures the failure to commute; the anticommutator $\{A,B\}=AB+BA$.
Together they split the product, $AB = \tfrac12(\{A,B\}+[A,B])$.

**WL** : the commutator is $AB-BA$ on two general operators.

```wl
With[{a = Array[\[FormalA], {2, 2}], b = Array[\[FormalB], {2, 2}]}, a . b - b . a]
```

The anticommutator is $AB+BA$.

```wl
With[{a = Array[\[FormalA], {2, 2}], b = Array[\[FormalB], {2, 2}]}, a . b + b . a]
```

The Wolfram Language also provides them as built-ins (passing `Dot` as the product). The commutator:

```wl
With[{a = Array[\[FormalA], {2, 2}], b = Array[\[FormalB], {2, 2}]}, Commutator[a, b, Dot]]
```

and the anticommutator:

```wl
With[{a = Array[\[FormalA], {2, 2}], b = Array[\[FormalB], {2, 2}]}, Anticommutator[a, b, Dot]]
```

**QF** : here the operands are *operators*, so the results are again operators, not arrays. The
built-in `Commutator` reduces directly on QuantumOperators (QF teaches it the operator product), and
returns a `QuantumOperator`.

```wl
With[{a = QuantumOperator[Array[\[FormalA], {2, 2}]], b = QuantumOperator[Array[\[FormalB], {2, 2}]]},
 Commutator[a, b]]
```

The built-in `Anticommutator` does *not* reduce on operators, so write it as `a @ b + b @ a`, with
`@` the operator product.

```wl
With[{a = QuantumOperator[Array[\[FormalA], {2, 2}]], b = QuantumOperator[Array[\[FormalB], {2, 2}]]},
 a @ b + b @ a]
```

An operator does far more than carry a matrix. Take the concrete Pauli case $[X,Y]=2iZ$: the
commutator is itself an operator object, shown here as its summary box.

```wl
c = Commutator[QuantumOperator["X"], QuantumOperator["Y"]]
```

It acts on a state, sending $|0\rangle$ to $2i|0\rangle$.

```wl
c[QuantumState["0"]]
```

It answers whether it is Hermitian (a commutator of Hermitian operators is anti-Hermitian, so `False`).

```wl
c["HermitianQ"]
```

Its eigenvalues are imaginary, $\pm 2i$.

```wl
c["Eigenvalues"]
```

And only last of all does it hand back its matrix, $2i\sigma_z$.

```wl
c["Matrix"] // Normal
```

The matrix is one view of a much richer object; that is the difference between QuantumFramework and
bare linear algebra.

### 2.2 [BSc] How do I build the Pauli operators and verify their algebra (Levi-Civita product rule)?

The Pauli operators close under commutation, $[\sigma_i,\sigma_j] = 2i\,\epsilon_{ijk}\,\sigma_k$.

**WL** : assemble the Pauli vector and check all nine commutators at once against the Levi-Civita
tensor.

```wl
With[{\[Sigma] = PauliMatrix[{1, 2, 3}]},
 Outer[#1 . #2 - #2 . #1 &, \[Sigma], \[Sigma], 1] ==
  2 I TensorContract[TensorProduct[LeviCivitaTensor[3], \[Sigma]], {3, 4}]]
```

**QF** : the same check staying at the *operator* level. Each `Commutator[#1, #2]` is an operator, and
every entry is compared directly to the operator $2i\,\epsilon_{ijk}\sigma_k$, with no matrices
extracted.

```wl
With[{\[Sigma] = QuantumOperator /@ {"X", "Y", "Z"}}, {com =
   Outer[Commutator[#1, #2] &, \[Sigma], \[Sigma], 1]},
 Table[com[[i, j]] == 2 I Sum[LeviCivitaTensor[3][[i, j, k]] \[Sigma][[k]], {k, 3}], {i, 3}, {j, 3}]]
```

The WL line returns `True` and the QF line returns a $3\times3$ table of `True`: every commutator
operator $[\sigma_i,\sigma_j]$ equals the operator $2i\,\epsilon_{ijk}\sigma_k$ (the diagonal
$[\sigma_i,\sigma_i]$ vanishing), so the whole $\mathfrak{su}(2)$ algebra is verified through operator
equality, never touching a matrix.

### 2.3a [BSc] How do I test that an operator is Hermitian?

An operator is Hermitian when it equals its own conjugate transpose, $A = A^\dagger$.

**WL** : test Hermiticity directly.

```wl
HermitianMatrixQ[PauliMatrix[1]]
```

**QF** : the operator answers `"HermitianQ"`.

```wl
QuantumOperator["X"]["HermitianQ"]
```

Both return `True`: $\sigma_x$ is Hermitian.

### 2.3b [BSc] Why must an observable be Hermitian?

Measurement outcomes are real numbers, and the possible outcomes of an observable are its eigenvalues.
Hermiticity is exactly the condition that forces those eigenvalues to be real, which is why every
observable is represented by a Hermitian operator.

**WL** : a generic Hermitian matrix has real eigenvalues.

```wl
FullSimplify[Eigenvalues[{{a1, a2 + I a3}, {a2 - I a3, a4}}] \[Element] Reals,
 (a1 | a2 | a3 | a4) \[Element] Reals]
```

**QF** : the same, read off the operator's `"Eigenvalues"`.

```wl
FullSimplify[QuantumOperator[{{a1, a2 + I a3}, {a2 - I a3, a4}}]["Eigenvalues"] \[Element] Reals,
 (a1 | a2 | a3 | a4) \[Element] Reals]
```

Both return `True`: every Hermitian operator has a real spectrum, so its measurement outcomes are real.
That real-spectrum guarantee is the physical reason observables must be Hermitian.

### 2.4 [BSc] How do I compute eigenvalues and eigenvectors and form the spectral decomposition?

Every observable is its own spectral sum, $A = \sum_a a\,|a\rangle\langle a|$. The built-in
`EigenvalueDecomposition` packages exactly this: it returns $\{S, D\}$ with the eigenvectors as the
columns of $S$, the eigenvalues on the diagonal of $D$, and $A = S\,D\,S^{-1}$.

**WL** : compute the decomposition once and bind it.

```wl
evd = EigenvalueDecomposition[{{a1, a2 + I a3}, {a2 - I a3, a4}}];
```

The eigenvalues are the diagonal of $D$.

```wl
Diagonal[evd[[2]]]
```

The eigenvectors are the columns of $S$ (rows of its transpose).

```wl
Transpose[evd[[1]]]
```

And the reconstruction $A = S\,D\,S^{-1}$ holds.

```wl
FullSimplify[{{a1, a2 + I a3}, {a2 - I a3, a4}} == evd[[1]] . evd[[2]] . Inverse[evd[[1]]],
 (a1 | a2 | a3 | a4) \[Element] Reals]
```

**QF** : keep it in the object model. The operator, shown as its summary box.

```wl
qo = QuantumOperator[{{a1, a2 + I a3}, {a2 - I a3, a4}}]
```

It hands back its eigenvalues,

```wl
qo["Eigenvalues"]
```

its eigenvectors,

```wl
qo["Eigenvectors"]
```

and `"Diagonalize"` returns the *operator* expressed in its own eigenbasis, a `QuantumOperator` whose
matrix is diagonal with the eigenvalues on it.

```wl
qo["Diagonalize"]
```

Both report the same two real eigenvalues $\tfrac12\bigl(a_1+a_4 \pm \sqrt{(a_1-a_4)^2+4(a_2^2+a_3^2)}\bigr)$.
On the WL side `EigenvalueDecomposition` reconstructs $A = S D S^{-1}$ (for a Hermitian operator the
eigenvectors are orthogonal, so normalizing $S$ to a unitary $U$ turns this into the spectral form
$A = U D U^\dagger = \sum_a a\,|a\rangle\langle a|$). On the QF side `"Diagonalize"` *is* that spectral
form, the operator already diagonal in its eigenbasis, still an object rather than a bare array.

### 2.5 [BSc] How do I evaluate a function of an operator $f(\hat A)$ by functional calculus?

Functional calculus defines $f(\hat A)$ by letting $f$ act on the eigenvalues in the eigenbasis,
$f(\hat A) = \sum_a f(a)\,|a\rangle\langle a|$. In WL the one-call tool that does exactly this is
`MatrixFunction`.

**WL** : apply a symbolic $f$ to $\sigma_x$ directly.

```wl
MatrixFunction[f, PauliMatrix[1]]
```

The result $f(\sigma_x) = \tfrac12\begin{pmatrix} f(1)+f(-1) & f(1)-f(-1) \\ f(1)-f(-1) & f(1)+f(-1)\end{pmatrix}$
is *not* $f$ taken entry by entry: the naive entrywise map would send the diagonal $0$s of $\sigma_x$
to $f(0)$, whereas functional calculus gives $\tfrac12\bigl(f(1)+f(-1)\bigr)$ there. The reason is the
spectral definition, $f$ acts on the eigenvalues $\pm1$ and the eigenprojectors of $\sigma_x$ spread
that action across the whole matrix.

To make that explicit, build the spectral sum $\sum_a f(a)\,|a\rangle\langle a|$ from the eigensystem.

```wl
spec = With[{sys = Eigensystem[PauliMatrix[1]]},
  Total[MapThread[f[#1] KroneckerProduct[#2, Conjugate[#2]] &, {sys[[1]], Normalize /@ sys[[2]]}]]]
```

It reproduces the built-in to the last entry.

```wl
spec == MatrixFunction[f, PauliMatrix[1]]
```

**QF** : the same spectral definition driven by the operator's own `"Eigensystem"`.

```wl
With[{sys = QuantumOperator["X"]["Eigensystem"]},
 Total[MapThread[f[#1] KroneckerProduct[#2, Conjugate[#2]] &, {sys[[1]], Normalize /@ sys[[2]]}]]]
```

So `MatrixFunction` is the direct answer and the spectral sum is what it evaluates: $f$ acts through
the eigenbasis, and the two agree exactly.

### 2.6 [BSc] How do I compute an expectation value and a variance on a pure state as $\langle\psi|\hat A|\psi\rangle$?

The mean is $\langle A\rangle = \langle\psi|A|\psi\rangle$ and the variance is
$\langle A^2\rangle - \langle A\rangle^2$.

**WL** : contract the bra, operator, and ket for the mean.

```wl
With[{a = PauliMatrix[3], psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 FullSimplify[Conjugate[psi] . a . psi, {\[Theta], \[Phi]} \[Element] Reals]]
```

The variance reuses $A^2$.

```wl
With[{a = PauliMatrix[3], psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 FullSimplify[Conjugate[psi] . a . a . psi - (Conjugate[psi] . a . psi)^2, {\[Theta], \[Phi]} \[Element] Reals]]
```

**QF** : the mean as a bra-operator-ket scalar.

```wl
With[{a = QuantumOperator["Z"], ps = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]},
 FullSimplify[(ps["Dagger"] @ a @ ps)["Scalar"], {\[Theta], \[Phi]} \[Element] Reals]]
```

The variance, again as scalars.

```wl
With[{a = QuantumOperator["Z"], ps = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]},
 FullSimplify[(ps["Dagger"] @ a @ a @ ps)["Scalar"] - ((ps["Dagger"] @ a @ ps)["Scalar"])^2,
  {\[Theta], \[Phi]} \[Element] Reals]]
```

Both give $\langle\sigma_z\rangle = \cos\theta$ and $\mathrm{Var}=\sin^2\theta$: the spin points along
$z$ with certainty at the poles and is maximally uncertain on the equator. The mixed-state generalization,
$\langle A\rangle = \mathrm{Tr}[A\rho]$ on any $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ (giving
$\langle\sigma_z\rangle = r_z$, $\mathrm{Var} = 1 - r_z^2$), is Part 3, 3.4.

### 2.7 [BSc] How do I write the resolution of the identity in an orthonormal basis?

Any orthonormal basis is complete, $\sum_i |i\rangle\langle i| = I$.

**WL** : sum the eigenprojectors of a generic Hermitian operator.

```wl
With[{vecs = Normalize /@ Eigenvectors[{{a1, a2 + I a3}, {a2 - I a3, a4}}]},
 FullSimplify[Total[KroneckerProduct[#, Conjugate[#]] & /@ vecs] == IdentityMatrix[2],
  (a1 | a2 | a3 | a4) \[Element] Reals]]
```

**QF** : the same with the operator's `"Eigenvectors"`.

```wl
With[{vecs = Normalize /@ QuantumOperator[{{a1, a2 + I a3}, {a2 - I a3, a4}}]["Eigenvectors"]},
 FullSimplify[Total[KroneckerProduct[#, Conjugate[#]] & /@ vecs] == IdentityMatrix[2],
  (a1 | a2 | a3 | a4) \[Element] Reals]]
```

Both return `True`: the eigenprojectors of any observable sum to the identity.

### 2.8 [BSc] How do I decide whether two observables are compatible (commute) and find a simultaneous eigenbasis?

Two observables are *compatible*, meaning they can hold definite values at the same time, exactly when
they commute, $[A,B]=0$. The reason commuting gives a shared eigenbasis is one line: if $A|a\rangle =
a|a\rangle$ and $AB = BA$, then $A\,(B|a\rangle) = B\,(A|a\rangle) = a\,(B|a\rangle)$, so $B$ keeps each
eigenspace of $A$ inside itself. When the eigenvalue $a$ is **non-degenerate** (a one-dimensional
eigenspace), $B|a\rangle$ can only be a multiple of $|a\rangle$, so $|a\rangle$ is automatically an
eigenvector of $B$ too. That is the whole job in the ordinary case: diagonalize one non-degenerate
observable, and its eigenvectors are already the simultaneous eigenbasis.

**Normal case.** Take a non-degenerate $A = X\otimes I + 2\,I\otimes Z$ and a compatible partner
$B = X\otimes Z$.

**WL** : compatibility is the vanishing commutator.

```wl
With[{a = KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]] + 2 KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]],
   b = KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]},
 a . b == b . a]
```

What makes this the easy case is that $A$ is non-degenerate: its four eigenvalues are distinct, so every
eigenspace is one-dimensional and the eigenbasis is unique.

```wl
Eigenvalues[
 KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]] + 2 KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]]
```

Diagonalize $A$ once: its eigenvectors also diagonalize $B$, so they *are* the simultaneous eigenbasis.
The check conjugates each operator into the candidate basis ($u$ holding the eigenvectors as rows) and
tests that the result is diagonal; the pair `{True, True}` says both are.

```wl
With[{a = KroneckerProduct[PauliMatrix[1], IdentityMatrix[2]] + 2 KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]],
   b = KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]},
 With[{u = Normalize /@ Eigenvectors[a]},
  DiagonalMatrixQ /@ (Chop[Conjugate[u] . # . Transpose[u]] & /@ {a, b})]]
```

**Special case: when neither observable is enough by itself.** The shortcut breaks when the observables
are degenerate, because then "diagonalize $A$" no longer names a unique basis. Take $A = X\otimes X$ and
$B = Z\otimes Z$: they commute, but each has $+1$ and $-1$ *doubled*.

```wl
{Eigenvalues[KroneckerProduct[PauliMatrix[1], PauliMatrix[1]]],
 Eigenvalues[KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]]}
```

A doubled eigenvalue is a two-dimensional eigenspace, and *any* orthonormal pair inside it is an equally
valid eigenbasis, so a choice that diagonalizes one operator need not diagonalize the other. Concretely,
the obvious eigenbasis of $X\otimes X$, each qubit in an $X$-eigenstate $|{\pm}\rangle$ (the rows of
$H\otimes H$), diagonalizes $X\otimes X$ but not $Z\otimes Z$: the pair comes back `{True, False}`.

```wl
With[{u = KroneckerProduct[HadamardMatrix[2], HadamardMatrix[2]]},
 DiagonalMatrixQ /@ (Chop[Conjugate[u] . # . Transpose[u]] & /@
   {KroneckerProduct[PauliMatrix[1], PauliMatrix[1]], KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]})]
```

The fix is to diagonalize a single *non-degenerate* combination $A + cB$ for a generic constant $c$. Its
eigenvalues are $a + c\,b$ over the joint values $(a,b)$, which a generic $c$ keeps distinct, so its
eigenbasis is unique; and because $A + cB$ commutes with both $A$ and $B$, that one basis diagonalizes
both. Take $c = \pi$, and the pair is `{True, True}`.

```wl
With[{a = KroneckerProduct[PauliMatrix[1], PauliMatrix[1]], b = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]},
 With[{u = Normalize /@ Eigenvectors[a + \[Pi] b]},
  DiagonalMatrixQ /@ (Chop[Conjugate[u] . # . Transpose[u]] & /@ {a, b})]]
```

**QF** : compatibility is the *operator* equality $AB = BA$, written `a @ b == b @ a`.

```wl
With[{a = QuantumOperator["XX"], b = QuantumOperator["ZZ"]},
 a @ b == b @ a]
```

The simultaneous eigenbasis is the eigenstates of `a + Pi b`. Each one holds a definite value of both $A$
and $B$, so each operator returns it to itself up to a phase, leaving its density matrix unchanged.

```wl
With[{a = QuantumOperator["XX"], b = QuantumOperator["ZZ"]},
 With[{basis = QuantumState /@ (a + \[Pi] b)["Eigenvectors"]},
  AllTrue[basis, (a[#])["DensityMatrix"] == #["DensityMatrix"] && (b[#])["DensityMatrix"] == #["DensityMatrix"] &]]]
```

So: compatible means commuting; in the ordinary case one non-degenerate observable already fixes the
shared eigenbasis, and when both observables are degenerate you restore uniqueness by diagonalizing a
generic combination $A + cB$ just once. The basis that $A + \pi B$ returns is the Bell basis; $\pi$ only
serves to break the degeneracy, and any value that does so yields the same simultaneous eigenbasis.

### 2.9a [MSc] How do I construct a complete set of commuting observables (CSCO)?

A complete set of commuting observables (CSCO) is a set of mutually commuting Hermitian operators whose
simultaneous eigenvalues label every basis state *uniquely*. Commuting is the easy half (handled in 2.8);
*complete* is the real content. A single observable is a CSCO only when it is non-degenerate: a degenerate
eigenvalue is a multi-dimensional eigenspace in which no basis is singled out, so the observable cannot
distinguish the states it contains. The construction is to adjoin further commuting observables, each
lifting some remaining degeneracy, until every simultaneous eigenspace is one-dimensional. Completeness is
then a computable test: once the observables commute, the number of *distinct* joint-eigenvalue tuples
must equal the Hilbert-space dimension, since only then does the label pin down a unique state.

**WL** : one observable is not enough. $Z\otimes I$ reads only the first qubit, so it is blind to the
second; its spectrum is doubly degenerate (two distinct eigenvalues for a four-dimensional space), and its
eigenbasis is not unique.

```wl
With[{a = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]]},
 {Eigenvalues[a], Length@DeleteDuplicates@Eigenvalues[a] < Length[a]}]
```

Adjoin the commuting observable $I\otimes Z$, which reads the second qubit. The joint-eigenvalue tuples are
now all distinct (four of them, equal to the dimension), so every state is uniquely labeled and the pair is
a CSCO.

```wl
With[{a = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], b = KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]},
 With[{joint = Transpose[{Diagonal[a], Diagonal[b]}]},
  {a . b == b . a, Length@DeleteDuplicates@joint == Length@joint}]]
```

**QF** : the same construction, with the operators named as strings. $Z\otimes I$ alone is degenerate,

```wl
With[{a = QuantumOperator["ZI"]},
 {a["Eigenvalues"], Length@DeleteDuplicates@a["Eigenvalues"] < Length@a["Eigenvalues"]}]
```

and adjoining $I\otimes Z$ completes the set.

```wl
With[{a = QuantumOperator["ZI"], b = QuantumOperator["IZ"]},
 With[{joint = Transpose[{Diagonal[Normal[a["Matrix"]]], Diagonal[Normal[b["Matrix"]]]}]},
  {a @ b == b @ a, Length@DeleteDuplicates@joint == Length@joint}]]
```

The single observable returns two distinct eigenvalues for a four-dimensional space (`True`: degenerate,
hence incomplete); the pair returns `{True, True}` (commuting, and its four joint tuples are distinct, so
complete). Here the joint eigenvalues are read straight off the diagonals because $Z\otimes I$ and
$I\otimes Z$ are already simultaneously diagonal; in general one first diagonalizes a non-degenerate
combination $A+cB$ (as in 2.8) and reads the tuples in that common eigenbasis. The unique labels this CSCO
assigns are the subject of 2.9b.

### 2.9b [MSc] How do I label basis states by the joint spectrum of a CSCO?

With $\{Z\otimes I, I\otimes Z\}$ established as a CSCO (2.9a), its joint eigenvalues are the *quantum
numbers* that name the basis: each state carries a tuple $(z_1, z_2)$, and because the set is complete
those tuples are all distinct. Reading them off (the operators are simultaneously diagonal) pairs every
computational basis state with its label.

**WL** : pair each basis state $|z_1 z_2\rangle$ with its joint eigenvalues.

```wl
With[{a = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], b = KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]},
 AssociationThread[Tuples[{0, 1}, 2], Transpose[{Diagonal[a], Diagonal[b]}]]]
```

**QF** : the same labels, with the operators named as strings.

```wl
With[{a = QuantumOperator["ZI"], b = QuantumOperator["IZ"]},
 AssociationThread[Tuples[{0, 1}, 2], Transpose[{Diagonal[Normal[a["Matrix"]]], Diagonal[Normal[b["Matrix"]]]}]]]
```

Both pair $|00\rangle\!\to\!(+,+)$, $|01\rangle\!\to\!(+,-)$, $|10\rangle\!\to\!(-,+)$,
$|11\rangle\!\to\!(-,-)$: the four eigenvalue tuples are distinct, so the joint spectrum is a complete,
unambiguous label for the basis. These tuples are the finite-dimensional analogue of the good quantum
numbers (like $n,\ell,m$) that index states in continuous problems.

### 2.10a [MSc] How do I define the Hilbert-Schmidt operator inner product $\langle A,B\rangle=\mathrm{Tr}[A^\dagger B]$ and show the Pauli operators are orthogonal?

Operators form a Hilbert space under their own inner product, the Hilbert-Schmidt product
$\langle A,B\rangle=\mathrm{Tr}[A^\dagger B]$. Under it the four operators $\{I,X,Y,Z\}$ are mutually
orthogonal, each of squared norm $\langle\sigma_j,\sigma_j\rangle = 2$.

**WL** : the Gram matrix $\langle\sigma_i,\sigma_j\rangle$ over $\{I,X,Y,Z\}$.

```wl
Table[Tr[ConjugateTranspose[PauliMatrix[i]] . PauliMatrix[j]], {i, 0, 3}, {j, 0, 3}]
```

**QF** : the same inner product on `QuantumOperator`s, $\mathrm{Tr}[A^\dagger B]$ being the trace of the
composed operator $A^\dagger B$.

```wl
With[{paulis = QuantumOperator /@ {"I", "X", "Y", "Z"}},
 Outer[Tr[Normal[(#1["Dagger"] @ #2)["Matrix"]]] &, paulis, paulis, 1]]
```

Both return $2\,\mathbb{1}_4$: the operators are mutually orthogonal ($\langle\sigma_i,\sigma_j\rangle = 0$
for $i\ne j$) and each has squared norm $2$, so $\{I,X,Y,Z\}/\sqrt2$ is an orthonormal basis of the
$2\times2$ operator space.

### 2.10b [MSc] How do I expand an operator in the Pauli basis?

Because $\{I,X,Y,Z\}$ is orthogonal under the Hilbert-Schmidt product, any $2\times2$ operator $A$ has a
unique Pauli expansion $A = \sum_j c_j\,\sigma_j$ with coefficients $c_j = \tfrac12\mathrm{Tr}[\sigma_j A]$
(the $\tfrac12$ is $1/\langle\sigma_j,\sigma_j\rangle$).

**WL** : compute the four coefficients directly from the definition.

```wl
With[{aa = {{a1, a2 + I a3}, {a2 - I a3, a4}}},
 AssociationThread[{"I", "X", "Y", "Z"}, FullSimplify @ Table[Tr[ConjugateTranspose[PauliMatrix[j]] . aa]/2, {j, 0, 3}]]]
```

**QF** : the operator decomposes itself.

```wl
FullSimplify /@ QuantumOperator[{{a1, a2 + I a3}, {a2 - I a3, a4}}]["PauliDecompose"]
```

Both return the coefficients $\{I:\tfrac{a_1+a_4}2,\ X:a_2,\ Y:-a_3,\ Z:\tfrac{a_1-a_4}2\}$: the
generic Hermitian operator is exactly $\tfrac{a_1+a_4}2 I + a_2\sigma_x - a_3\sigma_y + \tfrac{a_1-a_4}2\sigma_z$.

### 2.11 [BSc] How do I bound the ground-state energy from above using a trial state (the Rayleigh-Ritz variational method)?

The variational method estimates the ground-state energy $E_0$ from above using a guessed *trial state*
$|\psi\rangle$, without solving the eigenproblem. The key fact is that the **Rayleigh quotient** is an
upper bound on $E_0$,
$$
E[\psi] \;=\; \frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle} \;\ge\; E_0 ,
$$
with equality only when $|\psi\rangle$ is the ground state. (Expanding $|\psi\rangle = \sum_n c_n|n\rangle$
in the energy eigenbasis, $H|n\rangle = E_n|n\rangle$ with $E_0\le E_1\le\cdots$, makes $E[\psi]$ a weighted
average of the $E_n$, which cannot fall below the smallest.) *Rayleigh-Ritz* is the recipe that turns the
bound into an estimate: take a parametrized family of trial states and **minimize** $E[\psi]$ over the
parameters; the minimum is the lowest, hence best, upper bound that family can give. Below the trial state
is normalized, so $\langle\psi|\psi\rangle = 1$ and $E[\psi] = \langle\psi|H|\psi\rangle$.

**WL** : the trial energy over a one-parameter family.

```wl
exp = With[{h = {{1, 1}, {1, -1}}, psi = {Cos[t], Sin[t]}},
  FullSimplify[Conjugate[psi] . h . psi, t \[Element] Reals]]
```

The true ground energy is the lowest eigenvalue.

```wl
e0 = Min[Eigenvalues[{{1, 1}, {1, -1}}]]
```

Minimizing the trial energy reaches that bound.

```wl
Minimize[exp, t][[1]] == e0
```

**QF** : the Hamiltonian as a `QuantumOperator`, the trial energy as a bra-operator-ket scalar.

```wl
expQF = With[{h = QuantumOperator[{{1, 1}, {1, -1}}], ps = QuantumState[{Cos[t], Sin[t]}]},
  FullSimplify[(ps["Dagger"] @ h @ ps)["Scalar"], t \[Element] Reals]]
```

The ground energy from the operator's eigenvalues.

```wl
QuantumOperator[{{1, 1}, {1, -1}}]["Eigenvalues"] // Min
```

And the variational minimum reaches it.

```wl
Minimize[expQF, t][[1]] == Min[QuantumOperator[{{1, 1}, {1, -1}}]["Eigenvalues"]]
```

Both give the trial energy $\cos 2t + \sin 2t$, the true ground energy $-\sqrt2$, and `True`: the
minimum over the trial family reaches $E_0$, so the variational estimate is an upper bound that is
tight here.

### 2.12 [MSc] How do I compute first- and second-order energy shifts by non-degenerate time-independent perturbation theory?

Perturbation theory finds the spectrum of a Hamiltonian $H = H_0 + \lambda V$ that sits a small
perturbation $\lambda V$ away from a solvable one $H_0$, whose eigenstates $|n\rangle$ and energies $E_n$
are already known. Turning on $\lambda$ shifts each energy, and the shift is expanded in powers of
$\lambda$,
$$
E_n(\lambda) \;=\; E_n + \lambda\,E_n^{(1)} + \lambda^2 E_n^{(2)} + \cdots ,
$$
with the first- and second-order coefficients (the Rayleigh-Schrodinger shifts)
$$
E_n^{(1)} = \langle n|V|n\rangle , \qquad
E_n^{(2)} = \sum_{m\ne n} \frac{|\langle m|V|n\rangle|^2}{E_n - E_m} .
$$
The first order is just the diagonal matrix element of $V$ in the unperturbed state $|n\rangle$; the second
order sums the squared off-diagonal couplings $|\langle m|V|n\rangle|^2$ over the other states, each divided
by the energy gap $E_n - E_m$. *Non-degenerate* means the unperturbed energies $E_n$ are distinct, so those
gaps never vanish; degenerate levels require the separate treatment of 2.13.

**WL** : the first-order shifts are the diagonal matrix elements.

```wl
With[{V = {{w, v}, {Conjugate[v], -w}}, basis = IdentityMatrix[2]},
 Table[Conjugate[basis[[n]]] . V . basis[[n]], {n, 2}]]
```

The second-order shifts sum the off-diagonal couplings over the spectrum.

```wl
With[{en = {e1, e2}, V = {{w, v}, {Conjugate[v], -w}}, basis = IdentityMatrix[2]},
 Table[Sum[If[m == n, 0, Abs[Conjugate[basis[[m]]] . V . basis[[n]]]^2/(en[[n]] - en[[m]])], {m, 2}], {n, 2}]]
```

**QF** : the first-order shifts as bra-operator-ket scalars on the basis states.

```wl
With[{V = QuantumOperator[{{w, v}, {Conjugate[v], -w}}], basis = {QuantumState["0"], QuantumState["1"]}},
 Table[(basis[[n]]["Dagger"] @ V @ basis[[n]])["Scalar"], {n, 2}]]
```

The second-order shifts, same matrix elements summed.

```wl
With[{en = {e1, e2}, V = QuantumOperator[{{w, v}, {Conjugate[v], -w}}], basis = {QuantumState["0"], QuantumState["1"]}},
 Table[Sum[If[m == n, 0, Abs[(basis[[m]]["Dagger"] @ V @ basis[[n]])["Scalar"]]^2/(en[[n]] - en[[m]])], {m, 2}], {n, 2}]]
```

Both give first-order shifts $\{w,-w\}$ (the diagonal of $V$) and second-order shifts
$\left\{\dfrac{|v|^2}{e_1-e_2},\ \dfrac{|v|^2}{e_2-e_1}\right\}$: the off-diagonal coupling pushes the
levels apart, the lower one down and the upper one up.

### 2.13 [MSc] How do I resolve an avoided crossing by degenerate perturbation theory (level repulsion)?

Where two diagonal levels would cross, an off-diagonal coupling $\Delta$ opens a gap: the eigenvalues
of $\begin{pmatrix}\lambda & \Delta \\ \Delta & -\lambda\end{pmatrix}$ are
$\pm\sqrt{\Delta^2+\lambda^2}$. At the would-be crossing $\lambda = 0$ the two levels are degenerate, and
diagonalizing the coupling $\Delta$ within that degenerate subspace is exactly first-order degenerate
perturbation theory; in this $2\times2$ case the degenerate-subspace diagonalization and the exact
diagonalization coincide, giving the split $\pm\Delta$ at the crossing.

**WL** : the two eigenvalues of the two-level Hamiltonian.

```wl
With[{h = {{\[Lambda], \[CapitalDelta]}, {\[CapitalDelta], -\[Lambda]}}},
 FullSimplify[Eigenvalues[h], (\[Lambda] | \[CapitalDelta]) \[Element] Reals]]
```

Their difference is the gap.

```wl
With[{h = {{\[Lambda], \[CapitalDelta]}, {\[CapitalDelta], -\[Lambda]}}},
 FullSimplify[Abs[Subtract @@ Eigenvalues[h]], (\[Lambda] | \[CapitalDelta]) \[Element] Reals]]
```

**QF** : the eigenvalues through the operator.

```wl
With[{h = QuantumOperator[{{\[Lambda], \[CapitalDelta]}, {\[CapitalDelta], -\[Lambda]}}]},
 FullSimplify[h["Eigenvalues"], (\[Lambda] | \[CapitalDelta]) \[Element] Reals]]
```

And the gap.

```wl
With[{h = QuantumOperator[{{\[Lambda], \[CapitalDelta]}, {\[CapitalDelta], -\[Lambda]}}]},
 FullSimplify[Abs[Subtract @@ h["Eigenvalues"]], (\[Lambda] | \[CapitalDelta]) \[Element] Reals]]
```

Both give eigenvalues $\pm\sqrt{\Delta^2+\lambda^2}$ and a gap $2\sqrt{\Delta^2+\lambda^2}$. As
$\lambda$ sweeps through the would-be crossing at $\lambda=0$, the gap shrinks only to $2\Delta$, never
to zero: the coupling turns the crossing into an avoided crossing.
