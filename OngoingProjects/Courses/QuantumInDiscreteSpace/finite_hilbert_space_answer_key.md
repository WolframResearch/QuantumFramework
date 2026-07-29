# Answer Key: Quantum Questions in Finite-Dimensional Hilbert Space

Pedagogical answers to the 185 questions in
[`finite_hilbert_space_question_bank.md`](finite_hilbert_space_question_bank.md). Each answer
opens with the full question exactly as posed in the bank (in italics), then states the central
idea in physical terms, gives the closed-form result, explains the moving parts, and, where a
computation is the natural answer, demonstrates it with Wolfram Language or QuantumFramework
code.

## How to read this key

The code follows a computation-first rhythm: a sentence names what to compute, one cell computes
exactly one thing, and the following sentence says what the result means. Verifications are
written as an equality that the kernel reduces to `True`, so you watch the claim confirm itself.
Every cell has been executed in a Wolfram kernel (QuantumFramework 2.0.0). Within a section the
cells share one running session in reading order, so a symbol defined earlier (`plusx`, `psi`,
`U`) is still in scope later. Snippets follow the *minimal faithful idiom*: the lightest
construction that still shows the physics. Where a question is purely conceptual, the answer is
prose and no cell is forced onto it.

Notation: $\mathcal{H}$ is a finite-dimensional complex Hilbert space of dimension $n$;
$|\psi\rangle,|\phi\rangle$ are kets; $\langle\psi|$ is the dual bra; $\{|i\rangle\}$ is an
orthonormal basis; $A,B$ are operators; $\hat{n}$ is a unit vector; $\vec{\sigma}=(\sigma_x,\sigma_y,\sigma_z)$
are the Pauli matrices.

---

## A. Finite Hilbert Spaces, Bases, and Dimension

### 1. What makes a Hilbert space more than a vector space?

> *In an $n$-dimensional complex Hilbert space, what algebraic properties distinguish a vector space from a Hilbert space?*

A complex vector space gives you two operations: you can add two states and scale a state by a
complex number. That is enough to write superpositions, but not enough to speak of
probabilities, lengths, or angles. A **Hilbert space** is a complex vector space equipped with
an **inner product** $\langle\phi|\psi\rangle$, a map $\mathcal{H}\times\mathcal{H}\to\mathbb{C}$
that is linear in the ket, conjugate-linear in the bra, conjugate-symmetric
$\langle\phi|\psi\rangle=\overline{\langle\psi|\phi\rangle}$, and positive-definite
$\langle\psi|\psi\rangle>0$ for $|\psi\rangle\neq 0$. In other words, the inner product is the
one extra piece of structure that turns bare vectors into objects you can measure.

The inner product is what makes the space *physical*. It defines the norm
$\||\psi\rangle\|=\sqrt{\langle\psi|\psi\rangle}$ used to normalize states, it converts a ket
into the amplitude $\langle i|\psi\rangle$ whose modulus squared is a Born probability, and it
supplies the notion of orthogonality that lets distinct measurement outcomes be perfectly
distinguishable.

In finite dimension there is one further subtlety that quietly disappears. A general Hilbert
space is required to be **complete**: every Cauchy sequence of vectors converges to a vector in
the space. Completeness is a genuine constraint in infinite dimension (it is what separates the
square-summable sequences $\ell^2$ from a merely dense subspace), but every finite-dimensional
inner-product space is automatically complete. So in finite dimension "Hilbert space" reduces to
"complex vector space with an inner product," i.e. $\mathbb{C}^n$ with a fixed Hermitian
inner product. Completeness is the property you never have to check here and must always check
once the dimension is infinite.

### 2. When is a set of $n$ vectors an orthonormal basis?

> *What conditions must a set of $n$ vectors satisfy to be an orthonormal basis of an $n$-dimensional Hilbert space?*

For $n$ vectors $\{|e_i\rangle\}$ in an $n$-dimensional space, being an orthonormal basis is a
single condition: **orthonormality**,
$$\langle e_i|e_j\rangle=\delta_{ij}.$$
In other words, distinct basis vectors are perpendicular and each is a unit vector. The clean
operational test is that the **Gram matrix** $G_{ij}=\langle e_i|e_j\rangle$, the array of all
pairwise inner products, equals the identity.

Fix the two $S_x$ eigenvectors $|{+}x\rangle=\tfrac{1}{\sqrt2}(1,1)$ and
$|{-}x\rangle=\tfrac{1}{\sqrt2}(1,-1)$ and build their Gram matrix:

```wl
ClearAll[plusx, minusx, psi, coeffs, gram, U, newCoords, A, orig, active, passive, P, Aprime];
{plusx, minusx} = {{1, 1}, {1, -1}}/Sqrt[2];
gram = Outer[Conjugate[#1] . #2 &, {plusx, minusx}, {plusx, minusx}, 1]
```

As one can see, the off-diagonal inner products vanish and the diagonal ones are unit, so the
Gram matrix has collapsed to the identity. Confirm that in a single test:

```wl
gram == IdentityMatrix[2]
```

Orthonormality of $n$ vectors already forces linear independence and hence spanning, so this one
test certifies a basis; you never have to check completeness separately in finite dimension.
That is also the subtlety behind the word "conditions": orthonormality is *sufficient* only
because we fixed the count at $n$. Give me $n-1$ orthonormal vectors and the Gram test still
passes, but they do not span; give me $n+1$ and orthonormality is impossible. The Gram-identity
test decides orthonormality, and the dimension count decides basis-hood.

### 3. How do expansion coefficients become a column vector?

> *How does the expansion of a vector in one orthonormal basis determine its column-vector representation?*

Fix an orthonormal basis $\{|i\rangle\}$. Any state expands as
$$|\psi\rangle=\sum_i c_i\,|i\rangle,\qquad c_i=\langle i|\psi\rangle,$$
and the ordered list $(c_1,\dots,c_n)^\mathsf{T}$ *is* the column-vector representation of
$|\psi\rangle$ **in that basis**. In other words, the coefficient in front of $|i\rangle$ is
nothing more than the overlap of $|\psi\rangle$ with $|i\rangle$, which follows by taking
$\langle j|$ of both sides and using $\langle j|i\rangle=\delta_{ji}$.

Read off the coefficients of the state $|\psi\rangle=\tfrac15(3,4i)$ in the $S_x$ basis of the
previous cell, using $c_i=\langle f_i|\psi\rangle$:

```wl
psi = {3, 4 I}/5;
coeffs = Conjugate[{plusx, minusx}] . psi
```

These two complex numbers are the column vector of the same physical ket, now written in the
$S_x$ chart rather than the standard one. Confirm that reassembling the state from its
coefficients returns the original ket, so nothing was lost:

```wl
FullSimplify[coeffs[[1]] plusx + coeffs[[2]] minusx == psi]
```

As expected, the round trip is exact: passing to coefficients and back loses no information,
which is precisely what completeness of the basis buys us. Note that the column vector is not a
property of the state alone; it is the pair (state, basis). Change the basis and the same ket
gets a different column, which is exactly what the next question watches happen.

### 4. How do coordinates change under a change of basis?

> *How do the coordinates of a state vector change under a change of orthonormal basis?*

Let $\{|e_j\rangle\}$ and $\{|f_i\rangle\}$ be two orthonormal bases. The same physical ket has
coordinates $c^{(e)}_j=\langle e_j|\psi\rangle$ in the first and $c^{(f)}_i=\langle f_i|\psi\rangle$
in the second, related by the **transition matrix**
$$U_{ij}=\langle f_i|e_j\rangle,\qquad c^{(f)}=U\,c^{(e)},$$
which is unitary because both bases are orthonormal. In other words, $U$ just re-expresses the
old coordinates in the new basis: it moves the labels, not the vector. Operators transform along
with it by conjugation, $A^{(f)}=U A^{(e)} U^\dagger$.

Build the transition matrix from the standard basis to the $S_x$ basis. Because the $e_j$ are
standard, $\langle f_i|e_j\rangle$ is just the $j$-th component of $\langle f_i|$:

```wl
U = Conjugate[{plusx, minusx}]
```

As one can see, this is the (real) Hadamard matrix, the machine that rotates between the $z$ and
$x$ bases. Verify that it is unitary, as any map between two orthonormal bases must be:

```wl
UnitaryMatrixQ[U]
```

Now push the standard-basis coordinates of $|\psi\rangle$ through $U$ and confirm that
reassembling in the $f$-basis returns the same physical ket:

```wl
newCoords = U . psi;
FullSimplify[newCoords[[1]] plusx + newCoords[[2]] minusx == psi]
```

As expected, the coordinates changed but the arrow did not. That is the whole mental picture: a
change of basis is a passive relabeling of the *same* vector by new numbers, so any
basis-independent quantity (norm, expectation value) is untouched, as question 11 makes precise.

### 5. Active unitary transformation versus passive change of basis.

> *What is the distinction between an active unitary transformation of a state and a passive change of basis?*

Both operations are described by a unitary matrix, which is exactly why they are confused. The
distinction is *what actually changes*.

- **Active:** the physical state is transformed, $|\psi\rangle\mapsto U|\psi\rangle$, in a *fixed*
  basis. A different arrow now sits where the old one was, and expectation values generally
  change: $\langle\psi|U^\dagger A U|\psi\rangle$.
- **Passive:** the state is untouched; only the basis rotates. The *same* arrow gets new
  coordinates $c\mapsto Uc$, the operator becomes $UAU^\dagger$, and every physical prediction is
  unchanged.

In short, the same matrix $U$ means opposite things depending on whether it acts on the state or
on the basis, and the cleanest way to tell them apart is to watch an expectation value.

Record the expectation value $\langle\psi|A|\psi\rangle$ of a Hermitian observable in the
original state, to serve as a reference:

```wl
A = {{2, I}, {-I, 3}};
orig = Conjugate[psi] . A . psi
```

Now take the active picture: move the state to $U|\psi\rangle$ and test whether its expectation
value still equals the reference:

```wl
FullSimplify[Conjugate[U . psi] . A . (U . psi) == orig]
```

The equality fails, and that failure is the physics: an active transformation genuinely moves the
state, so what you measure changes. Contrast the passive picture, where we rotate both the
coordinates and the operator, $A\to UAU^\dagger$, and test again:

```wl
FullSimplify[Conjugate[newCoords] . (U . A . ConjugateTranspose[U]) . newCoords == orig]
```

This time the expectation is invariant: relabeling the basis changes no prediction. So the two
pictures are physically opposite, and which one you mean is a statement about the apparatus, not
about the algebra.

### 6. Why is $\dim=2s+1$ for spin-$s$?

> *Why is the dimension of a spin-$s$ Hilbert space equal to $2s+1$?*

The spin operators satisfy $[J_i,J_j]=i\,\epsilon_{ijk}J_k$. Diagonalize $J_z$ and $J^2$
together. The raising and lowering operators $J_\pm=J_x\pm iJ_y$ shift the $J_z$ eigenvalue $m$ by
$\pm1$, and $J_\pm|j,\pm j\rangle=0$ truncates the ladder at both ends. In other words, the
ladder has a top rung and a bottom rung with an integer number of steps between them, which
forces $j$ to be a nonnegative half-integer and the $m$ values to run
$$m=-j,\,-j+1,\,\dots,\,j-1,\,j,$$
that is $2j+1$ values, one basis vector each. Here $s\equiv j$.

Build $J_z$ for a few spins and read its eigenvalues, which are the allowed $m$ values:

```wl
Sort @ Eigenvalues[QuantumOperator["JZ"[#]]["Matrix"] // Normal] & /@ {1/2, 1, 3/2}
```

As one can see, each spectrum runs from $-s$ to $s$ in unit steps, so it holds exactly $2s+1$
values. The dimension is therefore not an assumption but a consequence of the angular-momentum
algebra together with the demand that the space be finite, which is why a spin-$\tfrac12$ is a
qubit and a spin-$1$ is a qutrit.

### 7. Completeness as a resolution of the identity.

> *In a finite-dimensional Hilbert space, how can completeness of an orthonormal basis be expressed as a resolution of the identity?*

That a basis is complete is captured by one operator equation, the **resolution of the identity**:
$$\sum_i |i\rangle\langle i| = I.$$
In other words, the rank-one projectors onto the basis vectors, added up, reconstruct the
do-nothing operator. This one line is the workhorse of Dirac notation: inserting
$I=\sum_i|i\rangle\langle i|$ between any two objects expands them in the basis.

Sum the outer products over the full standard basis of $\mathbb{C}^3$ and confirm the total is
the identity:

```wl
Sum[KroneckerProduct[UnitVector[3, i], UnitVector[3, i]], {i, 3}] == IdentityMatrix[3]
```

Completeness (the sum is the whole identity) is distinct from orthonormality (the cross terms
vanish); both are needed, and together they say the projectors $\{|i\rangle\langle i|\}$ form a
complete orthogonal measurement, as in question 53.

### 8. The dual vector and the inner product in bra-ket notation.

> *How does one construct the dual vector associated with a ket, and how does the inner product appear in bra-ket notation?*

To each ket $|\psi\rangle$ the inner product assigns a linear functional $\langle\psi|$, the
**bra**. Concretely, if $|\psi\rangle$ is the column vector $c$, then $\langle\psi|$ is the row
vector $c^\dagger$ (conjugate transpose), and
$$\langle\psi|\phi\rangle=\sum_i \overline{c_i}\,d_i .$$
In other words, forming the bra conjugates the amplitudes; it is not merely a transpose. The
conjugation is what makes $\langle\psi|\psi\rangle=\sum_i|c_i|^2\ge0$ a real, positive norm.

In QuantumFramework the dagger turns a ket into its bra; pair it with another ket to get the
overlap $\langle\psi|\phi\rangle$:

```wl
qpsi = QuantumState[{3, 4 I}/5]; qphi = QuantumState[{1, I}/Sqrt[2]];
(qpsi["Dagger"] @ qphi)["Scalar"]
```

That single complex number is the overlap of the two states. Confirm that the framework's dagger
reproduces exactly the hand-built dual pairing $\sum_i \overline{c_i} d_i$, so both routes name
the same object:

```wl
FullSimplify[(qpsi["Dagger"] @ qphi)["Scalar"] == Conjugate[{3, 4 I}/5] . ({1, I}/Sqrt[2])]
```

As expected, the two agree. Note that the correspondence ket $\leftrightarrow$ bra is
*antilinear*: the bra of $\alpha|\psi\rangle$ is $\overline{\alpha}\langle\psi|$. This is why a
global phase $e^{i\gamma}$ on a ket becomes $e^{-i\gamma}$ on the bra and cancels in
$\langle\psi|\psi\rangle$, the seed of the global-phase irrelevance of question 13.

### 9. Subspace versus invariant subspace versus eigenspace.

> *What is the difference between a subspace, an invariant subspace, and an eigenspace of an operator?*

These are three nested notions of "a smaller space inside $\mathcal{H}$," each stronger than the
last, and the distinction is always *relative to an operator*.

- A **subspace** $\mathcal{S}\subseteq\mathcal{H}$ is any set closed under addition and scalar
  multiplication. It refers to no operator: the $xy$-plane in $\mathbb{C}^3$ is a subspace.
- An **invariant subspace** of $A$ is a subspace that $A$ maps into itself,
  $A\mathcal{S}\subseteq\mathcal{S}$. Vectors may move *within* $\mathcal{S}$, they just cannot
  leave it. Invariant subspaces are how a big operator block-diagonalizes into independent pieces.
- An **eigenspace** $\mathcal{S}_\lambda=\{|v\rangle:A|v\rangle=\lambda|v\rangle\}$ is the special
  invariant subspace on which $A$ acts as a single scalar $\lambda$. Every eigenspace is
  invariant, but an invariant subspace need not be an eigenspace: $A$ may act nontrivially inside
  it.

The hierarchy is subspace $\supset$ invariant subspace $\supset$ eigenspace. In quantum mechanics
the eigenspaces of an observable are the outcome sectors: measuring $A$ projects the state into
one of them, and a **degenerate** eigenvalue is one whose eigenspace has dimension greater than
1 (questions 29, 30).

### 10. Projecting a vector onto a subspace.

> *How does one project an arbitrary vector onto a given subspace using an orthogonal projector?*

Let $\mathcal{S}$ have orthonormal basis $\{|b_k\rangle\}$. The **orthogonal projector** onto
$\mathcal{S}$ is
$$P_\mathcal{S}=\sum_k |b_k\rangle\langle b_k|,$$
and $P_\mathcal{S}|\psi\rangle$ keeps the components of $|\psi\rangle$ in $\mathcal{S}$ and
discards the rest. A projector is fixed by two properties: **idempotence** $P^2=P$ (projecting
twice is the same as projecting once) and **Hermiticity** $P^\dagger=P$ (so probabilities read
from it are real).

Build the projector onto $\mathrm{span}\{|{+}x\rangle\}$ and verify it is idempotent:

```wl
P = Outer[Times, plusx, Conjugate[plusx]];
FullSimplify[P . P == P]
```

As expected, applying the projector twice is the same as applying it once. Verify it is Hermitian
too, so that the probabilities read from it come out real:

```wl
FullSimplify[P == ConjugateTranspose[P]]
```

Finally, apply $P$ to $|\psi\rangle$ and confirm it returns exactly the $|{+}x\rangle$ component
$\langle{+}x|\psi\rangle\,|{+}x\rangle$:

```wl
FullSimplify[P . psi == (Conjugate[plusx] . psi) plusx]
```

Therefore $P$ keeps the part of $|\psi\rangle$ lying in the subspace and discards the rest. Two
facts make projectors the backbone of measurement theory: the complement $Q=I-P$ satisfies
$P+Q=I$, and the squared norm $\langle\psi|P|\psi\rangle=\|P|\psi\rangle\|^2$ is the Born
probability of finding the system in $\mathcal{S}$ (question 42).

### 11. What is basis-independent, and what is representation-dependent?

> *What information is basis-independent in a finite-dimensional quantum state, and what information is representation-dependent?*

The column vector of a state and the matrix of an operator are *representation-dependent*: they
are coordinates in a chosen basis and change under the transition matrix $U$ (question 4). What
survives a change of basis is the physics. In other words, coordinates are bookkeeping, and the
physics is whatever the relabeling leaves untouched: the **norm**, overlaps
$|\langle\phi|\psi\rangle|$, an operator's **spectrum** and **trace**, and **expectation values**.

Conjugate the observable $A$ from question 5 into the $S_x$ basis, $A'=UAU^\dagger$, and verify
the trace is unchanged:

```wl
Aprime = U . A . ConjugateTranspose[U];
FullSimplify[Tr[A] == Tr[Aprime]]
```

As expected, the trace is the same in both bases. Verify the spectrum is unchanged too:

```wl
FullSimplify[Sort[Eigenvalues[A]] == Sort[Eigenvalues[Aprime]]]
```

Therefore both the trace and the eigenvalues are basis-independent. The practical rule follows:
any number you could measure is basis-independent, and anything that references specific basis
labels (a particular amplitude $c_3$, an off-diagonal element $A_{12}$) is a bookkeeping artifact
of the chart you drew. This is the linear-algebra face of the statement that a pure state is a
ray, not a column of numbers (question 14).

### 12. Modeling a two-level system without wavefunctions.

> *How can a finite-dimensional Hilbert space model a two-level atom, a spin-$1/2$ particle, or a qubit without reference to position wavefunctions?*

Nothing in the quantum formalism requires a state to be a function $\psi(x)$ of position. The
axioms need only a Hilbert space, Hermitian observables, and the Born rule. A **two-level system**
takes the smallest nontrivial choice, $\mathcal{H}=\mathbb{C}^2$, and reads off all the physics
from $2\times2$ linear algebra. The two basis kets $|0\rangle,|1\rangle$ are whatever pair of
perfectly distinguishable alternatives the system offers:

- a **spin-$\tfrac12$**: $|0\rangle=|{\uparrow}\rangle$, $|1\rangle=|{\downarrow}\rangle$ along $z$;
- a **two-level atom**: $|0\rangle=$ ground, $|1\rangle=$ excited, with the rest of the atomic
  ladder far off-resonance and dropped (a finite truncation, question 148);
- a **qubit**: $|0\rangle,|1\rangle$ are abstract logical values with no spatial meaning at all.

All three share one Hilbert space, so they share one mathematics. A general pure state is
$|\psi\rangle=\cos\tfrac{\theta}{2}\,|0\rangle+e^{i\varphi}\sin\tfrac{\theta}{2}\,|1\rangle$
(question 17), observables are $2\times2$ Hermitian matrices $a_0 I+\vec a\cdot\vec\sigma$
(question 85), and dynamics is a rotation of the Bloch vector (question 64). The position
wavefunction is not a prerequisite for quantum mechanics; it is one particular
infinite-dimensional representation, and the finite case is where the conceptual structure is
cleanest and where quantum information lives.

---

*End of Section A (questions 1-12). Sections B-L follow the same format once this calibration
sample is approved.*
