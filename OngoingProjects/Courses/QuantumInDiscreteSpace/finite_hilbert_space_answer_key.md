# Answer Key: Quantum Questions in Finite-Dimensional Hilbert Space

Pedagogical answers to the 185 questions in
[`finite_hilbert_space_question_bank.md`](finite_hilbert_space_question_bank.md). Each answer
opens with the full question exactly as posed in the bank (in italics), then states the central
idea in physical terms, gives the closed-form result, explains the moving parts, and, where a
computation is the natural answer, demonstrates it with Wolfram Language or QuantumFramework
code.

## How to read this key

The code follows a computation-first rhythm: a sentence names what to compute, one cell computes
exactly one thing, and the next sentence reads the returned value back. Every cell has been
executed in a Wolfram kernel (QuantumFramework 2.0.0); the value quoted in the prose is the value
the cell returns, occasionally rewritten in a tidier but algebraically identical closed form. Within a section the cells share one running session in reading order, so a
symbol defined in an earlier cell (`plusx`, `psi`, `U`) is still in scope later. Snippets follow
the *minimal faithful idiom*: the lightest construction that still shows the physics. Where a
question is purely conceptual, the answer is prose and no cell is forced onto it.

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
$\langle\psi|\psi\rangle>0$ for $|\psi\rangle\neq 0$.

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
In dimension $n$, orthonormality of $n$ vectors *already forces* linear independence and hence
spanning, so you do not have to verify completeness separately. The clean operational test is
that the **Gram matrix** $G_{ij}=\langle e_i|e_j\rangle$ equals the identity. Equivalently, the
matrix whose columns are the $|e_i\rangle$ is unitary.

Take the two $S_x$ eigenvectors $|{+}x\rangle=\tfrac{1}{\sqrt2}(1,1)$ and
$|{-}x\rangle=\tfrac{1}{\sqrt2}(1,-1)$ and form their Gram matrix, the array of all pairwise
inner products:

```wl
{plusx, minusx} = {{1, 1}, {1, -1}}/Sqrt[2];
gram = Outer[Conjugate[#1] . #2 &, {plusx, minusx}, {plusx, minusx}, 1]
```

As one can see, the off-diagonal inner products vanish and the diagonal ones are 1. Confirm the
whole matrix is the identity in a single test:

```wl
gram === IdentityMatrix[2]
```

`True`: the pair is orthonormal, and because there are exactly $n=2$ of them, they form a basis.
The "conditions" language of the question is where the concept hides: orthonormality is
*sufficient* for a basis only because we already fixed the count at $n$. Give me $n-1$
orthonormal vectors and the Gram test still returns the identity, but they do not span; give me
$n+1$ and orthonormality is impossible. The Gram-equals-identity test decides orthonormality; the
dimension count decides basis-hood.

### 3. How do expansion coefficients become a column vector?

> *How does the expansion of a vector in one orthonormal basis determine its column-vector representation?*

Fix an orthonormal basis $\{|i\rangle\}$. Any state expands as
$$|\psi\rangle=\sum_i c_i\,|i\rangle,\qquad c_i=\langle i|\psi\rangle,$$
and the ordered list $(c_1,\dots,c_n)^\mathsf{T}$ *is* the column-vector representation of
$|\psi\rangle$ **in that basis**. The coefficient formula is just orthonormality applied to the
expansion: take $\langle j|$ of both sides and use $\langle j|i\rangle=\delta_{ji}$.

Read off the coefficients of the state $|\psi\rangle=\tfrac15(3,4i)$ in the $S_x$ basis of the
previous cell, using $c_i=\langle f_i|\psi\rangle$:

```wl
psi = {3, 4 I}/5;
coeffs = Conjugate[{plusx, minusx}] . psi
```

These two numbers, $\langle f_i|\psi\rangle=(3\pm 4i)/(5\sqrt2)$, are the column vector of the
same physical ket in the $S_x$ chart. Now reassemble the state from its coefficients and confirm
nothing was lost:

```wl
FullSimplify[coeffs[[1]] plusx + coeffs[[2]] minusx - psi]
```

`{0, 0}`: the reconstruction $\sum_i c_i|i\rangle=|\psi\rangle$ is exact. The column vector is not
a property of the state alone; it is the pair (state, basis). Change the basis and the same
physical ket gets a different column, which is exactly what the next question watches happen.

### 4. How do coordinates change under a change of basis?

> *How do the coordinates of a state vector change under a change of orthonormal basis?*

Let $\{|e_j\rangle\}$ and $\{|f_i\rangle\}$ be two orthonormal bases. The same physical ket has
coordinates $c^{(e)}_j=\langle e_j|\psi\rangle$ in the first and $c^{(f)}_i=\langle f_i|\psi\rangle$
in the second, related by the **transition matrix**
$$U_{ij}=\langle f_i|e_j\rangle,\qquad c^{(f)}=U\,c^{(e)},$$
which is unitary because both bases are orthonormal. Operators transform by conjugation,
$A^{(f)}=U A^{(e)} U^\dagger$, so the abstract relation $A|\psi\rangle$ holds in every basis.

Build the transition matrix from the standard basis to the $S_x$ basis. Because the $e_j$ are
standard, $\langle f_i|e_j\rangle$ is just the $j$-th component of $\langle f_i|$:

```wl
U = Conjugate[{plusx, minusx}]
```

This is the (real) Hadamard matrix. Confirm it is unitary, as any map between orthonormal bases
must be:

```wl
UnitaryMatrixQ[U]
```

`True`. Now push the standard-basis coordinates of $|\psi\rangle$ through $U$ and reassemble in
the $f$-basis to check the physical ket is unchanged:

```wl
newCoords = U . psi;
FullSimplify[newCoords[[1]] plusx + newCoords[[2]] minusx - psi]
```

`{0, 0}`: the coordinates changed, the arrow did not. That is the whole mental picture. A change
of basis is a passive relabeling of the *same* vector by new numbers, so any basis-independent
quantity (norm, expectation value) stays put, as question 11 makes precise.

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

Take a Hermitian observable and record its expectation in the original state:

```wl
A = {{2, I}, {-I, 3}};
orig = Conjugate[psi] . A . psi
```

The reference value is $\langle\psi|A|\psi\rangle=42/25$. Now apply the active picture: move the
state to $U|\psi\rangle$ and recompute the expectation:

```wl
active = Conjugate[U . psi] . A . (U . psi);
FullSimplify[active - orig]
```

A shift of $89/50$: the state genuinely moved, so what you measure changes. Contrast the passive
picture, where we rotate both the coordinates and the operator, $A\to UAU^\dagger$:

```wl
passive = Conjugate[newCoords] . (U . A . ConjugateTranspose[U]) . newCoords;
FullSimplify[passive - orig]
```

Exactly `0`: relabeling the basis changes no prediction. So the same matrix $U$ means two
physically opposite things depending on whether it acts on the state (active) or on the chart
(passive). Which one you intend is a physics statement, not a mathematical one.

### 6. Why is $\dim=2s+1$ for spin-$s$?

> *Why is the dimension of a spin-$s$ Hilbert space equal to $2s+1$?*

The spin operators satisfy $[J_i,J_j]=i\,\epsilon_{ijk}J_k$. Diagonalize $J_z$ and $J^2$
together. The raising and lowering operators $J_\pm=J_x\pm iJ_y$ shift the $J_z$ eigenvalue $m$ by
$\pm1$, and $J_\pm|j,\pm j\rangle=0$ truncates the ladder at both ends. Consistency forces $j$ to
be a nonnegative half-integer and the $m$ values to run
$$m=-j,\,-j+1,\,\dots,\,j-1,\,j,$$
which is $2j+1$ values, one basis vector each. Here $s\equiv j$.

Build $J_z$ for a few spins and read its eigenvalues, the allowed $m$ values:

```wl
Sort @ Eigenvalues[QuantumOperator["JZ"[#]]["Matrix"] // Normal] & /@ {1/2, 1, 3/2}
```

The three lists have lengths 2, 3, 4, that is $2s+1$, and their entries are exactly the ladders
$\{-\tfrac12,\tfrac12\}$, $\{-1,0,1\}$, $\{-\tfrac32,-\tfrac12,\tfrac12,\tfrac32\}$. The dimension
is not an assumption but a consequence of the angular-momentum algebra together with the demand
that the space be finite. This is why a spin-$\tfrac12$ is a qubit and a spin-$1$ is a qutrit.

### 7. Completeness as a resolution of the identity.

> *In a finite-dimensional Hilbert space, how can completeness of an orthonormal basis be expressed as a resolution of the identity?*

That a basis is complete is captured by one operator equation, the **resolution of the identity**:
$$\sum_i |i\rangle\langle i| = I.$$
Each $|i\rangle\langle i|$ is the rank-one projector onto $|i\rangle$; summing over a full
orthonormal basis reassembles the identity. This one line is the workhorse of Dirac notation:
inserting $I=\sum_i|i\rangle\langle i|$ between any two objects expands them in the basis.

Sum the outer products over the full standard basis of $\mathbb{C}^3$ and test the total against
the identity:

```wl
Sum[KroneckerProduct[UnitVector[3, i], UnitVector[3, i]], {i, 3}] === IdentityMatrix[3]
```

`True`: the basis is complete. Completeness (the sum is the whole identity) is distinct from
orthonormality (the cross terms vanish); both are needed, and together they say the projectors
$\{|i\rangle\langle i|\}$ form a complete orthogonal measurement, as in question 53.

### 8. The dual vector and the inner product in bra-ket notation.

> *How does one construct the dual vector associated with a ket, and how does the inner product appear in bra-ket notation?*

To each ket $|\psi\rangle$ the inner product assigns a linear functional $\langle\psi|$, the
**bra**. Concretely, if $|\psi\rangle$ is the column vector $c$, then $\langle\psi|$ is the row
vector $c^\dagger$ (conjugate transpose), and
$$\langle\psi|\phi\rangle=\sum_i \overline{c_i}\,d_i .$$
The conjugation is not cosmetic: it makes $\langle\psi|\psi\rangle=\sum_i|c_i|^2\ge0$ a real,
positive norm and enforces conjugate symmetry.

In QuantumFramework the dagger turns a ket into its bra; pair it with another ket to get the
overlap:

```wl
qpsi = QuantumState[{3, 4 I}/5]; qphi = QuantumState[{1, I}/Sqrt[2]];
(qpsi["Dagger"] @ qphi)["Scalar"]
```

The overlap is $\langle\psi|\phi\rangle=7/(5\sqrt2)$. Confirm this is exactly the hand-built dual
pairing $\sum_i \overline{c_i} d_i$, so the framework's dagger and the raw conjugate transpose
are the same object seen two ways:

```wl
(qpsi["Dagger"] @ qphi)["Scalar"] === Conjugate[{3, 4 I}/5] . ({1, I}/Sqrt[2])
```

`True`. The correspondence ket $\leftrightarrow$ bra is *antilinear*: the bra of
$\alpha|\psi\rangle$ is $\overline{\alpha}\langle\psi|$. This is why a global phase $e^{i\gamma}$
on a ket becomes $e^{-i\gamma}$ on the bra and cancels in $\langle\psi|\psi\rangle$, the seed of
the global-phase irrelevance of question 13.

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
twice is projecting once) and **Hermiticity** $P^\dagger=P$ (the projection is orthogonal).

Build the projector onto $\mathrm{span}\{|{+}x\rangle\}$ and check idempotence first:

```wl
P = Outer[Times, plusx, Conjugate[plusx]];
FullSimplify[P . P - P]
```

The zero matrix: $P^2=P$. Now check Hermiticity, so that probabilities read off from $P$ are real:

```wl
FullSimplify[P - ConjugateTranspose[P]]
```

Also zero: $P=P^\dagger$. Finally, apply $P$ to $|\psi\rangle$ and confirm it returns exactly the
$|{+}x\rangle$ component $\langle{+}x|\psi\rangle\,|{+}x\rangle$:

```wl
FullSimplify[P . psi - (Conjugate[plusx] . psi) plusx]
```

`{0, 0}`: the projection keeps only the $|{+}x\rangle$ part. Two facts make projectors the
backbone of measurement theory: the complement $Q=I-P$ satisfies $P+Q=I$, and the squared norm
$\langle\psi|P|\psi\rangle=\|P|\psi\rangle\|^2$ is the Born probability of finding the system in
$\mathcal{S}$ (question 42).

### 11. What is basis-independent, and what is representation-dependent?

> *What information is basis-independent in a finite-dimensional quantum state, and what information is representation-dependent?*

The column vector of a state and the matrix of an operator are *representation-dependent*: they
are coordinates in a chosen basis and change under the transition matrix $U$ (question 4). What
survives a change of basis is the physics: quantities built from the inner-product structure and
invariant under unitary conjugation, such as the **norm**, overlaps $|\langle\phi|\psi\rangle|$,
an operator's **spectrum** and **trace**, and **expectation values**.

Conjugate the observable $A$ from question 5 into the $S_x$ basis, $A'=UAU^\dagger$, and check the
trace is unchanged:

```wl
Aprime = U . A . ConjugateTranspose[U];
FullSimplify[Tr[A] - Tr[Aprime]]
```

`0`: the trace is basis-independent. Check the spectrum too:

```wl
FullSimplify[Sort[Eigenvalues[A]] - Sort[Eigenvalues[Aprime]]]
```

`{0, 0}`: identical eigenvalues in both bases. The practical rule follows: any number you could
measure is basis-independent, and anything that references specific basis labels (a particular
amplitude $c_3$, an off-diagonal element $A_{12}$) is a bookkeeping artifact of the chart you
drew. This is the linear-algebra face of the statement that a pure state is a ray, not a column
of numbers (question 14).

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
