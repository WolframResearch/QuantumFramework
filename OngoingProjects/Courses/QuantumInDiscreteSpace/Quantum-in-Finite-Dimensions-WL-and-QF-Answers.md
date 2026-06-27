---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework

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

## Part 1. Pure states and the Born rule

### 1.1 [BSc] How do I represent a qubit as a state vector $\{\alpha,\beta\}$ in $\mathbb{C}^2$?

A pure qubit is a two-component complex vector $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$.

**WL** : a qubit is literally a length-2 complex list, kept fully general.

```wl
{\[Alpha], \[Beta]}
```

**QF** : wrap the amplitudes in `QuantumState`. Evaluated without a trailing semicolon, the
construction displays as its **summary box**: a compact panel that *is* the live object, showing its
type, dimension, and basis, not a printout of a matrix. Properties are read off that object on demand.

```wl
qs = QuantumState[{\[Alpha], \[Beta]}]
```

Then `"StateVector"` reads one representation, the amplitudes, back out:

```wl
qs["StateVector"] // Normal
```

Both hold the amplitudes $\{\alpha,\beta\}$; the QF object additionally knows it is a single qubit in
the computational basis, and the summary box is that object itself, with the amplitude vector just one
field read from it.

### 1.2 [BSc] How do I read the Born-rule probabilities of a state?

The probability of outcome $i$ is $|\langle i|\psi\rangle|^2$, i.e. the squared modulus of the
amplitude (divided by the norm if the state is not yet normalized).

**WL** : square the moduli of the normalized amplitudes; `FullSimplify` reduces the nested modulus to
a single ratio (the numerator $\alpha\bar\alpha$ is $|\alpha|^2$).

```wl
Abs[Normalize[{\[Alpha], \[Beta]}]]^2 // FullSimplify
```

**QF** : ask the state for its probabilities; `"ProbabilitiesList"` normalizes internally and returns
them in basis order.

```wl
QuantumState[{\[Alpha], \[Beta]}]["ProbabilitiesList"]
```

Both give $\bigl\{|\alpha|^2,\,|\beta|^2\bigr\}/(|\alpha|^2+|\beta|^2)$: the Born rule weights the
outcomes by the squared amplitudes, unequally in general, and the weights sum to one.

### 1.3 [BSc] How do I check that a state is normalized, and normalize one that is not?

A physical state vector must have unit norm, $\langle\psi|\psi\rangle = 1$.

**WL** : `Norm` measures the length of the raw amplitude vector.

```wl
Norm[{\[Alpha], \[Beta]}]
```

`Normalize` rescales it to unit length, so the norm of the normalized vector is $1$.

```wl
Simplify[Norm @ Normalize[{\[Alpha], \[Beta]}]]
```

**QF** : the `"Norm"` property reports the same raw length.

```wl
QuantumState[{\[Alpha], \[Beta]}]["Norm"]
```

`"Normalize"` returns a new, unit-norm state object, whose `"Norm"` is then $1$.

```wl
QuantumState[{\[Alpha], \[Beta]}]["Normalize"]["Norm"]
```

Both report the raw norm $\sqrt{|\alpha|^2+|\beta|^2}$ and confirm the normalized state has norm $1$.

### 1.4 [BSc] How do I show that a global phase is physically unobservable?

Multiplying $|\psi\rangle$ by $e^{i\varphi}$ changes no measurable prediction, because every
probability depends only on $|c_i|^2$.

**WL** : compare the two probability vectors symbolically for all real phases.

```wl
FullSimplify[
 Abs[{\[Alpha], \[Beta]}]^2 == Abs[Exp[I \[CurlyPhi]] {\[Alpha], \[Beta]}]^2,
 \[CurlyPhi] \[Element] Reals]
```

**QF** : compare the two state objects directly. QuantumFramework's `==` is a *physical* comparison:
it equates states that differ only by a global phase (it uses an inner-product oracle, not the raw
amplitudes, so `===` would still see them as distinct).

```wl
QuantumState[{\[Alpha], \[Beta]}] == QuantumState[Exp[I \[CurlyPhi]] {\[Alpha], \[Beta]}]
```

Both return `True`, symbolically, for every phase $\varphi$: the WL test shows the Born
probabilities are phase-independent, and QuantumFramework declares the two state objects *equal*. The
phase labels no physical degree of freedom.

### 1.5 [BSc] How do I count the real parameters of a pure qubit and a pure qudit in $\mathbb{C}^d$?

A pure state in $\mathbb{C}^d$ has $2d$ real amplitude components; normalization removes one and the
global phase removes another, leaving $2d-2$. The rigorous count is the generic rank of the
parametrization's Jacobian: embed the family in $\mathbb{R}^4$ and count independent real directions.

**WL** : keep the global phase $g$, and the qubit family $e^{ig}\{\cos(\theta/2),\,e^{i\phi}\sin(\theta/2)\}$ fills the 3-sphere of normalized states, Jacobian rank $3$.

```wl
With[{psi = Exp[I g] {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 MatrixRank[D[Flatten[ReIm[ComplexExpand[psi]]], {{\[Theta], \[Phi], g}}]]]
```

Quotient the global phase (drop $g$) and the rank falls to $2$, the physical degrees of freedom.

```wl
With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 MatrixRank[D[Flatten[ReIm[ComplexExpand[psi]]], {{\[Theta], \[Phi]}}]]]
```

The general count for $\mathbb{C}^d$ is $2d-2$ ($2d$ components, less normalization, less global phase).

```wl
Table[2 d - 2, {d, 2, 5}]
```

**QF** : drive the general count from the object's own Hilbert-space dimension rather than a literal.

```wl
Table[2 (QuantumState["0", d]["Dimension"] - 1), {d, 2, 5}]
```

The Jacobian rank is $3$ with the phase and $2$ without it: a qubit has exactly two physical real
parameters (the two Bloch angles), and the difference of one is the unobservable global phase. The
general count is $\{2,4,6,8\}$ for $d = 2,3,4,5$.

### 1.6 [BSc] How do I compute an inner product $\langle\phi|\psi\rangle$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?

The overlap of two states, its absolute value, and its absolute square give the amplitude and the
probability to find one state in the other.

**WL** : the inner product is the *conjugated* dot product; the conjugation acts on the bra
$\langle\phi|$, which is why $\alpha_1,\beta_1$ appear conjugated. The cell returns the overlap, its
modulus, and its square together, as one inner product read three ways.

```wl
With[{phi = {\[Alpha]1, \[Beta]1}, psi = {\[Alpha]2, \[Beta]2}},
 With[{ip = Conjugate[phi] . psi}, {ip, Abs[ip], Abs[ip]^2}]]
```

**QF** : build the bra with `"Dagger"`, compose with the ket, and read the same three numbers off the
scalar.

```wl
With[{bra = QuantumState[{\[Alpha]1, \[Beta]1}]["Dagger"], ket = QuantumState[{\[Alpha]2, \[Beta]2}]},
 With[{ip = (bra @ ket)["Scalar"]}, {ip, Abs[ip], Abs[ip]^2}]]
```

Both give $\langle\phi|\psi\rangle = \alpha_2\bar\alpha_1 + \beta_2\bar\beta_1$, with its modulus and
its square. The conjugation on the bra is what makes this a sesquilinear inner product rather than a
plain bilinear dot product.

### 1.7 [BSc] How do I write a state both as a column vector of amplitudes and as a Dirac superposition $\sum_i c_i\,|i\rangle$?

The amplitude list and the ket-superposition are two views of one object.

**WL** : the amplitude vector is the list itself.

```wl
{\[Alpha], \[Beta]}
```

The Dirac form is the same numbers against the computational basis; the identity confirms the
expansion.

```wl
{\[Alpha], \[Beta]} == \[Alpha] {1, 0} + \[Beta] {0, 1}
```

**QF** : the state object, shown first as its summary box.

```wl
qs = QuantumState[{\[Alpha], \[Beta]}]
```

`"StateVector"` reads the amplitudes:

```wl
qs["StateVector"] // Normal
```

`"Formula"` typesets the same object as the Dirac superposition $\sum_i c_i|i\rangle$:

```wl
qs["Formula"]
```

The amplitude vector and the superposition $\alpha|0\rangle + \beta|1\rangle$ carry identical
information; both are representations read from the one object whose summary box is shown above.

### 1.8 [BSc] How do I get the Bloch vector of a qubit, and the Bloch-sphere angles?

The Bloch vector $\vec r = (\langle\sigma_x\rangle, \langle\sigma_y\rangle, \langle\sigma_z\rangle)$
places a pure qubit on the unit sphere, with $(\theta,\phi)$ its spherical angles.

**WL** : the Bloch vector is, *by definition*, the triple of Pauli expectation values
$\langle\psi|\sigma_j|\psi\rangle$. Evaluate that definition for a general $\{\alpha,\beta\}$ to get
the familiar component formula.

```wl
FullSimplify @ Table[Conjugate[{\[Alpha], \[Beta]}] . PauliMatrix[j] . {\[Alpha], \[Beta]}, {j, 3}]
```

For the angle-parametrized state the same definition gives the geometric vector on the sphere.

```wl
rvec = With[{\[Psi] = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
  FullSimplify[Table[Conjugate[\[Psi]] . PauliMatrix[j] . \[Psi], {j, 3}], {\[Theta], \[Phi]} \[Element] Reals]]
```

A built-in coordinate transform turns that Cartesian vector into spherical angles.

```wl
FullSimplify[ToSphericalCoordinates[rvec], 0 < \[Theta] < Pi && -Pi < \[Phi] < Pi]
```

**QF** : the state object, shown as its summary box.

```wl
qs = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]
```

`"BlochVector"` computes the same Pauli expectation values internally (as $\mathrm{Tr}[\rho\,\sigma_j]$).

```wl
FullSimplify[qs["BlochVector"], 0 < \[Theta] < Pi && -Pi < \[Phi] < Pi]
```

`"BlochSphericalCoordinates"` gives the polar and azimuthal angles.

```wl
FullSimplify[qs["BlochSphericalCoordinates"], 0 < \[Theta] < Pi && -Pi < \[Phi] < Pi]
```

Evaluating the definition returns $\{2\,\mathrm{Re}\,\bar\alpha\beta,\ 2\,\mathrm{Im}\,\bar\alpha\beta,\ |\alpha|^2-|\beta|^2\}$ for a general state and $(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$ for the parametrized one, with spherical coordinates $(1,\theta,\phi)$. So the component formula is a *result* of evaluating the Pauli expectation values, not the definition itself, and the parametrization angles are exactly the polar and azimuthal angles on the Bloch sphere.

### 1.9 [BSc] How do I prepare a qudit in $\mathbb{C}^d$ and a uniform superposition?

A qudit lives in $\mathbb{C}^d$; the computational state $|0\rangle$ and the equal-weight
superposition $\tfrac1{\sqrt d}\sum_k |k\rangle$ are the basic preparations.

**WL** : $|0\rangle$ is a unit vector.

```wl
With[{d = 3}, UnitVector[d, 1]]
```

The uniform state is a constant vector handed to `Normalize` (no hardcoded $1/\sqrt d$).

```wl
With[{d = 3}, Normalize @ ConstantArray[1, d]]
```

**QF** : name the basis state with its dimension.

```wl
With[{d = 3}, QuantumState["0", d]["StateVector"] // Normal]
```

Use the built-in `"UniformSuperposition"` state rather than constructing it by hand.

```wl
With[{d = 3}, QuantumState["UniformSuperposition", d]["StateVector"] // Normal]
```

Both give $|0\rangle = (1,0,0)$ and the uniform qutrit $\tfrac1{\sqrt3}(1,1,1)$; the construction is
written in $d$, so it works for any dimension.

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
With[{\[Sigma] = Table[PauliMatrix[j], {j, 3}]},
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

### 2.3 [BSc] How do I test that an operator is Hermitian, and why must observables be Hermitian?

An observable is Hermitian, $A=A^\dagger$, and that is exactly what forces its eigenvalues (the
possible measurement outcomes) to be real.

**WL** : test Hermiticity directly.

```wl
HermitianMatrixQ[PauliMatrix[1]]
```

A generic Hermitian matrix then has real eigenvalues.

```wl
FullSimplify[Eigenvalues[{{a1, a2 + I a3}, {a2 - I a3, a4}}] \[Element] Reals,
 (a1 | a2 | a3 | a4) \[Element] Reals]
```

**QF** : the operator answers `"HermitianQ"`.

```wl
QuantumOperator["X"]["HermitianQ"]
```

And its `"Eigenvalues"` are real.

```wl
FullSimplify[QuantumOperator[{{a1, a2 + I a3}, {a2 - I a3, a4}}]["Eigenvalues"] \[Element] Reals,
 (a1 | a2 | a3 | a4) \[Element] Reals]
```

Both report `True` twice: the operator is Hermitian, and Hermiticity guarantees a real spectrum.

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

By definition $f(A)$ acts on the eigenvalues in the eigenbasis, $f(A) = \sum_a f(a)\,|a\rangle\langle a|$.

**WL** : build that spectral sum with a symbolic $f$ from the operator's eigensystem.

```wl
spec = With[{sys = Eigensystem[PauliMatrix[1]]},
  Total[MapThread[f[#1] KroneckerProduct[#2, Conjugate[#2]] &, {sys[[1]], Normalize /@ sys[[2]]}]]]
```

Confirm it equals the built-in `MatrixFunction`.

```wl
spec == MatrixFunction[f, PauliMatrix[1]]
```

**QF** : the same definition driven by the operator's own `"Eigensystem"`.

```wl
With[{sys = QuantumOperator["X"]["Eigensystem"]},
 Total[MapThread[f[#1] KroneckerProduct[#2, Conjugate[#2]] &, {sys[[1]], Normalize /@ sys[[2]]}]]]
```

The result $f(\sigma_x) = \tfrac12\begin{pmatrix} f(1)+f(-1) & f(1)-f(-1) \\ f(1)-f(-1) & f(1)+f(-1)\end{pmatrix}$
shows $f$ acting through the eigenbasis, not entry by entry, and the WL check confirms it agrees with
`MatrixFunction`.

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
$z$ with certainty at the poles and is maximally uncertain on the equator.

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

Two observables are compatible exactly when they commute, and then they share an eigenbasis. Each of
$X\otimes X$ and $Z\otimes Z$ is degenerate (eigenvalues $\pm1$, each doubled), so neither alone pins
the basis down; a generic combination $A + cB$ is non-degenerate, and *its* eigenvectors are the
simultaneous eigenbasis. So we find the basis, not assume it.

**WL** : compatibility is the vanishing commutator.

```wl
With[{a = KroneckerProduct[PauliMatrix[1], PauliMatrix[1]], b = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]},
 a . b == b . a]
```

Diagonalizing $A + \pi B$ once gives the candidate common eigenbasis (its rows the simultaneous
eigenvectors); bind it.

```wl
u = With[{a = KroneckerProduct[PauliMatrix[1], PauliMatrix[1]], b = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]},
  Normalize /@ Eigenvectors[a + \[Pi] b]];
```

That basis diagonalizes $A = X\otimes X$.

```wl
With[{a = KroneckerProduct[PauliMatrix[1], PauliMatrix[1]]},
 DiagonalMatrixQ[Chop[Conjugate[u] . a . Transpose[u]]]]
```

And it diagonalizes $B = Z\otimes Z$, so the basis is simultaneous.

```wl
With[{b = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]]},
 DiagonalMatrixQ[Chop[Conjugate[u] . b . Transpose[u]]]]
```

**QF** : compatibility is the *operator* equality `a @ b == b @ a`.

```wl
With[{a = QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["X"]],
   b = QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["Z"]]},
 a @ b == b @ a]
```

The common eigenbasis is the eigenvectors of `a + Pi b`. Each found state is an eigenstate of
$A = X\otimes X$ (it returns unchanged up to phase), checked without leaving the object model.

```wl
With[{a = QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["X"]],
   b = QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["Z"]]},
 With[{basis = QuantumState /@ (a + \[Pi] b)["Eigenvectors"]},
  AllTrue[basis, a[#] == # &]]]
```

And each is equally an eigenstate of $B = Z\otimes Z$, confirming the basis is simultaneous.

```wl
With[{a = QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["X"]],
   b = QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["Z"]]},
 With[{basis = QuantumState /@ (a + \[Pi] b)["Eigenvectors"]},
  AllTrue[basis, b[#] == # &]]]
```

The commutator vanishes and the checks pass. The basis returned by $A+\pi B$ is the Bell basis,
*derived* rather than assumed: the coefficient $\pi$ only breaks the degeneracy, and any value that
does so gives the same simultaneous eigenbasis.

### 2.9 [MSc] How do I construct a complete set of commuting observables (CSCO) and label states by their joint spectrum?

A complete set of commuting observables assigns every basis state a unique tuple of eigenvalues.

**WL** : first, $Z\otimes I$ and $I\otimes Z$ commute.

```wl
With[{a = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], b = KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]},
 a . b == b . a]
```

Their joint eigenvalues then label the four states.

```wl
With[{a = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2]], b = KroneckerProduct[IdentityMatrix[2], PauliMatrix[3]]},
 Transpose[{Diagonal[a], Diagonal[b]}]]
```

**QF** : the same commutation check with `QuantumTensorProduct`.

```wl
With[{a = QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["I"]],
   b = QuantumTensorProduct[QuantumOperator["I"], QuantumOperator["Z"]]},
 a @ b == b @ a]
```

And the joint spectrum, read from the diagonals.

```wl
With[{a = QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["I"]],
   b = QuantumTensorProduct[QuantumOperator["I"], QuantumOperator["Z"]]},
 Transpose[{Diagonal[Normal[a["Matrix"]]], Diagonal[Normal[b["Matrix"]]]}]]
```

Both confirm the two observables commute and return the joint spectrum
$\{(+,+),(+,-),(-,+),(-,-)\}$: the four pairs are distinct, so the two observables form a complete
label for the basis.

### 2.10 [MSc] How do I define the Hilbert-Schmidt operator inner product $\langle A,B\rangle=\mathrm{Tr}[A^\dagger B]$ and expand an operator in the Pauli basis?

Operators carry their own inner product, $\langle A,B\rangle=\mathrm{Tr}[A^\dagger B]$, under which the
Paulis are orthogonal; the expansion coefficients are $c_j = \tfrac12\mathrm{Tr}[\sigma_j A]$.

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

### 2.11 [BSc] How do I get a Rayleigh-Ritz variational upper bound on the ground-state energy from a trial state?

For any normalized trial state, $\langle\psi|H|\psi\rangle \ge E_0$, with equality only at the ground
state; minimizing the trial energy is the variational method.

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

The Rayleigh-Schrodinger shifts are $E_n^{(1)} = \langle n|V|n\rangle$ and
$E_n^{(2)} = \sum_{m\ne n} \dfrac{|\langle m|V|n\rangle|^2}{E_n - E_m}$.

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
$\pm\sqrt{\Delta^2+\lambda^2}$.

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

---

## Part 3. States as density operators

### 3.1 [BSc] How do I build a density matrix from a classical ensemble and compute its purity $\mathrm{Tr}[\rho^2]$?

A classical mixture of pure states is $\rho = \sum_i p_i\,|\psi_i\rangle\langle\psi_i|$; its purity
$\mathrm{Tr}[\rho^2]$ runs from $1/d$ (maximally mixed) to $1$ (pure).

**WL** : sum the weighted projectors, each $|\psi_i\rangle\langle\psi_i| = $ `KroneckerProduct[psi, Conjugate[psi]]`.

```wl
rho = 3/4 KroneckerProduct[{1, 0}, Conjugate[{1, 0}]] +
   1/4 KroneckerProduct[{1, 1}/Sqrt[2], Conjugate[{1, 1}/Sqrt[2]]]
```

The purity is $\mathrm{Tr}[\rho^2]$.

```wl
Tr[rho . rho]
```

**QF** : a classical mixture is the sum of the component density matrices (NOT a sum of the state
vectors, which would build a coherent superposition); wrap that matrix as a state object.

```wl
rhoQF = QuantumState[3/4 QuantumState[{1, 0}]["DensityMatrix"] +
    1/4 QuantumState[{1, 1}/Sqrt[2]]["DensityMatrix"]]
```

`"Purity"` reads $\mathrm{Tr}[\rho^2]$ off the object.

```wl
rhoQF["Purity"]
```

Both give $\mathrm{Tr}[\rho^2] = 13/16 < 1$: the state is mixed. Summing the *state vectors* instead
would have built a coherent superposition (a pure state, of unit purity), a physically different
object from this classical mixture.

### 3.2 [BSc] How do I test whether a state is pure or mixed?

A state is pure iff any one of three equivalent conditions holds: $\mathrm{Tr}[\rho^2]=1$, $\rho$ has
rank one, or $\rho$ is idempotent ($\rho^2=\rho$, a projector).

**WL** : purity equals one only for a pure state.

```wl
Tr[rho . rho] == 1
```

A pure state's density matrix has rank one.

```wl
MatrixRank[rho] == 1
```

A pure state is idempotent, a projector $\rho^2=\rho$.

```wl
rho . rho == rho
```

**QF** : the object answers directly with `"PureStateQ"`.

```wl
rhoQF["PureStateQ"]
```

Contrast a genuine pure state, which returns `True`.

```wl
QuantumState[{1, I}/Sqrt[2]]["PureStateQ"]
```

The mixed state fails all three criteria (each cell returns `False`, the rank is $2$); the pure state
passes. The three tests are equivalent characterizations of the same rank-one property.

### 3.3 [BSc] How do I build the maximally mixed state $I/d$ and read its purity?

The maximally mixed state $I/d$ is the unique state of maximal entropy; its purity is the floor $1/d$.

**WL** : the normalized identity, and its purity.

```wl
With[{d = 3}, Tr[(IdentityMatrix[d]/d) . (IdentityMatrix[d]/d)]]
```

**QF** : the same state as an object, shown as its summary box.

```wl
QuantumState[IdentityMatrix[3]/3]
```

Its `"Purity"` is the floor $1/d$.

```wl
QuantumState[IdentityMatrix[3]/3]["Purity"]
```

Both give $1/3 = 1/d$, the minimum purity in dimension $d=3$: no state is less pure than $I/d$.

### 3.4 [BSc] How do I compute an expectation value and a variance of an observable on a general (mixed) state as $\mathrm{Tr}[\hat A\rho]$?

For any state, $\langle A\rangle = \mathrm{Tr}[A\rho]$ and $\mathrm{Var} = \mathrm{Tr}[A^2\rho] - \mathrm{Tr}[A\rho]^2$. Writing a general qubit state in Bloch form $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ makes the result transparent.

**WL** : build the general Bloch state and bind it ($\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

```wl
rhoB = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}]);
```

The mean of $\sigma_z$ is $\mathrm{Tr}[\sigma_z\rho]$.

```wl
Tr[PauliMatrix[3] . rhoB]
```

The variance reuses $\sigma_z^2 = I$.

```wl
Tr[PauliMatrix[3] . PauliMatrix[3] . rhoB] - Tr[PauliMatrix[3] . rhoB]^2
```

**QF** : the same general state as an object, shown as its summary box.

```wl
rhoBQF = QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])]
```

The mean, as a measurement-operator expectation on the state. The `"Mean"` runs through Born
probabilities $(1\pm r_z)/2$, so it is simplified under the physicality assumption $|r_z|\le 1$.

```wl
Simplify[QuantumMeasurementOperator[QuantumOperator["Z"]][rhoBQF]["Mean"], -1 <= rz <= 1]
```

The variance, from $\langle\sigma_z^2\rangle - \langle\sigma_z\rangle^2$, under the same assumption.

```wl
Simplify[
 QuantumMeasurementOperator[QuantumOperator["Z"] @ QuantumOperator["Z"]][rhoBQF]["Mean"] -
  QuantumMeasurementOperator[QuantumOperator["Z"]][rhoBQF]["Mean"]^2, -1 <= rz <= 1]
```

Both give $\langle\sigma_z\rangle = r_z$ (the expectation of a Pauli is its Bloch component) and
$\mathrm{Var} = 1 - r_z^2$: the spin is sharp ($\mathrm{Var}=0$) only at the poles $r_z=\pm1$. (The
WL trace $\mathrm{Tr}[\sigma_z\rho]=r_z$ holds for any $\rho$; the QF measurement form needs $\rho$ to
be a physical state, which is the $|r_z|\le1$ assumption.)

### 3.5 [BSc] How do I count the real parameters of a $d$-dimensional mixed state (density matrix)?

A $d\times d$ density matrix is Hermitian ($d^2$ real parameters) with unit trace (one constraint),
leaving $d^2-1$. As in 1.5, the rigorous count is a Jacobian rank.

**WL** : the general Hermitian unit-trace qubit $\{\{a, b+ic\},\{b-ic, 1-a\}\}$ has three real
parameters; its real-coordinate embedding has Jacobian rank $3 = 2^2-1$.

```wl
With[{rho = {{a, b + I c}, {b - I c, 1 - a}}},
 MatrixRank[D[Flatten[ReIm[Flatten[rho]]], {{a, b, c}}]]]
```

The general count for dimension $d$ is $d^2-1$.

```wl
Table[d^2 - 1, {d, 2, 4}]
```

**QF** : drive the general count from the object's own dimension.

```wl
Table[QuantumState["0", d]["Dimension"]^2 - 1, {d, 2, 4}]
```

The qubit rank is $3 = 2^2-1$; the general count is $\{3, 8, 15\}$ for $d=2,3,4$. (The pure-state
count $2d-2$ of 1.5 is smaller: a pure state is a measure-zero boundary of the full state space.)

### 3.6 [BSc] How do I map a qubit density matrix to its Bloch vector and read positivity and the pure-state boundary?

For a qubit, $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ with $r_j = \mathrm{Tr}[\sigma_j\rho]$;
positivity forces $|\vec r|\le 1$, and $|\vec r|=1$ exactly on the pure states.

**WL** : recover the Bloch vector from the general $\rho$ of 3.4 as the Pauli expectations.

```wl
Table[Tr[PauliMatrix[j] . rhoB], {j, 3}]
```

Purity and radius are tied by $\mathrm{Tr}[\rho^2] = \tfrac12(1+|\vec r|^2)$, so $|\vec r|=1$ iff the
state is pure.

```wl
FullSimplify[Tr[rhoB . rhoB]]
```

**QF** : the object computes the Bloch vector directly.

```wl
QuantumState[rhoB]["BlochVector"]
```

Both recover $\vec r = (r_x, r_y, r_z)$, and the purity is $\tfrac12(1 + r_x^2+r_y^2+r_z^2)$:
positivity ($\rho\succeq0$) is exactly $|\vec r|\le1$ (the Bloch ball), and the pure states are its
surface $|\vec r|=1$.

### 3.7 [MSc] How do I show that the state space is convex, with pure states as its extreme points?

The density operators form a convex set: any mixture $p\rho_1+(1-p)\rho_2$ is again a state. The pure
states are its extreme points, which for a qubit is the surface of the Bloch ball.

**WL** : every pure qubit sits on the sphere ($|\vec r|=1$), so pure states are extreme.

```wl
FullSimplify[
 Norm[Table[Tr[PauliMatrix[j] . KroneckerProduct[#, Conjugate[#]]], {j, 3}]] & @ {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]},
 {\[Theta], \[Phi]} \[Element] Reals]
```

A nontrivial mixture of two distinct pure states (here $|0\rangle$ and $|R\rangle=\{1,i\}/\sqrt2$, the
same pair used by the QF cells below) lands strictly inside (purity $<1$).

```wl
With[{p = 1/3},
 Tr[#.#] &[p KroneckerProduct[{1, 0}, Conjugate[{1, 0}]] + (1 - p) KroneckerProduct[{1, I}/Sqrt[2], Conjugate[{1, I}/Sqrt[2]]]]]
```

**QF** : build the same mixture as an object.

```wl
mix = QuantumState[1/3 QuantumState[{1, 0}]["DensityMatrix"] + 2/3 QuantumState[{1, I}/Sqrt[2]]["DensityMatrix"]]
```

Its purity is below one.

```wl
mix["Purity"]
```

Its Bloch radius is below one, so it lies in the interior of the ball.

```wl
Norm[mix["BlochVector"]]
```

The Bloch map is affine: the mixture's Bloch vector is the convex combination of the components'.

```wl
mix["BlochVector"] == 1/3 QuantumState[{1, 0}]["BlochVector"] + 2/3 QuantumState[{1, I}/Sqrt[2]]["BlochVector"]
```

Pure states have radius $1$ (extreme points, not expressible as a nontrivial mixture); every mixture
has radius and purity below $1$ and lies in the interior; and because the Bloch map is affine, mixing
states maps to mixing their Bloch vectors. That is the convex geometry of the Bloch ball.

### 3.8 [BSc] How do I compute the von Neumann entropy of a density matrix?

The von Neumann entropy $S(\rho) = -\mathrm{Tr}[\rho\log\rho] = -\sum_k \lambda_k\log\lambda_k$ is the
Shannon entropy of the eigenvalue spectrum: $0$ for a pure state, $\log_2 d$ for the maximally mixed.

**WL** : the Shannon entropy of the eigenvalues, in bits.

```wl
-Total[# Log2[#] & /@ Eigenvalues[{{3/4, 0}, {0, 1/4}}]]
```

**QF** : the object reports it directly (as a `Quantity` in bits).

```wl
QuantumState[{{3/4, 0}, {0, 1/4}}]["VonNeumannEntropy"]
```

Both give $-\tfrac34\log_2\tfrac34 - \tfrac14\log_2\tfrac14 \approx 0.811$ bits: the entropy is the
classical uncertainty of the eigenvalue distribution, vanishing on a pure state and reaching
$\log_2 d$ on the maximally mixed state.
