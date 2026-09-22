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
a single ratio (the numerator $\alpha\,\alpha^*$ is $|\alpha|^2$).

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
global phase another, so the naive count is $2d-2$. But that bookkeeping *assumes* each constraint
removes exactly one independent direction. To **measure** the count instead of assuming it, note what
"degrees of freedom" means geometrically: the states sweep out a surface, and the number of free real
parameters is the *dimension* of that surface. The dimension of a parametrized surface is the generic
**rank of its Jacobian**, the matrix of derivatives $\partial|\psi\rangle/\partial(\text{parameter})$:
the rank is how many linearly independent directions the state actually moves in as the parameters
vary, and a redundant parameter (whose motion duplicates another's) does not raise it. That is why the
Jacobian rank is the rigorous count, and why it exposes the global phase as a genuine redundancy.

**WL** : embed the qubit family (the two Bloch angles plus the global phase $g$) into $\mathbb{R}^4$ as
real and imaginary parts, and count its independent tangent directions. With the phase kept, the rank
is $3$: the normalized states fill a 3-sphere.

```wl
With[{psi = Exp[I g] {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 MatrixRank[D[Flatten[ReIm[ComplexExpand[psi]]], {{\[Theta], \[Phi], g}}]]]
```

Drop the phase $g$ and the rank falls to $2$, the physical state space. The drop of exactly one
direction is the global phase, a single redundancy.

```wl
With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
 MatrixRank[D[Flatten[ReIm[ComplexExpand[psi]]], {{\[Theta], \[Phi]}}]]]
```

For a general qudit, drop the hand-built chart and differentiate the one object that *already* has
both redundancies removed: the density matrix $\rho = |\psi\rangle\langle\psi| / \langle\psi|\psi\rangle$.
Build it from a raw, unnormalized vector $z \in \mathbb{C}^d$ and take the $2d$ real and imaginary
components of $z$ as the parameters. Because $\rho$ is unchanged by $z \to \lambda z$ for any complex
$\lambda$, it discards the normalization *and* the global phase at a single stroke, so its Jacobian
rank is the physical parameter count for any $d$, measured rather than assumed.

**WL** : form the projector $z z^\dagger / (z^\dagger z)$ by hand and rank its Jacobian for $d = 2, 3, 4$.

```wl
Table[
 With[{z = Array[a, d] + I Array[b, d]},
  MatrixRank[D[
    Flatten[ReIm[ComplexExpand[
      Flatten[Outer[Times, z, Conjugate[z]]/(Conjugate[z] . z)]]]],
    {Join[Array[a, d], Array[b, d]]}]]],
 {d, 2, 4}]
```

The ranks $\{2, 4, 6\}$ are exactly $2d-2$: each added Hilbert-space dimension contributes two new
real parameters.

**QF** : a `QuantumState` carries the same density matrix, so take it from the object instead of
building it. QF leaves a *symbolic* input unnormalized (its `"DensityMatrix"` is the bare
$z z^\dagger$, of trace $\langle\psi|\psi\rangle$), so divide by the trace and rank the identical map.

```wl
Table[
 With[{dm = QuantumState[Array[a, d] + I Array[b, d]]["DensityMatrix"]},
  MatrixRank[D[
    Flatten[ReIm[ComplexExpand[Flatten[dm/Tr[dm]]]]],
    {Join[Array[a, d], Array[b, d]]}]]],
 {d, 2, 4}]
```

So a qubit has exactly two real parameters, the two Bloch angles: the rank falls from $3$ (the
normalized vectors, a 3-sphere in $\mathbb{R}^4$) to $2$ (the physical states) precisely when the
global phase is removed. For any qudit the gauge freedom $z \to \lambda z$ is always two real
dimensions, one scale and one phase, so the density-matrix rank is $2d-2$ for every $d$. The measured
$\{2,4,6\}$ confirm it, extending to $\{2,4,6,8,10,\dots\}$ for $d = 2,3,4,5,6,\dots$. The count is now
measured at each step, not asserted.

### 1.6 [BSc] How do I compute an inner product $\langle\phi|\psi\rangle$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?

The overlap of two states, its absolute value, and its absolute square give the amplitude and the
probability to find one state in the other.

**WL** : the inner product is the *conjugated* dot product; the conjugation acts on the bra
$\langle\phi|$, which is why $\alpha_1,\beta_1$ appear conjugated. The cell returns the overlap, its
modulus, and its transition probability as a labeled Association, one inner product read three ways.

```wl
With[{phi = {\[Alpha]1, \[Beta]1}, psi = {\[Alpha]2, \[Beta]2}},
 With[{ip = Conjugate[phi] . psi},
  AssociationThread[
   {"\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]",
    "|\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]|",
    "|\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]\!\(\*SuperscriptBox[\(|\), \(2\)]\)"},
   {ip, Abs[ip], Abs[ip]^2}]]]
```

**QF** : build the bra with `"Dagger"`, compose with the ket, and read the same three quantities off
the scalar, under the same labels.

```wl
With[{bra = QuantumState[{\[Alpha]1, \[Beta]1}]["Dagger"], ket = QuantumState[{\[Alpha]2, \[Beta]2}]},
 With[{ip = (bra @ ket)["Scalar"]},
  AssociationThread[
   {"\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]",
    "|\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]|",
    "|\[LeftAngleBracket]\[Phi]|\[Psi]\[RightAngleBracket]\!\(\*SuperscriptBox[\(|\), \(2\)]\)"},
   {ip, Abs[ip], Abs[ip]^2}]]]
```

Both give $\langle\phi|\psi\rangle = \alpha_2\,\alpha_1^* + \beta_2\,\beta_1^*$, with its modulus and
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

Evaluating the definition returns $\{2\,\mathrm{Re}\,\alpha^*\beta,\ 2\,\mathrm{Im}\,\alpha^*\beta,\ |\alpha|^2-|\beta|^2\}$ for a general state and $(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$ for the parametrized one, with spherical coordinates $(1,\theta,\phi)$. So the component formula is a *result* of evaluating the Pauli expectation values, not the definition itself, and the parametrization angles are exactly the polar and azimuthal angles on the Bloch sphere.

### 1.9a [BSc] How do I prepare a qudit in the computational basis state $|0\rangle$?

A qudit lives in $\mathbb{C}^d$; its simplest preparation is the computational basis state $|0\rangle$,
the unit vector along the first axis.

**WL** : $|0\rangle$ is a unit vector in $\mathbb{C}^d$.

```wl
With[{d = 3}, UnitVector[d, 1]]
```

**QF** : name the basis state with its dimension.

```wl
With[{d = 3}, QuantumState["0", d]["StateVector"] // Normal]
```

Both give $|0\rangle = (1,0,0)$ in $d=3$; written in $d$, the preparation works for any dimension.

### 1.9b [BSc] How do I prepare a uniform superposition $\tfrac1{\sqrt d}\sum_k|k\rangle$ in $\mathbb{C}^d$?

The equal-weight superposition $\tfrac1{\sqrt d}\sum_k |k\rangle$ puts the same amplitude on every basis
state; it is the flat starting state of many algorithms.

**WL** : a constant vector handed to `Normalize` (no hardcoded $1/\sqrt d$).

```wl
With[{d = 3}, Normalize @ ConstantArray[1, d]]
```

**QF** : use the built-in `"UniformSuperposition"` state rather than building it by hand.

```wl
With[{d = 3}, QuantumState["UniformSuperposition", d]["StateVector"] // Normal]
```

Both give the uniform qutrit $\tfrac1{\sqrt3}(1,1,1)$; written in $d$, it works for any dimension.

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
$\langle\sigma_z\rangle = r_z$, $\mathrm{Var} = 1 - r_z^2$), is 3.4.

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
generic combination $A + cB$ just once. The basis that $A + \pi B$ returns is the Bell basis, *derived*
rather than assumed; $\pi$ only serves to break the degeneracy, and any value that does so yields the
same simultaneous eigenbasis.

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

---

## Part 3. States as density operators

### 3.1a [BSc] How do I build a density matrix from a classical ensemble?

A classical mixture of pure states, each $|\psi_i\rangle$ prepared with probability $p_i$, is the density
matrix $\rho = \sum_i p_i\,|\psi_i\rangle\langle\psi_i|$.

**WL** : sum the weighted projectors, each $|\psi_i\rangle\langle\psi_i| = $ `KroneckerProduct[psi, Conjugate[psi]]`.

```wl
rho = 3/4 KroneckerProduct[{1, 0}, Conjugate[{1, 0}]] +
   1/4 KroneckerProduct[{1, 1}/Sqrt[2], Conjugate[{1, 1}/Sqrt[2]]]
```

**QF** : a classical mixture is the sum of the component density matrices (NOT a sum of the state
vectors, which would build a coherent superposition); wrap that matrix as a state object.

```wl
rhoQF = QuantumState[3/4 QuantumState[{1, 0}]["DensityMatrix"] +
    1/4 QuantumState[{1, 1}/Sqrt[2]]["DensityMatrix"]]
```

Both build the same $\rho$. Summing the *state vectors* instead would have built a coherent superposition,
a pure state, a physically different object from this classical mixture.

### 3.1b [BSc] How do I compute the purity $\mathrm{Tr}[\rho^2]$ of a density matrix?

The purity $\mathrm{Tr}[\rho^2]$ measures how mixed a state is: it runs from $1/d$ (maximally mixed) up to
$1$ (pure), the value $1$ reached only when $\rho$ is a single projector.

**WL** : form the ensemble density matrix $\tfrac34|0\rangle\langle0| + \tfrac14|+\rangle\langle +|$ and
read its purity $\mathrm{Tr}[\rho^2]$.

```wl
With[{rho = 3/4 KroneckerProduct[{1, 0}, Conjugate[{1, 0}]] + 1/4 KroneckerProduct[{1, 1}/Sqrt[2], Conjugate[{1, 1}/Sqrt[2]]]},
 Tr[rho . rho]]
```

**QF** : `"Purity"` reads $\mathrm{Tr}[\rho^2]$ off the state object.

```wl
QuantumState[3/4 QuantumState[{1, 0}]["DensityMatrix"] + 1/4 QuantumState[{1, 1}/Sqrt[2]]["DensityMatrix"]]["Purity"]
```

Both give $\mathrm{Tr}[\rho^2] = 13/16 < 1$: the ensemble state is mixed, below the pure-state value $1$
and above the qubit floor $1/2$.

### 3.2 [BSc] How do I test whether a state is pure or mixed?

A state is pure iff any one of three equivalent conditions holds: $\mathrm{Tr}[\rho^2]=1$, $\rho$ has
rank one, or $\rho$ is idempotent ($\rho^2=\rho$, a projector).

**WL** : build a mixed state to test, the classical ensemble $\tfrac34|0\rangle\langle0| + \tfrac14|+\rangle\langle +|$.

```wl
rho = 3/4 KroneckerProduct[{1, 0}, Conjugate[{1, 0}]] + 1/4 KroneckerProduct[{1, 1}/Sqrt[2], Conjugate[{1, 1}/Sqrt[2]]];
```

Purity equals one only for a pure state.

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

**QF** : build the same mixture as a state object; it answers directly with `"PureStateQ"`.

```wl
rhoQF = QuantumState[3/4 QuantumState[{1, 0}]["DensityMatrix"] + 1/4 QuantumState[{1, 1}/Sqrt[2]]["DensityMatrix"]];
```

The mixed state is not pure.

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

**WL** : build the general qubit state in Bloch form ($\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

```wl
rhoB = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}]);
```

Recover the Bloch vector as the Pauli expectations $r_j = \mathrm{Tr}[\sigma_j\rho]$.

```wl
Table[Tr[PauliMatrix[j] . rhoB], {j, 3}]
```

Purity and radius are tied by $\mathrm{Tr}[\rho^2] = \tfrac12(1+|\vec r|^2)$, so $|\vec r|=1$ iff the
state is pure.

```wl
FullSimplify[Tr[rhoB . rhoB]]
```

**QF** : the object computes the Bloch vector directly from the same $\rho$.

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

### 3.9 [BSc] How do I test whether a matrix is a valid density operator?

A density operator must satisfy three conditions: Hermiticity ($\rho=\rho^\dagger$), unit trace
($\mathrm{Tr}\,\rho=1$), and positive semidefiniteness ($\rho\succeq0$, no negative eigenvalue). Together
they are exactly what makes $\rho$ a legitimate state.

**WL** : check the three conditions on a candidate matrix.

```wl
With[{rho = {{3/4, 0}, {0, 1/4}}},
 {HermitianMatrixQ[rho], Tr[rho] == 1, PositiveSemidefiniteMatrixQ[rho]}]
```

A Hermitian, unit-trace matrix that fails positivity is not a state.

```wl
With[{rho = {{3/4, 1}, {1, 1/4}}},
 {HermitianMatrixQ[rho], Tr[rho] == 1, PositiveSemidefiniteMatrixQ[rho]}]
```

**QF** : the object folds all three conditions into the single property `"PhysicalQ"`.

```wl
QuantumState[{{3/4, 0}, {0, 1/4}}]["PhysicalQ"]
```

The invalid matrix fails it.

```wl
QuantumState[{{3/4, 1}, {1, 1/4}}]["PhysicalQ"]
```

Both agree: the first matrix is a valid density operator (the three conditions all hold and `"PhysicalQ"`
is `True`); the second is Hermitian and unit-trace but not positive (one eigenvalue is negative), so it is
not a physical state.

---

## Part 4. Composite systems: tensor product and partial trace

### 4.1 [BSc] How do I form the tensor product of states and of operators?

Two subsystems combine by the tensor product: states as $|\psi\rangle\otimes|\phi\rangle$, operators as
$A\otimes B$, and the two are compatible, $(A\otimes B)(|\psi\rangle\otimes|\phi\rangle) =
(A|\psi\rangle)\otimes(B|\phi\rangle)$. In the computational basis the tensor product is the Kronecker
product, and the composite lives in the $d_1 d_2$-dimensional product space.

**WL** : the tensor product of two amplitude vectors is their (flattened) Kronecker product.

```wl
Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

The tensor product of two operators is the Kronecker product of their matrices, here $X\otimes Z$.

```wl
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]
```

**QF** : `QuantumTensorProduct` combines two states into one two-qubit state object.

```wl
QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]
```

Its `"StateVector"` is that same Kronecker product.

```wl
QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]["StateVector"] // Normal
```

It joins operators the same way; here $X\otimes Z$ kept as a `QuantumOperator`.

```wl
QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]
```

Both give the amplitude vector $\{ac,\,ad,\,bc,\,bd\}$ and the operator $X\otimes Z$: the composite state
factorizes into its parts, and the QF object additionally knows it is a two-qubit system (its factors,
order, and basis), a matrix being just one representation read from it.

### 4.2 [BSc] How do I take a partial trace of a two-party state and obtain a reduced density matrix?

The state of subsystem $A$ alone is the reduced density matrix $\rho_A = \mathrm{Tr}_B[\rho_{AB}]$, formed
by tracing out $B$. For a product state $\rho_A$ stays pure; for an *entangled* state it comes out mixed,
which is the operational signature of entanglement.

**WL** : trace qubit $B$ out of the Bell state $\rho = |\Phi^+\rangle\langle\Phi^+|$ by reshaping the
$4\times4$ matrix to a $2\times2\times2\times2$ tensor and contracting the two $B$ indices.

```wl
With[{bell = {1, 0, 0, 1}/Sqrt[2]},
 TensorContract[ArrayReshape[KroneckerProduct[bell, Conjugate[bell]], {2, 2, 2, 2}], {{2, 4}}]]
```

**QF** : `QuantumPartialTrace` traces out the listed qubit and returns the reduced state as an object.

```wl
QuantumPartialTrace[QuantumState["PhiPlus"], {2}]
```

Its `"Purity"` is below one, so the reduction of the entangled pair is mixed.

```wl
QuantumPartialTrace[QuantumState["PhiPlus"], {2}]["Purity"]
```

A *product* input, by contrast, reduces to a pure state (purity $1$).

```wl
QuantumPartialTrace[QuantumTensorProduct[QuantumState["0"], QuantumState["1"]], {2}]["Purity"]
```

Both give $\rho_A = I/2$ with purity $\tfrac12$: tracing a maximally entangled pair down to one qubit
yields the maximally mixed state, whereas a product state stays pure. The mixedness of the reduced state
is exactly the entanglement discarded along with $B$.

---

## Part 5. Projective measurement

### 5.1 [BSc] How do I perform a projective measurement and read outcomes, probabilities, and the mean?

A projective measurement of an observable $A$ returns one of $A$'s eigenvalues $a$ (the outcome) with
probability $\mathrm{Tr}[P_a\rho] = \langle a|\rho|a\rangle$ (the Born rule, $|\langle a|\psi\rangle|^2$ for a
pure state), and the mean is $\langle A\rangle = \mathrm{Tr}[A\rho]$. Take $A = Z$ on a general qubit
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$, $|\vec r|\le1$ ($\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

**WL** : the possible outcomes are the eigenvalues of the observable.

```wl
Eigenvalues[PauliMatrix[3]]
```

Each outcome's probability is the corresponding eigenvector sandwiched around $\rho$, $\langle a|\rho|a\rangle$.

```wl
With[{rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])},
 FullSimplify[Conjugate[#] . rho . # & /@ Eigenvectors[PauliMatrix[3]], -1 <= rz <= 1]]
```

The mean is $\mathrm{Tr}[Z\rho]$.

```wl
With[{rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])},
 FullSimplify[Tr[PauliMatrix[3] . rho], -1 <= rz <= 1]]
```

**QF** : applying the measurement operator to the state returns a `QuantumMeasurement`, whose summary box
lists the outcomes and their probabilities.

```wl
QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])]]
```

Its `"Mean"` is $\langle Z\rangle$.

```wl
FullSimplify[QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])]]["Mean"], -1 <= rz <= 1]
```

Its `"ProbabilitiesList"` gives the outcome probabilities.

```wl
FullSimplify[QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])]]["ProbabilitiesList"], -1 <= rz <= 1]
```

Both give outcomes $\pm1$ with Born probabilities $\tfrac{1+r_z}2$ (for $+1$) and $\tfrac{1-r_z}2$ (for
$-1$) and mean $\langle Z\rangle = r_z$: a $Z$ measurement reads only the populations (the $r_z$ component),
so the coherences $r_x, r_y$ drop out, and the mean is the difference of the two outcome probabilities. For a
pure state ($|\vec r| = 1$) these are the familiar $\cos^2\tfrac\theta2, \sin^2\tfrac\theta2, \cos\theta$.

### 5.2 [BSc] How do I measure in a non-computational basis and get the post-measurement (collapsed) states?

Measuring an observable whose eigenbasis is not the computational one, here $X$ with eigenstates
$|\pm\rangle = (|0\rangle\pm|1\rangle)/\sqrt2$, collapses a general pure qubit
$|\psi\rangle = \{\cos\tfrac\theta2, e^{i\varphi}\sin\tfrac\theta2\}$ onto one of those eigenstates with
probability $|\langle\pm|\psi\rangle|^2$.

**WL** : the collapse targets are the normalized eigenvectors of $X$.

```wl
Normalize /@ Eigenvectors[PauliMatrix[1]]
```

Their probabilities are the squared overlaps with the state.

```wl
FullSimplify[ComplexExpand[With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}}, Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[1]])]], {\[Theta], \[Phi]} \[Element] Reals]
```

**QF** : the $X$-basis measurement object.

```wl
QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]]
```

Its `"States"` are the post-measurement states, each collapsed onto an $X$ eigenstate and weighted by its
measurement amplitude.

```wl
QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}]]["States"]
```

Both give the collapse basis $|\pm\rangle$ with probabilities $\tfrac12(1\pm\sin\theta\cos\varphi)$: unlike
a $Z$ measurement, the $X$ outcome depends on the relative phase $\varphi$, so a real state would hide it. A
measurement in the $X$ basis leaves the qubit in an $X$ eigenstate. The post-measurement objects are those
eigenstates, each scaled by its outcome amplitude (its squared norm is the probability).

### 5.3 [BSc] How do I simulate finite-shot statistics and watch the empirical frequencies approach the Born rule?

Real measurements return a finite sample of outcomes whose frequencies converge to the Born probabilities
as the shot count grows. This is the one place the two answers do not agree cell-for-cell: each draws its
own random sample, and what they share is the limit. Measure $Z$ on the state
$\{\tfrac{\sqrt3}2, \tfrac12\}$, whose Born probabilities are the unequal $\tfrac34, \tfrac14$.

**WL** : draw $2000$ outcomes with the Born weights and tabulate their frequencies.

```wl
SeedRandom[1]; N[Counts[RandomChoice[{3/4, 1/4} -> {1, -1}, 2000]]/2000]
```

**QF** : `"SimulatedCounts"` draws the same-size sample from the measurement and returns the raw counts.

```wl
SeedRandom[1]; QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[{Sqrt[3]/2, 1/2}]]["SimulatedCounts", 2000]
```

Both cluster near the Born prediction $\{1500, 500\}$ (frequencies $\approx \{0.75, 0.25\}$): the two samples
differ in their individual draws but each approaches $\tfrac34, \tfrac14$, and the gap shrinks like
$1/\sqrt N$. The unequal weights make the convergence discriminating, the frequencies track two distinct
probabilities rather than a symmetric coin flip. The Born rule is the large-sample limit of the measured
frequencies.

### 5.4 [BSc] How do I apply a non-selective projective measurement (the Lüders channel)?

A non-selective measurement records that a measurement happened but not which outcome, replacing $\rho$ by
$\sum_k P_k\,\rho\,P_k$ (the Lüders channel). For a $Z$ measurement this erases the off-diagonal coherences
and leaves the diagonal: pure dephasing. Apply it to a general qubit state in Bloch form
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ (pure or mixed; $\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

**WL** : sandwich $\rho$ between the two $Z$ projectors and sum.

```wl
With[{rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])},
 Sum[p . rho . p, {p, {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}}]]
```

**QF** : the measurement's `"PostMeasurementState"` is that non-selective mixture, as a state object.

```wl
QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])]]["PostMeasurementState"]
```

Both give $\mathrm{diag}\big(\tfrac{1+r_z}2, \tfrac{1-r_z}2\big)$: the coherences $r_x, r_y$ are erased while
the populations $\tfrac{1\pm r_z}2$ are untouched, for any state, pure or mixed. This is decoherence in the
measured basis, the same map a dephasing environment applies, and it is what separates a measurement that is
*read* from one that is merely coupled.

---

## Part 6. Uncertainty and incompatibility

### 6.1 [MSc] How do I verify the Robertson-Schrodinger uncertainty relation for a given state?

For observables $A, B$ on *any* state $\rho$ (pure or mixed), with standard deviations
$\sigma_A = \sqrt{\langle A^2\rangle - \langle A\rangle^2}$ and $\sigma_B$ likewise, the Robertson-Schrodinger
relation bounds the product of variances,
$$
\sigma_A^2\,\sigma_B^2 \;\ge\; \Big(\tfrac12\langle\{A,B\}\rangle - \langle A\rangle\langle B\rangle\Big)^2
  + \Big(\tfrac1{2i}\langle[A,B]\rangle\Big)^2 ,
$$
the first term a covariance, the second the Robertson (commutator) piece, with $\langle\cdot\rangle =
\mathrm{Tr}[\cdot\,\rho]$. Take $A=X, B=Z$ on a general qubit $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$,
$|\vec r|\le1$; the telling quantity is the slack, how far the state sits above the bound.

**WL** : with $\langle\cdot\rangle$ the trace against $\rho$, compute the left side minus the right side.

```wl
With[{rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}]), x = PauliMatrix[1], z = PauliMatrix[3]},
 With[{ev = Tr[# . rho] &},
  FullSimplify[(ev[x . x] - ev[x]^2) (ev[z . z] - ev[z]^2) -
    ((ev[(x . z + z . x)/2] - ev[x] ev[z])^2 + (ev[x . z - z . x]/(2 I))^2)]]]
```

**QF** : the same, the general state as a `QuantumState` and its density matrix, the operators as
`QuantumOperator`s and the commutator via `Commutator`.

```wl
With[{qs = QuantumState[1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])], x = QuantumOperator["X"], z = QuantumOperator["Z"]},
 With[{ev = Tr[Normal[#["Matrix"]] . Normal[qs["DensityMatrix"]]] &},
  FullSimplify[(ev[x @ x] - ev[x]^2) (ev[z @ z] - ev[z]^2) -
    ((ev[(x @ z + z @ x)/2] - ev[x] ev[z])^2 + (ev[Commutator[x, z]]/(2 I))^2)]]]
```

Both return $1 - |\vec r|^2 \ge 0$: the relation holds for every state, and the slack is exactly the
mixedness, $1 - |\vec r|^2 = 2(1 - \mathrm{Tr}[\rho^2])$. It vanishes only on the pure states ($|\vec r|=1$),
which therefore saturate Robertson-Schrodinger; a mixed state's extra uncertainty lifts it strictly above
the bound. The commutator term alone is the weaker Robertson relation $\sigma_A\sigma_B \ge
\tfrac12|\langle[A,B]\rangle|$; Schrodinger's covariance term is what tightens it.

### 6.2 [MSc] How do I find a minimum-uncertainty state that saturates the Robertson bound?

The Robertson bound $\sigma_A\sigma_B \ge \tfrac12|\langle[A,B]\rangle|$ is saturated exactly where the
*Robertson gap* $g = \sigma_A^2\sigma_B^2 - \big(\tfrac1{2i}\langle[A,B]\rangle\big)^2$ vanishes. So do not
guess a state: compute $g$ over the whole pure-qubit family and solve $g = 0$. Take $A=X, B=Z$.

**WL** : the gap as a function of the Bloch angles.

```wl
gap = With[{psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}, x = PauliMatrix[1], z = PauliMatrix[3]},
  With[{ev = Conjugate[psi] . # . psi &},
   FullSimplify[(ev[x . x] - ev[x]^2) (ev[z . z] - ev[z]^2) - (ev[x . z - z . x]/(2 I))^2, {\[Theta], \[Phi]} \[Element] Reals]]]
```

Solving $g = 0$ finds the minimum-uncertainty states.

```wl
Reduce[gap == 0 && 0 < \[Theta] < Pi && -Pi < \[Phi] <= Pi, {\[Theta], \[Phi]}]
```

One such state sits where the two branches cross, $\theta = \phi = \tfrac\pi2$.

```wl
{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]} /. {\[Theta] -> Pi/2, \[Phi] -> Pi/2} // FullSimplify
```

**QF** : the same gap from the `QuantumState` and `QuantumOperator`s (the solve is identical).

```wl
With[{psi = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}], x = QuantumOperator["X"], z = QuantumOperator["Z"]},
 With[{ev = (psi["Dagger"] @ # @ psi)["Scalar"] &},
  FullSimplify[(ev[x @ x] - ev[x]^2) (ev[z @ z] - ev[z]^2) - (ev[Commutator[x, z]]/(2 I))^2, {\[Theta], \[Phi]} \[Element] Reals]]]
```

The gap is $g = \langle X\rangle^2\langle Z\rangle^2 = \cos^2\theta\,\cos^2\phi\,\sin^2\theta \ge 0$, so it
vanishes exactly where $\langle X\rangle = 0$ or $\langle Z\rangle = 0$: the equator $\theta = \tfrac\pi2$
($\langle Z\rangle = 0$) and the meridians $\phi = \pm\tfrac\pi2$ ($\langle X\rangle = 0$). At their crossing
both vanish, giving the $Y$ eigenstate $\{1, i\}/\sqrt2$, the most symmetric minimum-uncertainty state, which
fell out of the solve rather than being assumed. There the covariance $\langle X\rangle\langle Z\rangle$ is
zero, which is what lets the plain Robertson bound be reached.

### 6.3 [BSc] How do I compute the classical Shannon entropy $H(p) = -\sum_i p_i\log_2 p_i$ of a measurement distribution?

A measurement yields a probability distribution; its Shannon entropy $H(p) = -\sum_i p_i\log_2 p_i$ (in
bits) measures the outcome uncertainty. Take the $Z$ distribution of a state with populations
$\{\tfrac34, \tfrac14\}$.

**WL** : the entropy is the negative weighted sum of $\log_2$ probabilities.

```wl
-Total[# Log2[#] & /@ {3/4, 1/4}]
```

**QF** : a `QuantumMeasurement` reports its `"Entropy"` directly, as a `Quantity` in bits.

```wl
QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[{Sqrt[3]/2, 1/2}]]["Entropy"]
```

Both give $-\tfrac34\log_2\tfrac34 - \tfrac14\log_2\tfrac14 \approx 0.811$ bits: a biased two-outcome
distribution carries less than one bit, reaching exactly one bit only for the uniform $\{\tfrac12,\tfrac12\}$.

### 6.4 [MSc] How do I construct a maximal set of mutually unbiased bases for a qubit?

Two bases are mutually unbiased when every cross-overlap is flat, $|\langle e|f\rangle|^2 = 1/d$: a definite
state in one is maximally uncertain in the other. Rather than name the bases, construct them: start from the
computational basis, parametrize a basis unbiased to it, then solve for the phase that makes two such
families mutually unbiased. A qubit ($d=2$) admits $d+1 = 3$ mutually unbiased bases (MUBs).

**WL** : a basis unbiased to the computational one has equal-modulus amplitudes, $\{1, e^{i\alpha}\}/\sqrt2$
(with its orthogonal partner); its overlaps with $|0\rangle, |1\rangle$ are both $1/2$ for any phase.

```wl
With[{u = {1, Exp[I \[Alpha]]}/Sqrt[2]}, FullSimplify[Abs[Conjugate[#] . u]^2 & /@ IdentityMatrix[2], \[Alpha] \[Element] Reals]]
```

Two such families are mutually unbiased only at a special phase separation; solve for it, fixing the first at
$\alpha = 0$.

```wl
With[{ov = FullSimplify[Abs[Conjugate[{1, 1}/Sqrt[2]] . ({1, Exp[I \[Beta]]}/Sqrt[2])]^2, \[Beta] \[Element] Reals]},
 Reduce[ov == 1/2 && 0 <= \[Beta] < 2 Pi, \[Beta]]]
```

The phase must be $\beta = \tfrac\pi2$, so the maximal set is the computational basis together with the
$\alpha = 0$ and $\alpha = \tfrac\pi2$ families.

```wl
With[{mub = Function[a, {{1, Exp[I a]}, {1, -Exp[I a]}}/Sqrt[2]]},
 Prepend[mub /@ {0, Pi/2}, IdentityMatrix[2]]]
```

**QF** : those three constructed bases are (up to phase) the eigenbases of $Z, X, Y$; read them off the
operators.

```wl
Function[A, Normal /@ QuantumOperator[A]["Eigenvectors"]] /@ {"Z", "X", "Y"}
```

Solving the unbiasedness conditions delivers three bases and no more: the computational ($Z$) basis and the
phase families at $\alpha = 0$ and $\tfrac\pi2$, which are the $X$ and $Y$ eigenbases. Both the count
$d+1 = 3$ and the identity of the bases (the Pauli eigenbases) emerge from the construction rather than being
assumed.

### 6.4b [BSc] How do I verify that a set of bases is mutually unbiased?

Given candidate bases, the test is direct: every cross-overlap between *different* bases must equal
$|\langle e|f\rangle|^2 = 1/d$. Check it on the three qubit MUBs, the eigenbases of $X, Y, Z$.

**WL** : form the three eigenbases and compute all three pairwise overlap matrices.

```wl
With[{bx = Normalize /@ Eigenvectors[PauliMatrix[1]], by = Normalize /@ Eigenvectors[PauliMatrix[2]], bz = Normalize /@ Eigenvectors[PauliMatrix[3]]},
 Outer[Abs[Conjugate[#1] . #2]^2 &, #[[1]], #[[2]], 1] & /@ {{bx, by}, {bx, bz}, {by, bz}}]
```

**QF** : the same from the operators' `"Eigenvectors"`, over all pairs.

```wl
With[{b = Function[A, Normal /@ QuantumOperator[A]["Eigenvectors"]] /@ {"X", "Y", "Z"}},
 Outer[Abs[Conjugate[#1] . #2]^2 &, #[[1]], #[[2]], 1] & /@ Subsets[b, {2}]]
```

Every entry of all three pairwise matrices is $\tfrac12 = 1/d$: the $X, Y, Z$ eigenbases are pairwise mutually
unbiased, confirming they are the three MUBs of a qubit. A single entry different from $\tfrac12$ would
disqualify the set, which is what makes this a genuine test rather than a restatement.

### 6.5 [MSc] How do I verify the Maassen-Uffink entropic uncertainty relation $H(A)+H(B)\ge-\log_2 c^2$ for two observables?

Let $H(A)$ denote the Shannon entropy (6.3) of the outcome distribution of measuring $A$. The Maassen-Uffink
relation bounds the *summed* entropies of two measurements, $H(A) + H(B) \ge -\log_2 c^2$, where
$c = \max_{a,b}|\langle a|b\rangle|$ is the largest overlap between an $A$ eigenstate and a $B$ eigenstate.
Take $A = X, B = Z$ and a complex off-axis state $\{\cos\tfrac\pi6, e^{2\pi i/5}\sin\tfrac\pi6\}$, whose
relative phase keeps both measurements genuinely uncertain; the overlap constant $c$, and hence the bound,
is computed below rather than assumed.

**WL** : first compute the overlap constant $c = \max_{a,b}|\langle a|b\rangle|$ from the two eigenbases,
which sets the bound $-\log_2 c^2$.

```wl
With[{eX = Normalize /@ Eigenvectors[PauliMatrix[1]], eZ = Normalize /@ Eigenvectors[PauliMatrix[3]]},
 With[{c = Max[Flatten[Outer[Abs[Conjugate[#1] . #2] &, eX, eZ, 1]]]}, {c, -Log2[c^2]}]]
```

This returns $\{1/\sqrt2, 1\}$: for the mutually unbiased $X, Z$ every cross-overlap is $1/\sqrt2$, so the
bound $-\log_2 c^2$ is exactly $1$ bit. Then sum the Shannon entropies of the $X$ and $Z$ measurement
distributions and compare to that bound.

```wl
With[{psi = {Cos[Pi/6], Exp[I 2 Pi/5] Sin[Pi/6]}, h = -Total[If[# == 0, 0, # Log2[#]] &[#] & /@ #] &},
 {h[Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[1]])] +
    h[Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[3]])], 1} // N]
```

**QF** : sum the two measurements' `"Entropy"` (stripping the bit unit) and compare to $1$.

```wl
With[{psi = QuantumState[{Cos[Pi/6], Exp[I 2 Pi/5] Sin[Pi/6]}]},
 {QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["X"]][psi]["Entropy"]] +
    QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["Z"]][psi]["Entropy"]], 1} // N]
```

Both give $H(X) + H(Z) \approx 1.76 \ge 1$: the summed uncertainty of two complementary measurements cannot
fall below the computed bound, so a state cannot be sharp in both $X$ and $Z$ at once.

### 6.6 [MSc] How do I find a state that saturates the Maassen-Uffink entropic bound (maximal complementarity)?

The Maassen-Uffink bound of 6.5, $H(X) + H(Z) \ge 1$ bit for the complementary $X, Z$ (with $H(A)$ the
Shannon entropy of measuring $A$, from 6.3 and 6.5), is achievable with equality: a $Z$ eigenstate
$|0\rangle$ has $H(Z) = 0$ (certain) and $H(X) = 1$ (maximally uncertain), so the sum equals the bound
exactly.

**WL** : the summed entropy for $|0\rangle$.

```wl
With[{psi = {1, 0}, h = -Total[If[# == 0, 0, # Log2[#]] &[#] & /@ #] &},
 h[Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[1]])] +
  h[Abs[Conjugate[#] . psi]^2 & /@ (Normalize /@ Eigenvectors[PauliMatrix[3]])]]
```

**QF** : the same from the two measurements' `"Entropy"`.

```wl
With[{psi = QuantumState["0"]},
 QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["X"]][psi]["Entropy"]] +
  QuantityMagnitude[QuantumMeasurementOperator[QuantumOperator["Z"]][psi]["Entropy"]]]
```

Both give exactly $1$: knowing the qubit is $|0\rangle$ pins $Z$ completely ($H(Z) = 0$) at the cost of
total ignorance of $X$ ($H(X) = 1$). Mutually unbiased bases realize maximal complementarity, the entropic
bound is met with equality.
