---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 1)

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
rank is the physical parameter count for any $d$.

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
$\{2,4,6\}$ confirm it, extending to $\{2,4,6,8,10,\dots\}$ for $d = 2,3,4,5,6,\dots$.

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
