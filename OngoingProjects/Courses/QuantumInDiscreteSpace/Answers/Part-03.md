---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 3)

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
leaving $d^2-1$. As in Part 1, 1.5, the rigorous count is a Jacobian rank.

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
