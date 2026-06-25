# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework

A companion answer key to `Question-List.md`. For each question it gives **two** worked answers:

1. **WL** : native Wolfram Language only (plain vectors and matrices, `PauliMatrix`, `Conjugate`,
   `Normalize`, `Eigenvectors`, `MatrixExp`, and friends). Nothing beyond the built-in language.
2. **QF** : the same task done through the **QuantumFramework** paclet, leaning on its objects and
   property downvalues (`QuantumState[...]["..."]`) rather than rebuilding matrices by hand.

Three habits run through the answers. They stay **symbolic** wherever possible (general amplitudes
$\alpha,\beta$, general angles $\theta,\phi$), so each operation is exercised in full generality
rather than on a lucky special case. They write the computation **explicitly** rather than hiding it
in a helper function, because the code is itself part of the explanation. And they prefer **built-in
features** (a `Normalize`, a named state) over hand-rolled constructions, so nothing is hardcoded:
every cell computes its result. Run the WL answers in a bare kernel; the QF answers need the paclet
loaded once.

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

**QF** : wrap the amplitudes in `QuantumState`; the object carries the basis with it, and
`["StateVector"]` returns the amplitudes back.

```wl
qs = QuantumState[{\[Alpha], \[Beta]}];
qs["StateVector"] // Normal
```

Both hold the amplitudes $\{\alpha,\beta\}$; the QF object additionally knows it is a single qubit in
the computational basis.

### 1.2 [BSc] How do I read the Born-rule probabilities of a state?

The probability of outcome $i$ is $|\langle i|\psi\rangle|^2$, i.e. the squared modulus of the
amplitude (divided by the norm if the state is not yet normalized).

**WL** : square the moduli of the normalized amplitudes.

```wl
Abs[Normalize[{\[Alpha], \[Beta]}]]^2
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

**WL** : `Norm` measures it, `Normalize` fixes it.

```wl
{Norm[{\[Alpha], \[Beta]}], Simplify[Norm @ Normalize[{\[Alpha], \[Beta]}]]}
```

**QF** : the `"Norm"` and `"Normalize"` properties do the same; `"Normalize"` returns a new state.

```wl
{QuantumState[{\[Alpha], \[Beta]}]["Norm"], QuantumState[{\[Alpha], \[Beta]}]["Normalize"]["Norm"]}
```

Both report the raw norm $\sqrt{|\alpha|^2+|\beta|^2}$ and confirm the normalized state has norm $1$.

### 1.4 [BSc] How do I show that a global phase is physically unobservable?

Multiplying $|\psi\rangle$ by $e^{i\varphi}$ changes no measurable prediction, because every
probability depends only on $|c_i|^2$.

**WL** : compare the probability vectors symbolically for all real phases.

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
global phase removes another, leaving $2d-2$.

**WL** : the unit-norm family $\{\cos(\theta/2),\, e^{i\phi}\sin(\theta/2)\}$ shows the qubit needs
only two real angles; the general count is $2d-2$.

```wl
{FullSimplify[
  Norm[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}] == 1,
  {\[Theta], \[Phi]} \[Element] Reals],
 Table[2 d - 2, {d, 2, 5}]}
```

**QF** : drive the same count from the object's own Hilbert-space dimension rather than a literal.

```wl
Table[2 (QuantumState["0", d]["Dimension"] - 1), {d, 2, 5}]
```

The qubit family is normalized for all real angles (so two parameters suffice), and the count is
$\{2,4,6,8\}$ for $d = 2,3,4,5$.

### 1.6 [BSc] How do I compute an inner product $\langle\phi|\psi\rangle$, its modulus, and the transition probability $|\langle\phi|\psi\rangle|^2$?

The overlap of two states and its absolute square give the amplitude and probability to find one in
the other.

**WL** : the inner product is the *conjugated* dot product; the conjugation acts on the bra
$\langle\phi|$, which is why the $\alpha_1,\beta_1$ appear conjugated.

```wl
With[{phi = {\[Alpha]1, \[Beta]1}, psi = {\[Alpha]2, \[Beta]2}},
 With[{ip = Conjugate[phi] . psi}, {ip, Abs[ip], Abs[ip]^2}]]
```

**QF** : build the bra with `"Dagger"`, compose with the ket, and read the scalar.

```wl
With[{bra = QuantumState[{\[Alpha]1, \[Beta]1}]["Dagger"], ket = QuantumState[{\[Alpha]2, \[Beta]2}]},
 With[{ip = (bra @ ket)["Scalar"]}, {ip, Abs[ip], Abs[ip]^2}]]
```

Both give $\langle\phi|\psi\rangle = \alpha_2\bar\alpha_1 + \beta_2\bar\beta_1$, with its modulus and
its square. The conjugation on the bra is what makes this a sesquilinear inner product rather than a
plain bilinear dot product.

### 1.7 [BSc] How do I write a state both as a column vector of amplitudes and as a Dirac superposition $\sum_i c_i\,|i\rangle$?

The amplitude list and the ket-superposition are two views of one object.

**WL** : the amplitude vector is the list itself; the Dirac form is the same numbers against the
computational basis.

```wl
{{\[Alpha], \[Beta]}, {\[Alpha], \[Beta]} == \[Alpha] {1, 0} + \[Beta] {0, 1}}
```

**QF** : `"StateVector"` gives the amplitudes, `"Formula"` renders the $\sum_i c_i|i\rangle$ form.

```wl
qs = QuantumState[{\[Alpha], \[Beta]}];
{qs["StateVector"] // Normal, qs["Formula"]}
```

The amplitude vector and the superposition $\alpha|0\rangle + \beta|1\rangle$ carry identical
information; the `True` confirms the basis expansion, and QF's `"Formula"` typesets the superposition.

### 1.8 [BSc] How do I get the Bloch vector of a qubit, and the Bloch-sphere angles?

The Bloch vector $\vec r = (\langle\sigma_x\rangle, \langle\sigma_y\rangle, \langle\sigma_z\rangle)$
places a pure qubit on the unit sphere, with $(\theta,\phi)$ its spherical angles.

**WL** : the Bloch vector is, *by definition*, the triple of Pauli expectation values
$\langle\psi|\sigma_j|\psi\rangle$. Compute them directly with `PauliMatrix`. For a general
$\{\alpha,\beta\}$ they return the familiar component formula; for the angle-parametrized state they
give the geometric vector, which the built-in then turns into spherical angles.

```wl
{
 FullSimplify @ Table[Conjugate[{\[Alpha], \[Beta]}] . PauliMatrix[j] . {\[Alpha], \[Beta]}, {j, 3}],
 With[{\[Psi] = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}},
  With[{r = FullSimplify[Table[Conjugate[\[Psi]] . PauliMatrix[j] . \[Psi], {j, 3}], {\[Theta], \[Phi]} \[Element] Reals]},
   {r, FullSimplify[ToSphericalCoordinates[r], 0 < \[Theta] < Pi && -Pi < \[Phi] < Pi]}]]
}
```

**QF** : `"BlochVector"` computes the same Pauli expectation values internally (as $\mathrm{Tr}[\rho\,\sigma_j]$); `"BlochSphericalCoordinates"` gives the angles.

```wl
qs = QuantumState[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}];
FullSimplify[{qs["BlochVector"], qs["BlochSphericalCoordinates"]},
 0 < \[Theta] < Pi && -Pi < \[Phi] < Pi]
```

Evaluating the definition returns $\{2\,\mathrm{Re}\,\bar\alpha\beta,\ 2\,\mathrm{Im}\,\bar\alpha\beta,\ |\alpha|^2-|\beta|^2\}$ for a general state and $(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$ for the parametrized one, with spherical coordinates $(1,\theta,\phi)$. So the component formula is a *result* of evaluating the Pauli expectation values, not the definition itself, and the parametrization angles are exactly the polar and azimuthal angles on the Bloch sphere.

### 1.9 [BSc] How do I prepare a qudit in $\mathbb{C}^d$ and a uniform superposition?

A qudit lives in $\mathbb{C}^d$; the computational state $|0\rangle$ and the equal-weight
superposition $\tfrac1{\sqrt d}\sum_k |k\rangle$ are the basic preparations.

**WL** : $|0\rangle$ is a unit vector; the uniform state is a constant vector handed to `Normalize`
(no hardcoded $1/\sqrt d$).

```wl
With[{d = 3}, {UnitVector[d, 1], Normalize @ ConstantArray[1, d]}]
```

**QF** : name the basis state with its dimension, and use the built-in `"UniformSuperposition"`
state rather than constructing it by hand.

```wl
With[{d = 3},
 {QuantumState["0", d]["StateVector"] // Normal,
  QuantumState["UniformSuperposition", d]["StateVector"] // Normal}]
```

Both give $|0\rangle = (1,0,0)$ and the uniform qutrit $\tfrac1{\sqrt3}(1,1,1)$; the construction is
written in $d$, so it works for any dimension.
