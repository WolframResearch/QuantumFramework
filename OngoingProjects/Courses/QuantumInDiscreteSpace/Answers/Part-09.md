## Part 9. Single-qubit operations and SU(2)

### 9.1a [BSc] How do I build the standard named single-qubit gates $X, Y, Z, H, S, T$?

The named single-qubit gates are the workhorses of quantum computing: the three Paulis $X, Y, Z$, the
Hadamard $H = \tfrac1{\sqrt2}\begin{psmallmatrix}1&1\\1&-1\end{psmallmatrix}$, and the phase gates
$S = \mathrm{diag}(1,i)$ and $T = \mathrm{diag}(1,e^{i\pi/4})$ (with $S = T^2$).

**WL** : each is a small explicit matrix (the Paulis from `PauliMatrix`, the Hadamard from
`HadamardMatrix`, the phase gates diagonal).

```wl
MatrixForm /@ <|"X" -> PauliMatrix[1], "Y" -> PauliMatrix[2], "Z" -> PauliMatrix[3], "H" -> HadamardMatrix[2], "S" -> DiagonalMatrix[{1, I}], "T" -> DiagonalMatrix[{1, Exp[I Pi/4]}]|>
```

**QF** : each is a named `QuantumOperator`; here their matrices, read off the objects.

```wl
AssociationMap[MatrixForm[QuantumOperator[#]["Matrix"]] &, {"X", "Y", "Z", "H", "S", "T"}]
```

Both give the same six matrices. The QF names carry more than a matrix: each `QuantumOperator["H"]` knows
its arity, its qubit order, that it is unitary and Hermitian (for the Paulis and $H$), and composes with
the others by `@`, so the named gate is a reusable building block, not a frozen array.

### 9.1b [BSc] How do I build the axis-angle rotation gates $R_{\hat n}(\theta)=e^{-i\theta\,\hat n\cdot\vec\sigma/2}$?

A rotation gate turns the qubit by angle $\theta$ about a Bloch axis $\hat n$: $R_{\hat n}(\theta) =
e^{-i\theta\,\hat n\cdot\vec\sigma/2}$. The three coordinate rotations are $R_x, R_y, R_z$; here is the
general one and the $R_z$ special case.

**WL** : the matrix exponential of $-i\tfrac\theta2\,\hat n\cdot\vec\sigma$ for a general unit axis
$\hat n$ from its spherical angles.

```wl
With[{n = FromSphericalCoordinates[{1, \[Alpha], \[Beta]}]}, FullSimplify[MatrixExp[-I (\[Theta]/2) n . PauliMatrix[{1, 2, 3}]], {\[Alpha], \[Beta], \[Theta]} \[Element] Reals]]
```

Sending the polar angle to zero points the axis along $\hat z$, so for $R_z$ we get the substitution below.
It is made in the resulting matrix rather than in the axis, since at the pole the azimuth $\beta$ is
degenerate and `FromSphericalCoordinates` rejects it there.

```wl
% /. {\[Alpha] -> 0}
```

**QF** : the named rotation gates `"RX"`, `"RY"`, `"RZ"` (call-form with the angle); here $R_z(\theta)$.

```wl
QuantumOperator["RZ"[\[Theta]]]
```

The object prints as itself; its matrix is the same diagonal the WL side arrived at.

```wl
%["Matrix"] // MatrixForm
```

Both produce a unit-determinant unitary that rotates the Bloch sphere by $\theta$ about $\hat n$; for
$\hat n = \hat z$ it is $R_z(\theta) = \mathrm{diag}(e^{-i\theta/2}, e^{+i\theta/2})$. The QF object also
carries the angle as a parameter it can substitute or compose.

### 9.2 [BSc] How do I derive the closed-form axis-angle exponential $e^{-i\alpha\,\hat n\cdot\vec\sigma/2} = \cos(\alpha/2)I - i\sin(\alpha/2)\,\hat n\cdot\vec\sigma$?

Because $(\hat n\cdot\vec\sigma)^2 = I$ for a unit $\hat n$, the exponential series splits into even and
odd terms and closes in a single line: $e^{-i\alpha\,\hat n\cdot\vec\sigma/2} = \cos\tfrac\alpha2\,I -
i\sin\tfrac\alpha2\,\hat n\cdot\vec\sigma$. Verify this identity for a general axis and angle.

**WL** : compare the matrix exponential with the proposed closed form; the difference should vanish.

```wl
With[{n = FromSphericalCoordinates[{1, \[Theta], \[Phi]}]}, Simplify[MatrixExp[-I (\[Alpha]/2) n . PauliMatrix[{1, 2, 3}]] == Cos[\[Alpha]/2] IdentityMatrix[2] - I Sin[\[Alpha]/2] n . PauliMatrix[{1, 2, 3}], {\[Theta], \[Phi], \[Alpha]} \[Element] Reals]]
```

**QF** : the same rotation built as a `QuantumOperator` from the axis, whose matrix equals the closed form.

```wl
With[{n = FromSphericalCoordinates[{1, \[Theta], \[Phi]}]}, Simplify[QuantumOperator[MatrixExp[-I (\[Alpha]/2) n . PauliMatrix[{1, 2, 3}]]] == Cos[\[Alpha]/2] QuantumOperator["I"] - I Sin[\[Alpha]/2] QuantumOperator[n . PauliMatrix[{1, 2, 3}]], {\[Theta], \[Phi], \[Alpha]} \[Element] Reals]]
```

Both return True: the entire one-parameter rotation is captured by two scalars, $\cos\tfrac\alpha2$ and
$\sin\tfrac\alpha2$, because the generator squares to the identity. This is the single most useful
identity for single-qubit gates, and it makes $R_{\hat n}(2\pi) = -I$ manifest (9.5).

### 9.3 [BSc] How do I verify the unitarity and unit-determinant conditions that define SU(2)?

The single-qubit rotations form the group $SU(2)$: unitary ($U^\dagger U = I$, preserving probability)
and unit determinant ($\det U = 1$, distinguishing $SU(2)$ from the larger $U(2)$). Check both for a
general rotation $R_{\hat n}(\theta)$.

**WL** : test $U^\dagger U = I$ and $\det U = 1$ on the general axis-angle matrix.

```wl
With[{n = FromSphericalCoordinates[{1, \[Alpha], \[Beta]}]}, With[{u = Cos[\[Theta]/2] IdentityMatrix[2] - I Sin[\[Theta]/2] n . PauliMatrix[{1, 2, 3}]}, Simplify[{u . ConjugateTranspose[u] == IdentityMatrix[2], Det[u] == 1}, {\[Alpha], \[Beta], \[Theta]} \[Element] Reals]]]
```

**WL (numerical)** : the same two conditions on a random $SU(2)$ element, drawn as a Haar-random
$U(2)$ matrix and rescaled by $1/\sqrt{\det U}$ so its determinant becomes $1$.

```wl
su = With[{u = RandomVariate@CircularUnitaryMatrixDistribution[2]}, u/Sqrt[Det[u]]];
UnitaryMatrixQ[su]
Chop@Det[su] == 1
```

**QF** : compose the operator with its own `"Dagger"` and simplify that product to the identity, then take
the determinant of its matrix. The `"UnitaryQ"` property is no help with a symbolic angle: having no way to
know $\theta$ is real, it answers `False` rather than complaining, so the reality assumption has to be
supplied explicitly to `FullSimplify` instead.

```wl
With[{u = QuantumOperator["RZ"[\[Theta]]]}, {FullSimplify[Normal[(u["Dagger"] @ u)["Matrix"]] == IdentityMatrix[2], \[Theta] \[Element] Reals], Simplify[Det[Normal[u["Matrix"]]] == 1]}]
```

Both confirm $\{$unitary, $\det = 1\}$: the unit determinant follows because the generator
$\hat n\cdot\vec\sigma$ is traceless, so $\det e^{-i\theta\,\hat n\cdot\vec\sigma/2} = e^{-i\theta\,
\mathrm{Tr}[\hat n\cdot\vec\sigma]/2} = e^0 = 1$. Every single-qubit rotation lives in $SU(2)$.

### 9.4 [BSc] How do I decompose an arbitrary single-qubit unitary into Euler ($Z$-$Y$-$Z$) angles?

Any single-qubit unitary factors into three rotations about alternating axes, $U = e^{i\phi}\,
R_z(\alpha)\,R_y(\beta)\,R_z(\gamma)$ (the $Z$-$Y$-$Z$ Euler decomposition), so three real angles plus a
global phase $\phi$ parametrize all of $U(2)$. The triple $\{\alpha,\beta,\gamma\}$ is exactly the set of
Euler angles of the $SO(3)$ rotation that $U$ performs on the Bloch sphere, so it can be read either from
that geometric rotation or straight from the operator. Take a random $SU(2)$ element and extract its
$Z$-$Y$-$Z$ angles both ways.

**WL** : send the $SU(2)$ matrix to its $SO(3)$ Bloch rotation $R_{ij} = \tfrac12\mathrm{Tr}[\sigma_i\,U\,
\sigma_j\,U^\dagger]$, then read that rotation's Euler angles with the built-in `EulerAngles`, whose default
axis order is $Z$-$Y$-$Z$. The built-in reports the two outer angles in $(-\pi,\pi]$, so `Mod` by $2\pi$
brings them into $[0,2\pi)$.

```wl
su2so3[u_] := Chop@Table[Tr[PauliMatrix[i] . u . PauliMatrix[j] . ConjugateTranspose[u]]/2, {i, 3}, {j, 3}];
su2 = With[{u = RandomVariate@CircularUnitaryMatrixDistribution[2]}, u/Sqrt[Det[u]]];
Mod[EulerAngles@Chop@su2so3[su2], 2 Pi]
```

**QF** : the `QuantumOperator` answers `"EulerAngles"` about itself directly, in the same $Z$-$Y$-$Z$
convention, but reports the two outer angles already in $[0,2\pi)$.

```wl
Mod[QuantumOperator[su2]["EulerAngles"], 2 Pi]
```

Both routes return the same triple $\{\alpha,\beta,\gamma\}$ for the same random $U$, with the middle angle
$\beta\in[0,\pi]$ either way. The only difference is bookkeeping of the two outer angles: the built-in
`EulerAngles` reports them in $(-\pi,\pi]$ while QF reports them in $[0,2\pi)$, and $\mathrm{Mod}[\cdot,2\pi]$
folds both into the common window $[0,2\pi)$, where they coincide to machine precision. The agreement makes
the double cover concrete: the three angles of the $SU(2)$ operator are precisely the three Euler angles of
the $SO(3)$ rotation it induces, so every one-qubit gate is three turns about the $z$ and $y$ axes, the fact
that makes an arbitrary single-qubit gate implementable from just $R_z$ and $R_y$ hardware.

### 9.5 [MSc] How do I show that SU(2) double-covers SO(3) and that a $2\pi$ rotation is not the identity?

$SU(2)$ is a *double* cover of the rotation group $SO(3)$: the map sending a qubit rotation $U$ to the
$SO(3)$ rotation of its Bloch vector, $R_{ij} = \tfrac12\mathrm{Tr}[\sigma_i\,U\,\sigma_j\,U^\dagger]$, is
two-to-one, because $U$ and $-U$ give the same $R$. The signature is that a full $2\pi$ qubit rotation is
$R_{\hat n}(2\pi) = -I$, not the identity; only $4\pi$ returns to $I$.

**WL** : the covering map itself, $U \mapsto R_{ij} = \tfrac12\mathrm{Tr}[\sigma_i\,U\,\sigma_j\,U^\dagger]$,
which reads off the $SO(3)$ rotation a qubit gate performs on the Bloch vector.

```wl
su2so3[u_] := Chop@Table[Tr[PauliMatrix[i] . u . PauliMatrix[j] . ConjugateTranspose[u]]/2, {i, 3}, {j, 3}];
```

A random unit vector for the rotation axis, so nothing below is an accident of a privileged direction.

```wl
n = RandomPoint@Sphere[];
```

Turning about that axis: a full $2\pi$ rotation returns $-I$ rather than the identity, and only the double
turn $4\pi$ brings the qubit back to $+I$.

```wl
MatrixExp[-I (2 \[Pi])/2 n . PauliMatrix[{1, 2, 3}]] // Chop
MatrixExp[-I (4 \[Pi])/2 n . PauliMatrix[{1, 2, 3}]] // Chop
```

On the $SO(3)$ side the sign is gone: `su2so3` sends the axis operator $\hat n\cdot\vec\sigma$ to its real
Bloch-space rotation (eigenvalues $\pm1$), and exponentiating that rotation over a full $2\pi$ collapses to
the $3\times3$ identity.

```wl
MatrixExp[-I 2 \[Pi] su2so3[n . PauliMatrix[{1, 2, 3}]]] // Chop
```

**QF** : the same, with $R_z(2\pi)$ a `QuantumOperator`; its matrix is $-I$ while its adjoint action on the
Paulis is the identity rotation.

```wl
With[{u = QuantumOperator["RZ"[2 Pi]]}, {Normal[u["Matrix"]], Table[Tr[PauliMatrix[i] . Normal[u["Matrix"]] . PauliMatrix[j] . ConjugateTranspose[Normal[u["Matrix"]]]]/2, {i, 3}, {j, 3}]}]
```

Both routes make the same point: a $2\pi$ rotation equals $-I$ (a global sign, physically a state $\to
-$state) while the $SO(3)$ rotation it induces on the Bloch vector is the identity, the same $SO(3)$ element
that $+I$ induces, and only a $4\pi$ turn brings the qubit itself back to $+I$. Two distinct $SU(2)$ elements
$\pm I$ sit above one $SO(3)$ rotation. The unobservable $2\pi$ sign is real: it shows up in interference
between a rotated and an unrotated path (the neutron $4\pi$ symmetry).

### 9.6 [MSc] How do I realize a single-qubit unitary to within a target error $\epsilon$ as a sequence of $\{H, T\}$ gates (the Solovay-Kitaev idea)?

The gates $\{H, T\}$ generate a dense subgroup of $SU(2)$, so any single-qubit unitary can be approximated
to arbitrary accuracy $\epsilon$ by a finite word in $H$ and $T$; the Solovay-Kitaev theorem says the word
length grows only polylogarithmically in $1/\epsilon$. Here we demonstrate the *density* directly by brute
force (not the efficient recursion): search words in $\{H,T\}$ for the best approximation to a target
rotation $U = R_{\hat n}(0.8)$ about the tilted axis $\hat n = (\hat x + \hat z)/\sqrt2$, measured by the
phase-invariant error $1 - \tfrac12|\mathrm{Tr}[W^\dagger U]|$.

**WL** : the two generators, kept numeric because the search multiplies out thousands of products.

```wl
gates = N[<|"H" -> HadamardMatrix[2], "T" -> DiagonalMatrix[{1, Exp[I Pi/4]}]|>];
```

The target, a rotation by $0.8$ about the tilted axis $\hat n = (\hat x + \hat z)/\sqrt2$, tilted precisely
so that no short word lands on it exactly.

```wl
utarget = N[MatrixExp[-I (0.4) (PauliMatrix[1] + PauliMatrix[3])/Sqrt[2]]];
```

The cost of a word: multiply its gates out with `Fold` and score it by
$1 - \tfrac12|\mathrm{Tr}[W^\dagger U]|$, which is blind to global phase and so measures only the physically
meaningful distance to the target.

```wl
err[w_] := 1 - Abs[Tr[ConjugateTranspose[Fold[Dot, IdentityMatrix[2], gates /@ w]] . utarget]]/2;
```

Enumerate every word up to length $L$ and keep the best error, for three increasing length budgets.

```wl
Table[Min[err /@ Flatten[Table[Tuples[{"H", "T"}, len], {len, 1, L}], 1]], {L, {4, 8, 12}}]
```

Take the best word found and confirm it as an explicit gate sequence.

```wl
best = MinimalBy[Flatten[Table[Tuples[{"H", "T"}, len], {len, 1, 12}], 1], err, 1][[1]]
```

**QF** : build that winning word as a `QuantumCircuitOperator` of named $H$ and $T$ gates; its overall
unitary approximates the target to the same error.

```wl
With[{qc = QuantumCircuitOperator[Thread[best -> 1]]}, 1 - Abs[Tr[ConjugateTranspose[Normal[qc["QuantumOperator"]["Matrix"]]] . utarget]]/2]
```

The best-approximation error never increases as longer $\{H,T\}$ words are allowed and drops sharply once
enough length is available (here from $\sim2\times10^{-2}$ to $\sim3\times10^{-5}$), and the QF circuit
built from the winning word reproduces that error: a finite, discrete gate set suffices to reach any
single-qubit unitary as closely as desired, the theoretical basis for compiling arbitrary rotations onto
fault-tolerant hardware whose native gates are exactly $\{H, T\}$ (and the two-qubit CNOT of Part 10).
