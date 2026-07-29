## Part 5. Projective measurement

### 5.1 [BSc] How do I perform a projective measurement and read outcomes, probabilities, and the mean?

A projective measurement of an observable $A$ returns one of $A$'s eigenvalues $a$ (the outcome) with
probability $\mathrm{Tr}[P_a\rho] = \langle a|\rho|a\rangle$ (the Born rule, $|\langle a|\psi\rangle|^2$ for a
pure state), and the mean is $\langle A\rangle = \mathrm{Tr}[A\rho]$. In short, the eigenvalues are what you
can see and $\rho$ weights them. We measure $A = Z$ on a general qubit, building the state and the
measurement each once and reusing them.

The general mixed qubit is the Bloch state $\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ with $|\vec r|\le1$
($\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`); define it once:

```wl
rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])
```

**WL** : the possible outcomes are the eigenvalues of the observable, here in the order $\{-1, +1\}$:

```wl
Eigenvalues[PauliMatrix[3]]
```

Each outcome's probability is the corresponding eigenvector sandwiched around $\rho$,
$\langle a|\rho|a\rangle$, listed in that same eigenvalue order:

```wl
FullSimplify[Conjugate[#] . rho . # & /@ Eigenvectors[PauliMatrix[3]], -1 <= rz <= 1]
```

The mean is $\mathrm{Tr}[Z\rho]$:

```wl
FullSimplify[Tr[PauliMatrix[3] . rho], -1 <= rz <= 1]
```

**QF** : applying the measurement operator to the state returns a `QuantumMeasurement` object, whose summary
box lists the outcomes and their probabilities. Build it once from the same Bloch state, using the
`"BlochVector"` constructor:

```wl
meas = QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState["BlochVector"[{rx, ry, rz}]]]
```

Its `"Mean"` is $\langle Z\rangle$:

```wl
FullSimplify[meas["Mean"], -1 <= rz <= 1]
```

Its `"Probabilities"` gives the same two weights as an association keyed by each outcome, so the pairing
needs no ordering convention:

```wl
FullSimplify[#, -1 <= rz <= 1] & /@ meas["Probabilities"]
```

Both assign $+1 \mapsto \tfrac{1+r_z}2$ and $-1 \mapsto \tfrac{1-r_z}2$ (the WL side an ordered list in
eigenvalue order, the QF side an association keyed by outcome) and give the mean $\langle Z\rangle = r_z$: a $Z$ measurement reads only
the populations (the $r_z$ component), so the coherences $r_x, r_y$ drop out, and the mean is the difference
of the two outcome probabilities. For a pure state ($|\vec r| = 1$) these are the familiar
$\cos^2\tfrac\theta2, \sin^2\tfrac\theta2, \cos\theta$.

The same three reads carry to any dimension unchanged. For a qutrit, take $J_z = \mathrm{diag}(1,0,-1)$ and
build it together with its eigensystem: the eigenvalues are the outcomes, the (normalized) eigenvectors the
outcome states.

```wl
jz3 = DiagonalMatrix[{1, 0, -1}]; {jzEigenVal, jzEigenVec} = MapAt[Normalize, Eigensystem[jz3], {2, All}]
```

Measure it on a concrete normalized state:

```wl
psi3 = Normalize[{1, I, 2}]
```

The Born rule pairs each eigenvalue with the squared overlap $|\langle a|\psi_3\rangle|^2$ of its
eigenvector; `AssociationThread` keys those weights by outcome, so the reading needs no ordering convention:

```wl
AssociationThread[jzEigenVal, Abs[Conjugate[#] . psi3]^2 & /@ jzEigenVec]
```

and the mean is $\langle J_z\rangle = \langle\psi_3|J_z|\psi_3\rangle$:

```wl
Conjugate[psi3] . jz3 . psi3
```

The QF measurement object carries the same reads. Applying it to the state returns a `QuantumMeasurement`
whose summary box already lists the three outcomes and their weights:

```wl
meas3 = QuantumMeasurementOperator[QuantumOperator[jz3]][QuantumState[psi3]]
```

Its `"Mean"` is $\langle J_z\rangle$:

```wl
meas3["Mean"]
```

and its `"Probabilities"` are the same outcome-keyed weights as the `AssociationThread` above:

```wl
meas3["Probabilities"]
```

Both give outcomes $\{1, 0, -1\}$ with $P(1) = P(0) = \tfrac16$ and $P(-1) = \tfrac23$ and mean
$\langle J_z\rangle = -\tfrac12$, all in exact arithmetic: the qubit machinery is dimension-agnostic, only
the number of outcomes grows.

### 5.2 [BSc] How do I measure in a non-computational basis and get the post-measurement (collapsed) states?

Measuring an observable whose eigenbasis is not the computational one is a change of basis: express the
state in that eigenbasis, and the new amplitudes give the outcome probabilities by the Born rule, while the
measurement collapses the qubit onto one of the eigenstates. We measure $X$, with eigenstates
$|\pm\rangle = (|0\rangle\pm|1\rangle)/\sqrt2$, on a general pure qubit $|\psi\rangle$, defined once:

```wl
psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]};
```

**WL** : the outcomes live in the $X$ eigenbasis, so build it first (normalized), in eigenvalue order
$\{-1, +1\}$:

```wl
eigenX = Normalize /@ Eigenvectors[PauliMatrix[1]]
```

The amplitudes of $|\psi\rangle$ in the new basis are its overlaps with those eigenvectors,
$\langle\pm|\psi\rangle$; name them `amplX`:

```wl
amplX = FullSimplify@ComplexExpand[Conjugate[#] . psi & /@ eigenX]
```

Squaring those amplitudes gives the Born probabilities $|\langle\pm|\psi\rangle|^2$, the first of two ways to
read them here; the QF measurement below is the second:

```wl
FullSimplify /@ ComplexExpand[Abs[amplX]^2]
```

**QF** : applying the `"X"` measurement operator to the state returns a `QuantumMeasurement` object:

```wl
measX = QuantumMeasurementOperator["X"][QuantumState[psi]]
```

Its `"Probabilities"` are exactly those squared amplitudes, an association keyed by each $X$ outcome $x_\pm$
so the pairing needs no ordering convention:

```wl
FullSimplify[#, {\[Theta], \[Phi]} \[Element] Reals] & /@ measX["Probabilities"]
```

Both give the outcome probabilities $\tfrac12(1\pm\sin\theta\cos\varphi) = \tfrac12(1\pm\langle X\rangle)$
(the WL list in eigenvalue order $\{-1,+1\}$, the QF association labeling each outcome $x_\pm$), and the
measurement leaves the qubit in $|\pm\rangle$. Unlike a $Z$ measurement (5.1), the $X$ outcome depends on the
relative phase $\varphi$ through $\cos\varphi$, so a real state would hide it: the choice of measurement
basis decides which feature of the state becomes visible.

The change of basis works in any dimension. For a qutrit, measure in the Fourier basis, the qudit
generalization of the $X$ basis, whose amplitudes are the discrete Fourier transform of the state:

```wl
psi3 = Normalize[{1, 1 + I, 1}]
```

WL: apply the $3\times3$ discrete Fourier transform to get the amplitudes,

```wl
FourierMatrix[3] . psi3 // FullSimplify
```

and square them for the Born probabilities:

```wl
Abs[FourierMatrix[3] . psi3]^2 // FullSimplify
```

QF: that Fourier basis is exactly the named `"X"[3]`, the qudit-$d{=}3$ generalization of the $X$ basis, so
re-expressing the qutrit in it needs no explicit matrix; its `"StateVector"` reproduces those amplitudes:

```wl
FullSimplify /@ Normal@QuantumState[QuantumState[psi3], "X"[3]]["StateVector"]
```

and its `"Probabilities"` are the same weights, keyed by outcome:

```wl
FullSimplify /@ QuantumState[QuantumState[psi3], "X"[3]]["Probabilities"]
```

Both give the Fourier-basis probabilities $\{\tfrac56, \tfrac1{12}, \tfrac1{12}\}$, now exact: the
change-of-basis measurement is the identical construction at $d=3$, with the $3\times3$ discrete Fourier
transform playing the role the $2\times2$ Hadamard played for the qubit.

### 5.3 [BSc] How do I simulate finite-shot statistics and watch the empirical frequencies approach the Born rule?

Real measurements return a finite sample of outcomes whose frequencies converge to the Born probabilities
as the shot count grows. This is the one place the two answers do not agree cell-for-cell: each draws its
own random sample, and what they share is the limit. Measure $Z$ on the state
$\{\tfrac{\sqrt3}2, \tfrac12\}$, whose Born probabilities are the unequal $\tfrac34, \tfrac14$.

**WL** : draw $2000$ outcomes with the Born weights and tabulate their frequencies.

```wl
SeedRandom[1]; N[Counts[RandomChoice[{3/4, 1/4} -> {1, -1}, 2000]]/2000]
```

**QF** : `"SimulatedMeasurement"` draws the same-size sample of outcomes; `Counts` tabulates it and dividing
by the shot count gives the empirical frequencies, keyed by outcome.

```wl
SeedRandom[1]; Counts[QuantumMeasurementOperator["Z"][QuantumState[{Sqrt[3]/2, 1/2}]]["SimulatedMeasurement", 2000]]/2000 // N
```

Both cluster near the Born prediction, frequencies $\approx \{0.75, 0.25\}$: the two samples differ in their
individual draws but each approaches $\tfrac34, \tfrac14$, and the gap shrinks like $1/\sqrt N$. The unequal
weights make the convergence discriminating, the frequencies track two distinct probabilities rather than a
symmetric coin flip. The Born rule is the large-sample limit of the measured frequencies.

### 5.4 [BSc] How do I apply a non-selective projective measurement (the Lüders channel)?

A non-selective measurement records that a measurement happened but not which outcome, replacing $\rho$ by
$\sum_k P_k\,\rho\,P_k$ (the Lüders channel). For a $Z$ measurement this erases the off-diagonal coherences
and leaves the diagonal: pure dephasing. Apply it to a general qubit state in Bloch form
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ (pure or mixed; $\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`); define it
once:

```wl
rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])
```

**WL** : sandwich $\rho$ between the two $Z$ projectors and sum.

```wl
Sum[p . rho . p, {p, {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}}]
```

**QF** : the measurement's `"PostMeasurementState"` is that non-selective mixture, as a state object.

```wl
QuantumMeasurementOperator["Z"][QuantumState[rho]]["PostMeasurementState"]
```

Reading its density matrix off that state object gives the same diagonal the WL sandwich produced:

```wl
%["DensityMatrix"] // Normal // FullSimplify
```

Both give $\mathrm{diag}\big(\tfrac{1+r_z}2, \tfrac{1-r_z}2\big)$: the coherences $r_x, r_y$ are erased while
the populations $\tfrac{1\pm r_z}2$ are untouched, for any state, pure or mixed. This is decoherence in the
measured basis, the same map a dephasing environment applies, and it is what separates a measurement that is
*read* from one that is merely coupled.
