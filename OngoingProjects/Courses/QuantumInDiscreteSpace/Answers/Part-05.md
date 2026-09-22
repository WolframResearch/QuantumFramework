## Part 5. Projective measurement

### 5.1 [BSc] How do I perform a projective measurement and read outcomes, probabilities, and the mean?

A projective measurement of an observable $A$ returns one of $A$'s eigenvalues $a$ (the outcome) with
probability $\mathrm{Tr}[P_a\rho] = \langle a|\rho|a\rangle$ (the Born rule, $|\langle a|\psi\rangle|^2$ for a
pure state), and the mean is $\langle A\rangle = \mathrm{Tr}[A\rho]$: the eigenvalues are what can be seen
and $\rho$ weights them. Measure $A = Z$ on the general mixed qubit, the Bloch state
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ with $|\vec r|\le1$ ($\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

**WL** : build the state once.

```wl
rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])
```

The possible outcomes are the eigenvalues of the observable, here in the order $\{-1, +1\}$.

```wl
Eigenvalues[PauliMatrix[3]]
```

Each outcome's probability is the corresponding eigenvector sandwiched around $\rho$,
$\langle a|\rho|a\rangle$, listed in that same eigenvalue order.

```wl
probs = FullSimplify[Conjugate[#] . rho . # & /@ Eigenvectors[PauliMatrix[3]], -1 <= rz <= 1]
```

The mean is $\mathrm{Tr}[Z\rho]$.

```wl
FullSimplify[Tr[PauliMatrix[3] . rho], -1 <= rz <= 1]
```

For a pure state ($|\vec r| = 1$), substituting the Bloch vector of $\{\cos\tfrac\theta2, e^{i\varphi}\sin\tfrac\theta2\}$
returns the familiar half-angle populations and the mean $\cos\theta$.

```wl
FullSimplify[{probs, Tr[PauliMatrix[3] . rho]} /. Thread[{rx, ry, rz} -> FromSphericalCoordinates[{1, \[Theta], \[Phi]}]], 0 <= \[Theta] <= Pi]
```

**QF** : applying the measurement operator to the state returns a `QuantumMeasurement` object, whose
summary box lists the outcomes and their probabilities. Build it once from the same Bloch state, using
the `"BlochVector"` constructor.

```wl
meas = QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState["BlochVector"[{rx, ry, rz}]]]
```

Its `"Mean"` is $\langle Z\rangle$.

```wl
FullSimplify[meas["Mean"], -1 <= rz <= 1]
```

Its `"Probabilities"` gives the two weights as an association whose keys are the framework's outcome
labels, each displaying its eigenvalue, so the pairing is carried by the key rather than by list position.

```wl
FullSimplify[#, -1 <= rz <= 1] & /@ meas["Probabilities"]
```

The two routes agree: the mean matches $\mathrm{Tr}[Z\rho]$, and pairing each eigenvalue with its
probability (`"EigenvalueVectors"` against `"ProbabilitiesList"`) reproduces the WL eigenvalue-to-weight
map key for key, not merely as an unordered set.

```wl
{FullSimplify[meas["Mean"] == Tr[PauliMatrix[3] . rho], -1 <= rz <= 1], FullSimplify[Values[KeySort[AssociationThread[Flatten[meas["EigenvalueVectors"]], meas["ProbabilitiesList"]]]] == Values[KeySort[AssociationThread[Eigenvalues[PauliMatrix[3]], probs]]], -1 <= rz <= 1]}
```

Both assign $+1 \mapsto \tfrac{1+r_z}2$ and $-1 \mapsto \tfrac{1-r_z}2$ and give the mean
$\langle Z\rangle = r_z$: a $Z$ measurement reads only the populations (the $r_z$ component), so the
coherences $r_x, r_y$ drop out, and the mean is the difference of the two outcome probabilities. On a pure
state the populations collapse to $\cos^2\tfrac\theta2, \sin^2\tfrac\theta2$ and the mean to $\cos\theta$.

The same three reads carry to any dimension unchanged. For a qutrit, take $J_z = \mathrm{diag}(1,0,-1)$
on a concrete normalized state.

**WL** : build $J_z$ together with its eigensystem: the eigenvalues are the outcomes, the (normalized)
eigenvectors the outcome states.

```wl
jz3 = DiagonalMatrix[{1, 0, -1}]; {jzEigenVal, jzEigenVec} = MapAt[Normalize, Eigensystem[jz3], {2, All}]
```

The state to measure, chosen with three distinct Born weights so the eigenvalue-to-probability pairing
has something to test.

```wl
psi3 = Normalize[{1, 2 I, 3}]
```

The Born rule pairs each eigenvalue with the squared overlap $|\langle a|\psi_3\rangle|^2$ of its
eigenvector; `AssociationThread` keys those weights by eigenvalue.

```wl
wl3 = AssociationThread[jzEigenVal, Abs[Conjugate[#] . psi3]^2 & /@ jzEigenVec]
```

The mean is $\langle J_z\rangle = \langle\psi_3|J_z|\psi_3\rangle$.

```wl
Conjugate[psi3] . jz3 . psi3
```

**QF** : the measurement object carries the same reads; applying it to the state returns a
`QuantumMeasurement` whose summary box already lists the three outcomes and their weights.

```wl
meas3 = QuantumMeasurementOperator[QuantumOperator[jz3]][QuantumState[psi3]]
```

Its `"Mean"` is $\langle J_z\rangle$.

```wl
meas3["Mean"]
```

Its `"Probabilities"`, keyed by outcome label, carry the same weights; the two routes agree on the mean
and, pairing eigenvalue with probability, on the whole eigenvalue-to-weight map.

```wl
{meas3["Mean"] == Conjugate[psi3] . jz3 . psi3, Values[KeySort[AssociationThread[Flatten[meas3["EigenvalueVectors"]], meas3["ProbabilitiesList"]]]] == Values[KeySort[wl3]]}
```

Both give outcomes $\{1, 0, -1\}$ with the three distinct weights $\tfrac1{14}, \tfrac27, \tfrac9{14}$ and
mean $\langle J_z\rangle = -\tfrac47$, all in exact arithmetic, and the eigenvalue-to-probability map agrees
key for key: the qubit machinery is dimension-agnostic, only the number of outcomes grows.

### 5.2 [BSc] How do I measure in a non-computational basis and get the post-measurement (collapsed) states?

Measuring an observable whose eigenbasis is not the computational one is a change of basis: express the
state in that eigenbasis, and the new amplitudes give the outcome probabilities by the Born rule, while
the measurement collapses the qubit onto the eigenstate of the outcome found. Conditioned on outcome $a$
the update is the Lüders rule, $\rho \mapsto P_a\rho P_a/\mathrm{Tr}[P_a\rho]$ with
$P_a = |a\rangle\langle a|$; for a rank-one $P_a$ the result is $|a\rangle\langle a|$ itself, whatever the
state was, and the rule is undefined only for an outcome of probability zero, where
$\mathrm{Tr}[P_a\rho] = 0$. Measure $X$, with eigenstates $|\pm\rangle = (|0\rangle\pm|1\rangle)/\sqrt2$, on
a general pure qubit $|\psi\rangle$.

**WL** : the state, defined once.

```wl
psi = {Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]};
```

The outcomes live in the $X$ eigenbasis, so build it first (normalized), in eigenvalue order $\{-1, +1\}$.

```wl
eigenX = Normalize /@ Eigenvectors[PauliMatrix[1]]
```

The amplitudes of $|\psi\rangle$ in the new basis are its overlaps with those eigenvectors,
$\langle\pm|\psi\rangle$.

```wl
amplX = FullSimplify@ComplexExpand[Conjugate[#] . psi & /@ eigenX]
```

Squaring those amplitudes gives the Born probabilities $|\langle\pm|\psi\rangle|^2$.

```wl
FullSimplify /@ ComplexExpand[Abs[amplX]^2]
```

For the collapse, the state as a density matrix.

```wl
rhoPsi = KroneckerProduct[psi, Conjugate[psi]];
```

The Lüders update for each outcome, in the same eigenvalue order, returns the two collapsed states as
density matrices.

```wl
collapsed = FullSimplify[Table[With[{P = KroneckerProduct[e, Conjugate[e]]}, P . rhoPsi . P/Tr[P . rhoPsi]], {e, eigenX}], {\[Theta], \[Phi]} \[Element] Reals]
```

**QF** : applying the $X$ measurement operator to the state returns a `QuantumMeasurement` object.

```wl
measX = QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState[psi]]
```

Its `"Probabilities"` are exactly those squared amplitudes, an association keyed by the outcome labels,
the $+1$ outcome listed first.

```wl
FullSimplify[#, {\[Theta], \[Phi]} \[Element] Reals] & /@ measX["Probabilities"]
```

Its `"StateAssociation"` holds the post-measurement state of each outcome under the same labels, here
listed with the $-1$ outcome first; each is sub-normalized, its squared norm being that outcome's
probability.

```wl
Simplify[ComplexExpand[#["Norm"]^2], {\[Theta], \[Phi]} \[Element] Reals] & /@ Values[measX["StateAssociation"]]
```

Normalizing each and reading its density matrix gives the collapsed states.

```wl
collapsedQF = FullSimplify[ComplexExpand[Normal[#["Normalized"]["DensityMatrix"]]], {\[Theta], \[Phi]} \[Element] Reals] & /@ Values[measX["StateAssociation"]]
```

Taken in that listed order, $-1$ then $+1$, they are the projectors the Lüders update produced.

```wl
collapsedQF == collapsed
```

At the edge where an outcome is impossible, the update is undefined. Measure $X$ on the $X$ eigenstate
$|+\rangle$: the $-1$ outcome has probability zero, so its Lüders quotient is $0/0$, and the framework
keeps only the one outcome that can occur.

```wl
QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState["Plus"]]["Probabilities"]
```

```wl
Length[QuantumMeasurementOperator[QuantumOperator["X"]][QuantumState["Plus"]]["StateAssociation"]]
```

Both give the outcome probabilities $\tfrac12(1\pm\sin\theta\cos\varphi) = \tfrac12(1\pm\langle X\rangle)$
and, for either outcome, the collapsed state $|\pm\rangle\langle\pm|$, which keeps no memory of $\theta$ or
$\varphi$ beyond the outcome itself: the pre-measurement state decides only how likely each outcome is.
On the $X$ eigenstate the opposite outcome drops out, one probability going to zero and the other to one,
which is where the update rule has nothing to normalize. Unlike a $Z$ measurement (5.1), the $X$ outcome
depends on the relative phase $\varphi$ through $\cos\varphi$, so a real state would hide it: the choice of
measurement basis decides which feature of the state becomes visible.

The change of basis works in any dimension. For a qutrit, measure in the Fourier basis, the qudit
generalization of the $X$ basis, whose amplitudes are the discrete Fourier transform of the state. Choose
a state whose three Fourier weights are all distinct, so the direction of the transform is not hidden by a
coincidence.

**WL** : the state to measure.

```wl
psi3 = Normalize[{1, 2 + I, 3 I}]
```

Apply the $3\times3$ discrete Fourier transform to get the amplitudes in the Fourier basis.

```wl
Simplify@ComplexExpand[FourierMatrix[3] . psi3]
```

Square them for the Born probabilities.

```wl
FullSimplify[Abs[FourierMatrix[3] . psi3]^2]
```

The collapse targets are the Fourier basis states, the rank-one projectors onto the rows of the transform;
the Lüders update sends the qutrit onto one of them.

```wl
collapsed3 = RootReduce[Table[With[{P = KroneckerProduct[e, Conjugate[e]]}, P . KroneckerProduct[psi3, Conjugate[psi3]] . P/Tr[P . KroneckerProduct[psi3, Conjugate[psi3]]]], {e, Conjugate[FourierMatrix[3]]}]];
```

Each is a pure projector: rank one and idempotent.

```wl
{MatrixRank[#], # . # == #} & /@ collapsed3
```

**QF** : that Fourier basis is exactly the named `"X"[3]`, the qudit-$d{=}3$ generalization of the $X$
basis, so re-expressing the qutrit in it needs no explicit matrix; its `"StateVector"` reproduces those
amplitudes.

```wl
Simplify[RootReduce[Normal@QuantumState[QuantumState[psi3], "X"[3]]["StateVector"]] == RootReduce[FourierMatrix[3] . psi3]]
```

Measuring the shift operator `"X"[3]` (whose eigenbasis is the Fourier basis) returns the measurement
object; its `"Probabilities"` are the same weights.

```wl
measF = QuantumMeasurementOperator[QuantumOperator["X"[3]]][QuantumState[psi3]];
RootReduce /@ measF["Probabilities"]
```

Its post-measurement states, normalized, are the Fourier-basis projectors: the same set the Lüders update
produced.

```wl
Sort[RootReduce[Normal[#["DensityMatrix"]]/#["Norm"]^2] & /@ Values[measF["StateAssociation"]]] == Sort[collapsed3]
```

Both give the Fourier-basis probabilities $\{\tfrac59, \tfrac{2(5-2\sqrt3)}{45}, \tfrac{2(5+2\sqrt3)}{45}\}$,
exact and all distinct, and collapse the qutrit onto a Fourier basis state: the change-of-basis measurement
is the identical construction at $d=3$, with the $3\times3$ discrete Fourier transform playing the role the
$2\times2$ Hadamard played for the qubit. The distinct weights are what make the transform's direction
testable, since the inverse transform would return the last two weights swapped.

### 5.3 [BSc] How do I simulate finite-shot statistics and watch the empirical frequencies approach the Born rule?

Real measurements return a finite sample of outcomes whose frequencies converge to the Born probabilities
as the shot count grows. This is the one place the two answers do not agree cell-for-cell: each draws its
own random sample, and what they share is the limit. Measure $Z$ on the state
$\{\tfrac{\sqrt3}2, \tfrac12\}$, whose Born probabilities are the unequal $\tfrac34, \tfrac14$ computed
from the amplitudes rather than typed in.

**WL** : draw $2000$ outcomes weighted by the Born probabilities $|\langle a|\psi\rangle|^2$ and tabulate
their frequencies, keyed in sorted order.

```wl
SeedRandom[1]; N[KeySort[Counts[RandomChoice[Abs[{Sqrt[3]/2, 1/2}]^2 -> {1, -1}, 2000]]]/2000]
```

The sampling error shrinks like $1/\sqrt N$: averaged over $100$ runs at each shot count, the mean absolute
deviation of the $+1$ frequency, scaled by $\sqrt N$, holds near the constant $\sqrt{p(1-p)}\sqrt{2/\pi}$.

```wl
SeedRandom[3]; N[Table[{n, Sqrt[n] Mean[Table[Abs[Count[RandomChoice[Abs[{Sqrt[3]/2, 1/2}]^2 -> {1, -1}, n], 1]/n - First[Abs[{Sqrt[3]/2, 1/2}]^2]], 100]]}, {n, {100, 1000, 10000}}]]
```

**QF** : `"SimulatedCounts"` draws a same-size sample under its own seed and returns the outcome tally
directly; dividing by the shot count gives the empirical frequencies.

```wl
SeedRandom[2]; N[QuantumMeasurementOperator[QuantumOperator["Z"]][QuantumState[{Sqrt[3]/2, 1/2}]]["SimulatedCounts", 2000]/2000]
```

Both cluster near the Born prediction, three quarters for $+1$ and one quarter for $-1$: the two samples
differ in their individual draws but each approaches $\tfrac34, \tfrac14$, and the scaled deviation staying
flat across three decades of $N$ is the $1/\sqrt N$ law. The unequal weights make the convergence
discriminating, the frequencies track two distinct probabilities rather than a symmetric coin flip. The
Born rule is the large-sample limit of the measured frequencies.

### 5.4 [BSc] How do I apply a non-selective projective measurement (the Lüders channel)?

A non-selective measurement records that a measurement happened but not which outcome, replacing $\rho$ by
$\sum_k P_k\,\rho\,P_k$ (the Lüders channel). For a $Z$ measurement this erases the off-diagonal coherences
and leaves the diagonal: pure dephasing. Apply it to a general qubit state in Bloch form
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ (pure or mixed; $\vec\sigma = $ `PauliMatrix[{1, 2, 3}]`).

**WL** : the state, defined once.

```wl
rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])
```

The $Z$ projectors are the outer products of the computational basis, so build them from the identity
rather than typing them out, then sandwich $\rho$ between each and sum.

```wl
proj = KroneckerProduct[#, Conjugate[#]] & /@ IdentityMatrix[2];
luders = Total[# . rho . # & /@ proj]
```

The map is a channel and a projection: it preserves the trace, it is idempotent (measuring twice is
measuring once), and on the Bloch vector it keeps $r_z$ and kills $r_x, r_y$.

```wl
{Tr[luders] == 1, Total[# . luders . # & /@ proj] == luders, Simplify[Tr[luders . #] & /@ PauliMatrix[{1, 2, 3}]]}
```

**QF** : the measurement's `"PostMeasurementState"` is that non-selective mixture, as a state object.

```wl
post = QuantumMeasurementOperator["Z"][QuantumState[rho]]["PostMeasurementState"]
```

Reading its density matrix off that state object gives the same diagonal the WL sandwich produced.

```wl
FullSimplify[Normal[post["DensityMatrix"]] == luders]
```

Both give $\mathrm{diag}\big(\tfrac{1+r_z}2, \tfrac{1-r_z}2\big)$, with Bloch vector $(0, 0, r_z)$: the
coherences $r_x, r_y$ are erased while the populations $\tfrac{1\pm r_z}2$ are untouched, for any state,
pure or mixed, and applying the channel again changes nothing. This is decoherence in the measured basis,
the same map a dephasing environment applies, and it is what separates a measurement that is *read* from
one that is merely coupled.
