---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 5)

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
