## Part 12. Mixed states: distinguishability and thermal states

A density matrix carries two kinds of uncertainty at once, the superposition inside a state vector and the
classical ignorance of which state vector one holds, and this Part works with that combination directly. It
starts with the noisy entangled family whose parameter crosses a threshold from entangled to separable, shows
that every mixed state is a pure state of a larger system seen in part, and that the ensemble behind a density
matrix is never unique. It then builds the thermal state a finite system settles into, and turns to telling two
states apart: the trace distance, the fidelity, and the measurement that achieves the smallest possible error.
It closes with the question of what a probability rule on a Hilbert space is allowed to be at all.

### 12.1 [BSc] How do I build a Werner state $\rho_W(\lambda) = \lambda\,|\Psi^-\rangle\langle\Psi^-| + (1-\lambda)\,\tfrac{I}{4}$ and locate its entanglement threshold?

A Werner state is a singlet diluted by white noise, and what picks the family out is a symmetry: these are the
two-qubit states left unchanged when the *same* unitary acts on both qubits. That invariance leaves exactly one
free parameter, so every statement below depends on which parameter is meant, and the convention has to be
fixed before a threshold can be quoted. Take

$$\rho_W(\lambda) = \lambda\,|\Psi^-\rangle\langle\Psi^-| + (1-\lambda)\,\frac{I}{4}, \qquad |\Psi^-\rangle = \frac{|01\rangle - |10\rangle}{\sqrt2},$$

so $\lambda$ is the weight on the singlet and $1-\lambda$ the weight on the uniform mixture. Read as a mixing
probability $\lambda$ runs over $0 \le \lambda \le 1$, from white noise to a pure singlet, but positivity of
$\rho_W$ is the only real constraint and it permits slightly more, which the spectrum below settles. As
$\lambda$ falls the noise washes the entanglement out, and the value of $\lambda$ where that happens is the
threshold. Since the partial transpose decides separability exactly for two qubits (Part 11, 11.6 and 11.12b),
it sits wherever the smallest eigenvalue of $\rho_W^{T_B}$ changes sign.

**WL** : the singlet, and the family it seeds.

```wl
werner = With[{singlet = Normalize@Flatten[KroneckerProduct[{1, 0}, {0, 1}] - KroneckerProduct[{0, 1}, {1, 0}]]}, \[FormalLambda] KroneckerProduct[singlet, singlet] + (1 - \[FormalLambda]) IdentityMatrix[4]/4]
```

The defining symmetry, checked against a generic single-qubit unitary applied to both factors at once.

```wl
With[{u = {{Cos[\[FormalTheta]/2] Exp[I \[FormalPhi]], Sin[\[FormalTheta]/2] Exp[I \[FormalPsi]]}, {-Sin[\[FormalTheta]/2] Exp[-I \[FormalPsi]], Cos[\[FormalTheta]/2] Exp[-I \[FormalPhi]]}}}, Simplify[{ConjugateTranspose[u] . u == IdentityMatrix[2], KroneckerProduct[u, u] . werner . ConjugateTranspose[KroneckerProduct[u, u]] == werner}, Element[{\[FormalTheta], \[FormalPhi], \[FormalPsi]}, Reals]]]
```

Its spectrum separates the singlet from the threefold triplet, and demanding both be nonnegative fixes how far
$\lambda$ may run.

```wl
{Simplify@Eigenvalues[werner], Reduce[And @@ Thread[Eigenvalues[werner] >= 0], \[FormalLambda], Reals]}
```

Transposing the second qubit alone moves one eigenvalue against the other three.

```wl
Simplify@Eigenvalues@ArrayReshape[Transpose[ArrayReshape[werner, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]
```

The threshold is where that eigenvalue turns negative.

```wl
Reduce[Min[Eigenvalues@ArrayReshape[Transpose[ArrayReshape[werner, {2, 2, 2, 2}], {1, 4, 3, 2}], {4, 4}]] < 0 && -1/3 <= \[FormalLambda] <= 1, \[FormalLambda], Reals]
```

The concurrence of Part 11, 11.7 is an exact measure rather than a criterion, so it locates the same point and
also says how much entanglement survives above it.

```wl
With[{y = KroneckerProduct[PauliMatrix[2], PauliMatrix[2]]}, Simplify[Max[0, 2 Max[#] - Total[#]] &@Sqrt@Simplify@Eigenvalues[werner . y . Conjugate[werner] . y], -1/3 <= \[FormalLambda] <= 1]]
```

**QF** : the family is a named state whose second argument is the local dimension. It carries its own
convention, and comparing thresholds across the two is only meaningful once that convention is written down:
the framework's parameter is the weight on the *triplet*,

$$\rho_W^{\text{QF}}(p) = (1-p)\,|\Psi^-\rangle\langle\Psi^-| + \frac{p}{3}\big(I - |\Psi^-\rangle\langle\Psi^-|\big), \qquad 0 \le p \le 1,$$

running from the pure singlet at $p = 0$ to the maximally mixed triplet at $p = 1$, and passing through the
maximally mixed state at $p = \tfrac34$.

```wl
QuantumState["Werner"[\[FormalP], 2]]
```

```wl
With[{singlet = Normalize@Flatten[KroneckerProduct[{1, 0}, {0, 1}] - KroneckerProduct[{0, 1}, {1, 0}]]}, Simplify[Normal@QuantumState["Werner"[\[FormalP], 2]]["DensityMatrix"] == (1 - \[FormalP]) KroneckerProduct[singlet, singlet] + (\[FormalP]/3) (IdentityMatrix[4] - KroneckerProduct[singlet, singlet])]]
```

The spectrum shows the same split: $1-p$ on the singlet, the rest shared equally by the three triplet states.

```wl
Simplify@QuantumState["Werner"[\[FormalP], 2]]["Eigenvalues"]
```

One family in two coordinates, so matching the density matrices converts between them and the conversion holds
identically, not just at sampled values.

```wl
{Solve[Normal[QuantumState["Werner"[\[FormalP], 2]]["DensityMatrix"]] == werner, \[FormalLambda]], Simplify[Normal@QuantumState["Werner"[\[FormalP], 2]]["DensityMatrix"] == (werner /. \[FormalLambda] -> (3 - 4 \[FormalP])/3)]}
```

The partial transpose is a named operation, and the threshold comes out in the framework's own parameter.

```wl
{Simplify@QuantumPartialTranspose[QuantumState["Werner"[\[FormalP], 2]], {2}]["Eigenvalues"], Reduce[Min[QuantumPartialTranspose[QuantumState["Werner"[\[FormalP], 2]], {2}]["Eigenvalues"]] < 0 && 0 <= \[FormalP] <= 1, \[FormalP], Reals]}
```

Both monotones are named too, and each vanishes on exactly one side of that point.

```wl
Simplify[{QuantumEntanglementMonotone[QuantumState["Werner"[\[FormalP], 2]], "Concurrence"], QuantumEntanglementMonotone[QuantumState["Werner"[\[FormalP], 2]], "Negativity"]}, 0 <= \[FormalP] <= 1]
```

The conversion $\lambda = (3-4p)/3$ turns either threshold into the other.

```wl
Reduce[(3 - 4 \[FormalP])/3 > 1/3 && 0 <= \[FormalP] <= 1, \[FormalP], Reals]
```

One family, two coordinates, one threshold. In the singlet-weight convention $\rho_W(\lambda)$ the spectrum is
$(1+3\lambda)/4$ on the singlet against $(1-\lambda)/4$ three times, so positivity gives
$-\tfrac13 \le \lambda \le 1$: the parameter reaches slightly below $0$, where the mixing reading stops making
sense and the balance tips toward the triplet. The partial transpose holds three eigenvalues at
$(1+\lambda)/4$ and sends the fourth to $(1-3\lambda)/4$, negative exactly on $\tfrac13 < \lambda \le 1$, so a
singlet survives dilution until it holds a third of the weight. The concurrence $\max\{0,(3\lambda-1)/2\}$
vanishes at that same $\lambda = \tfrac13$ and rises to $1$ at $\lambda = 1$, and criterion and measure cannot
disagree here because the partial transpose is exact at $2\times2$ (Part 11, 11.12b). The framework's
triplet-weight convention $\rho_W^{\text{QF}}(p)$ runs the opposite way over $0 \le p \le 1$, spectrum
$\{1-p, p/3, p/3, p/3\}$, and matching the density matrices gives $\lambda = (3-4p)/3$ as an identity. Its
partial transpose has eigenvalues $(3-2p)/6$ three times and $p-\tfrac12$, so entanglement holds on
$0 \le p < \tfrac12$, with concurrence $1-2p$ and negativity $\max\{0,(1-2p)/2\}$ reaching zero there. The two
thresholds $\lambda = \tfrac13$ and $p = \tfrac12$ are the same physical point, since the conversion carries
one inequality onto the other: a singlet tolerates white noise up to two thirds of the total weight, and no
further.

### 12.2 [BSc] How do I purify a mixed state $\rho$ into a pure $|\Psi\rangle$ on a larger system, so that $\mathrm{Tr}_B|\Psi\rangle\langle\Psi| = \rho$?

A mixed state is what a pure state looks like when part of the world is out of view, and purification makes
that literal: adjoin an ancilla $B$ and produce a pure $|\Psi\rangle$ on $AB$ whose reduced state on $A$ is the
$\rho$ one started with. The spectral decomposition supplies it. Writing
$\rho = \sum_i p_i |e_i\rangle\langle e_i|$, the vector

$$|\Psi\rangle = \sum_i \sqrt{p_i}\,|e_i\rangle_A \otimes |i\rangle_B$$

is pure, and tracing out $B$ against the orthonormal $|i\rangle_B$ returns each $p_i |e_i\rangle\langle e_i|$
and nothing else, so the ancilla needs only as many dimensions as $\rho$ has nonzero eigenvalues. Take the
generic mixed qubit in Bloch form, $\rho = \tfrac12(I + r\,\hat n\cdot\vec\sigma)$ with $0 \le r \le 1$ and
$\hat n$ the unit vector at angles $(\theta,\varphi)$: its eigenvalues are $(1\pm r)/2$, so one ancilla qubit
is enough.

**WL** : the state, with the Bloch vector split into its length and its direction.

```wl
rho = With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, 1/2 (IdentityMatrix[2] + \[FormalR] nvec . PauliMatrix[{1, 2, 3}])]
```

Its eigenvectors are the spin-up and spin-down states along $\hat n$, orthonormal, and weighting their
projectors by $(1\pm r)/2$ rebuilds $\rho$.

```wl
{up, down} = {{Cos[\[FormalTheta]/2], Exp[I \[FormalPhi]] Sin[\[FormalTheta]/2]}, {-Exp[-I \[FormalPhi]] Sin[\[FormalTheta]/2], Cos[\[FormalTheta]/2]}}; Simplify[{Conjugate[up] . up, Conjugate[down] . down, Conjugate[up] . down, rho == (1 + \[FormalR])/2 KroneckerProduct[up, Conjugate[up]] + (1 - \[FormalR])/2 KroneckerProduct[down, Conjugate[down]]}, Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

Pairing each eigenvector with its own ancilla ket, amplitudes the square roots of the eigenvalues, gives the
purification.

```wl
psi = Simplify[Sqrt[(1 + \[FormalR])/2] Flatten@KroneckerProduct[up, {1, 0}] + Sqrt[(1 - \[FormalR])/2] Flatten@KroneckerProduct[down, {0, 1}], Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

It is a unit vector, and contracting the ancilla index against itself returns the original state exactly.

```wl
Simplify[{Conjugate[psi] . psi, ArrayReshape[TensorContract[ArrayReshape[KroneckerProduct[psi, Conjugate[psi]], {2, 2, 2, 2}], {{2, 4}}], {2, 2}] == rho}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

Tracing the other way leaves the ancilla in a state carrying the same spectrum, so whatever mixedness $\rho$
has is shared symmetrically across the cut.

```wl
Simplify[Eigenvalues@Simplify[ArrayReshape[TensorContract[ArrayReshape[KroneckerProduct[psi, Conjugate[psi]], {2, 2, 2, 2}], {{1, 3}}], {2, 2}], Element[{\[FormalTheta], \[FormalPhi]}, Reals]], 0 <= \[FormalR] <= 1]
```

The choice of ancilla basis was free, so a unitary acting on $B$ alone sends one purification to another with
the same reduced state.

```wl
v = {{Cos[\[FormalA]/2] Exp[I \[FormalB]], Sin[\[FormalA]/2] Exp[I \[FormalC]]}, {-Sin[\[FormalA]/2] Exp[-I \[FormalC]], Cos[\[FormalA]/2] Exp[-I \[FormalB]]}}; Simplify[{ConjugateTranspose[v] . v == IdentityMatrix[2], ComplexExpand[ArrayReshape[TensorContract[ArrayReshape[KroneckerProduct[#, Conjugate[#]] &@(KroneckerProduct[IdentityMatrix[2], v] . psi), {2, 2, 2, 2}], {{2, 4}}], {2, 2}]] == rho}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi], \[FormalA], \[FormalB], \[FormalC]}, Reals]]
```

**QF** : the Bloch vector names the state directly.

```wl
qs = QuantumState["BlochVector"[\[FormalR] {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}]]
```

`"Purify"` adjoins the ancilla and returns a state on two qubits that reports itself pure.

```wl
{qs["Dimensions"], qs["Purify"]["Dimensions"], qs["Purify"]["PureStateQ"]}
```

`QuantumPartialTrace` over the ancilla recovers the original, and `"Unpurify"` is the same step named as the
inverse of `"Purify"`.

```wl
Simplify[{Normal@QuantumPartialTrace[qs["Purify"], {2}]["DensityMatrix"] == Normal@qs["DensityMatrix"], Normal@qs["Purify"]["Unpurify"]["DensityMatrix"] == Normal@qs["DensityMatrix"]}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

The generic mixed qubit has eigenvalues $(1\pm r)/2$ along $\pm\hat n$, so one ancilla qubit carries its
purification, and the amplitudes $\sqrt{(1\pm r)/2}$ are exactly the Schmidt coefficients of Part 11, 11.3.
Both routes return $\rho$ on the nose: the explicit vector traces back to the Bloch matrix and QF's `"Purify"`
inverts under `QuantumPartialTrace` and under `"Unpurify"`, all as symbolic identities in $r$, $\theta$ and
$\varphi$. Tracing the other way leaves the ancilla with the spectrum $\{(1-r)/2, (1+r)/2\}$, the same one
$\rho$ has, which is the pure-state statement that the two sides of a cut are equally mixed. The two
purifications are not the same vector, because the ancilla basis is arbitrary: a unitary on $B$ alone maps one
to the other and leaves the reduced state fixed, which the generic $V$ confirms. At $r = 1$ the smaller
eigenvalue vanishes, $|\Psi\rangle$ collapses to a product, and no ancilla is needed at all: purification buys
nothing for a state that was already pure, and buys the most at $r = 0$, where it returns a maximally entangled
pair.

### 12.3 [MSc] How do I show that different ensembles give the same density matrix, two of them doing so exactly when $|\tilde\varphi_j\rangle = \sum_i U_{ji}|\tilde\psi_i\rangle$ for an isometry $U$?

A density matrix records $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ and forgets the list that built it, so
the preparer's ensemble is not recoverable from the state. Absorbing each weight into its vector,
$|\tilde\psi_i\rangle = \sqrt{p_i}\,|\psi_i\rangle$, turns the sum into
$\rho = \sum_i |\tilde\psi_i\rangle\langle\tilde\psi_i|$, and in that form the freedom is a single line of
algebra. Given any $U$ with $U^\dagger U = I$, set $|\tilde\varphi_j\rangle = \sum_i U_{ji}|\tilde\psi_i\rangle$;
then

$$\sum_j |\tilde\varphi_j\rangle\langle\tilde\varphi_j| = \sum_{i,k} \Big(\sum_j U^*_{jk} U_{ji}\Big) |\tilde\psi_i\rangle\langle\tilde\psi_k| = \sum_i |\tilde\psi_i\rangle\langle\tilde\psi_i| = \rho,$$

and the Hughston-Jozsa-Wootters theorem says the converse too: every ensemble for $\rho$ arises this way, from
the eigen-ensemble, with $U$ an isometry rather than a square unitary when the new list is longer. Since $U$
need not be square, neither the weights, nor the states, nor even how many of them there are is fixed by
$\rho$. Take the generic mixed qubit $\rho = \tfrac12(I + r\,\hat n\cdot\vec\sigma)$ again.

**WL** : the state, and the ensemble its own spectrum provides, weights folded into the vectors.

```wl
rho = With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, 1/2 (IdentityMatrix[2] + \[FormalR] nvec . PauliMatrix[{1, 2, 3}])]; ens = {Sqrt[(1 + \[FormalR])/2] {Cos[\[FormalTheta]/2], Exp[I \[FormalPhi]] Sin[\[FormalTheta]/2]}, Sqrt[(1 - \[FormalR])/2] {-Exp[-I \[FormalPhi]] Sin[\[FormalTheta]/2], Cos[\[FormalTheta]/2]}}
```

Its members square-sum to $\rho$, and their norms are the eigenvalues, so this is the orthogonal ensemble the
spectral decomposition names.

```wl
Simplify[{Total[KroneckerProduct[#, Conjugate[#]] & /@ ens] == rho, Conjugate[#] . # & /@ ens}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

Mixing the two members with a generic $\mathrm{SU}(2)$ matrix produces another ensemble for the same $\rho$,
with its own weights and its own mutual overlap.

```wl
u = {{Cos[\[FormalA]/2] Exp[I \[FormalB]], Sin[\[FormalA]/2] Exp[I \[FormalC]]}, {-Sin[\[FormalA]/2] Exp[-I \[FormalC]], Cos[\[FormalA]/2] Exp[-I \[FormalB]]}}; Simplify[ComplexExpand[{Total[KroneckerProduct[#, Conjugate[#]] & /@ (u . ens)] == rho, Conjugate[#] . # & /@ (u . ens), Conjugate[First[u . ens]] . Last[u . ens]}], 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi], \[FormalA], \[FormalB], \[FormalC]}, Reals]]
```

The isometry need not be square, so the member count is not fixed either: four columns of length two, orthonormal
as columns, turn two states into four.

```wl
v = 1/2 {{1, 1}, {1, -1}, {1, I}, {1, -I}}; Simplify[ComplexExpand[{ConjugateTranspose[v] . v == IdentityMatrix[2], Total[KroneckerProduct[#, Conjugate[#]] & /@ (v . ens)] == rho, Conjugate[#] . # & /@ (v . ens)}], 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

**QF** : the freedom is not just algebraic. Stacking the ensemble members against orthonormal ancilla kets is
exactly the purification of Part 12, 12.2, so an observer holding the ancilla can *choose* which ensemble the
system is in, by choosing what to measure.

```wl
qpsi = QuantumState[Total[MapThread[Flatten@KroneckerProduct[#1, #2] &, {ens, IdentityMatrix[2]}]], {2, 2}]
```

Measuring the ancilla in the computational basis returns the eigen-ensemble weights; measuring it along $X$
returns an unbiased pair instead.

```wl
Simplify[{Values@QuantumMeasurementOperator[{2}][qpsi]["Probabilities"], Values@QuantumMeasurementOperator[QuantumOperator["X", {2}]][qpsi]["Probabilities"]}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

Either way the conditional states of the system, weighted by the outcome probabilities they come with, rebuild
the same $\rho$, which is what stops the choice from being detectable on the system alone.

```wl
Simplify[Total[Normal[QuantumPartialTrace[#, {2}]["DensityMatrix"]] & /@ Values[#["StateAssociation"]]] == rho & /@ {QuantumMeasurementOperator[{2}][qpsi], QuantumMeasurementOperator[QuantumOperator["X", {2}]][qpsi]}, 0 <= \[FormalR] <= 1 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

The spectral ensemble carries weights $(1\pm r)/2$ on the two orthogonal eigenvectors, and rotating it by a
generic $\mathrm{SU}(2)$ leaves $\rho$ untouched while moving the weights to $(1 \pm r\cos a)/2$: one free angle
already gives a continuum of ensembles, every one of them a legitimate answer to "how was this prepared". Their
two members overlap by $-\tfrac{r}{2}\sin a\,e^{-i(b+c)}$, nonzero away from $a = 0$, so a generic ensemble for
a mixed state is not orthogonal and its members are not distinguishable states. The rectangular isometry pushes
further, turning two states into four of equal weight $\tfrac14$ for the same $\rho$, so the number of
preparations is not fixed by the state either. What the framework adds is the mechanism: the purification of
12.2 makes the ancilla a control on the ensemble, and measuring it in the computational basis gives
probabilities $(1\pm r)/2$ while measuring along $X$ gives $\tfrac12$ and $\tfrac12$, two genuinely different
ensembles conjured on the system by a choice made elsewhere. Both rebuild $\rho$ exactly, which is no-signalling
(Part 11, 11.5) seen from the preparer's side: the remote choice steers which ensemble the system is in, and
that is precisely the information $\rho$ does not carry.

### 12.4 [BSc] How do I build the Gibbs (thermal) state $\rho = e^{-\beta \hat H}/Z$ with $Z = \mathrm{Tr}\,e^{-\beta \hat H}$, and compute thermal averages $\langle A\rangle = \mathrm{Tr}[A\rho]$?

A system in equilibrium with a bath at inverse temperature $\beta = 1/k_BT$ sits in the Gibbs state, the
exponential of the Hamiltonian normalized to unit trace. Everything thermodynamic follows from the one scalar
$Z$: the mean energy is $-\partial_\beta \log Z$, the free energy is $F = -\beta^{-1}\log Z$, and the entropy
closes the triangle as $S = \beta(\langle H\rangle - F)$. Take the general two-level system, splitting $\omega$
about an arbitrary axis, $\hat H = \tfrac{\omega}{2}\,\hat n\cdot\vec\sigma$ with $\hat n$ at angles
$(\theta,\varphi)$, so nothing depends on aligning the field with a coordinate axis.

**WL** : the Hamiltonian, written in the same Bloch form the states of 12.2 and 12.3 use.

```wl
ham = With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, \[FormalW]/2 nvec . PauliMatrix[{1, 2, 3}]]
```

Its two levels sit at $\pm\omega/2$, and the partition function is the trace of the matrix exponential.

```wl
{Simplify[Eigenvalues[ham], \[FormalW] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]], Simplify@ExpToTrig@Simplify[Tr@MatrixExp[-\[FormalB] ham], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]}
```

Dividing the exponential by that trace gives the state, and it is again a Bloch vector: pointing along
$-\hat n$, of length set by the temperature alone.

```wl
gibbs = Simplify[MatrixExp[-\[FormalB] ham]/Tr@MatrixExp[-\[FormalB] ham], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]; Simplify[gibbs == With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, 1/2 (IdentityMatrix[2] - Tanh[\[FormalB] \[FormalW]/2] nvec . PauliMatrix[{1, 2, 3}])], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

A thermal average is a trace against that state. Doing it for $\hat H$ itself, and differentiating $\log Z$,
are two routes to the same number.

```wl
Simplify@ExpToTrig@Simplify[{Tr[ham . gibbs], -D[Log[Tr@MatrixExp[-\[FormalB] ham]], \[FormalB]]}, \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

The entropy comes from the same spectrum, and the thermodynamic identity relating it to the energy and the
free energy is an algebraic consequence rather than a further assumption.

```wl
FullSimplify[{-Total[# Log[#] & /@ Eigenvalues[gibbs]] == Log[2 Cosh[\[FormalB] \[FormalW]/2]] - (\[FormalB] \[FormalW]/2) Tanh[\[FormalB] \[FormalW]/2], -Total[# Log[#] & /@ Eigenvalues[gibbs]] == \[FormalB] (Tr[ham . gibbs] + Log[Tr@MatrixExp[-\[FormalB] ham]]/\[FormalB])}, \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

Both temperature limits are limits of the single number $\tanh(\beta\omega/2)$, and the cold end lands on the
eigenprojector for the lower level.

```wl
With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, {Limit[Tanh[\[FormalB] \[FormalW]/2], \[FormalB] -> 0], Limit[Tanh[\[FormalB] \[FormalW]/2], \[FormalB] -> Infinity, Assumptions -> \[FormalW] > 0], Simplify[ham . (1/2 (IdentityMatrix[2] - nvec . PauliMatrix[{1, 2, 3}])) == -(\[FormalW]/2) (1/2 (IdentityMatrix[2] - nvec . PauliMatrix[{1, 2, 3}])), \[FormalW] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]}]
```

**QF** : the Hamiltonian is a combination of named Pauli operators rather than a matrix.

```wl
hq = With[{nvec = {Sin[\[FormalTheta]] Cos[\[FormalPhi]], Sin[\[FormalTheta]] Sin[\[FormalPhi]], Cos[\[FormalTheta]]}}, \[FormalW]/2 (nvec[[1]] QuantumOperator["X"] + nvec[[2]] QuantumOperator["Y"] + nvec[[3]] QuantumOperator["Z"])]
```

`MatrixExp` acts on the operator directly and returns another operator, so the Gibbs state is its matrix handed
to `QuantumState` and normalized. The state reports its own Bloch vector, which is where the temperature shows
up.

```wl
gq = QuantumState[MatrixExp[-\[FormalB] hq]["Matrix"], 2]["Normalized"]; Simplify@ExpToTrig@Simplify[gq["BlochVector"], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

The mean energy is the trace of the operator against the state.

```wl
Simplify@ExpToTrig@Simplify[Tr[Normal[hq["Matrix"]] . Normal[gq["DensityMatrix"]]], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]
```

`"VonNeumannEntropy"` is a `Quantity` carrying its unit, and it counts in bits, so it is the WL entropy divided
by $\log 2$.

```wl
{QuantityUnit[gq["VonNeumannEntropy"]], FullSimplify[QuantityMagnitude[gq["VonNeumannEntropy"]] Log[2] == Log[2 Cosh[\[FormalB] \[FormalW]/2]] - (\[FormalB] \[FormalW]/2) Tanh[\[FormalB] \[FormalW]/2], \[FormalW] > 0 && \[FormalB] > 0 && Element[{\[FormalTheta], \[FormalPhi]}, Reals]]}
```

The two-level Gibbs state is $\tfrac12\big(I - \tanh(\beta\omega/2)\,\hat n\cdot\vec\sigma\big)$, verified as an
identity for every axis and temperature, so the whole thermal family is the segment of the Bloch line along
$-\hat n$ of length $\tanh(\beta\omega/2)$, running from the centre at infinite temperature out toward the
surface as it cools. The minus sign is Boltzmann weighting: population accumulates on the lower level, opposite
the field. The partition function is $Z = 2\cosh(\beta\omega/2)$, and the mean energy
$\langle H\rangle = -\tfrac{\omega}{2}\tanh(\beta\omega/2)$ arrives identically from the trace $\mathrm{Tr}[H\rho]$
and from $-\partial_\beta \log Z$, which is what makes $Z$ a generating function rather than a normalization.
The entropy is $\log 2\cosh(\beta\omega/2) - \tfrac{\beta\omega}{2}\tanh(\beta\omega/2)$, running from
$\log 2$ at $\beta = 0$ down to $0$ as $\beta \to \infty$, and $S = \beta(\langle H\rangle - F)$ comes back
`True` rather than being imposed. The cold limit is the eigenprojector of the lower level, confirmed by
$H P = -\tfrac{\omega}{2}P$: a system at zero temperature is in its ground state, here as a computed
consequence of $\tanh \to 1$. The framework agrees on every count, with one unit convention to watch, since
`"VonNeumannEntropy"` returns `Quantity[..., "Bits"]` and matches the natural-log entropy only after
multiplication by $\log 2$.

### 12.5 [BSc] How do I compare two states by the trace distance $D(\rho,\sigma) = \tfrac12\|\rho-\sigma\|_1$ and, for pure states, the fidelity $F = |\langle\psi|\varphi\rangle|^2$?

Two states can be close in two different senses, and the two carry different meanings. The trace distance
$D = \tfrac12\|\rho-\sigma\|_1 = \tfrac12\sum_k|\lambda_k(\rho-\sigma)|$ is a metric, and it is the largest gap
in outcome probability that any single measurement can be made to show, which is why it governs the
discrimination problem of 12.7. The fidelity is an overlap rather than a distance, equal to $1$ for identical
states and $0$ for orthogonal ones. Both conventions need naming before any formula is quoted: here $F$ is the
*squared* overlap $|\langle\psi|\varphi\rangle|^2$, and the trace norm carries the factor $\tfrac12$ that puts
$D$ in $[0,1]$. For a qubit both quantities become geometry in the Bloch ball, where $D$ is half the Euclidean
separation of the Bloch vectors, so the general mixed state is the right thing to start from.

**WL** : two generic mixed qubits, each carrying its own Bloch vector.

```wl
{rho1, rho2} = 1/2 (IdentityMatrix[2] + # . PauliMatrix[{1, 2, 3}]) & /@ {{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}}
```

The difference of two states is traceless, so its two eigenvalues are equal and opposite and the trace norm
collapses to the length of the Bloch separation.

```wl
Simplify[Total[Abs[Eigenvalues[rho1 - rho2]]]/2 == EuclideanDistance[{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}]/2, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals]]
```

The fidelity needs purity, so drop to two generic pure qubits. Their squared overlap is fixed by the angle
between their Bloch vectors alone, the phases cancelling.

```wl
FullSimplify[ComplexExpand[Abs[Conjugate[{Cos[\[FormalT]/2], Exp[I \[FormalP]] Sin[\[FormalT]/2]}] . {Cos[\[FormalA]/2], Exp[I \[FormalB]] Sin[\[FormalA]/2]}]^2] == (1 + {Sin[\[FormalT]] Cos[\[FormalP]], Sin[\[FormalT]] Sin[\[FormalP]], Cos[\[FormalT]]} . {Sin[\[FormalA]] Cos[\[FormalB]], Sin[\[FormalA]] Sin[\[FormalB]], Cos[\[FormalA]]})/2, Element[{\[FormalT], \[FormalP], \[FormalA], \[FormalB]}, Reals]]
```

On the Bloch sphere the two measures are then locked together, since a chord and the angle it subtends carry
the same information.

```wl
Simplify[EuclideanDistance[{Sin[\[FormalT]] Cos[\[FormalP]], Sin[\[FormalT]] Sin[\[FormalP]], Cos[\[FormalT]]}, {Sin[\[FormalA]] Cos[\[FormalB]], Sin[\[FormalA]] Sin[\[FormalB]], Cos[\[FormalA]]}]/2 == Sqrt[1 - (1 + {Sin[\[FormalT]] Cos[\[FormalP]], Sin[\[FormalT]] Sin[\[FormalP]], Cos[\[FormalT]]} . {Sin[\[FormalA]] Cos[\[FormalB]], Sin[\[FormalA]] Sin[\[FormalB]], Cos[\[FormalA]]})/2], Element[{\[FormalT], \[FormalP], \[FormalA], \[FormalB]}, Reals]]
```

**QF** : `QuantumDistance` takes the measure as a third argument, where the trace distance is named `"Trace"`
and `"Bloch"` is the Bloch-vector separation. For a qubit those are the same number.

```wl
{q1, q2} = QuantumState["BlochVector"[#]] & /@ {{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}}; Simplify[{QuantumDistance[q1, q2, "Trace"] == QuantumDistance[q1, q2, "Bloch"], QuantumDistance[q1, q2, "Trace"] == EuclideanDistance[{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}]/2}, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals]]
```

For two pure qubits a common unitary rotates one to the north pole and the other into the $x$-$z$ plane, so a
single angle $\theta$ between the Bloch vectors carries everything. Three numbers, three different conventions.

```wl
{p1, p2} = QuantumState /@ {{1, 0}, {Cos[\[FormalT]/2], Sin[\[FormalT]/2]}}; FullSimplify[{QuantumSimilarity[p1, p2, "Fidelity"], QuantumDistance[p1, p2], QuantumDistance[p1, p2, "Trace"]}, 0 <= \[FormalT] <= Pi]
```

The first is $\sqrt F$ and the second is a distance built from it, not the fidelity, so the two are
complementary rather than equal, and only the third is the $D$ defined above.

```wl
FullSimplify[{QuantumDistance[p1, p2] + QuantumSimilarity[p1, p2, "Fidelity"] == 1, QuantumDistance[p1, p2, "Trace"] == Sqrt[1 - QuantumSimilarity[p1, p2, "Fidelity"]^2]}, 0 <= \[FormalT] <= Pi]
```

For qubits the trace distance is exactly half the Euclidean separation of Bloch vectors, proved here as an
identity in all six components rather than sampled, so the abstract $\tfrac12\|\rho-\sigma\|_1$ is a ruler laid
across the Bloch ball: $0$ at the centre against itself, $1$ between antipodal pure states, and never more,
since the ball has diameter $2$. The pure-state fidelity is $(1 + \hat n_1\cdot\hat n_2)/2 = \cos^2(\Theta/2)$
with $\Theta$ the angle between the Bloch vectors, so it depends on nothing else, and the two phases drop out
of the closed form. The two measures are then one quantity in two coordinates, $D = \sqrt{1-F}$, orthogonal
states sitting at $D = 1$, $F = 0$ and identical ones at $D = 0$, $F = 1$. In the framework three conventions
sit side by side and must not be confused: `QuantumDistance[q1, q2, "Trace"]` is $D = \sin(\theta/2)$ and
agrees with `"Bloch"` identically; `QuantumSimilarity[q1, q2, "Fidelity"]` is $\cos(\theta/2) = \sqrt F$, the
overlap and not its square; and the bare `QuantumDistance[q1, q2]` defaults to the `"Fidelity"` measure, which
returns $2\sin^2(\theta/4) = 1 - \sqrt F$, a *distance* derived from the fidelity rather than the fidelity
itself, the two summing to $1$. Reading that default as a fidelity inverts every comparison it appears in.

### 12.6 [MSc] How do I compute the Uhlmann fidelity $F(\rho,\sigma) = \big(\mathrm{Tr}\sqrt{\sqrt{\rho}\,\sigma\sqrt{\rho}}\big)^2$ between two mixed states?

The pure-state overlap of 12.5 has no meaning for mixed states, and the replacement is Uhlmann's: sandwich
$\sigma$ between square roots of $\rho$, take the square root of that, trace, and square. The nesting looks
forbidding but two observations remove it. First, $\sqrt\rho\,\sigma\sqrt\rho$ is similar to $\rho\sigma$, so
they share a spectrum and $\mathrm{Tr}\sqrt{\sqrt\rho\sigma\sqrt\rho} = \sum_k\sqrt{\lambda_k(\rho\sigma)}$,
with no matrix square root left. Second, for a qubit there are only two such eigenvalues, and a sum of two
square roots squares to something built from their sum and product alone, both of which are polynomial in the
states. The convention here is the *squared* trace, so $F = 1$ for identical states and $F$ reduces to
$|\langle\psi|\varphi\rangle|^2$ on pure ones.

**WL** : two generic mixed qubits again, each with its own Bloch vector.

```wl
{rho1, rho2} = 1/2 (IdentityMatrix[2] + # . PauliMatrix[{1, 2, 3}]) & /@ {{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}}
```

The two eigenvalues of the product are never needed individually: only their sum and their product enter, and
those are a trace and a pair of determinants.

```wl
Simplify[{Total[Eigenvalues[rho1 . rho2]] == Tr[rho1 . rho2], Times @@ Eigenvalues[rho1 . rho2] == Det[rho1] Det[rho2]}, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals]]
```

Squaring a sum of two square roots returns the sum plus twice the root of the product, which is the step that
turns the definition into a formula.

```wl
Simplify[(Sqrt[\[FormalL]] + Sqrt[\[FormalM]])^2 == \[FormalL] + \[FormalM] + 2 Sqrt[\[FormalL] \[FormalM]], \[FormalL] >= 0 && \[FormalM] >= 0]
```

So $F = \mathrm{Tr}[\rho\sigma] + 2\sqrt{\det\rho\det\sigma}$, and in Bloch coordinates every piece is
elementary: the trace is an inner product of Bloch vectors and each determinant measures how far its state sits
from the surface.

```wl
Simplify[Tr[rho1 . rho2] + 2 Sqrt[Det[rho1] Det[rho2]] == (1 + {\[FormalX], \[FormalY], \[FormalZ]} . {\[FormalU], \[FormalV], \[FormalW]} + Sqrt[(1 - {\[FormalX], \[FormalY], \[FormalZ]} . {\[FormalX], \[FormalY], \[FormalZ]}) (1 - {\[FormalU], \[FormalV], \[FormalW]} . {\[FormalU], \[FormalV], \[FormalW]})])/2, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals]]
```

That chain still has to be held against the definition itself, matrix square roots and all. Fidelity is
invariant under a common unitary, so put the first Bloch vector on the $z$ axis and rotate the second into the
$x$-$z$ plane: two lengths $a$, $b$ and one angle $c$ between them carry every pair of qubit states.

```wl
{sig1, sig2} = {1/2 (IdentityMatrix[2] + {0, 0, \[FormalA]} . PauliMatrix[{1, 2, 3}]), 1/2 (IdentityMatrix[2] + {\[FormalB] Sin[\[FormalC]], 0, \[FormalB] Cos[\[FormalC]]} . PauliMatrix[{1, 2, 3}])}; FullSimplify[{Tr[MatrixPower[MatrixPower[sig1, 1/2] . sig2 . MatrixPower[sig1, 1/2], 1/2]]^2 == (1 + \[FormalA] \[FormalB] Cos[\[FormalC]] + Sqrt[(1 - \[FormalA]^2) (1 - \[FormalB]^2)])/2, Tr[MatrixPower[MatrixPower[sig1, 1/2] . sig2 . MatrixPower[sig1, 1/2], 1/2]] == Tr[MatrixPower[sig1 . sig2, 1/2]]}, 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi]
```

Only $\rho$ carries the square roots, so the definition is lopsided in a way the quantity is not.

```wl
FullSimplify[Tr[MatrixPower[MatrixPower[sig1, 1/2] . sig2 . MatrixPower[sig1, 1/2], 1/2]] == Tr[MatrixPower[MatrixPower[sig2, 1/2] . sig1 . MatrixPower[sig2, 1/2], 1/2]], 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi]
```

Two limits pin the formula down: coincident states, and the pure-state case that has to reproduce 12.5.

```wl
{FullSimplify[(1 + \[FormalA] \[FormalB] Cos[\[FormalC]] + Sqrt[(1 - \[FormalA]^2) (1 - \[FormalB]^2)])/2 /. {\[FormalB] -> \[FormalA], \[FormalC] -> 0}, 0 <= \[FormalA] <= 1], FullSimplify[(1 + \[FormalA] \[FormalB] Cos[\[FormalC]] + Sqrt[(1 - \[FormalA]^2) (1 - \[FormalB]^2)])/2 /. {\[FormalA] -> 1, \[FormalB] -> 1}, 0 <= \[FormalC] <= Pi]}
```

**QF** : `QuantumSimilarity` with the `"Fidelity"` measure computes $\mathrm{Tr}[(\rho\sigma)^{1/2}]$, which is
$\sqrt F$ and not $F$, so squaring it is what matches the definition above.

```wl
{qs1, qs2} = QuantumState["BlochVector"[#]] & /@ {{0, 0, \[FormalA]}, {\[FormalB] Sin[\[FormalC]], 0, \[FormalB] Cos[\[FormalC]]}}; FullSimplify[QuantumSimilarity[qs1, qs2, "Fidelity"]^2 == (1 + \[FormalA] \[FormalB] Cos[\[FormalC]] + Sqrt[(1 - \[FormalA]^2) (1 - \[FormalB]^2)])/2, 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi]
```

The framework computes the product form, so holding it against the sandwiched definition is the check that the
two really are one quantity.

```wl
FullSimplify[QuantumSimilarity[qs1, qs2, "Fidelity"] == Tr[MatrixPower[MatrixPower[sig1, 1/2] . sig2 . MatrixPower[sig1, 1/2], 1/2]], 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi]
```

The nested square roots come apart because $\sqrt\rho\,\sigma\sqrt\rho$ and $\rho\sigma$ are similar, confirmed
here as an identity rather than quoted, so the trace of the root is the sum of the roots of the eigenvalues of
the plain product. For a qubit that sum has two terms, its square is
$\mathrm{Tr}[\rho\sigma] + 2\sqrt{\det\rho\det\sigma}$, and in Bloch coordinates

$$F = \tfrac12\Big(1 + \vec r_1\cdot\vec r_2 + \sqrt{(1-|\vec r_1|^2)(1-|\vec r_2|^2)}\Big),$$

an identity in all six components with no matrix functions anywhere. The definition computed directly, square
roots included, returns exactly that, which is the check that the shortcut is the same quantity. The second
term is what mixed states add: it vanishes as soon as either state is pure, and then $F$ collapses to
$(1 + \vec r_1\cdot\vec r_2)/2$, the $\cos^2(\Theta/2)$ of 12.5. Coincident states give $F = 1$ for any purity,
including two copies of the maximally mixed state, where the Bloch term vanishes and the determinant term
supplies the whole of it. Exchanging the two states leaves the sandwiched expression alone, so
$F(\rho,\sigma) = F(\sigma,\rho)$ despite a definition that privileges $\rho$. The framework agrees once the
convention is respected: `QuantumSimilarity` with `"Fidelity"` is $\sqrt F$, so its square is the Uhlmann
fidelity, and since the framework evaluates the product form $\mathrm{Tr}[(\rho\sigma)^{1/2}]$ rather than the
sandwich, matching it against $\mathrm{Tr}\sqrt{\sqrt\rho\sigma\sqrt\rho}$ is what certifies the two as one
quantity rather than two that happen to agree on the closed form.

### 12.7 [MSc] How do I find the optimal state-discrimination measurement and the Helstrom error probability $P_{\mathrm{err}} = \tfrac12\big(1 - \|p_1\rho_1 - p_2\rho_2\|_1\big)$?

One of two known states arrives, $\rho_1$ with prior $p$ and $\rho_2$ with prior $1-p$, and a single measurement
must name which. The apparatus is a two-outcome POVM $\{E_1, E_2\}$, outcome $E_2$ meaning "it was $\rho_2$", so

$$P_{\mathrm{err}} = p\,\mathrm{Tr}[\rho_1 E_2] + (1-p)\,\mathrm{Tr}[\rho_2 E_1] = (1-p) + \mathrm{Tr}[\Lambda E_2], \qquad \Lambda = p\rho_1 - (1-p)\rho_2 .$$

Only $\mathrm{Tr}[\Lambda E_2]$ is free, and over $0 \le E_2 \le I$ it is made as negative as possible by putting
$E_2$ on the negative eigenspace of the Helstrom matrix $\Lambda$ and nowhere else. So the optimal measurement
*is* the pair of spectral projectors of $\Lambda$, and what it leaves behind is the sum of $\Lambda$'s negative
eigenvalues, which is $\tfrac12(1 - \|\Lambda\|_1)$ once $\mathrm{Tr}\,\Lambda = 2p-1$ is used. For a qubit
$\Lambda$ has one eigenvalue of each sign, so the measurement is a single projector and the entire answer is one
direction in the Bloch ball.

**WL** : two generic mixed qubits, an arbitrary prior, and the matrix whose spectrum settles it.

```wl
{rho1, rho2} = 1/2 (IdentityMatrix[2] + # . PauliMatrix[{1, 2, 3}]) & /@ {{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}}; lam = \[FormalP] rho1 - (1 - \[FormalP]) rho2
```

$\Lambda$ is a difference of states, so its Bloch part is the weighted difference
$\vec d = p\,\vec r_1 - (1-p)\,\vec r_2$, and the two projectors along $\pm\hat d$ are exactly its
eigenprojectors. That is the measurement: complementary, and diagonalizing $\Lambda$.

```wl
dv = \[FormalP] {\[FormalX], \[FormalY], \[FormalZ]} - (1 - \[FormalP]) {\[FormalU], \[FormalV], \[FormalW]}; {ep, em} = 1/2 (IdentityMatrix[2] + # Normalize[dv] . PauliMatrix[{1, 2, 3}]) & /@ {1, -1}; Simplify[{ep . ep == ep, em . em == em, ep + em == IdentityMatrix[2], lam . ep == ((2 \[FormalP] - 1 + Sqrt[dv . dv])/2) ep, lam . em == ((2 \[FormalP] - 1 - Sqrt[dv . dv])/2) em}, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals] && 0 <= \[FormalP] <= 1]
```

Guessing $\rho_2$ on the negative projector and $\rho_1$ on the positive one, the error is the sum of
$\Lambda$'s negative part, and the trace norm turns that into the Helstrom value.

```wl
{Simplify[\[FormalP] Tr[rho1 . em] + (1 - \[FormalP]) Tr[rho2 . ep] == (1 - Sqrt[dv . dv])/2, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals] && 0 <= \[FormalP] <= 1], Resolve[ForAll[{\[FormalC], \[FormalD]}, \[FormalD] >= Abs[\[FormalC]], Abs[(\[FormalC] + \[FormalD])/2] + Abs[(\[FormalC] - \[FormalD])/2] == \[FormalD]], Reals]}
```

Nothing else does better, and the reason is visible rather than argued: a projector along any direction $\hat m$
errs by an amount depending on $\hat m$ through a single inner product, largest when $\hat m$ is aligned with
$\vec d$.

```wl
Simplify[\[FormalP] Tr[rho1 . (IdentityMatrix[2] - 1/2 (IdentityMatrix[2] + {\[FormalM], \[FormalN], \[FormalO]} . PauliMatrix[{1, 2, 3}]))] + (1 - \[FormalP]) Tr[rho2 . (1/2 (IdentityMatrix[2] + {\[FormalM], \[FormalN], \[FormalO]} . PauliMatrix[{1, 2, 3}]))] == (1 - dv . {\[FormalM], \[FormalN], \[FormalO]})/2, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW], \[FormalM], \[FormalN], \[FormalO]}, Reals] && 0 <= \[FormalP] <= 1]
```

```wl
Reduce[ForAll[\[FormalT], 0 <= \[FormalT] <= Pi, (1 - \[FormalD] Cos[\[FormalT]])/2 >= (1 - \[FormalD])/2] && 0 <= \[FormalD] <= 1, \[FormalD], Reals]
```

At equal priors the length $|\vec d|$ is the trace distance of 12.5, which is what gives that metric its
operational meaning.

```wl
Simplify[(Sqrt[dv . dv] /. \[FormalP] -> 1/2) == EuclideanDistance[{\[FormalX], \[FormalY], \[FormalZ]}, {\[FormalU], \[FormalV], \[FormalW]}]/2, Element[{\[FormalX], \[FormalY], \[FormalZ], \[FormalU], \[FormalV], \[FormalW]}, Reals]]
```

**QF** : the same optimal effect, scored against framework states in the gauge of 12.6, reproduces the bound
written with `QuantumDistance` instead of a norm.

```wl
{qs1, qs2} = QuantumState["BlochVector"[#]] & /@ {{0, 0, \[FormalA]}, {\[FormalB] Sin[\[FormalC]], 0, \[FormalB] Cos[\[FormalC]]}}; eff = QuantumOperator[1/2 (IdentityMatrix[2] + Simplify[Normalize[{0, 0, \[FormalA]} - {\[FormalB] Sin[\[FormalC]], 0, \[FormalB] Cos[\[FormalC]]}], 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi] . PauliMatrix[{1, 2, 3}])]; FullSimplify[(Tr[Normal[qs1["DensityMatrix"]] . (IdentityMatrix[2] - Normal[eff["Matrix"]])] + Tr[Normal[qs2["DensityMatrix"]] . Normal[eff["Matrix"]]])/2 == (1 - EuclideanDistance[{0, 0, \[FormalA]}, {\[FormalB] Sin[\[FormalC]], 0, \[FormalB] Cos[\[FormalC]]}]/2)/2, 0 <= \[FormalA] <= 1 && 0 <= \[FormalB] <= 1 && 0 <= \[FormalC] <= Pi]
```

Run as a measurement on a definite pair, the two outcome probabilities and the error they imply come out against
the bound.

```wl
With[{n1 = QuantumState["BlochVector"[{0, 0, 1}]], n2 = QuantumState["BlochVector"[{Sin[Pi/3], 0, Cos[Pi/3]}]]}, With[{nm = QuantumMeasurementOperator[QuantumOperator[Normalize[{0, 0, 1} - {Sin[Pi/3], 0, Cos[Pi/3]}] . PauliMatrix[{1, 2, 3}]]]}, Chop@N@{Values@nm[n1]["Probabilities"], (Last@Values@nm[n1]["Probabilities"] + First@Values@nm[n2]["Probabilities"])/2, (1 - QuantumDistance[n1, n2, "Trace"])/2}]]
```

The optimal apparatus is the spectral measurement of $\Lambda = p\rho_1 - (1-p)\rho_2$, and for a qubit that is
the projector pair along $\pm\hat d$ with $\vec d = p\vec r_1 - (1-p)\vec r_2$: idempotent, complementary, and
returning $\Lambda$'s eigenvalues $\big((2p-1)\pm|\vec d|\big)/2$, all five checks `True`. Scoring it gives
$P_{\mathrm{err}} = \tfrac12(1 - |\vec d|)$, and since the two eigenvalues carry opposite signs their absolute
values sum to $|\vec d| = \|\Lambda\|_1$, which is the Helstrom formula reached without quoting it. A projector
along any other direction errs by $\big(1 - \vec d\cdot\hat m\big)/2$, an identity in all nine components, so
alignment with $\vec d$ is optimal and `Reduce` returns the whole range rather than an isolated case. At equal
priors $|\vec d|$ is exactly the trace distance $D$ of 12.5, so $P_{\mathrm{err}} = (1-D)/2$: the trace distance
is the margin by which the best measurement beats a coin flip, $D = 1$ giving certainty and $D = 0$ leaving pure
guesswork. The framework reproduces the bound in its own vocabulary, and on $|0\rangle$ against a state
$60^\circ$ away the measurement returns $\tfrac34$ and $\tfrac14$ with error $\tfrac14$, matching
$(1 - D)/2$ exactly.

### 12.8a [MSc] How do I verify Gleason's frame-function constraint for a qutrit, that every frame function $f$ with $\sum_i f(e_i) = 1$ over each orthonormal basis is $f(\psi) = \langle\psi|\rho|\psi\rangle$ for some density operator?

A frame function assigns a number $f(\psi) \in [0,1]$ to every unit vector, subject to one rule: over any
orthonormal basis the values sum to $1$. No linearity is assumed, only that constraint, repeated for every
basis at once. Gleason's theorem says that for $d \ge 3$ this forces $f(\psi) = \langle\psi|\rho|\psi\rangle$
for a unique density operator, so the Born rule is not an extra postulate but the only consistent way to put
probabilities on subspaces. One direction is a short calculation and is done below; the converse is the theorem
and is quoted. What can be settled here besides the easy direction is that the correspondence is one-to-one and
already saturated by finitely many bases: in $d = 3$ four mutually unbiased bases give twelve numbers, and those
pin down all nine real parameters of $\rho$.

**WL** : a generic qutrit density operator, Hermitian with unit trace by construction.

```wl
rho = {{\[FormalA], \[FormalB] + I \[FormalC], \[FormalD] + I \[FormalE]}, {\[FormalB] - I \[FormalC], \[FormalF], \[FormalG] + I \[FormalH]}, {\[FormalD] - I \[FormalE], \[FormalG] - I \[FormalH], 1 - \[FormalA] - \[FormalF]}}; Simplify[{rho == ConjugateTranspose[rho], Tr[rho]}, Element[{\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}, Reals]]
```

The Born values sum over a basis to the trace of $\rho$ against $\sum_i|e_i\rangle\langle e_i|$, for *any* three
vectors, orthonormal or not. So the frame condition is nothing but completeness of the basis, and that is the
whole of the easy direction.

```wl
With[{ev = Array[Subscript[\[FormalW], ##] &, {3, 3}]}, Simplify[Total[Table[Conjugate[ev[[i]]] . rho . ev[[i]], {i, 3}]] == Tr[rho . Total[Table[KroneckerProduct[ev[[i]], Conjugate[ev[[i]]]], {i, 3}]]]]]
```

Four mutually unbiased bases exist in dimension three, built from powers of a cube root of unity alongside the
computational basis. Each is orthonormal and complete, and vectors drawn from different bases overlap with the
same modulus, which is what makes the four of them together informationally complete.

```wl
mubs = With[{w = Exp[2 Pi I/3]}, Join[{IdentityMatrix[3]}, Table[Table[Table[w^(k m^2 + j m)/Sqrt[3], {m, 0, 2}], {j, 0, 2}], {k, 0, 2}]]]; Simplify[{Table[Conjugate[b] . Transpose[b] == IdentityMatrix[3], {b, mubs}], Table[Total[KroneckerProduct[#, Conjugate[#]] & /@ b] == IdentityMatrix[3], {b, mubs}], Union@Flatten@Table[Abs[Conjugate[mubs[[i, r]]] . mubs[[j, s]]]^2, {i, 4}, {j, i + 1, 4}, {r, 3}, {s, 3}]}]
```

The frame function on those twelve vectors is affine in the eight parameters, and the coefficient matrix has
full column rank, so the values determine the state and no two density operators can share a frame function.

```wl
With[{pars = {\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}, fv = Table[Conjugate[v] . rho . v, {v, Flatten[mubs, 1]}]}, With[{mat = Simplify[Outer[D, fv, pars], Element[pars, Reals]]}, {Dimensions[mat], MatrixRank[mat], Simplify[fv == mat . pars + (fv /. Thread[pars -> 0]), Element[pars, Reals]], Simplify[LinearSolve[mat, fv - (fv /. Thread[pars -> 0])] == pars, Element[pars, Reals]]}]]
```

That settles injectivity. The Gleason direction, that nothing *but* a density operator produces a frame
function, is a counting statement once the vectors are a fixed finite set. On the twelve vectors above the
constraint is four equations, one per basis, so the assignments satisfying it form an affine space of dimension
$12 - 4 = 8$. The Born assignments already live in that space and their linear part has rank $8$, so the two
coincide and there is no room left for anything else.

```wl
With[{vecs = Flatten[mubs, 1], pars = {\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}}, With[{cons = Table[Boole[MemberQ[b, v]], {b, mubs}, {v, vecs}], fv = Table[Conjugate[v] . rho . v, {v, vecs}]}, With[{mat = Simplify[Outer[D, fv, pars], Element[pars, Reals]]}, {Dimensions[cons], MatrixRank[cons], Length@NullSpace[cons], MatrixRank[mat], Simplify[cons . mat == ConstantArray[0, {4, 8}]], Simplify[cons . (fv /. Thread[pars -> 0])]}]]]
```

**QF** : the same generic state as an object, carrying the same eight parameters.

```wl
qt = QuantumState[rho, 3]
```

Each basis is a `QuantumBasis` built from its vectors, and re-expressing the state there is the whole
computation: the frame function's values on that basis sit on the diagonal, and the frame condition is `Tr` of
the re-expressed state.

```wl
Simplify[Table[{Diagonal[Normal@QuantumState[qt, QuantumBasis[b]]["DensityMatrix"]], Tr[QuantumState[qt, QuantumBasis[b]]]}, {b, mubs}], Element[{\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}, Reals]]
```

The diagonal is what to read for a mixed state, and $|\cdot|^2$ on `"AmplitudesList"` is its counterpart for a
pure one. `"Probabilities"` is a different quantity here, since it divides by the total absolute value and so
reports entries summing to one whatever it is handed, which on a parametrized state means dividing by
$|a| + |f| + |1 - a - f|$ rather than by the trace.

Collecting those twelve values, the framework route reaches the same correspondence on its own: full column rank
recovers the state from them, and against the four basis constraints the counting closes as before.

```wl
fvQF = Simplify[Flatten@Table[Diagonal[Normal@QuantumState[qt, QuantumBasis[b]]["DensityMatrix"]], {b, mubs}], Element[{\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}, Reals]]; With[{pars = {\[FormalA], \[FormalB], \[FormalC], \[FormalD], \[FormalE], \[FormalF], \[FormalG], \[FormalH]}, cons = Table[Boole[MemberQ[b, v]], {b, mubs}, {v, Flatten[mubs, 1]}]}, With[{matQF = Simplify[Outer[D, fvQF, pars], Element[pars, Reals]]}, {MatrixRank[matQF], Simplify[LinearSolve[matQF, fvQF - (fvQF /. Thread[pars -> 0])] == pars, Element[pars, Reals]], Length@NullSpace[cons], Simplify[cons . matQF == ConstantArray[0, {4, 8}]]}]]
```

The frame condition costs nothing to satisfy from a density operator: the sum over a basis equals
$\mathrm{Tr}[\rho\sum_i|e_i\rangle\langle e_i|]$ as an identity in nine arbitrary vector components, so once the
basis is complete the sum is $\mathrm{Tr}\rho = 1$ and orthonormality enters at exactly one point. The reverse
direction is Gleason's theorem, but its uniqueness half is a finite computation: the four mutually unbiased
bases come back orthonormal, complete, and pairwise overlapping at $\tfrac13$, and the twelve frame values they
produce are affine in the eight parameters with a $12\times8$ coefficient matrix of rank $8$, so `LinearSolve`
inverts it and returns the parameters exactly: the map from density operators to frame functions is injective.
The other direction closes on that configuration by counting. The frame condition there is a $4\times12$ system
of rank $4$, leaving an affine solution space of dimension $8$; the Born assignments satisfy it, since
`cons . mat` vanishes and the offset gives $\{1,1,1,1\}$, and their linear part has rank $8$. An $8$-dimensional
family inside an $8$-dimensional space is the whole of it, so on these twelve vectors every frame function is
$\langle v|\rho|v\rangle$ and the correspondence is exact in both directions. Two things are left to the theorem
rather than the computation: the count shows every solution is $\langle v|H|v\rangle$ for a Hermitian $H$ of
unit trace, and it is the requirement $f \ge 0$ that makes $H$ positive; and Gleason quantifies over every basis
in the continuum, not a chosen twelve. The framework needs no separate calculation for any of this: carrying the
same generic state into each `QuantumBasis` puts the twelve frame values on the diagonal, and `Tr` of the
re-expressed state is the frame condition, returning $1$ in all four bases, with the diagonal rather than
`"Probabilities"` carrying the values a parametrized state has. Those twelve values then carry the same
correspondence without borrowing the WL construction: rank $8$ again, `LinearSolve` again returning the
parameters, and the same $8$-dimensional match against the four constraints.

### 12.8b [MSc] Why does the qubit ($d = 2$) fail Gleason, and what does the POVM (Busch) version restore?

Gleason's rigidity is geometric, and the Bloch picture shows where it comes from. A measurement's outcomes are
points on the surface of the generalized Bloch ball, and the frame condition says the probabilities attached to
them sum to one. For a qubit a projective measurement is a pair of *antipodal* points, so a rule depending on
the state only through the scalar product $\vec r\cdot\hat n$ faces one equation, relating $\hat n$ to $-\hat n$
and nothing else. Every odd deviation from linearity survives that, continuity included, so there are
infinitely many non-Born rules. From $d = 3$ upward the eigenprojectors of a measurement are the vertices of a
regular simplex inscribed in the ball, $N$ of them at once, and requiring $N$ probabilities to sum to one at
every simplex is a far stronger functional constraint: it forces linearity in $\vec r\cdot\hat n$, which is the
Born rule. Busch's theorem repairs the qubit by admitting POVMs, whose outcomes need not come in antipodal
pairs, so a qubit inherits the same three-outcome constraint that a qutrit has for free. Below, a rule is a
function $g$ applied to outcome probabilities, and the question is which $g$ survives additivity.

**WL** : take the smoothstep $g(s) = 3s^2 - 2s^3$. It respects the two-outcome constraint exactly, and it is a
legitimate probability assignment: inside $[0,1]$ and increasing, so nothing pathological is being smuggled in.

```wl
g[\[FormalS]_] := 3 \[FormalS]^2 - 2 \[FormalS]^3; {Simplify[g[\[FormalP]] + g[1 - \[FormalP]] == 1], Reduce[ForAll[\[FormalP], 0 <= \[FormalP] <= 1, 0 <= g[\[FormalP]] <= 1], Reals], Reduce[ForAll[\[FormalP], 0 <= \[FormalP] <= 1, (D[g[\[FormalS]], \[FormalS]] /. \[FormalS] -> \[FormalP]) >= 0], Reals]}
```

It is nevertheless not the Born rule, agreeing with it at three isolated points and nowhere else.

```wl
Reduce[g[\[FormalP]] == \[FormalP] && 0 <= \[FormalP] <= 1, \[FormalP], Reals]
```

Three outcomes destroy it. Additivity over a triple summing to one holds only where one of the probabilities
vanishes, that is on the boundary, and at the centre it returns $\tfrac79$ instead of $1$.

```wl
{Simplify[g[\[FormalU]] + g[\[FormalV]] + g[1 - \[FormalU] - \[FormalV]] == 1], 3 g[1/3]}
```

What three outcomes leave standing is exactly the affine rules, and an impossible outcome must be assigned zero,
which fixes the remaining freedom and returns the Born rule.

```wl
{Reduce[ForAll[{\[FormalU], \[FormalV]}, (\[FormalA] + \[FormalB] \[FormalU]) + (\[FormalA] + \[FormalB] \[FormalV]) + (\[FormalA] + \[FormalB] (1 - \[FormalU] - \[FormalV])) == 1], {\[FormalA], \[FormalB]}, Reals], Reduce[\[FormalB] == 1 - 3 \[FormalA] && \[FormalA] == 0, {\[FormalA], \[FormalB]}, Reals]}
```

**QF** : a qubit can be given three outcomes without leaving dimension two, by measuring a POVM instead of a
basis. Three coplanar directions at $120^\circ$ furnish one, and the framework recognizes it as such.

```wl
trine = Table[{Sin[2 Pi k/3], 0, Cos[2 Pi k/3]}, {k, 0, 2}]; eff = (2/3) (1/2 (IdentityMatrix[2] + # . PauliMatrix[{1, 2, 3}])) & /@ trine; {QuantumMeasurementOperator[eff]["POVMQ"], Simplify@Total@Values@QuantumMeasurementOperator[eff][QuantumState["BlochVector"[{\[FormalX], \[FormalY], \[FormalZ]}]]]["Probabilities"]}
```

On the maximally mixed qubit its three outcomes are equally likely, the Born values sum to one, and the
smoothstep rule does not.

```wl
With[{mm = QuantumMeasurementOperator[eff][QuantumState[IdentityMatrix[2]/2, 2]]["Probabilities"]}, Simplify@{Values[mm], Total@Values[mm], Total[g /@ Values[mm]]}]
```

A qutrit needs no POVM for this: an ordinary basis already gives three outcomes, and re-expressing the
maximally mixed state in one produces the same three probabilities and the same failure.

```wl
With[{b = With[{w = Exp[2 Pi I/3]}, Table[Table[w^(m^2 + j m)/Sqrt[3], {m, 0, 2}], {j, 0, 2}]]}, With[{diag = Diagonal[Normal@QuantumState[QuantumState[IdentityMatrix[3]/3, 3], QuantumBasis[b]]["DensityMatrix"]]}, Simplify@{diag, Total[diag], Total[g /@ diag]}]]
```

The smoothstep is a genuine frame rule for a qubit measured projectively: $g(p) + g(1-p) = 1$ identically, with
values in $[0,1]$ and a nonnegative derivative throughout, so continuity and monotonicity are no obstacle to
being non-Born. It coincides with $p$ only at $0$, $\tfrac12$ and $1$, so it is a different rule, and at
$d = 2$ nothing rules it out: an antipodal pair supplies one equation, and one equation cannot pin down a
function. Three outcomes can. Additivity over a triple summing to one holds only on the boundary
$uv(u+v-1) = 0$, and at the centre the rule returns $\tfrac79$, missing by $\tfrac29$; among affine rules
$a + bs$ the triple constraint forces $b = 1-3a$, and assigning zero to an impossible outcome sets $a = 0$ and
$b = 1$, leaving $g(s) = s$ alone. That is Gleason's mechanism in one line of algebra, and it explains the
dimensional divide: a qutrit basis has three vertices, so the constraint applies already to projective
measurements, while a qubit basis has two and needs a POVM to supply a third. Both routes exhibit the same
arithmetic, the trine on a maximally mixed qubit and a mutually unbiased basis on a maximally mixed qutrit each
giving probabilities $\{\tfrac13,\tfrac13,\tfrac13\}$, Born total $1$, smoothstep total $\tfrac79$. Busch's
theorem is the statement that widening the constraint this way leaves only $\mathrm{Tr}[\rho E]$, so the Born
rule returns in dimension two from POVMs rather than from bases.
