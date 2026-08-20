## Part 14. Quantum information and entropy

Entropy is one functional used twice. Given a probability distribution it is Shannon's, and given a state it
is von Neumann's, and the second is the first applied to the spectrum, so the whole Part is an exercise in
asking which distribution a number refers to. From that single quantity the rest follows: the relative entropy
that measures how badly one state impersonates another, the joint and conditional and mutual entropies whose
classical inequalities survive only in part (a quantum conditional entropy can go negative), the coherent
information that a channel leaves intact, the data-processing inequality that forbids processing from creating
information, the Holevo bound on what an ensemble can carry, the capacities assembled from those pieces, and
the coherence of a state read against a fixed basis.

### 14.1 [MSc] How do I compute and compare the Shannon and von Neumann entropies?

There is only one formula. For a probability distribution $p$ and for a state $\rho$,

$$H(p) = -\sum_i p_i \log_2 p_i, \qquad S(\rho) = -\mathrm{Tr}\,\rho \log_2 \rho = H(\mathrm{spec}\,\rho),$$

with the convention $0\log 0 = 0$, which is not a technicality: it is what makes a pure state have entropy
exactly $0$ rather than an indeterminate expression. Reading $S$ as $H$ of the eigenvalue list is the whole
content of the second equality, since a matrix function acts on eigenvalues and leaves eigenvectors alone.

So "comparing" the two is only a question once a distribution is named, and the honest comparison is against
the one an experiment actually produces. Measuring $\rho$ in an orthonormal basis $\{|i\rangle\}$ gives
outcome probabilities $p_i = \langle i|\rho|i\rangle$, the diagonal of $\rho$ in that basis, and

$$H(p) \ge S(\rho),$$

with equality exactly when the basis is an eigenbasis of $\rho$. The diagonal is a doubly stochastic image of
the spectrum, so it is majorized by it, and $H$ is Schur concave: measuring in the wrong basis discards the
coherences and can only spread the outcomes further. Written on the standing mixed qubit
$\rho = \tfrac12(I + \vec r\cdot\vec\sigma)$ the two sides become one function of one variable each. The
spectrum is $(1 \pm |\vec r|)/2$ and the $Z$ outcomes are $(1 \pm r_z)/2$, so with
$h(x) = -x\log_2 x - (1-x)\log_2(1-x)$,

$$S(\rho) = h\!\left(\frac{1+|\vec r|}{2}\right), \qquad H(p_Z) = h\!\left(\frac{1+r_z}{2}\right),$$

and since $h(\tfrac{1+r}{2})$ decreases as the radius $r$ grows while $|r_z| \le |\vec r|$ always, the
inequality is exactly the statement that the measured radius cannot exceed the true one.

**WL** : the state, and the one functional both entropies are.

```wl
rho = 1/2 (IdentityMatrix[2] + {rx, ry, rz} . PauliMatrix[{1, 2, 3}])
```

```wl
shannon[p_] := -# . Log[2, #] & @ DeleteCases[p, 0]
```

The spectrum depends on the Bloch vector only through its length.

```wl
Simplify[Eigenvalues[rho], Element[{rx, ry, rz}, Reals]]
```

The von Neumann entropy computed from its definition as a matrix function, and the claim that it is the
Shannon entropy of that spectrum.

```wl
Simplify[-Tr[rho . MatrixLog[rho]]/Log[2], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

```wl
Simplify[-Tr[rho . MatrixLog[rho]]/Log[2] == shannon[Eigenvalues[rho]], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

A measurement in the computational basis reads the diagonal, and its Shannon entropy is a different number
built from the same formula.

```wl
{Diagonal[rho], shannon[Diagonal[rho]]}
```

The comparison needs no search. The binary entropy falls as the radius grows, and the measured radius $|r_z|$
never exceeds the true radius $|\vec r|$, so the measured distribution is the more spread of the two.

```wl
{Simplify[D[shannon[{(1 + r)/2, (1 - r)/2}], r] <= 0, 0 <= r < 1], Reduce[ForAll[{rx, ry, rz}, Element[{rx, ry, rz}, Reals], Abs[rz] <= Sqrt[rx^2 + ry^2 + rz^2]]]}
```

Equality is not generic: it holds exactly on the $z$ axis, where the computational basis diagonalizes $\rho$.

```wl
Simplify[(shannon[Diagonal[rho]] - shannon[Eigenvalues[rho]]) /. {rx -> 0, ry -> 0}, -1 < rz < 1]
```

Off that axis the gap is a positive number of bits, the coherence the measurement threw away.

```wl
N[(shannon[Diagonal[rho]] - shannon[Eigenvalues[rho]]) /. {rx -> 1/3, ry -> 1/4, rz -> 1/5}]
```

Two pure states show that the entropies answer different questions. Both are noiseless, so $S = 0$; the one
aligned with the measurement axis is certain, the one across it is a fair coin. Substituting into the state
before applying the functional is what lets the $0\log 0 = 0$ convention act.

```wl
{shannon[Eigenvalues[rho] /. {rx -> 1, ry -> 0, rz -> 0}], shannon[Diagonal[rho] /. {rx -> 1, ry -> 0, rz -> 0}], shannon[Diagonal[rho] /. {rx -> 0, ry -> 0, rz -> 1}]}
```

**QF** : both entropies are properties, and both carry a unit. The bare property returns a `Quantity` in bits,
and a numeric second argument selects the base, so `2` gives bits and `E` gives nats.

```wl
QuantumState[rho]["VonNeumannEntropy"]
```

```wl
Simplify[QuantumState[rho]["VonNeumannEntropy", 2] == shannon[Eigenvalues[rho]], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

The Shannon side belongs to the measurement rather than to the state: a `QuantumMeasurement` reports the
entropy of its own outcome distribution, which is the diagonal the $Z$ observable exposes.

```wl
Simplify[QuantumMeasurementOperator[PauliMatrix[3]][QuantumState[rho]]["Entropy", 2] == shannon[Diagonal[rho]], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

On the $z$ axis the observable and the state share an eigenbasis, and the two entropies coincide as objects,
not merely numerically.

```wl
With[{diagonal = rho /. {rx -> 0, ry -> 0}}, Simplify[QuantumState[diagonal]["VonNeumannEntropy", 2] == QuantumMeasurementOperator[PauliMatrix[3]][QuantumState[diagonal]]["Entropy", 2], -1 < rz < 1]]
```

The pure equatorial state is the sharpest case: no mixedness at all, one full bit of measurement randomness.

```wl
{QuantumState[{1, 1}/Sqrt[2]]["VonNeumannEntropy", 2], QuantumMeasurementOperator[PauliMatrix[3]][QuantumState[{1, 1}/Sqrt[2]]]["Entropy", 2]}
```

The maximum is set by dimension, $S \le \log_2 d$, reached only by the maximally mixed state.

```wl
{QuantumState[IdentityMatrix[3]/3, QuditBasis[3]]["VonNeumannEntropy", 2], Log[2, 3]}
```

One functional, two arguments. Applied to the spectrum of $\rho$ it is the von Neumann entropy and measures
how mixed the state is, a basis independent number that vanishes on every pure state. Applied to the outcome
probabilities of a measurement it is the Shannon entropy and measures how unpredictable that particular
experiment is, a basis dependent number that vanishes only when one outcome is certain. The two meet,
$H(p) = S(\rho)$, exactly when the measurement is diagonal in the state, which is why the eigenbasis is the
one place a quantum state behaves like a classical distribution; everywhere else $H(p) - S(\rho) > 0$ is the
coherence that the chosen measurement failed to see. For the qubit above that gap is
$h(\tfrac{1+r_z}{2}) - h(\tfrac{1+|\vec r|}{2})$, which is $0.131$ bits at $\vec r = (\tfrac13, \tfrac14, \tfrac15)$
and a full bit for a pure state measured across its own axis.

### 14.2 [MSc] How do I compute the quantum relative entropy and check its nonnegativity?

The relative entropy asks how badly one state impersonates another,

$$S(\rho\,\|\,\sigma) = \mathrm{Tr}\,\rho(\log_2\rho - \log_2\sigma) = -S(\rho) - \mathrm{Tr}\,\rho\log_2\sigma,$$

and Klein's inequality says the answer is never negative, $S(\rho\,\|\,\sigma) \ge 0$, with equality only at
$\rho = \sigma$. That single fact is the engine of the rest of the Part: the entropy bound of 14.1, the
data-processing inequality of 14.5 and the Holevo bound of 14.6 are all this inequality applied to a
particular pair.

Two properties keep it from being a distance. It is not symmetric, so $S(\rho\,\|\,\sigma)$ and
$S(\sigma\,\|\,\rho)$ are different numbers and one must say which is meant, and it is not even finite in
general: it is finite exactly when $\mathrm{supp}\,\rho \subseteq \mathrm{supp}\,\sigma$, and $+\infty$
otherwise. The infinity is physical rather than pathological. If $\sigma$ assigns probability zero to an
outcome that $\rho$ allows, then a single measurement can rule $\sigma$ out with certainty, so no finite
amount of evidence is needed to separate them.

Two special cases anchor the general statement. Against the maximally mixed state the definition collapses to
the entropy bound,

$$S\!\left(\rho\,\middle\|\,\frac{I}{d}\right) = \log_2 d - S(\rho),$$

so nonnegativity there is exactly the statement $S(\rho) \le \log_2 d$ from 14.1. And when $\rho$ and $\sigma$
commute, both are diagonal in one basis and the quantity reduces to the classical Kullback-Leibler divergence
$\sum_i p_i \log_2(p_i/q_i)$ of their spectra, so the quantum statement contains the classical one rather than
replacing it.

**WL** : a Bloch state builder, and the definition read straight off the formula.

```wl
bloch[v_] := 1/2 (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])
```

```wl
relativeEntropy[a_, b_] := Tr[a . (MatrixLog[a] - MatrixLog[b])]/Log[2]
```

It vanishes when the two states agree, identically in the Bloch vector rather than at sampled points.

```wl
Simplify[relativeEntropy[bloch[{rx, ry, rz}], bloch[{rx, ry, rz}]], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

Against the maximally mixed state it is the entropy bound in disguise, which is the one case where
nonnegativity is provable in closed form rather than sampled.

```wl
Simplify[relativeEntropy[bloch[{rx, ry, rz}], IdentityMatrix[2]/2] == 1 + Tr[bloch[{rx, ry, rz}] . MatrixLog[bloch[{rx, ry, rz}]]]/Log[2], rx^2 + ry^2 + rz^2 <= 1 && Element[{rx, ry, rz}, Reals]]
```

Written in the Bloch radius, no radius makes it negative, so `Reduce` returns `False` for the whole search.

```wl
Reduce[1 + ((1 + r)/2 Log[2, (1 + r)/2] + (1 - r)/2 Log[2, (1 - r)/2]) < 0 && 0 <= r < 1, r, Reals]
```

For a general pair the symbolic form is unwieldy, so the check is a sampled one: two hundred independent pairs
drawn uniformly from the Bloch ball, none of them negative.

```wl
Min @ Table[Re[relativeEntropy @@ (bloch /@ RandomPoint[Ball[{0, 0, 0}], 2])], 200]
```

Equality is sharp. Moving the second state by one part in a hundred already costs a positive number of bits.

```wl
Chop @ N @ {relativeEntropy[bloch[{1/5, 0, 1/2}], bloch[{1/5, 0, 1/2}]], relativeEntropy[bloch[{1/5, 0, 1/2}], bloch[{1/5, 0, 51/100}]]}
```

The asymmetry is not a small effect. Comparing the maximally mixed state with a nearly pure one gives
different answers in the two directions, and mistaking one for the other is a factor of nearly two here.

```wl
Chop @ N @ {relativeEntropy[bloch[{0, 0, 0}], bloch[{0, 0, 9/10}]], relativeEntropy[bloch[{0, 0, 9/10}], bloch[{0, 0, 0}]]}
```

Driving $\sigma$ to a pure state closes its support and the quantity diverges, which is the limit rather than
an error.

```wl
Limit[relativeEntropy[IdentityMatrix[2]/2, bloch[{0, 0, 1 - eps}]], eps -> 0, Direction -> "FromAbove"]
```

For commuting states it is the classical divergence, whose minimum over the whole square is zero, attained
only on the diagonal $p = q$.

```wl
NMinimize[{p Log[2, p/q] + (1 - p) Log[2, (1 - p)/(1 - q)], 0 < p < 1, 0 < q < 1}, {p, q}]
```

**QF** : the framework carries it as a named distance, returning a `Quantity` in bits like the entropies of
14.1. On a pair of well separated mixed states it agrees with the definition.

```wl
Chop @ N @ {QuantityMagnitude @ QuantumDistance[QuantumState[bloch[{1/5, 0, 1/2}]], QuantumState[bloch[{0, 1/4, -1/3}]], "RelativeEntropy"], relativeEntropy[bloch[{1/5, 0, 1/2}], bloch[{0, 1/4, -1/3}]]}
```

Two branches of that implementation are worth knowing before trusting it. When the second argument is pure the
function returns the entropy of the first argument instead of the divergence, so the answer that should be
$+\infty$ comes back finite and equal to $S(\rho)$.

```wl
{QuantumDistance[QuantumState[bloch[{1/5, 0, 1/2}]], QuantumState[{1, 0}], "RelativeEntropy"] === QuantumState[bloch[{1/5, 0, 1/2}]]["VonNeumannEntropy"], Limit[relativeEntropy[bloch[{1/5, 0, 1/2}], bloch[{0, 0, 1 - eps}]], eps -> 0, Direction -> "FromAbove"]}
```

And when the cross term $|\mathrm{Tr}\,\rho\log\sigma|$ exceeds $2\log d$, the implementation falls back to the
*reversed* divergence. The asymmetric pair above sits past that threshold, so what comes back is
$S(\sigma\,\|\,\rho)$, the smaller of the two numbers.

```wl
N @ {Abs[Tr[bloch[{0, 0, 0}] . MatrixLog[bloch[{0, 0, 9/10}]]]], 2 Log[2]}
```

```wl
Chop @ N @ {QuantityMagnitude @ QuantumDistance[QuantumState[bloch[{0, 0, 0}]], QuantumState[bloch[{0, 0, 9/10}]], "RelativeEntropy"], relativeEntropy[bloch[{0, 0, 9/10}], bloch[{0, 0, 0}]], relativeEntropy[bloch[{0, 0, 0}], bloch[{0, 0, 9/10}]]}
```

Since the quantity is unbounded, the paired similarity maps it through $2^{-S}$ rather than $1 - S$, so it
lands in $(0, 1]$ with $1$ at $\rho = \sigma$.

```wl
Chop @ N @ QuantumSimilarity[QuantumState[bloch[{1/5, 0, 1/2}]], QuantumState[bloch[{0, 1/4, -1/3}]], "RelativeEntropy"]
```

Nonnegativity is the whole content, and it is sharp: the quantity is zero only when the two states are equal
and positive otherwise, which is what makes it usable as evidence. It is not a metric, being asymmetric and
unbounded, so it is a divergence, and the direction has to be stated every time. Against $I/d$ it reproduces
$\log_2 d - S(\rho)$, on commuting states it reduces to the Kullback-Leibler divergence, and it goes to
infinity exactly when the second state assigns zero probability where the first does not. In QF the named
distance is faithful for well separated mixed pairs but substitutes $S(\rho)$ when $\sigma$ is pure and
silently reverses the arguments once $|\mathrm{Tr}\,\rho\log\sigma| > 2\log d$, so for a divergence near either
boundary the definition above is the reliable route.

### 14.3 [MSc] How do I compute joint, conditional, and mutual quantum entropies (and find a negative conditional entropy)?

For a bipartite state $\rho_{AB}$ with marginals $\rho_A = \mathrm{Tr}_B\,\rho_{AB}$ and
$\rho_B = \mathrm{Tr}_A\,\rho_{AB}$, three numbers come out of the single functional of 14.1,

$$S(A,B) = S(\rho_{AB}), \qquad S(A|B) = S(A,B) - S(B), \qquad I(A{:}B) = S(A) + S(B) - S(A,B),$$

the joint entropy being that functional applied to the whole state and the other two differences of it. The
whole question is which of the classical relations among them survive.

That the conditional entropy is *defined* as a difference is the crux. Classically
$H(A|B) = \sum_b p_b H(A|B{=}b)$ averages genuine entropies over the outcomes of $B$ and so cannot be negative.
Here there is no joint distribution over a pair of values to condition on, only two states and a subtraction,
and nothing makes the result positive.

Two classical facts survive. Subadditivity $S(A,B) \le S(A) + S(B)$ is exactly $I(A{:}B) \ge 0$, and it is
Klein's inequality of 14.2 in disguise: since $\log(\rho_A\otimes\rho_B) = \log\rho_A\otimes I + I\otimes\log\rho_B$,
the divergence of the state from its own product of marginals expands to

$$S(\rho_{AB} \,\|\, \rho_A \otimes \rho_B) = I(A{:}B),$$

so the mutual information is nonnegative for the same reason a relative entropy is. What fails is the lower
bound $H(A,B) \ge \max\{H(A), H(B)\}$, that a whole is at least as uncertain as either of its parts. Its
quantum replacement is the Araki-Lieb triangle inequality

$$|S(A) - S(B)| \le S(A,B) \le S(A) + S(B), \qquad\text{that is}\qquad S(A|B) \ge -S(A),$$

so the conditional entropy is bounded below, by $-S(A)$ rather than by zero.

A pure entangled state is where the classical bound fails hardest, since $S(A,B) = 0$ leaves the parts free to
carry entropy the whole does not. Local unitaries change no entropy and bring any bipartite pure state to its
Schmidt form (Part 11, 11.3), so for everything below

$$|\psi(t)\rangle = \cos t\,|00\rangle + \sin t\,|11\rangle, \qquad 0 < t < \frac{\pi}{2},$$

is the general entangled pure state of two qubits, with $t = \pi/4$ the maximally entangled point and the two
excluded endpoints the products.

**WL** : the von Neumann entropy read off the spectrum, so that the $0\log 0 = 0$ convention acts on a state
with rank-deficient support.

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

Each marginal is a contraction of the reshaped density matrix over the other factor's pair of indices.

```wl
reduceA[m_] := TensorContract[ArrayReshape[m, {2, 2, 2, 2}], {{2, 4}}]
reduceB[m_] := TensorContract[ArrayReshape[m, {2, 2, 2, 2}], {{1, 3}}]
```

The Schmidt family as a density matrix.

```wl
rhoAB = With[{psi = {Cos[t], 0, 0, Sin[t]}}, Outer[Times, psi, Conjugate[psi]]]
```

One functional, three arguments: the two parts and the whole.

```wl
{sA, sB, sAB} = Simplify[vonNeumann /@ {reduceA[rhoAB], reduceB[rhoAB], rhoAB}, 0 < t < Pi/2]
```

The joint, conditional and mutual entropies are assembled from those three by the definitions above.

```wl
Simplify[{sAB, sAB - sB, sA + sB - sAB}, 0 < t < Pi/2]
```

Measured in units of the marginal entropy, which is positive across the open interval, they collapse to
constants: the whole family has one and the same entropic structure, only rescaled.

```wl
Simplify[{sAB, sAB - sB, sA + sB - sAB}/sA, 0 < t < Pi/2]
```

The extreme values over the family are a search, not a quotation. The lowest conditional entropy and the
highest mutual information are reached at the same angle.

```wl
FullSimplify @ {Minimize[{sAB - sB, 0 <= t <= Pi/2}, t], Maximize[{sA + sB - sAB, 0 <= t <= Pi/2}, t]}
```

Erasing the one coherence between the two Schmidt terms leaves a state with the same populations, correlated
but only classically.

```wl
dephased = DiagonalMatrix[Diagonal[rhoAB]]
```

Dephasing in the Schmidt basis leaves both marginals untouched, so $S(A)$ and $S(B)$ are the same numbers as
before and whatever changes is the whole's doing.

```wl
Simplify[{reduceA[dephased], reduceB[dephased]} == {reduceA[rhoAB], reduceB[rhoAB]}, 0 < t < Pi/2]
```

The same three quantities in the same units, for that state.

```wl
Simplify[With[{d = vonNeumann[dephased], a = vonNeumann[reduceA[dephased]], b = vonNeumann[reduceB[dephased]]}, {d, d - b, a + b - d}]/sA, 0 < t < Pi/2]
```

The two surviving inequalities are claims about every state, and the pure family is too narrow to test them,
so draw mixed states at random from a Ginibre matrix normalized to unit trace.

```wl
randomState := #/Tr[#] & @ (# . ConjugateTranspose[#] &) @ RandomComplex[{-1 - I, 1 + I}, {4, 4}]
```

Subadditivity needs no sampling once the mutual information is recognized as the relative entropy against the
product of marginals, an identity a single generic state already exposes.

```wl
Chop @ Re @ With[{r = randomState}, Tr[r . (MatrixLog[r] - MatrixLog[KroneckerProduct[reduceA[r], reduceB[r]]])]/Log[2] - (vonNeumann[reduceA[r]] + vonNeumann[reduceB[r]] - vonNeumann[r])]
```

The Araki-Lieb bound is the one that limits how negative a conditional entropy can get, and its slack over two
hundred random mixed states stays on one side of zero.

```wl
Chop @ Re @ Min @ Table[With[{r = randomState}, vonNeumann[r] - Abs[vonNeumann[reduceA[r]] - vonNeumann[reduceB[r]]]], 200]
```

**QF** : the four amplitudes factor into two qudits on their own, so the family is a two-qubit object without
being told so.

```wl
qs = QuantumState[{Cos[t], 0, 0, Sin[t]}]
```

The marginals are partial traces, and they come back as states, not matrices.

```wl
{qa, qb} = QuantumPartialTrace[qs, #] & /@ {{2}, {1}}
```

The same three quantities, each entropy read off its own object with base $2$ selected by the second argument.

```wl
qtriple = Simplify[{qs["VonNeumannEntropy", 2], qs["VonNeumannEntropy", 2] - qb["VonNeumannEntropy", 2], qa["VonNeumannEntropy", 2] + qb["VonNeumannEntropy", 2] - qs["VonNeumannEntropy", 2]}, 0 < t < Pi/2]
```

```wl
qtriple === Simplify[{sAB, sAB - sB, sA + sB - sAB}, 0 < t < Pi/2]
```

Two of the three are named quantities. `"MutualInformationI"` assembles $I(A{:}B)$ from the same three
entropies, and `"EntanglementEntropy"` is the marginal entropy, which on a pure state is minus the conditional
entropy; both carry a unit, as the entropies of 14.1 do.

```wl
Simplify[QuantumEntanglementMonotone[qs, #] & /@ {"MutualInformationI", "EntanglementEntropy"}, 0 < t < Pi/2]
```

```wl
Simplify[QuantumEntanglementMonotone[qs, #] & /@ {"MutualInformationI", "EntanglementEntropy"}, 0 < t < Pi/2] === (Quantity[#, "Bits"] & /@ {sA + sB - sAB, sB - sAB})
```

One caveat before the named quantity is used as a test. Mutual information counts total correlation, classical
included, so it is positive on a separable state; take the maximally entangled state with its coherence erased.

```wl
qdephased = QuantumState[DiagonalMatrix[Diagonal[QuantumState["PhiPlus"]["DensityMatrix"]]], QuantumBasis[2, 2]]
```

Its mutual information is a full bit, while the default entanglement test, which uses the realignment
criterion, reports the state separable.

```wl
{QuantumEntanglementMonotone[qdephased, "MutualInformationI"], QuantumEntangledQ[qdephased]}
```

So the monotone measures correlation, not entanglement. `QuantumEntangledQ` will not convert it into an
entanglement test either: asked for the mutual-information method, which it takes as its *third* argument after
the bipartition, it returns `Indeterminate` for a maximally entangled state as readily as for this separable
one, so the answer carries no information about the state.

```wl
QuantumEntangledQ[#, Automatic, "MutualInformationI"] & /@ {qdephased, QuantumState["PhiPlus"]}
```

Three numbers from one functional, and only the joint entropy is an entropy at all; the other two are
differences, which is exactly why one of them can be negative. On the pure family $S(A,B) = 0$, so
$S(A|B) = -S(A)$ and $I(A{:}B) = 2S(A)$: the conditional entropy is negative wherever the state is entangled,
saturating the Araki-Lieb bound at every angle, and it reaches $-1$ bit at $t = \pi/4$, where the mutual
information reaches $2$ bits, twice the most a pair of classical bits can share. Dephasing that state in its
Schmidt basis leaves both marginals exactly where they were and moves the pair to $S(A|B) = 0$ and
$I(A{:}B) = S(A)$, so half the correlation lived in the single coherence, and the negative part of the
conditional entropy is precisely the part with no classical counterpart. That part is not a bookkeeping
artifact: $-S(A|B)$ is the coherent information, and in the state-merging protocol $S(A|B)$ is the quantum
communication cost of handing $A$ to whoever holds $B$, so a negative value says the transfer costs nothing and
leaves $|S(A|B)|$ ebits over.

### 14.4 [MSc] How do I compute the coherent information $I_c(\rho,\mathcal{N})=S(\mathcal{N}(\rho))-S((\mathcal{N}\otimes\mathrm{id})(|\psi_\rho\rangle\langle\psi_\rho|))$ of a state through a channel?

The entropy of the output alone cannot say how much *quantum* information survived a channel. A channel that
discards its input and emits a fixed pure state has a zero-entropy output and has destroyed everything, so the
quantity has to be told what went in. Purify $\rho$ into $|\psi_\rho\rangle_{RA}$, so that an untouched
reference $R$ holds a perfect record of $A$, send $A$ alone through $\mathcal{N}$, and weigh what the output
keeps of that record against what the pair keeps:

$$I_c(\rho, \mathcal{N}) = S(B) - S(RB), \qquad \rho_{RB} = (\mathrm{id}_R \otimes \mathcal{N})(|\psi_\rho\rangle\langle\psi_\rho|).$$

Read against 14.3 this is one of the same three quantities with its sign reversed, $I_c = -S(R|B)$, so positive
coherent information is exactly the regime where the conditional entropy of the reference given the output is
negative, and the state-merging reading carries over: the channel leaves ebits rather than costing them.

Nothing depends on which purification is used. Any two purifications of $\rho$ differ by an isometry acting on
$R$ alone, which commutes with $\mathcal{N}\otimes\mathrm{id}$ and changes no entropy, so $I_c$ is a function of
$\rho$ and $\mathcal{N}$ only.

The second form comes from Stinespring. Writing $\mathcal{N}(\rho) = \mathrm{Tr}_E[V\rho V^\dagger]$ for an
isometry $V$ into system and environment, the triple $RBE$ is pure whenever $RA$ was, so $S(RB) = S(E)$ and

$$I_c = S(B) - S(E),$$

the output's entropy minus the environment's. The channel is good exactly when it keeps more than it leaks, and
$\mathcal{N}^c$, the map carrying $\rho$ to that environment state, is the complementary channel.

For a concrete channel take amplitude damping, an excited state decaying to the ground state with probability
$\gamma$, with Kraus operators

$$K_0 = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-\gamma}\end{pmatrix}, \qquad K_1 = \begin{pmatrix} 0 & \sqrt{\gamma} \\ 0 & 0\end{pmatrix},$$

on the family $\rho_p = (1-p)\,|0\rangle\langle 0| + p\,|1\rangle\langle 1|$. That family is the one that decides
this channel: damping commutes with rotations about $z$, so a coherence between the two levels contributes
nothing the diagonal does not already carry, and the last cell below prices what an off-axis state costs.

**WL** : the entropy functional, read off the spectrum so the $0\log 0 = 0$ convention acts.

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

Every state below turns out diagonal, so the binary entropy is the closed form to compare against.

```wl
h2[x_] := -x Log[2, x] - (1 - x) Log[2, 1 - x]
```

The channel as its Kraus operators.

```wl
kraus = {{{1, 0}, {0, Sqrt[1 - g]}}, {{0, Sqrt[g]}, {0, 0}}}
```

A Kraus set is a channel only if it preserves the trace, which is a condition on the operators alone.

```wl
Simplify[Total[ConjugateTranspose[#] . # & /@ kraus] == IdentityMatrix[2], 0 < g < 1]
```

The purification of $\rho_p$ in Schmidt form, the reference carrying the second factor.

```wl
psi = {Sqrt[1 - p], 0, 0, Sqrt[p]}
```

It is a purification because tracing the reference away returns the state it was built from.

```wl
Simplify[TensorContract[ArrayReshape[Outer[Times, psi, Conjugate[psi]], {2, 2, 2, 2}], {{2, 4}}] == DiagonalMatrix[{1 - p, p}], 0 < p < 1]
```

The channel acts on the first factor only, the reference passing through untouched.

```wl
rhoRB = Total[With[{k = KroneckerProduct[#, IdentityMatrix[2]]}, k . Outer[Times, psi, Conjugate[psi]] . ConjugateTranspose[k]] & /@ kraus]
```

The output on its own is that state with the reference traced away.

```wl
rhoB = TensorContract[ArrayReshape[rhoRB, {2, 2, 2, 2}], {{2, 4}}]
```

The two entropies the definition asks for.

```wl
{entB, entRB} = FullSimplify[vonNeumann /@ {rhoB, rhoRB}, 0 < p < 1 && 0 < g < 1]
```

Both are binary entropies of a damped population, one of what survives and one of what leaks.

```wl
FullSimplify[{entB, entRB} == {h2[(1 - g) p], h2[g p]}, 0 < p < 1 && 0 < g < 1]
```

The complementary channel sends the input to the environment, its matrix elements the overlaps of the Kraus
operators in that state.

```wl
rhoE = Outer[Tr[#1 . DiagonalMatrix[{1 - p, p}] . ConjugateTranspose[#2]] &, kraus, kraus, 1]
```

The environment's entropy is the pair's, which is the Stinespring identity and the second form of $I_c$.

```wl
FullSimplify[vonNeumann[rhoE] == entRB, 0 < p < 1 && 0 < g < 1]
```

Which of the two binary entropies wins is settled by monotonicity, since the binary entropy rises across the
lower half of its range.

```wl
Simplify[D[h2[x], x] > 0, 0 < x < 1/2]
```

The best any input can do, on either side of $\gamma = 1/2$ and at it.

```wl
FullSimplify @ Table[Maximize[{h2[(1 - g) x] - h2[g x], 0 <= x <= 1}, x], {g, {1/4, 1/2, 3/4}}]
```

**QF** : the channel is named, and an order fixes which qudit it acts on.

```wl
channel = QuantumChannel["AmplitudeDamping"[g], {1}]
```

A state purifies itself, with the reference supplied as a second qudit.

```wl
qpsi = QuantumState[DiagonalMatrix[{1 - p, p}]]["Purify"]
```

What it built is the Schmidt form the WL side wrote out by hand.

```wl
FullSimplify[Normal @ qpsi["StateVector"], 0 < p < 1] === psi
```

Applying the channel to the first qudit leaves the reference alone, and the result is a state, not a matrix.

```wl
qout = channel[qpsi]
```

The same two entropies, each read off its own object.

```wl
qpair = FullSimplify[{QuantumPartialTrace[qout, {2}]["VonNeumannEntropy", 2], qout["VonNeumannEntropy", 2]}, 0 < p < 1 && 0 < g < 1]
```

```wl
qpair === {entB, entRB}
```

A channel carries its own dilation: `channel["Operator"]` keeps the Kraus index as an extra qudit instead of
tracing it, so applying it gives the pure triple $RBE$ whose existence the Stinespring form assumes.

```wl
qtri = channel["Operator"][qpsi]
```

Its three marginals are the environment, the output and the reference, in that order.

```wl
FullSimplify[Normal @ QuantumPartialTrace[qtri, #]["DensityMatrix"] & /@ {{2, 3}, {1, 3}, {1, 2}}, 0 < p < 1 && 0 < g < 1]
```

```wl
FullSimplify[QuantumPartialTrace[qtri, {2, 3}]["VonNeumannEntropy", 2] == qpair[[2]], 0 < p < 1 && 0 < g < 1]
```

Rotating the reference is the freedom in the choice of purification, and it moves neither entropy.

```wl
With[{rotated = channel[QuantumOperator["H", {2}][qpsi]]}, FullSimplify[{QuantumPartialTrace[rotated, {2}]["VonNeumannEntropy", 2], rotated["VonNeumannEntropy", 2]} == qpair, 0 < p < 1 && 0 < g < 1]]
```

Nothing in the recipe used the input's diagonality, so an off-axis mixed state goes through unchanged, and the
comparison is against the diagonal state carrying the same population.

```wl
With[{r = {1/3, 1/4, 1/5}}, Chop @ N @ {With[{out = QuantumChannel["AmplitudeDamping"[2/5], {1}][QuantumState[(IdentityMatrix[2] + r . PauliMatrix[{1, 2, 3}])/2]["Purify"]]}, QuantumPartialTrace[out, {2}]["VonNeumannEntropy", 2] - out["VonNeumannEntropy", 2]], h2[(1 - 2/5) (1 - r[[3]])/2] - h2[2/5 (1 - r[[3]])/2]}]
```

The coherent information is a difference of two entropies, and which two is the whole content. Against the
output it is the pair $S(B) - S(RB)$, and against the environment the same number is $S(B) - S(E)$, so a channel
transmits quantum information exactly when its output holds more of the reference's record than its environment
does. On the damping family both terms collapse to binary entropies of a damped population, $I_c = h_2((1-\gamma)p) - h_2(\gamma p)$,
and the competition is between the fraction $(1-\gamma)p$ that survives and the fraction $\gamma p$ that leaks.
At $\gamma = 1/2$ those are equal, the channel is its own complement, and $I_c$ vanishes for every input, not
merely at its optimum. Past that point the environment gets the larger share, the best any input can manage is
zero, and the channel is antidegradable: the environment could reconstruct the output, so no quantum information
crosses at all and the quantum capacity is zero. Below it the maximum is strictly positive, $\log_2(4/3)$ bits
at $\gamma = 1/4$, reached at $p = 4/9$ rather than at the maximally mixed input, which is why the capacity
question is an optimization and not a substitution. Coherence in the input is not rewarded here: the off-axis
state scores below the diagonal state with the same population, since damping is covariant under phase rotations
and the coherence it carries survives only to be traced away.

### 14.5 [MSc] How do I verify the data-processing inequality under a channel?

The relative entropy of 14.2 measures how well one state can be told from another, so the one thing a noisy
process must never do is raise it. For every channel $\mathcal{N}$ and every pair of states,

$$S\big(\mathcal{N}(\rho)\,\big\|\,\mathcal{N}(\sigma)\big) \le S(\rho\,\|\,\sigma).$$

Processing cannot manufacture distinguishability. Anything that could be learned about which of $\rho$ or
$\sigma$ was prepared, by measuring after the channel, could already have been learned by measuring before it,
since the channel itself is one of the things a measurement could have done first. Equality is the exceptional
case: it holds exactly when the channel is reversible on that particular pair, in the sense that a recovery map
returns both states unchanged.

This is the master inequality of the Part. The Holevo bound of 14.6 and the entropy bounds already used are all
this statement applied to a chosen pair, and the nearest one is the mutual information of 14.3. Reading it as a
relative entropy,

$$I(A{:}B) = S\big(\rho_{AB}\,\big\|\,\rho_A \otimes \rho_B\big),$$

a channel acting on $B$ alone carries $\rho_{AB}$ to the processed state and carries $\rho_A \otimes \rho_B$ to
the product of the *new* marginals, so both arguments move together and the inequality applies verbatim. No
separate proof is needed, which is worth more than a second verification: the corollary is the theorem with
different letters.

For a channel take depolarizing noise, which mixes a state with the maximally mixed one and so shrinks its Bloch
vector, $\vec r \mapsto (1-p)\vec r$, with $p = 0$ the identity and $p = 1$ total erasure. The two states are
taken non-commuting, so that what is being contracted is genuinely quantum distinguishability rather than the
classical divergence of two diagonal distributions.

**WL** : the state builder and the divergence, read straight off the definition.

```wl
bloch[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2
```

```wl
relativeEntropy[a_, b_] := Tr[a . (MatrixLog[a] - MatrixLog[b])]/Log[2]
```

The channel as its four Kraus operators, the identity kept with weight $1 - 3p/4$ and each Pauli error with
weight $p/4$.

```wl
krausOps = {Sqrt[1 - 3 p/4] IdentityMatrix[2], Sqrt[p/4] PauliMatrix[1], Sqrt[p/4] PauliMatrix[2], Sqrt[p/4] PauliMatrix[3]}
```

```wl
depolarize[m_] := Total[# . m . ConjugateTranspose[#] & /@ krausOps]
```

What that does to a state is visible only in the Bloch picture, where it is a single contraction toward the
centre.

```wl
FullSimplify[depolarize[bloch[{rx, ry, rz}]] == bloch[(1 - p) {rx, ry, rz}], 0 < p < 1 && Element[{rx, ry, rz}, Reals]]
```

A pair whose Bloch vectors point in different directions, so the two states do not share an eigenbasis.

```wl
{rho, sigma} = bloch /@ {{0, 0, 3/5}, 4/5 {Sin[Pi/3], 0, Cos[Pi/3]}}
```

```wl
Simplify[rho . sigma != sigma . rho]
```

The divergence before the channel is a number, and after it a function of the noise.

```wl
before = FullSimplify @ relativeEntropy[rho, sigma]
```

```wl
after = FullSimplify[relativeEntropy[depolarize[rho], depolarize[sigma]], 0 < p < 1]
```

The inequality is then a question with no free parameters left, and the empty violating set is the answer.

```wl
Reduce[after > before && 0 < p < 1, p, Reals]
```

Both ends of the noise range are forced: at $p = 0$ nothing has happened, and at $p = 1$ both states have become
the same maximally mixed state and nothing can tell them apart.

```wl
{Simplify[after /. p -> 0] == before, Limit[after, p -> 1, Direction -> "FromBelow"]}
```

One pair proves nothing about all pairs. Relative entropy is invariant under a common unitary and depolarizing
noise is unitarily covariant, so a single rotation puts the first Bloch vector on $z$ and the second in the
$xz$ plane, leaving two lengths, one angle and the noise to search over.

```wl
slack[al_?NumericQ, be_?NumericQ, th_?NumericQ, q_?NumericQ] := Re[relativeEntropy[bloch[{0, 0, al}], bloch[be {Sin[th], 0, Cos[th]}]] - relativeEntropy[bloch[{0, 0, (1 - q) al}], bloch[(1 - q) be {Sin[th], 0, Cos[th]}]]]
```

```wl
NMinimize[{slack[al, be, th, q], 0 < al < 1 && 0 < be < 1 && 0 <= th <= Pi && 0 < q < 1}, {al, be, th, q}]
```

For the corollary the entropy is read off the spectrum, so that a pure state costs nothing.

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

```wl
reduceA[m_] := TensorContract[ArrayReshape[m, {2, 2, 2, 2}], {{2, 4}}]
reduceB[m_] := TensorContract[ArrayReshape[m, {2, 2, 2, 2}], {{1, 3}}]
```

```wl
mutual[m_] := vonNeumann[reduceA[m]] + vonNeumann[reduceB[m]] - vonNeumann[m]
```

An entangled pure pair, with the noise acting on the second qubit only.

```wl
rhoAB = With[{psi = {Cos[Pi/5], 0, 0, Sin[Pi/5]}}, Outer[Times, psi, Conjugate[psi]]]
```

```wl
localB[m_] := Total[With[{k = KroneckerProduct[IdentityMatrix[2], #]}, k . m . ConjugateTranspose[k]] & /@ krausOps]
```

The correlation the two qubits share, as the noise on one of them is turned up. Substituting into the state
before applying the functional is what lets the $0\log 0 = 0$ convention act at the two pure endpoints.

```wl
N @ Table[mutual[localB[rhoAB] /. p -> pv], {pv, {0, 1/4, 1/2, 3/4, 1}}]
```

That decay needs no separate argument once the mutual information is recognized as a divergence against the
product of the marginals.

```wl
With[{m = localB[rhoAB] /. p -> 2/5}, Chop @ N @ {relativeEntropy[m, KroneckerProduct[reduceA[m], reduceB[m]]], mutual[m]}]
```

And once the channel is seen to carry that product to the product of the processed marginals, which is what
makes the pair move together and the corollary immediate.

```wl
FullSimplify[localB[KroneckerProduct[reduceA[rhoAB], reduceB[rhoAB]]] == KroneckerProduct[reduceA[localB[rhoAB]], reduceB[localB[rhoAB]]], 0 < p < 1]
```

**QF** : the state is named by its Bloch vector, the channel is named, and the divergence is the named distance
of 14.2, so none of the three helpers above is needed and the noise stays symbolic throughout.

```wl
{qrho, qsigma} = QuantumState["BlochVector"[#]] & /@ {{0, 0, 3/5}, 4/5 {Sin[Pi/3], 0, Cos[Pi/3]}}
channel = QuantumChannel["Depolarizing"[p]]
```

The contraction toward the centre is returned rather than asserted, since the Bloch vector is a property of the
processed state.

```wl
FullSimplify[channel[QuantumState["BlochVector"[{rx, ry, rz}]]]["BlochVector"], 0 < p < 1 && Element[{rx, ry, rz}, Reals]]
```

The named distance does carry a symbolic parameter, but only if the domain is supplied: its raw output is a
2322-leaf expression of nested conjugated square roots, and the assumption $0<p<1$ is what collapses it to
sixty-odd leaves. Without it the result is unreadable and nothing downstream will decide.

```wl
rel[a_, b_] := FullSimplify[QuantityMagnitude @ QuantumDistance[a, b, "RelativeEntropy"], 0 < p < 1]
```

```wl
{qbefore, qafter} = {rel[qrho, qsigma], rel[channel[qrho], channel[qsigma]]}
```

The same two quantities the definition gave, which is the check that the named distance is the divergence and
not a variant of it.

```wl
FullSimplify[{qbefore - before, qafter - after}, 0 < p < 1]
```

So the inequality is again a question with no free parameters, and the empty violating set is again the answer,
now reached without a Kraus operator or a matrix logarithm written out.

```wl
Reduce[qafter > qbefore && 0 < p < 1, p, Reals]
```

The closed form also settles more than the inequality: both endpoints, and that the fall is strict everywhere
between them rather than merely non-increasing.

```wl
{Simplify[qafter /. p -> 0] == qbefore, Limit[qafter, p -> 1, Direction -> "FromBelow"], Reduce[D[qafter, p] > 0 && 0 < p < 1, p, Reals]}
```

One pair and one channel family still prove nothing about all channels. Here the named library is the argument:
the inequality is not about depolarizing noise, and every channel in it can be asked the same question in the
same line, where each would otherwise need its Kraus set written out.

```wl
Table[With[{c = QuantumChannel[ch]}, {Head[ch], rel[c[qrho], c[qsigma]] <= rel[qrho, qsigma]}],
  {ch, {"BitFlip"[1/3], "PhaseFlip"[1/3], "BitPhaseFlip"[1/3], "Depolarizing"[2/5],
        "AmplitudeDamping"[1/3], "PhaseDamping"[1/3], "ResetError"[1/5, 1/5]}}]
```

The corollary is a named monotone, with the noise on the second qubit turned up through the same range.

```wl
qpsi = QuantumState[{Cos[Pi/5], 0, 0, Sin[Pi/5]}]
```

```wl
N @ Table[QuantityMagnitude @ QuantumEntanglementMonotone[
   QuantumChannel["Depolarizing"[pv], {2}][qpsi], "MutualInformationI"], {pv, {0, 1/4, 1/2, 3/4, 1}}]
```

One inequality, and everything else in the Part is a corollary of it with the arguments chosen well. Its content
is that a channel is a blur: the divergence between two states can only fall, and it falls strictly unless the
channel happens to be reversible on that pair. On the depolarizing family the fall is total, from the initial
divergence down to zero as $p \to 1$, because at full strength every state becomes the same maximally mixed
state and no measurement whatsoever can distinguish inputs that have become identical. The mutual information
inherits this without a new argument, since it is the divergence of a joint state from the product of its
marginals and a local channel moves both of those together, so correlation with an untouched reference decays
under noise for exactly the reason distinguishability does. The equality case is the useful edge: a channel that
loses no relative entropy on a pair has destroyed nothing recoverable about that pair, which is the entropic
statement of reversibility, and it is why the coherent information of 14.4 turns negative precisely when a
channel has leaked to its environment what it can no longer be asked to give back.

### 14.6 [MSc] How do I compute the Holevo bound $\chi=S\!\big(\sum_i p_i\rho_i\big)-\sum_i p_i\,S(\rho_i)$ on the accessible information of an ensemble?

An ensemble $\{p_i, \rho_i\}$ is a way of carrying a classical choice on a quantum system: the sender draws an
index $i$ with probability $p_i$ and hands over the state $\rho_i$, and the receiver, holding one copy and free
to choose any measurement whatsoever, tries to work out which $i$ was drawn. Write $X$ for the index the sender
drew and $Y$ for the outcome the measurement returns; then the mutual information $I(X{:}Y)$ of 14.3 is how much
of the choice survived the trip. The Holevo quantity $\chi$ bounds that above for every measurement at once,
which is what makes it worth computing: it is a property of the encoding alone, fixed before any receiver is
specified.

Its form is the failure of the entropy to be linear. The receiver who is told nothing holds the average state
$\bar\rho = \sum_i p_i \rho_i$ and faces $S(\bar\rho)$ worth of uncertainty; the receiver who is told $i$ still
faces $S(\rho_i)$, on average $\sum_i p_i S(\rho_i)$. The difference is what being told the index was worth,
and since $S$ is concave the difference cannot be negative. So $\chi$ is a concavity gap, and it collapses to
zero exactly when learning $i$ tells the receiver nothing about the state, which happens only when every
$\rho_i$ carrying weight is the same state.

Two rewritings make the bound a corollary of what the Part already has rather than a new theorem. Averaging the
divergence of each member against the ensemble average,

$$\chi = \sum_i p_i \, S(\rho_i \,\|\, \bar\rho),$$

turns it into a relative entropy, so Klein's inequality of 14.2 delivers $\chi \ge 0$ and the data-processing
inequality of 14.5 delivers its decrease under any channel applied to the states. Alternatively, recording the
index itself in an orthonormal register gives the classical-quantum state

$$\rho_{XB} = \sum_i p_i \, |i\rangle\langle i| \otimes \rho_i,$$

whose mutual information $I(X{:}B)$ is $\chi$ on the nose, so the quantity is the mutual information of 14.3
read on a state built to carry a classical variable in one half. The upper bound $\chi \le S(\bar\rho) \le
\log_2 d$ follows from the first form, and it is the sharp statement of how much a channel of a given dimension
can carry: a qubit is capped at one bit regardless of how many states the ensemble crowds into it.

**WL** : the generic mixed qubit, the entropy read off the spectrum, and the definition assembled from the two.

```wl
bloch[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2
```

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

```wl
holevo[ps_, rhos_] := vonNeumann[ps . rhos] - ps . (vonNeumann /@ rhos)
```

Two pure states whose Bloch vectors are separated by $\theta$, so the overlap is $\cos(\theta/2)$ and the
Hilbert-space angle is half the Bloch angle. Each member is pure, so the subtracted average costs nothing and
what is left is the entropy of the mixture alone.

```wl
chiPure = FullSimplify[holevo[{1/2, 1/2}, bloch /@ {{0, 0, 1}, {Sin[th], 0, Cos[th]}}], 0 < th < Pi]
```

The two ends of the range are the orthogonal pair and the coincident pair.

```wl
{FullSimplify[chiPure /. th -> Pi], Limit[chiPure, th -> 0, Direction -> "FromAbove"]}
```

Mixing the members costs the subtracted term, and the bound drops strictly below the entropy of the average
that a pure ensemble would have reached.

```wl
ens = bloch /@ {{0, 0, 1/2}, {1/2 Sin[th], 0, 1/2 Cos[th]}}
```

```wl
FullSimplify[{vonNeumann[Total[ens]/2], holevo[{1/2, 1/2}, ens]} /. th -> Pi/2]
```

The average-divergence form is an identity, not an approximation, and it holds in the angle rather than at
sampled points.

```wl
relativeEntropy[a_, b_] := Tr[a . (MatrixLog[a] - MatrixLog[b])]/Log[2]
```

```wl
FullSimplify[holevo[{1/2, 1/2}, ens] - Total[relativeEntropy[#, Total[ens]/2] & /@ ens]/2, 0 < th < Pi]
```

Equality in the concavity is the degenerate ensemble, identically in the Bloch radius rather than at a chosen
state.

```wl
FullSimplify[holevo[{1/2, 1/2}, {bloch[{0, 0, r}], bloch[{0, 0, r}]}], 0 < r < 1]
```

The same quantity for three pure states at equal angles around a great circle, maximized over the priors, so
the search reports both the best the ensemble can do and the weights that reach it.

```wl
NMaximize[{holevo[{p1, p2, 1 - p1 - p2}, bloch /@ {{0, 0, 1}, {Sin[2 Pi/3], 0, Cos[2 Pi/3]}, {-Sin[2 Pi/3], 0, Cos[2 Pi/3]}}], 0 < p1 < 1 && 0 < p2 < 1 && p1 + p2 < 1}, {p1, p2}]
```

**QF** : the entropy is a property of the state and the weighted sum of the members is the average state, so the
definition assembles from the objects with no matrix anywhere in it.

```wl
qens = QuantumState["BlochVector"[#]] & /@ {{0, 0, 1/2}, {1/2 Sin[th], 0, 1/2 Cos[th]}}
```

```wl
qholevo[ps_, qs_] := Total[ps qs]["VonNeumannEntropy"] - ps . (#["VonNeumannEntropy"] & /@ qs)
```

```wl
qchi = FullSimplify[QuantityMagnitude @ qholevo[{1/2, 1/2}, qens], 0 < th < Pi]
```

The register form is the one QF reaches without being given the definition at all. The classical-quantum state
is the weighted sum of each register basis state tensored against the state it labels, and a weighted sum of states
whose type is `"Matrix"` is a mixture rather than a superposition, which is what the encoding calls for.

```wl
qcq = Total[MapThread[#1 QuantumTensorProduct[#2, #3] &, {{1/2, 1/2}, QuantumState /@ {"0", "1"}, qens}]]
```

Asking that state for the mutual information of 14.3.

```wl
qmi = FullSimplify[QuantityMagnitude @ QuantumEntanglementMonotone[qcq, "MutualInformationI"], 0 < th < Pi]
```

Both readings against the definition assembled by hand, as identities in the angle and not numerical
coincidences.

```wl
FullSimplify[{qchi, qmi} - holevo[{1/2, 1/2}, ens], 0 < th < Pi]
```

The quantity is a concavity gap, so it is nonnegative for the same reason the entropy is concave, and it is zero
only for an ensemble that encodes nothing because all of its members coincide. Read as an average divergence
against the ensemble mean it inherits everything proved in 14.2 and 14.5 without a new argument, which is why it
can only fall when the states are pushed through a channel; read as the mutual information of a register with
the system it is the quantity of 14.3 evaluated on a state built for the purpose, and QF returns it from that
reading alone. The ceiling $\log_2 d$ is the part with operational teeth: no choice of ensemble and no cleverness
in the priors extracts more than one bit per qubit, and non-orthogonal members do strictly worse, which is the
entropic reason a quantum channel of a given dimension has a finite classical capacity and the starting point
for 14.7.

### 14.7 [MSc] How do I assemble a channel-capacity estimate for a simple channel, the classical (Holevo) capacity or the quantum (coherent-information) capacity, from the entropic quantities above?

Neither capacity is a new quantity. A channel is handed a state and returns one, so every entropic quantity in
the Part can be evaluated on its output, and a capacity is what remains after the best input has been chosen.
Sending each member of an ensemble through $\mathcal{N}$ and taking the Holevo quantity of 14.6 on what comes
out, then maximizing over the encoding, gives the classical quantity

$$\chi(\mathcal{N}) = \max_{\{p_i, \rho_i\}} \left[ S\!\left(\mathcal{N}\!\left(\sum_i p_i \rho_i\right)\right) - \sum_i p_i \, S(\mathcal{N}(\rho_i)) \right],$$

and maximizing the coherent information of 14.4 over inputs gives the quantum one,

$$Q^{(1)}(\mathcal{N}) = \max_{\rho} I_c(\rho, \mathcal{N}).$$

The word *estimate* in the question is doing real work, and it is worth being exact about why. Both expressions
are single-letter: they ask what one use of the channel achieves. The capacities are the limits
$C = \lim_n \chi(\mathcal{N}^{\otimes n})/n$ and $Q = \lim_n Q^{(1)}(\mathcal{N}^{\otimes n})/n$, and neither
quantity is additive in general, so the single-letter value is a lower bound rather than the answer. Two classes
of channel escape the regularization: $\chi$ is additive for entanglement-breaking channels, and
$Q^{(1)} = Q$ for degradable ones, where the environment can be reconstructed from the output. The dephasing
channel below is degradable, so its coherent-information estimate is its quantum capacity exactly, which is what
makes it worth using as the worked case.

Take dephasing, which leaves the computational basis alone and randomizes the relative phase with probability
$p$, so $\mathcal{N}(\rho) = (1-p)\rho + p\,Z\rho Z$. Its classical side needs no search at all. The ceiling of
14.6 is $\chi \le \log_2 d$, one bit for a qubit, so an ensemble that reaches one bit is optimal by that bound
alone and no optimization is required to certify it. The quantum side is a genuine maximization, and the two
answers are different functions of $p$, which is the point of asking both.

**WL** : the entropy off the spectrum, and the binary entropy the diagonal states will reduce to.

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

```wl
h2[x_] := -x Log[2, x] - (1 - x) Log[2, 1 - x]
```

The channel as its Kraus pair, which is a channel only if the pair preserves the trace.

```wl
kraus = {Sqrt[1 - p] IdentityMatrix[2], Sqrt[p] PauliMatrix[3]}
```

```wl
Simplify[Total[ConjugateTranspose[#] . # & /@ kraus] == IdentityMatrix[2], 0 < p < 1]
```

```wl
chan[m_] := Total[# . m . ConjugateTranspose[#] & /@ kraus]
```

```wl
bloch[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2
```

```wl
holevo[ps_, rhos_] := vonNeumann[ps . rhos] - ps . (vonNeumann /@ rhos)
```

A code is a choice of ensemble, so the Holevo quantity of the output is a function of the basis it is written in.

```wl
wlChi[ps_, vs_] := holevo[ps, chan /@ (bloch /@ vs)]
```

For a code written in the basis the channel does not touch.

```wl
FullSimplify[wlChi[{1/2, 1/2}, {{0, 0, 1}, {0, 0, -1}}], 0 < p < 1]
```

And for one written in the basis the channel randomizes, which is what makes the choice of encoding a real one
rather than a formality.

```wl
FullSimplify[wlChi[{1/2, 1/2}, {{1, 0, 0}, {-1, 0, 0}}] == 1 - h2[p], 0 < p < 1]
```

For the quantum side the input is purified and only the first factor is sent, exactly as in 14.4, with the
population left free so the maximization has something to range over.

```wl
psi = {Sqrt[(1 + r)/2], 0, 0, Sqrt[(1 - r)/2]}
```

```wl
rhoRB = Total[With[{k = KroneckerProduct[#, IdentityMatrix[2]]}, k . Outer[Times, psi, Conjugate[psi]] . ConjugateTranspose[k]] & /@ kraus]
```

```wl
rhoB = TensorContract[ArrayReshape[rhoRB, {2, 2, 2, 2}], {{2, 4}}]
```

```wl
ic = FullSimplify[vonNeumann[rhoB] - vonNeumann[rhoRB], 0 < p < 1 && -1 < r < 1]
```

The best input at three noise strengths, the optimizer reporting where it sits as well as what it is worth.

```wl
FullSimplify @ Table[Maximize[{ic /. p -> pv, -1 <= r <= 1}, r], {pv, {1/8, 1/4, 1/2}}]
```

The three numbers the two estimates are made of, collected as functions of the noise rather than at sampled
points, so the QF side has something to be checked against.

```wl
wlCaps = FullSimplify[{wlChi[{1/2, 1/2}, {{0, 0, 1}, {0, 0, -1}}], wlChi[{1/2, 1/2}, {{1, 0, 0}, {-1, 0, 0}}], Limit[ic, r -> 0]}, 0 < p < 1]
```

**QF** : the channel is named, and the classical estimate is the register state of 14.6 with the channel acting
on the half that travels.

```wl
qens[vs_] := QuantumState["BlochVector"[#]] & /@ vs
```

The register carries one orthogonal state per ensemble member, so its dimension is set by the code rather than
fixed at a qubit.

```wl
reg[n_] := QuantumState["BasisState"[{#}], n] & /@ Range[0, n - 1]
```

Applying a state to a state is the tensor product with the argument placed first, so `Construct` pairs each
member with its label and leaves the register as the first qudit, and contracting against the priors is the
mixture.

```wl
qcq[ps_, qs_] := ps . MapThread[Construct, {qs, reg[Length[ps]]}]
```

```wl
qchi[ps_, qs_] := QuantityMagnitude @ QuantumEntanglementMonotone[QuantumChannel["PhaseFlip"[p], {2}][qcq[ps, qs]], "MutualInformationI"]
```

For the quantum estimate the optimal input is the one the maximization above lands on, and a state purifies
itself.

```wl
qout = QuantumChannel["PhaseFlip"[p], {1}][QuantumState[IdentityMatrix[2]/2]["Purify"]]
```

The two entropies the coherent information is built from, each read off its own object.

```wl
FullSimplify[{QuantumPartialTrace[qout, {2}]["VonNeumannEntropy", 2], qout["VonNeumannEntropy", 2]}, 0 < p < 1]
```

The same three numbers, every one of them read off an object rather than assembled from a matrix.

```wl
qfCaps = FullSimplify[{qchi[{1/2, 1/2}, qens[{{0, 0, 1}, {0, 0, -1}}]], qchi[{1/2, 1/2}, qens[{{1, 0, 0}, {-1, 0, 0}}]], QuantumPartialTrace[qout, {2}]["VonNeumannEntropy", 2] - qout["VonNeumannEntropy", 2]}, 0 < p < 1]
```

The two routes agree on both capacities as identities in the noise, which the printed closed forms do not make
obvious on their own since `FullSimplify` leaves them in different shapes.

```wl
FullSimplify[qfCaps - wlCaps, 0 < p < 1]
```

A capacity is an optimization of a quantity the Part already had, and the two optimizations answer different
questions about the same channel. Dephasing leaves the computational basis untouched, so a code written there
loses nothing and reaches the dimensional ceiling of 14.6 for every noise strength, while a code written in the
conjugate basis pays $h_2(p)$ for the phase it cannot keep: the classical estimate is a property of the best
encoding, not of the channel acting on any particular state. The quantum estimate is the one that degrades,
because the reference is entangled with what is sent and the phase the channel randomizes is exactly what that
entanglement is made of, so what survives is one bit less the entropy of the noise. At half dephasing the two
part company completely, a full classical bit crossing a channel through which no quantum information passes at
all, which is why the two capacities are separate resources and not two scales for one. Since dephasing is
degradable the quantum number is the capacity and not merely a bound on it, and where degradability fails, as for
the damping channel of 14.4 past $\gamma = 1/2$, the estimate has to be floored at zero and the single-letter
value stops being the last word.

### 14.8 [MSc] How do I quantify the coherence of a state in a fixed basis by the relative-entropy $C_{\mathrm{rel}}(\rho)=S(\rho_{\mathrm{diag}})-S(\rho)$ or $\ell_1$ measure $C_{\ell_1}(\rho)=\sum_{i\neq j}|\rho_{ij}|$?

Coherence is not a property of a state alone. Fixing a basis splits the states into the diagonal ones, which
carry none of it, and the rest, and both measures above are ways of asking how far $\rho$ sits from its own
diagonal part $\rho_{\mathrm{diag}} = \Delta(\rho)$, where $\Delta$ deletes the off-diagonal entries. That map is
a channel, the one that measures in the chosen basis and forgets the outcome, so the question is how much is lost
when it acts.

The relative-entropy measure is not a new quantity at all. Since $\Delta(\rho)$ and $\rho$ share their diagonal,
the two terms of $S(\Delta\rho) - S(\rho)$ reassemble into the divergence of 14.2 against the diagonal state,

$$C_{\mathrm{rel}}(\rho) = S(\rho \,\|\, \Delta(\rho)),$$

so Klein's inequality already gives $C_{\mathrm{rel}} \ge 0$ with equality only on the diagonal states, and no
separate positivity argument is needed. The $\ell_1$ measure has no such reading: it is a sum of moduli of matrix
entries, fixed by the basis rather than derived from an entropy, which makes it the easier of the two to compute
and the harder to connect to anything else in the Part.

For a qubit write the Bloch vector in polar form, $\vec r = \rho(\sin\theta, 0, \cos\theta)$, so that $\rho$ is
the length and $\theta$ the tilt away from the axis the basis picks out. The two measures then depend on those
two numbers differently, which is what makes them inequivalent rather than two scales for one thing.

**WL** : the generic mixed qubit, the entropy off its spectrum, and the binary entropy the diagonal states
reduce to.

```wl
bloch[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2
```

```wl
vonNeumann[m_] := -# . Log[2, #] & @ DeleteCases[Eigenvalues[m], 0]
```

```wl
h2[x_] := -x Log[2, x] - (1 - x) Log[2, 1 - x]
```

The incoherent projection is the diagonal, and the state under study is tilted off the axis so that it carries
coherence to measure.

```wl
delta[m_] := DiagonalMatrix[Diagonal[m]]
```

```wl
rho = bloch[rad {Sin[th], 0, Cos[th]}]
```

The relative-entropy measure from its definition, as a difference of two binary entropies of the two lengths
that matter, the tilted projection and the full radius.

```wl
FullSimplify[vonNeumann[delta[rho]] - vonNeumann[rho] == h2[(1 + rad Cos[th])/2] - h2[(1 + rad)/2], 0 < rad < 1 && 0 < th < Pi]
```

That it is the divergence against the diagonal state is an identity, so the nonnegativity of 14.2 carries over
without a new argument.

```wl
relativeEntropy[a_, b_] := Tr[a . (MatrixLog[a] - MatrixLog[b])]/Log[2]
```

```wl
FullSimplify[vonNeumann[delta[rho]] - vonNeumann[rho] - relativeEntropy[rho, delta[rho]], 0 < rad < 1 && 0 < th < Pi]
```

The $\ell_1$ measure reads off the equator alone, seeing the radius and the tilt only through their product.

```wl
FullSimplify[Total[Abs /@ Flatten[rho - delta[rho]]] == rad Sin[th], 0 < rad < 1 && 0 < th < Pi]
```

Both vanish on the axis, where the state is already diagonal, identically in the radius rather than at a chosen
length.

```wl
FullSimplify[{h2[(1 + rad Cos[th])/2] - h2[(1 + rad)/2], rad Sin[th]} /. th -> 0, 0 < rad < 1]
```

Two states carrying the same $\ell_1$ coherence, which the other measure separates, so neither quantity is a
function of the other.

```wl
FullSimplify[{{1/2, h2[(1 + 1/2 Cos[Pi/2])/2] - h2[(1 + 1/2)/2]}, {1/Sqrt[2] Sin[Pi/4], h2[(1 + Cos[Pi/4]/Sqrt[2])/2] - h2[(1 + 1/Sqrt[2])/2]}}]
```

**QF** : the incoherent projection is a named channel, full phase damping, and the entropies are properties of
the states it returns.

```wl
qrho = QuantumState[bloch[rad {Sin[th], 0, Cos[th]}]]
```

```wl
qdelta = QuantumChannel["PhaseDamping"[1]]
```

```wl
FullSimplify[Normal @ qdelta[qrho]["DensityMatrix"] == delta[rho], 0 < rad < 1 && 0 < th < Pi]
```

The measure is then the difference of two entropies, each read off its own object.

```wl
qcrel = FullSimplify[qdelta[qrho]["VonNeumannEntropy", 2] - qrho["VonNeumannEntropy", 2], 0 < rad < 1 && 0 < th < Pi]
```

```wl
FullSimplify[qcrel - (vonNeumann[delta[rho]] - vonNeumann[rho]), 0 < rad < 1 && 0 < th < Pi]
```

The $\ell_1$ measure is a sum over matrix entries, but the state hands over its own entries and its trace, and
the diagonal of a density matrix is real and nonnegative, so taking moduli leaves the diagonal alone and
subtracting the trace removes exactly it. No reference to the basis projection is needed and the form does not
care about the dimension.

```wl
cl1[qs_] := Total[Abs[qs["AmplitudesList"]]] - Tr[qs]
```

```wl
FullSimplify[cl1[qrho] == rad Sin[th], 0 < rad < 1 && 0 < th < Pi]
```

Coherence is a statement about a state and a basis together, and the diagonal states are the ones the basis
declares to have none. The relative-entropy measure turns out to be the divergence of 14.2 against the state's
own diagonal, so it is nonnegative for the reason Klein's inequality is true and vanishes exactly on the
incoherent states, and it is the quantity that behaves well under the operations that cannot create coherence.
The $\ell_1$ measure agrees with it on those two boundary facts and disagrees everywhere in between: it sees the
Bloch vector only through the length of its equatorial part, so it cannot tell a nearly pure tilted state from a
short one lying flat, while the entropic measure prices the purity as well as the tilt. Two states with the same
$\ell_1$ coherence and different relative-entropy coherence settle that neither is a function of the other, which
is why the coherence of a state is a family of answers rather than one number, and why a resource theory has to
name its measure along with its basis.
