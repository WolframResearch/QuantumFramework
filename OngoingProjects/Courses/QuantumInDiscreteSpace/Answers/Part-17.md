## Part 17. Tomography, estimation, and metrology

### 17.1 [MSc] How do I reconstruct an unknown density matrix by state tomography?

A qubit state has three real parameters and a projective measurement returns one bit, so a single setting cannot
determine $\rho$ no matter how many copies are spent on it. Measuring $Z$ alone fixes $r_z$ and leaves the other
two components untouched, which is why tomography is a set of experiments. What is needed is a collection of
POVM elements spanning the Hermitian operators, so that the map from state to outcome probabilities can be
inverted; such a set is informationally complete. Three Pauli bases supply six projectors for a qubit and are
more than enough.

Once the set is complete the inversion is linear. Each outcome probability is $p_k = \mathrm{Tr}[\rho E_k]$,
linear in the unknown, so stacking the $E_k$ and applying the pseudo-inverse to the observed frequencies returns
$\rho$. For a qubit measured in the three Pauli bases that machinery collapses to something recognizable, since
$\langle \sigma_k \rangle = p_+ - p_-$ in each setting and

$$\rho = \frac{1}{2}\left(I + \sum_k \langle \sigma_k \rangle \, \sigma_k\right),$$

so the Pauli expansion *is* the tomographic inversion, with the expectation values read off as count
differences.

Two things separate that identity from a working protocol. Frequencies are not probabilities: with $N$ copies per
setting each $\langle\sigma_k\rangle$ carries a statistical error of order $1/\sqrt{N}$, so the reconstruction
improves only as the square root of the effort. And nothing in the linear solve knows that $\rho$ must be
positive. The estimated Bloch vector is built from three independently noisy numbers, and for a state near the
boundary the fluctuations can push $|\vec r_{\mathrm{est}}| > 1$, which is not a state at all. That is the common
case for nearly pure states and modest $N$, and it is why the linear estimate must be followed by a projection
back onto the state space or replaced by a likelihood maximization over physical states.

**WL** : the Bloch parametrization, and the inversion identity.

```wl
bloch[v_] := (IdentityMatrix[2] + v . PauliMatrix[{1, 2, 3}])/2
```

The three expectation values determine the state identically in the Bloch vector, which is the statement that
the three Pauli bases are informationally complete for a qubit.

```wl
FullSimplify[bloch[Table[Tr[bloch[{rx, ry, rz}] . PauliMatrix[k]], {k, 3}]] == bloch[{rx, ry, rz}], Element[{rx, ry, rz}, Reals]]
```

An unknown to be reconstructed, off every axis so all three settings carry signal, and with $|\vec r|$ just
under $1$ so the positivity constraint binds.

```wl
rTrue = {2/5, 3/10, 17/20}
```

```wl
truth = bloch[rTrue]
```

The two projectors of each setting, with the $-1$ outcome first, which is the order the framework's POVM
elements use.

```wl
proj[k_] := {(IdentityMatrix[2] - PauliMatrix[k])/2, (IdentityMatrix[2] + PauliMatrix[k])/2}
```

```wl
probs = Table[Tr[truth . #] & /@ proj[k], {k, 3}]
```

Twenty copies per setting.

```wl
SeedRandom[1]
counts = RandomVariate[MultinomialDistribution[20, N[#]]] & /@ probs
```

Each expectation value is the count difference in its setting, and the three together are the estimated Bloch
vector.

```wl
rEst = (#[[2]] - #[[1]])/Total[#] & /@ counts
```

Whether that is a state at all is a question about its length.

```wl
N @ Norm[rEst]
```

```wl
Norm[rEst] <= 1
```

Rescaling to the boundary is the cheapest repair.

```wl
rFix = rEst/Max[1, Norm[rEst]]
```

The statistical error is the distance from the true Bloch vector, and reading its trend needs an average over
many runs at each number of copies.

```wl
error[n_] := Norm[(#[[2]] - #[[1]])/Total[#] & /@ (RandomVariate[MultinomialDistribution[n, N[#]]] & /@ probs) - rTrue]
```

```wl
SeedRandom[2]
N @ Table[{n, Mean @ Table[error[n], 200]}, {n, {20, 200, 2000, 20000}}]
```

**QF** : the measurements are named operators, one per Pauli, and the framework carries the whole estimation
pipeline.

```wl
pauli = QuantumMeasurementOperator[QuantumOperator[#]] & /@ {"X", "Y", "Z"}
```

```wl
qtruth = QuantumState["BlochVector"[rTrue]]
```

Two properties of the same operator list its outcomes in opposite orders: the probability list of the $Z$
measurement begins with the probability of $|0\rangle$, while its POVM elements begin with the projector onto
$|1\rangle$.

```wl
Table[q[qtruth]["ProbabilitiesList"], {q, pauli}]
```

Reversed onto the POVM order, they are the WL side's probabilities, which is the statement that the two sides
describe one experiment.

```wl
Reverse /@ Table[q[qtruth]["ProbabilitiesList"], {q, pauli}] === probs
```

`QuantumMeasurementSimulation` draws the counts, taking the state and the operators and returning the
association the estimator consumes.

```wl
SeedRandom[1]
data = QuantumMeasurementSimulation[qtruth, pauli, 20]
```

```wl
est = QuantumStateEstimate[data]
```

The estimator reports whether the linear solve was invertible and whether what it returned was a state.

```wl
{est["Invertible"], est["PhysicalInversion"], est["TotalCounts"], est["CountsPerMeasurement"]}
```

Two estimates come back: the linear inversion made physical, and the state maximizing the likelihood of the
observed counts.

```wl
{est["InversionState"]["BlochVector"], est["MaximumLikelihoodState"]["BlochVector"]} // N // Chop
```

Handing the estimator the counts the WL side drew separates the inversion map from the sample it was fed, and
that map is the same rescaling.

```wl
Chop @ N[rFix - QuantumStateEstimate[AssociationThread[pauli, counts]]["InversionState"]["BlochVector"]]
```

`QuantumSimilarity` with its default `"Fidelity"` measure returns the root fidelity
$\mathrm{Tr}\sqrt{\sqrt{\rho}\,\sigma\sqrt{\rho}}$, which is $|\langle a|b\rangle|$ on pure states, so the two
numbers below are $\sqrt{F}$ and not $F$.

```wl
N @ {QuantumSimilarity[qtruth, est["InversionState"]], QuantumSimilarity[qtruth, est["MaximumLikelihoodState"]]}
```

Raising the number of copies drives the statistical error down, and on a large sample the same pipeline returns
the state it started from.

```wl
SeedRandom[3]
With[{big = QuantumStateEstimate[QuantumMeasurementSimulation[qtruth, pauli, 200000]]}, {big["InversionState"]["BlochVector"] // N // Chop, N @ QuantumSimilarity[qtruth, big["InversionState"]]}]
```

Tomography is a linear solve wrapped in two problems the solve does not see. Informational completeness is what
makes the solve possible, and for a qubit the three Pauli bases deliver it with the inversion collapsing to the
Pauli expansion, so the reconstruction is three expectation values read as count differences. Finite statistics
then make each of those a noisy number whose error falls only as $1/\sqrt{N}$, and because the three are
estimated independently the resulting vector need not lie in the Bloch ball at all. Near the boundary the length
test and QF's `"PhysicalInversion"` report the same failure, and the repair is either a projection back onto the
state space or a likelihood maximization restricted to it. For a qubit the framework's projection is the
rescaling to the boundary, while the maximum-likelihood state is a different point of the ball, so the cost of
the positivity constraint shows up as the gap between the two fidelities.

### 17.2 [MSc] How do I reconstruct an unknown channel by process tomography?

A channel is linear, so it is fixed by its action on a spanning set of operators, and the Choi state packages
that action as a single state on $d^2$ dimensions (Part 13, 13.3),

$$\rho_{\mathcal{E}} = \frac{1}{d}\sum_{i,j}|i\rangle\langle j| \otimes \mathcal{E}(|i\rangle\langle j|)
= (\mathrm{id}\otimes\mathcal{E})\big(|\Phi^+\rangle\langle\Phi^+|\big).$$

That fixes the whole procedure. Probe the channel with a spanning set of inputs, reconstruct each output by the
state tomography of Part 17, 17.1, assemble the Choi state from those outputs, and turn it back into a channel by
reading off Kraus operators from its spectrum. The last step is the one that makes this a reconstruction and not
just a table of measurements: what comes out is an operator you can apply to a state the channel was never probed
with.

For a qubit four inputs suffice, and they have to span the Hermitian operators rather than merely be distinct
states. A spanning set of *states* does not contain $|0\rangle\langle 1|$, which is not a state at all, so its
image is never measured; linearity supplies it from
$|0\rangle\langle 1| = P_+ + iP_{+i} - \tfrac{1+i}{2}(P_0 + P_1)$.

Two things go wrong at finite counts, and they pull in opposite directions. The assembled Choi keeps unit trace,
because each reconstructed output is a state, so the map is trace preserving. It need not be positive, and a
negative eigenvalue means the reconstruction is not completely positive. Discarding the negative part is the
repair, the analogue of rescaling to the Bloch ball, but it removes weight and so breaks trace preservation,
which then has to be restored. Complete positivity and trace preservation have to be imposed together.

**WL** : a channel whose damping is drawn at random and never inspected.

```wl
SeedRandom[11]
gammaTrue = RandomReal[];
kraus = {{{1, 0}, {0, Sqrt[1 - gammaTrue]}}, {{0, Sqrt[gammaTrue]}, {0, 0}}};
map[s_] := Total[# . s . ConjugateTranspose[#] & /@ kraus]
```

The four probes.

```wl
{p0, p1, pplus, pi} = {{{1, 0}, {0, 0}}, {{0, 0}, {0, 1}},
   (IdentityMatrix[2] + PauliMatrix[1])/2, (IdentityMatrix[2] + PauliMatrix[2])/2}
```

One output tomographed: three Pauli settings, count differences, and the rescaling of Part 17, 17.1.

```wl
proj[k_] := {(IdentityMatrix[2] - PauliMatrix[k])/2, (IdentityMatrix[2] + PauliMatrix[k])/2};
tomo[s_, n_] := Module[{out = map[s], probs, counts, r},
   probs = Table[Re[Tr[out . #]] & /@ proj[k], {k, 3}];
   counts = Table[RandomVariate[MultinomialDistribution[n, N[probs[[k]]]]], {k, 3}];
   r = N[(#[[2]] - #[[1]])/Total[#] & /@ counts];
   (IdentityMatrix[2] + (r/Max[1, Norm[r]]) . PauliMatrix[{1, 2, 3}])/2]
```

```wl
est = tomo[#, 2000] & /@ {p0, p1, pplus, pi}
```

The four outputs assembled into the Choi state, with the unmeasured image supplied by linearity. Since
$|i\rangle\langle j| \otimes A$ places $A$ in block $(i,j)$, the assembly is just the images stacked as blocks.

```wl
e01 = est[[3]] + I est[[4]] - ((1 + I)/2) (est[[1]] + est[[2]]);
J = ArrayFlatten[{{est[[1]], e01}, {ConjugateTranspose[e01], est[[2]]}}]/2;
Chop @ {Tr[J], Eigenvalues[J]}
```

Kraus operators are the eigenvectors of the Choi state, unstacked into matrices and weighted by the square root
of the eigenvalues. Negative eigenvalues are dropped, and the surviving set is renormalized so that it closes the
completeness relation again.

```wl
krausFrom[lam_, vec_] := Module[{kr, m},
   kr = Select[MapThread[Sqrt[2 Max[Re[#1], 0]] Transpose[ArrayReshape[Normalize[#2], {2, 2}]] &, {lam, vec}],
     Max[Abs[#]] > 10^-6 &];
   m = MatrixPower[Total[ConjugateTranspose[#] . # & /@ kr], -1/2];
   # . m & /@ kr]
```

```wl
kEst = krausFrom @@ Eigensystem[J];
{Length[kEst], Chop @ Max @ Abs[Total[ConjugateTranspose[#] . # & /@ kEst] - IdentityMatrix[2]]}
```

The reconstruction is now a map, and the test is a state it was never probed with.

```wl
mapEst[s_] := Total[# . s . ConjugateTranspose[#] & /@ kEst];
With[{fresh = (IdentityMatrix[2] + N[{2/5, -3/5, 1/2}] . PauliMatrix[{1, 2, 3}])/2},
  Chop @ Max @ Abs[mapEst[fresh] - map[fresh]]]
```

**QF** : the same unknown channel as an object, probed with named states.

```wl
ch = QuantumChannel[kraus]
```

```wl
inputs = QuantumState /@ {"0", "1", "Plus", "Right"}
```

Each output is tomographed by the pipeline of Part 17, 17.1, which returns a physical state for each.

```wl
pauli = QuantumMeasurementOperator[QuantumOperator[#]] & /@ {"X", "Y", "Z"};
outs = QuantumStateEstimate[QuantumMeasurementSimulation[ch[#], pauli, 2000]] & /@ inputs
```

The four matrix units are the ordered pairs of computational kets against bras, each pairing with the image of
that unit. Read as a matrix state, the assembly is the Choi state.

```wl
units = QuantumOperator[#1 @ #2["Dagger"]] & @@@ Tuples[QuantumState /@ {"0", "1"}, 2]
```

```wl
rho = #["MaximumLikelihoodState"]["Operator"] & /@ outs;
q01 = rho[[3]] + I rho[[4]] - ((1 + I)/2) (rho[[1]] + rho[[2]]);
choiEst = (Inner[QuantumTensorProduct, units, {rho[[1]], q01, q01["Dagger"], rho[[2]]}, Plus]/2)["MatrixQuantumState"]
```

```wl
Chop @ N @ choiEst["Eigenvalues"]
```

Its own spectrum feeds the same extraction, and the result is a channel.

```wl
chEst = QuantumChannel[krausFrom[choiEst["Eigenvalues"], choiEst["Eigenvectors"]]]
```

```wl
{chEst["TracePreservingQ"], Chop @ N @ (QuantumChannel[chEst, {2}] @ QuantumState["PhiPlus"])["Eigenvalues"]}
```

The reconstruction is compared with the channel it came from on states that were never probes.

```wl
N @ (QuantumSimilarity[chEst[#], ch[#]] & /@ {QuantumState["BlochVector"[{2/5, -3/5, 1/2}]], QuantumState["Minus"], QuantumState["1"]})
```

Process tomography is a pipeline with four stages: probe with a spanning set, tomograph each output, assemble the
Choi state, and read the Kraus operators off its spectrum. The middle stage is Part 17, 17.1 unchanged, and the
outer two are what make it a channel reconstruction. The statistics spoil the third stage in a way the state
problem does not prepare you for. Unit trace survives, since each output was a state, so trace preservation comes
free. Positivity does not survive, and a negative eigenvalue means the assembled object is not a channel at all,
even though every output it was built from was individually repaired into a legitimate state. Dropping the
negative part restores complete positivity but removes weight, which breaks the trace preservation that had been
free, so the Kraus set has to be renormalized afterwards. That is the honest content of the last stage: the two
constraints are not independent, and a serious implementation fits the Choi state under both at once instead of
repairing one and then the other. What comes out here is nonetheless a genuine channel, and it agrees with the
unknown one on states that played no part in building it.

*A note on the outcome ordering.* For these eigen-decomposed Pauli operators `"ProbabilitiesList"` runs opposite
to `"POVMElements"`, as the two cells above show. Both `QuantumMeasurementSimulation` and `QuantumStateEstimate`
work from `"POVMElements"`, so the shipped pipeline is self-consistent and the ordering matters only when counts
are built by hand: index them by `"ProbabilitiesList"` and the reconstruction returns the antipodal Bloch vector,
$-\vec r$ in place of $\vec r$, with no warning.

### 17.3 [MSc] How do I compute the quantum Fisher information and the Cramer-Rao bound for a phase-estimation problem?

Tomography reconstructs an object. Metrology asks a narrower question: one number $\phi$ is imprinted on a probe
by $U_\phi = e^{-i\phi G}$, and we want to know how well it can be estimated, and by which measurement.

Fix a measurement with outcomes $k$ and probabilities $p_k(\phi)$. Its Fisher information is

$$F(\phi) = \sum_k \frac{\left(\partial_\phi p_k\right)^2}{p_k},$$

and the Cramer-Rao bound says any unbiased estimator built from $N$ repetitions obeys
$\mathrm{Var}(\hat\phi) \ge 1/(N F)$. Fisher information is the curvature of the outcome distribution against the
parameter: if the probabilities barely move as $\phi$ changes, $F$ is small and no estimator can be sharp.

That much is classical and depends on the measurement chosen. The quantum Fisher information is its maximum over
every POVM, so it is a property of the state family alone. For a pure family it is

$$F_Q = 4\left(\langle\partial_\phi\psi|\partial_\phi\psi\rangle - \left|\langle\psi|\partial_\phi\psi\rangle\right|^2\right)
= 4\,\mathrm{Var}_\psi(G),$$

the second equality holding when the parameter enters through $U_\phi = e^{-i\phi G}$. The bracket is exactly the
Fubini-Study metric, so $F_Q = 4g$ and the quantum Cramer-Rao bound reads $\mathrm{Var}(\hat\phi)\ge 1/(N F_Q)$.
The generator being a variance is the whole of metrology in one line: a probe helps only insofar as it is spread
out in the quantity conjugate to the parameter.

**WL** : a qubit probe in $|+\rangle$ with the phase written on it by a rotation about $z$, so $G = Z/2$.

```wl
psi = MatrixExp[-I \[Phi] PauliMatrix[3]/2] . {1, 1}/Sqrt[2]
```

Born probabilities of a Pauli measurement, taken from the expectation value.

```wl
born[m_] := Simplify[ComplexExpand[{1 - #, 1 + #}/2] &[Conjugate[psi] . m . psi], Element[{\[Phi], \[Alpha]}, Reals]]
```

```wl
fisher[p_] := Simplify[Total[D[p, \[Phi]]^2/p], Element[{\[Phi], \[Alpha]}, Reals]]
```

The three Pauli measurements do not carry the same information about $\phi$.

```wl
{born[#], fisher[born[#]]} & /@ PauliMatrix[{1, 2, 3}]
```

The quantum Fisher information needs no measurement at all, only the generator.

```wl
With[{G = PauliMatrix[3]/2},
  Simplify[4 (Conjugate[psi] . G . G . psi - (Conjugate[psi] . G . psi)^2), Element[\[Phi], Reals]]]
```

Both numbers came out free of $\phi$, for two different reasons, and only one of them is general. Tilting the
measurement axis out of the equatorial plane by $\alpha$ interpolates between the two Pauli extremes and shows
what a generic measurement does.

```wl
fisher[born[Sin[\[Alpha]] PauliMatrix[1] + Cos[\[Alpha]] PauliMatrix[3]]]
```

```wl
{#, Simplify[fisher[born[Sin[\[Alpha]] PauliMatrix[1] + Cos[\[Alpha]] PauliMatrix[3]]] /. \[Alpha] -> #]} & /@ {0, Pi/4, Pi/2}
```

**QF** : the same family, with the parameter declared so the framework can differentiate it.

```wl
probe[n_, st_] := QuantumState[QuantumCircuitOperator[Table["RZ"[\[Phi]] -> k, {k, n}]][st], "Parameters" -> {\[Phi]}]
```

`FubiniStudyMetricTensor` is that bracket, so four times it is the quantum Fisher information.

```wl
Simplify[4 FubiniStudyMetricTensor[probe[1, QuantumState["Plus"]], "Matrix"][[1, 1]], Element[\[Phi], Reals]]
```

Spending $n$ probes independently, against spending them on one entangled probe.

```wl
Table[{n,
   Simplify[4 FubiniStudyMetricTensor[probe[n, QuantumState[StringRepeat["+", n]]], "Matrix"][[1, 1]], Element[\[Phi], Reals]],
   Simplify[4 FubiniStudyMetricTensor[probe[n, QuantumState["GHZ"[n]]], "Matrix"][[1, 1]], Element[\[Phi], Reals]]},
  {n, 4}]
```

The single qubit fixes the scale: $F_Q = 1$, so one probe buys a variance floor of one radian squared and $N$
repetitions bring it to $1/N$. Two of the three Pauli measurements reach that floor and the third learns nothing,
which is the sharpest way to see that Fisher information belongs to the pairing of state with measurement while
the quantum Fisher information belongs to the state alone. Measuring along $z$ is blind because the phase rotates
the Bloch vector about $z$ and never moves its $z$ component, so the outcome distribution is flat in $\phi$.

The absence of $\phi$ from those first answers is worth separating into its two causes. $F_Q$ is constant for a
reason that always applies to a phase: the generator commutes with the evolution it generates, so
$\mathrm{Var}(G)$ is the same on every member of the family and the bound cannot depend on where in the orbit the
probe sits. The constancy of the equatorial measurements is not general at all. The tilted family gives

$$F(\alpha,\phi) = \frac{\sin^2\alpha\,\sin^2\phi}{1 - \sin^2\alpha\cos^2\phi},$$

which does depend on $\phi$ everywhere except at the endpoints $\alpha = 0$ and $\alpha = \pi/2$. At
$\alpha = \pi/4$ it vanishes at $\phi = 0$ and peaks at $1/2$, half the bound, so this measurement is both
suboptimal and worst exactly where the estimate would be needed. A measurement that saturates the quantum
Cramer-Rao bound at every $\phi$ is the exception, and finding one is the design problem metrology is about.

The scaling is the reason the question is worth asking. Independent probes give $F_Q = n$, so the precision
improves as $1/\sqrt n$, the standard quantum limit, and it is exactly the improvement of averaging $n$
independent trials. A GHZ probe gives $F_Q = n^2$ and precision $1/n$, the Heisenberg limit. The generator on the
entangled probe is the same sum $\tfrac12\sum_i Z_i$, but the GHZ state is spread across its two extreme
eigenvalues rather than across each qubit separately, and the variance of that sum grows as $n^2$ instead of $n$.
Entanglement buys a factor $n$ in variance, and the $4\,\mathrm{Var}(G)$ formula says so directly.

### 17.4 [MSc] How do I estimate average gate fidelity by randomized benchmarking, fitting the survival probability to $A p^m + B$ and reading the depolarizing rate from the decay $p$?

Process tomography reconstructs a channel completely and pays for it: Part 17, 17.2 needed four state
tomographies, and every error in preparing and measuring the probes went straight into the answer. Calibrating a
gate set does not need the whole channel, only one number, and randomized benchmarking gets that number without
caring about preparation or measurement.

The protocol: apply $m$ Clifford gates drawn uniformly at random, then the one Clifford that inverts their
product, and record whether the system came back to where it started. Noiseless gates return survival $1$ at
every $m$. Real gates decay, and the claim is that the decay is a single exponential

$$P_{\mathrm{survive}}(m) = A\,p^m + B .$$

A single exponential appears whatever the noise is, because of the Clifford twirl. Averaging any channel over the
Clifford group, $\overline{\mathcal{E}} = |\mathcal{C}|^{-1}\sum_C \mathcal{C}^\dagger \circ \mathcal{E} \circ \mathcal{C}$,
produces a depolarizing channel, and depolarizing shrinks the Bloch vector by a fixed factor per step.
Randomizing the sequence is what carries out that average, so $p$ is the shrink factor of the twirled noise, that
is the surviving depolarizing parameter is $1 - p$. From it the error per Clifford and the average gate fidelity
follow for dimension $d$:

$$r = \frac{(d-1)(1-p)}{d} = \frac{1-p}{2}, \qquad F_{\mathrm{avg}} = 1 - r = \frac{1+p}{2} \quad (d=2).$$

Everything else about the noise is invisible, and $A$ and $B$ absorb the state preparation and measurement errors.
That is the trade: one number instead of a channel, in exchange for immunity to the errors that dominate a
tomography.

**WL** : the single-qubit Clifford group, obtained by closing the two generators $H$ and $S$ under multiplication,
with global phases quotiented out so that the group has its 24 elements and not more.

```wl
canonical[m_] := Module[{z = First @ Select[Flatten[m], Abs[#] > 10^-9 &]}, Simplify[m z^-1 Abs[z]]];
cliff = FixedPoint[
   DeleteDuplicates[canonical /@ Flatten[{#, # . {{1, 1}, {1, -1}}/Sqrt[2] & /@ #, # . {{1, 0}, {0, I}} & /@ #}, 1]] &,
   {IdentityMatrix[2]}];
Length[cliff]
```

The noise to be benchmarked, applied after every gate. Amplitude damping is chosen because it is neither
depolarizing nor even unital, so the twirl has something to do.

```wl
damping = 1/5;
damp[s_] := Total[# . s . ConjugateTranspose[#] & /@ {{{1, 0}, {0, Sqrt[1 - damping]}}, {{0, Sqrt[damping]}, {0, 0}}}]
```

Any qubit channel acts on the Bloch vector as $\vec r \mapsto M\vec r + \vec c$. This reads off the pair by
sending in the three axis states and the centre.

```wl
blochOf[f_] := Module[{c = Table[Tr[f[IdentityMatrix[2]/2] . PauliMatrix[k]], {k, 3}]},
   {Transpose @ Table[Table[Tr[f[(IdentityMatrix[2] + UnitVector[3, j] . PauliMatrix[{1, 2, 3}])/2] . PauliMatrix[k]], {k, 3}] - c, {j, 3}], c}]
```

For the damping the contraction is unequal along the three axes and the centre is pulled toward $|0\rangle$, so
this noise is a long way from depolarizing.

```wl
Simplify @ blochOf[damp]
```

Averaging the same channel over the group is the twirl. What comes back has no preferred axis and no shift: a
depolarizing channel, which is the theorem the whole protocol rests on.

```wl
Simplify @ blochOf[Mean[Table[ConjugateTranspose[c] . damp[c . # . ConjugateTranspose[c]] . c, {c, cliff}]] &]
```

One benchmarking sequence: start in $|0\rangle$, apply $m$ random Cliffords with the noise after each, then the
inverting Clifford followed by the noise once more, and report the probability of finding $|0\rangle$ again.

```wl
step[rho_, c_] := damp[c . rho . ConjugateTranspose[c]];
survival[m_] := Module[{seq = RandomChoice[N[cliff], m]},
   Re @ step[Fold[step, {{1., 0}, {0, 0}}, seq], ConjugateTranspose[Dot @@ Reverse[seq]]][[1, 1]]]
```

The decay curve is the survival averaged over many sequences at each length, which is what performs the twirl.

```wl
SeedRandom[7]
rbData = Table[{m, Mean[Table[survival[m], 1000]]}, {m, Range[1, 41, 4]}]
```

Fitting the three-parameter model separates the part that carries the gate error from the parts that do not.

```wl
rbFit = NonlinearModelFit[rbData, a p^m + b, {{a, .4}, {p, .86}, {b, .6}}, m];
{a, p, b} /. rbFit["BestFitParameters"]
```

Everything asked for follows from the decay constant alone.

```wl
pFit = p /. rbFit["BestFitParameters"];
<|"decay p" -> pFit, "depolarizing parameter" -> 1 - pFit,
  "error per Clifford" -> (1 - pFit)/2, "average gate fidelity" -> (1 + pFit)/2|>
```

**QF** : the same noise as a channel object, so the value the benchmark is chasing can be computed exactly and
compared against it.

```wl
noise = QuantumChannel["AmplitudeDamping"[damping]]
```

The average gate fidelity of any channel is fixed by the overlap of its Choi state (Part 17, 17.2) with the
maximally entangled state, the entanglement fidelity $F_{\mathrm{pro}}$, through
$F_{\mathrm{avg}} = (d F_{\mathrm{pro}} + 1)/(d+1)$.

```wl
avgFidelity[ch_] := Simplify[(2 QuantumSimilarity[QuantumChannel[ch, {2}] @ QuantumState["PhiPlus"], QuantumState["PhiPlus"]]^2 + 1)/3]
```

```wl
avgFidelity[noise]
```

The decay constant the benchmark should return is that fidelity read backwards.

```wl
Simplify[2 avgFidelity[noise] - 1]
```

The twirl does not merely resemble a depolarizing channel, it is one: the named channel carrying the fitted
depolarizing parameter has exactly the average fidelity of the damping it replaced.

```wl
{avgFidelity[noise], avgFidelity[QuantumChannel["Depolarizing"[1 - 4 (1 + Sqrt[5])/15]]]}
```

Measured against computed, the benchmark on one side and the channel object on the other.

```wl
N @ {pFit, Simplify[2 avgFidelity[noise] - 1]}
```

Randomized benchmarking is one scalar pulled out of a curve, and the curve is exponential only because the
randomization twirls the noise into a depolarizing channel. The raw damping contracts the Bloch ball unequally
and moves its centre; the twirled channel is a bare contraction, and the single factor left is what the fit
returns. From it the depolarizing parameter is $1-p$, the error per Clifford is half that, and the average gate
fidelity is one minus the error. Nothing else about the channel survives, which is both the limitation and the
point: the same averaging that discards the anisotropy discards the preparation and measurement errors too.

The fitted offsets deserve a look rather than a shrug. For unital noise $B$ would be $1/d$, the survival of a
maximally mixed state; here it lands near $0.6$, because damping pulls toward $|0\rangle$ and the last noise acts
after the final twirled step. That bias is absorbed by $B$ and leaves $p$ untouched. Reading a fidelity off $B$,
or off the survival at any single sequence length, would report non-unitality and readout error as if they were
gate error. Only $p$ is the gate.
