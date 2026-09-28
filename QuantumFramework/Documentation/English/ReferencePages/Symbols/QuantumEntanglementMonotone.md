---
Template: Symbol
Name: QuantumEntanglementMonotone
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumEntanglementMonotone
Keywords: [entanglement monotone, concurrence, negativity, log negativity, entanglement entropy, Renyi entropy, quantum discord, mutual information, classical correlation, bipartition, Schmidt decomposition]
SeeAlso: [QuantumState, QuantumPartialTrace, QuantumEntangledQ, QuantumMeasurementOperator]
RelatedTutorials: [[Wolfram Quantum Framework Tutorial](Tutorial)]
RelatedGuides: [WolframQuantumComputationFramework]
---

<!--
md source for Documentation/English/ReferencePages/Symbols/QuantumEntanglementMonotone.nb,
recovered from the legacy .nb (which had no .md twin) via NotebookToMarkdown, then
modernized and re-verified cell by cell against the working-tree kernel carrying the
entanglement-measures audit (commits db97269d, c9de71be on branch
claude/elegant-spence-d8ebb3). Deliberate deviations from the legacy .nb, each verified
against that kernel:
  1. The named-state list forms {"Name", param} parse as literal amplitude arrays on
     kernel 2.1.1 (silently, giving a wrong one-qudit state), which produced the bulk of
     the legacy page's error cascade. They are rewritten to the function-call form
     "Name"[param]: "RandomPure"[3] (three qubits), "Werner"[p], "RandomMixed"[2],
     "UniformMixture"[2]. The integer is the subsystem (qubit) count.
  2. The mixed-state EntanglementEntropy example now returns Indeterminate and issues
     QuantumEntanglementMonotone::mixedentropy: the audit made the reduced von Neumann
     and Renyi entropies pure-state-only measures, since the reduced entropy of a mixed
     state counts classical mixing, not entanglement. This is the one intended message on
     the page; it documents the guard and reinforces the caption.
  3. "Discord", "MutualInformationI", and "MutualInformationJ" (new public measures) are
     documented: three rows in the measures table, a Details note, and a Scope example
     showing that the discord depends on the chosen measurement.
  4. % / %% output references are replaced with named variables so each example unit reads
     self-contained. The two Applications plots wrap each measure in a NumericQ-guarded
     helper, so Plot's symbolic trial evaluates nothing and the build is message-free.
  5. The legacy "Generalizations and Extensions" section (the concurrence-vector example)
     is folded into Scope, and the empty "Options" section dropped: MarkdownToNotebook
     silently drops section titles outside its taxonomy, and this symbol has no options.
-->

## Usage

<code>[QuantumEntanglementMonotone]()[*qs*, *bipart*, *t*]</code> computes the entanglement monotone using the measure or metric *t* on the quantum state *qs* between the subsystems of the bipartition list *bipart*.

<code>[QuantumEntanglementMonotone]()[*qs*, *m*, *bipart*, *t*]</code> uses the [QuantumMeasurementOperator]() *m* on the second subsystem, for the measurement-dependent measures `"Discord"` and `"MutualInformationJ"`.

## Details & Options

- The following measures or metrics are supported:

|   |   |
|---|---|
| `"Concurrence"` | the generalized concurrence, exact for two qubits and for pure states of any dimension |
| `"Negativity"` | the negativity (Peres-Horodecki criterion) |
| `"LogNegativity"` | the logarithmic negativity (Peres-Horodecki criterion) |
| `"EntanglementEntropy"` | the von Neumann entropy of a reduced subsystem, for a pure state |
| `{"RenyiEntropy", α}` | the Renyi-*α* entropy of a reduced subsystem, for a pure state |
| `"Realignment"` | the realignment metric: the summed singular values of the realigned density matrix, minus one |
| `"ConcurrenceVector"` | the generalized concurrence vector, whose norm is the concurrence |
| `"MutualInformationI"` | the quantum mutual information $I(A:B) = S(\rho_A) + S(\rho_B) - S(\rho_{AB})$ |
| `"MutualInformationJ"` | the one-sided classical correlation $J$ for a given measurement on the second subsystem |
| `"Discord"` | the quantum discord relative to a given measurement, $I - J$ |

- The bipartition list *bipart* has the form <code>{{$i_{1}$, $i_{2}$, …}, {$j_{1}$, $j_{2}$, …}}</code>. For example, <code>{{1, 2}, {3}}</code> splits the state into the subsystem <code>{1, 2}</code> and the subsystem <code>{3}</code>. If no bipartition list *bipart* is given, <code>[qs]()["Bipartition"]</code> is performed first, which splits the system into its first subsystem and the rest (that is, <code>{{1}, {2, 3, …}}</code>).

- The measure or metric *t* is `"Concurrence"` by default, unless specified otherwise.

- For `"Concurrence"`, [QuantumEntanglementMonotone]() uses the Rungta-Buzek-Caves-Hillery-Milburn generalization of the Hill-Wootters concurrence (Phys. Rev. A 64, 042315 (2001)), which applies in any dimension. It is exact for two qubits and for pure states of any dimension, and a lower bound on the convex-roof concurrence for a higher-dimensional mixed state, where (unlike the two-qubit case) it is basis-dependent and not invariant under local unitaries. The value is the norm of the concurrence vector, whose components are $\max(0,\ \sigma_1 - \sum_{i \geq 2} \sigma_i)$ with $\sigma_i$ the singular values of $\omega = \sqrt{\rho}\ \sqrt{\tilde{\rho}}$, where $\tilde{\rho} = (Y_\alpha \otimes Y_\beta)\ \rho^{*}\ (Y_\alpha \otimes Y_\beta)$ is the spin-flipped state built from a pair of antisymmetric $\mathfrak{su}(d)$ generators. For a bipartition of dimensions $d_1 \times d_2$ the concurrence vector has $\frac{d_1 (d_1 - 1)}{2} \cdot \frac{d_2 (d_2 - 1)}{2}$ components. For a pure state its norm reduces to $\sqrt{2\,(1 - \operatorname{Tr}[\rho_A^{2}])}$, with $\rho_A$ the reduced density matrix of one part.

- If no *α* is given for `"RenyiEntropy"`, it is set to $\alpha = 1/2$. The alias `"RenyiEntanglementEntropy"` behaves identically.

- The `"EntanglementEntropy"` and `{"RenyiEntropy", α}` measures quantify entanglement only for a pure global state. For a genuinely mixed input they return [Indeterminate]() and issue a `QuantumEntanglementMonotone::mixedentropy` message, since the reduced von Neumann or Renyi entropy of a mixed state counts classical mixing rather than entanglement. Use `"Negativity"`, `"LogNegativity"`, or `"Concurrence"` for a mixed state.

- `"Discord"` returns $I(A:B) - J$ for the measurement supplied as the second argument (a [QuantumMeasurementOperator]()), or a projective computational-basis measurement on the second subsystem when none is given. It is therefore the discord for that fixed measurement, an upper bound on the Ollivier-Zurek quantum discord, which minimizes over all local measurements (equivalently maximizes $J$). The value depends on the chosen measurement, and the computational-basis default carries no special physical status. Likewise `"MutualInformationJ"` is the classical correlation for the given measurement, whose maximum over measurements is the Ollivier-Zurek classical correlation. (References: Ollivier & Zurek, Phys. Rev. Lett. 88, 017901 (2001); Henderson & Vedral, J. Phys. A 34, 6899 (2001).)

## Basic Examples

Compute the concurrence of a quantum state:

```wl
QuantumEntanglementMonotone[QuantumState[{0.6, 0, 0, 0.8}], "Concurrence"]
```

<!-- => 0.96 -->

---

Compute the entanglement entropy:

```wl
QuantumEntanglementMonotone[QuantumState[{0.6, 0, 0, 0.8}], "EntanglementEntropy"]
```

<!-- => Quantity[0.9426831892554921, "Bits"] -->

---

Compute the negativity:

```wl
QuantumEntanglementMonotone[QuantumState[{0.6, 0, 0, 0.8}], {{1}, {2}}, "Negativity"]
```

<!-- => 0.48 -->

---

Compute the logarithmic negativity:

```wl
QuantumEntanglementMonotone[QuantumState[{0.6, 0, 0, 0.8}], {{1}, {2}}, "LogNegativity"]
```

<!-- => 0.9708536543404835 -->

---

Compute the Renyi *α*-entropy with $\alpha = 0.25$:

```wl
QuantumEntanglementMonotone[QuantumState[{0.6, 0, 0, 0.8}], {{1}, {2}}, {"RenyiEntropy", .25}]
```

<!-- => 0.9853394393457473 -->

## Scope

The bipartition list need not span the whole state:

```wl
st = QuantumState["RandomPure"[3]];
logNeg13 = QuantumEntanglementMonotone[st, {{1}, {3}}, "LogNegativity"]
```

<!-- => a nonnegative real (random pure state), e.g. 0.31 -->

Tracing out the excluded subsystem first gives the same value:

```wl
logNeg13 == QuantumEntanglementMonotone[QuantumPartialTrace[st, {2}], "LogNegativity"]
```

<!-- => True -->

---

The parts of a bipartition may consist of several subsystems:

```wl
QuantumEntanglementMonotone[QuantumState["GHZ"], {{1, 2}, {3}}, "EntanglementEntropy"]
```

<!-- => Quantity[1, "Bits"] -->

---

[QuantumEntanglementMonotone]() works for a state of any dimension. Consider the maximally entangled two-qutrit state:

```wl
qutrit = QuantumState[{1/Sqrt[3], 0, 0, 0, 1/Sqrt[3], 0, 0, 0, 1/Sqrt[3]}, 3];
qutrit["Formula"]
```

<!-- => (|00> + |11> + |22>)/Sqrt[3] -->

```wl
QuantumEntanglementMonotone[qutrit, {{1}, {2}}]
```

<!-- => 2/Sqrt[3] -->

---

The `"ConcurrenceVector"` gives the per-generator-pair components whose norm is the concurrence. Generate a $3 \times 2 \times 2$ random pure state:

```wl
SeedRandom[1234];
st = QuantumState["RandomPure", {3, 2, 2}]
```

<!-- => QuantumState (12-dimensional, dimensions {3, 2, 2}) -->

Find its concurrence vector across the first subsystem and the rest:

```wl
cv = QuantumEntanglementMonotone[st, "ConcurrenceVector"]
```

<!-- => a length-18 list of reals -->

The vector has $\frac{d_1 (d_1 - 1)}{2} \cdot \frac{d_2 (d_2 - 1)}{2}$ components, with $d_1 = 3$ and $d_2 = 2 \times 2 = 4$ the dimensions of the two parts:

```wl
(3 (3 - 1))/2 (4 (4 - 1))/2 == Length[cv]
```

<!-- => True -->

Its norm is the concurrence:

```wl
Norm[cv]
```

<!-- => 0.9075219094957584 -->

Since the state is pure, the norm equals $\sqrt{2\,(1 - \operatorname{Tr}[\rho_A^{2}])}$, with $\rho_A$ the reduced density matrix of one part:

```wl
Chop[With[{r = QuantumPartialTrace[st["Bipartition"], {1}]["DensityMatrix"]},
   Sqrt[2 (1 - Tr[r . r])]] - Norm[cv], 10^-6] == 0
```

<!-- => True -->

---

The measures `"MutualInformationI"`, `"MutualInformationJ"`, and `"Discord"` separate total, classical, and quantum correlations. Take the classically correlated two-qubit state with density matrix $\mathrm{diag}(1/2, 0, 0, 1/2)$:

```wl
classical = QuantumState[DiagonalMatrix[{1/2, 0, 0, 1/2}], 2]
```

<!-- => QuantumState (two qubits, mixed) -->

Its total correlation is one bit of mutual information:

```wl
QuantumEntanglementMonotone[classical, "MutualInformationI"]
```

<!-- => Quantity[1, "Bits"] -->

The default measurement is a projective measurement in the computational basis of the second qubit. The state is perfectly correlated in that basis, so all of its correlation is classical and the discord vanishes:

```wl
QuantumEntanglementMonotone[classical, "Discord"]
```

<!-- => Quantity[0, "Bits"] -->

Measuring the second qubit in a different basis leaves correlations that no single measurement outcome accounts for, so the discord for that measurement is nonzero. The value depends on the measurement, and the computational basis is not privileged:

```wl
QuantumEntanglementMonotone[classical, QuantumMeasurementOperator["X"], "Discord"]
```

<!-- => Quantity[1, "Bits"] -->

## Applications

The entanglement of a `"W"` state is more robust than `"GHZ"` entanglement, because a `"W"` state is not bi-separable:

```wl
QuantumEntanglementMonotone[QuantumState["W"]]
```

<!-- => (2 Sqrt[2])/3 -->

---

Plot the concurrence and the logarithmic negativity of a Werner state, and read off the range of the parameter for which the state is separable:

```wl
wernerLogNeg[p_?NumericQ] := QuantumEntanglementMonotone[QuantumState["Werner"[p]], {{1}, {2}}, "LogNegativity"];
wernerConc[p_?NumericQ] := QuantumEntanglementMonotone[QuantumState["Werner"[p]], {{1}, {2}}, "Concurrence"];
Plot[{wernerLogNeg[p], wernerConc[p]}, {p, 0, 1},
 PlotLegends -> {"LogNegativity", "Concurrence"}]
```

<!-- => Graphics; both measures fall to 0 as the Werner state becomes separable -->

---

Plot several entanglement monotones for the state $\beta |00\rangle + \sqrt{1 - \beta^2}\, |11\rangle$:

```wl
betaConc[b_?NumericQ] := QuantumEntanglementMonotone[QuantumState[{b, 0, 0, Sqrt[1 - b^2]}], {{1}, {2}}, "Concurrence"];
betaLogNeg[b_?NumericQ] := QuantumEntanglementMonotone[QuantumState[{b, 0, 0, Sqrt[1 - b^2]}], {{1}, {2}}, "LogNegativity"];
betaEntEntropy[b_?NumericQ] := QuantumEntanglementMonotone[QuantumState[{b, 0, 0, Sqrt[1 - b^2]}], {{1}, {2}}, "EntanglementEntropy"];
Plot[{betaConc[b], betaLogNeg[b], betaEntEntropy[b]}, {b, 0, 1},
 PlotLegends -> {"Concurrence", "LogNegativity", "Entanglement Entropy"}]
```

<!-- => Graphics; three curves, all vanishing at the separable endpoints beta = 0, 1 -->

## Properties and Relations

The concurrence is defined for a mixed state:

```wl
QuantumEntanglementMonotone[QuantumState["RandomMixed"[2]]]
```

<!-- => a small nonnegative real (random mixed state) -->

The entanglement entropy is not defined for a mixed state: on a genuinely mixed input, `"EntanglementEntropy"` returns [Indeterminate]() and issues a `QuantumEntanglementMonotone::mixedentropy` message, since the reduced entropy of a mixed state measures classical mixing rather than entanglement:

```wl
QuantumEntanglementMonotone[QuantumState["UniformMixture"[2]], "EntanglementEntropy"]
```

<!-- => QuantumEntanglementMonotone::mixedentropy message, then Indeterminate -->

---

For a qubit pair, [QuantumEntanglementMonotone]() gives 1 for a maximally entangled state:

```wl
QuantumEntanglementMonotone[QuantumState["PsiPlus"]]
```

<!-- => 1 -->

It gives 0 for a separable state:

```wl
QuantumEntanglementMonotone[QuantumState["01"]]
```

<!-- => 0 -->
