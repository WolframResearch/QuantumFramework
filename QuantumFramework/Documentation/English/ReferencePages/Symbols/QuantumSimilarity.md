---
Template: Symbol
Name: QuantumSimilarity
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumSimilarity
Keywords: [quantum similarity, fidelity, state comparison, similarity measure, complement of distance]
SeeAlso: [QuantumDistance, QuantumState]
RelatedTutorials: [[Wolfram Quantum Framework Tutorial](Tutorial)]
RelatedGuides: [WolframQuantumComputationFramework]
---

<!--
Source of QuantumSimilarity.nb, built with MarkdownToNotebook. Every example state outside
Possible Issues is a density matrix (positive semidefinite, unit trace) or a normalized state
vector, the inputs for which the similarity lies in [0, 1]; the one pair that is not, and the
QuantumDistance::notphysical message it issues, is under Possible Issues.
Tests/CleanupRegression.wlt replays these examples ("B1-Examples-Match" and
"B1-Examples-Match-NonPSD"), so a change to one needs the same change to the other.
-->

## Usage

<code>[QuantumSimilarity]()[$qs_{1}$,$qs_{2}$,*t*]</code> returns the similarity between two quantum discrete states using measure *t*, a number in [0,1] for density matrices and normalized state vectors.

## Details & Options

- <code>[QuantumSimilarity]()[…]</code> returns the similarity between two input states, a number in [0,1] when both are numeric and each is a density matrix or a normalized state vector. In <code>[QuantumSimilarity]()[$qs_{1}$,$qs_{2}$,*t*]</code> the following similarity measures *t* are supported:

|   |   |
|---|---|
| `"Fidelity"` | Similarity based on the fidelity between two states |
| `"Trace"` | Similarity based on the trace distance between two states |
| `"Bures"` | Similarity based on the Bures distance between two states |
| `"BuresAngle"` | Similarity based on the Bures angle between two states |
| `"HilbertSchmidt"` | Similarity based on the Hilbert-Schmidt distance |
| `"Bloch"` | Similarity based on the Euclidean distance between two Bloch vectors (single qubit only) |
| `"RelativeEntropy"` | Similarity based on the relative von Neumann entropy between two states |
| `"RelativePurity"` | $Tr[\rho_{1}.\rho_{2}]$ |

- If no measure *t* is given, `"Fidelity"` is used.

- Each similarity is constructed from the corresponding [QuantumDistance]() measure *d* by the rules below, so that identical states give similarity `1` and states with orthogonal supports give similarity `0`, with two exceptions: the `"RelativePurity"` similarity of a mixed state with itself is its purity $Tr[\rho^{2}]$, less than `1`, and the `"HilbertSchmidt"` similarity of two states with orthogonal supports is $1-\sqrt{(Tr[\rho_{1}^{2}]+Tr[\rho_{2}^{2}])/2}$, above `0` unless both states are pure. For `"Fidelity"`, `"Trace"`, `"Bloch"`, and `"RelativePurity"`, the similarity is $1-d$. For `"Bures"` and `"HilbertSchmidt"` (whose distance ranges in $[0,\sqrt{2}]$), the similarity is $1-\frac{d}{\sqrt{2}}$. For `"BuresAngle"` (whose distance ranges in $[0,\frac{\pi }{2}]$), the similarity is $1-\frac{d}{\frac{\pi }{2}}$. The `"RelativeEntropy"` case is unbounded and uses an exponential decay: similarity is $2^{-d}$. See [QuantumDistance]() for the underlying distance formulas.

- The range [0,1] holds for density matrices, which are positive semidefinite with unit trace; a normalized state vector qualifies. The inputs are used as given: a state whose density matrix does not have unit trace, such as an unnormalized state vector, is not rescaled first, except by `"Bloch"`, which reads the Bloch vector of the rescaled state; <code>*qs*["Normalized"]</code> gives the rescaled state. An input whose density matrix is numeric and has an eigenvalue below $-10^{-8}$ issues a `QuantumDistance::notphysical` message; a symbolic input is not checked.

## Basic Examples

Find the similarity between two pure states in terms of fidelity:

```wl
QuantumSimilarity[QuantumState["0"], QuantumState["1"]]
```
<!-- => 0 -->

---

Find the similarity between two mixed states in terms of fidelity:

```wl
QuantumSimilarity[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{{1/2, 0}, {0, 1/2}}]]
```
<!-- => Sqrt[3/2]/2 + 1/(2 Sqrt[2]) -->

---

Find the similarity between a pure state and a mixed state:

```wl
QuantumSimilarity[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{1, 0}]]
```
<!-- => 1/2 -->

---

Find the similarity between two quantum states using the trace measure:

```wl
QuantumSimilarity[QuantumState[{{1/4, 1/4}, {1/4, 3/4}}], 
 QuantumState["+"], "Trace"]
```
<!-- => 1 - 1/(2 Sqrt[2]) -->

## Scope

Find the similarity between multiqubit states:

```wl
QuantumSimilarity[QuantumState["GHZ"], 
 QuantumState["W"], "HilbertSchmidt"]
```
<!-- => 0 -->

---

Find the similarity between qudit states:

```wl
QuantumSimilarity[QuantumState[{1, 0, 0}, 3], 
 QuantumState[{1/Sqrt[3], 0, Sqrt[2/3]}, 3]]
```
<!-- => 1/Sqrt[3] -->

## Properties and Relations

The maximally mixed state and $|+\rangle$ have the same populations in the computational basis, $\{1/2, 1/2\}$, but different fidelity similarities to $\mathrm{diag}(1/4, 3/4)$. The maximally mixed state commutes with it, so its similarity is the overlap $\sum_i \sqrt{p_i q_i}$ of the populations $p_i$ and $q_i$ of the two states; $|+\rangle$ does not commute with it, and its similarity is smaller:

```wl
QuantumSimilarity[QuantumState[{{1/4, 0}, {0, 3/4}}], #] & /@ 
 {QuantumState[{{1/2, 0}, {0, 1/2}}], QuantumState["+"]}
```
<!-- => {Sqrt[3/2]/2 + 1/(2 Sqrt[2]), 1/Sqrt[2]} -->

## Possible Issues

These two Hermitian matrices have unit trace but each has a negative eigenvalue, so neither is a density matrix; the similarity falls outside [0,1], and a message is issued:

```wl
QuantumSimilarity[QuantumState[{{1/4, 1}, {1, 3/4}}], 
 QuantumState[{{1/2, 2}, {2, 1/2}}], "Trace"]
```
<!-- => QuantumDistance::notphysical message, then 1 - Sqrt[17]/4 -->
