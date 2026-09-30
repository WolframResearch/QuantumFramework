---
Template: Symbol
Name: QuantumDistance
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumDistance
Keywords: [quantum distance, fidelity, trace distance, Bures distance, Bures angle, Hilbert-Schmidt distance, Bloch distance, relative entropy, distance measure]
SeeAlso: [QuantumState, QuantumSimilarity]
RelatedTutorials: [[Wolfram Quantum Framework Tutorial](Tutorial)]
RelatedGuides: [WolframQuantumComputationFramework]
---

<!--
Source of QuantumDistance.nb, built with MarkdownToNotebook. Every example state outside
Possible Issues is a density matrix (positive semidefinite, unit trace) or a normalized
state vector, and none issues a message; the examples under Possible Issues show the
QuantumDistance::notphysical and QuantumDistance::notnormalized messages and the Failure
that "Bloch" returns for a qutrit. The first four Basic Examples, the first two Scope
examples, and the first two Possible Issues examples use the same inputs as the
QuantumSimilarity page, so a change to them belongs on both pages.
-->

## Usage

<code>[QuantumDistance]()[$qs_{1}$,$qs_{2}$,*t*]</code> returns the distance between two quantum discrete states using measure *t*.

## Details & Options

- <code>[QuantumDistance]()[…]</code> returns the distance between two input states, a number for numeric input, a quantity in bits for `"RelativeEntropy"`, and a <code>[Failure]()</code> object for `"Bloch"` beyond a single qubit. In <code>[QuantumDistance]()[$qs_{1}$,$qs_{2}$,*t*]</code> the following distance measures *t* are supported:

|   |   |
|---|---|
| `"Fidelity"` | Distance between two states in terms of fidelity |
| `"Trace"` | Trace distance between two states |
| `"Bures"` | Bures distance between two states |
| `"BuresAngle"` | Bures angle between two states |
| `"HilbertSchmidt"` | Hilbert-Schmidt distance |
| `"Bloch"` | Half the Euclidean distance between two Bloch vectors (single qubit only) |
| `"RelativeEntropy"` | Relative von Neumann entropy between two states |
| `"RelativePurity"` | $1-Tr[\rho_{1}.\rho_{2}]$ |

- If no measure *t* is given, `"Fidelity"` is used.

- The distance based on fidelity is calculated as $1-Tr[\sqrt{\rho_{1}^{1/2}.\rho_{2}.\rho_{1}^{1/2}}]$ with $\rho_{1,2}$ the corresponding density matrices and $Tr[\sqrt{\rho_{1}^{1/2}.\rho_{2}.\rho_{1}^{1/2}}]$ the fidelity. The trace distance is $\frac{1}{2}Tr[\sqrt{(\rho_{1}-\rho_{2})^{2}}]$. The Bures distance is $\sqrt{2\, (1\, -\, \mathit{fidelity})}$, and Bures angle is <code>[ArcCos]()[*fidelity*]</code>. The Hilbert-Schmidt distance is the [Hilbert-Schmidt norm](https://mathworld.wolfram.com/Hilbert-SchmidtNorm.html) of the difference between density matrices $\sqrt{Tr[(\rho_{1}-\rho_{2})^{2}]}$. The Bloch distance is half the [EuclideanDistance]() between two Bloch vectors. And the relative entropy is <code>[Tr]()[$\rho _{1}$.[MatrixLog]()[$\rho _{1}$]-$\rho _{1}$.[MatrixLog]()[$\rho _{2}$]]/[Log]()[2]</code>, returned as a quantity in bits.

- The measures are defined for density matrices, which are Hermitian and positive semidefinite with unit trace. An input whose trace is a positive number other than 1, such as a state vector whose squared norm is not 1, is divided by its trace first. An exact positive trace other than 1 issues a `QuantumDistance::notnormalized` message; a machine-precision trace is divided in any case, and issues the message when it is more than $10^{-8}$ from 1. A symbolic input is rescaled only when its trace evaluates to a number, and is otherwise used as given; <code>*qs*["Normalized"]</code> gives the state it stands for. An input whose trace is 0, negative, or not real, or a numeric input that, once divided by its trace, is not Hermitian (tested exactly when its entries, once divided, are exact rational or complex rational numbers, and otherwise by [HermitianMatrixQ]() with a tolerance of $10^{-8}$) or has a machine-precision eigenvalue below $-10^{-8}$, issues a `QuantumDistance::notphysical` message. The entries of an exact input multiplied by an irrational number, as $\rho/\sqrt{2}$ for a matrix $\rho$ of rational entries, need not be rational once divided by the trace, and such an input is then tested with the tolerance. A machine-precision trace counts as 0 when its absolute value is at most $10^{-8}$ times the Frobenius norm of the numeric entries of the matrix, as positive when its real part exceeds, and its imaginary part is at most, $10^{-8}$ times the larger of that norm and the absolute value of the trace, and otherwise as negative or not real. A machine-precision state vector whose squared norm underflows to 0 counts as the zero vector. An exact trace whose sign cannot be decided, numerically or, for an algebraic trace, with [RootReduce](), is used as given, and a matrix with a symbolic entry is checked only through its trace.

- The `"Bloch"` distance is defined for a single qubit, whose Bloch vectors give the trace distance; for any other dimension it returns a <code>[Failure]()</code> object.

## Basic Examples

Find the distance between two pure states in terms of fidelity:

```wl
QuantumDistance[QuantumState["0"], QuantumState["1"]]
```
<!-- => 1 -->

---

Find the distance between two mixed states in terms of fidelity:

```wl
QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{{1/2, 0}, {0, 1/2}}]]
```
<!-- => 1 - Sqrt[3/2]/2 - 1/(2 Sqrt[2]) -->

---

Find the distance between a pure state and a mixed state:

```wl
QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{1, 0}]]
```
<!-- => 1/2 -->

---

Find the distance between two quantum states using the trace metric:

```wl
QuantumDistance[QuantumState[{{1/4, 1/4}, {1/4, 3/4}}], 
 QuantumState["+"], "Trace"]
```
<!-- => 1/(2 Sqrt[2]) -->

---

Find the distance between a pure and mixed state in terms of the Bures angle:

```wl
QuantumDistance[QuantumState[{{1/3, 0}, {0, 2/3}}], 
 QuantumState[{1, 0}], "BuresAngle"]
```
<!-- => ArcCos[1/Sqrt[3]] -->

## Scope

Find the distance between multiqubit states:

```wl
QuantumDistance[QuantumState["GHZ"], 
 QuantumState["W"], "HilbertSchmidt"]
```
<!-- => Sqrt[2] -->

---

Find the distance between qudit states:

```wl
QuantumDistance[QuantumState[{1, 0, 0}, 3], 
 QuantumState[{1/Sqrt[3], 0, Sqrt[2/3]}, 3]]
```
<!-- => 1 - 1/Sqrt[3] -->

---

Find the distance between a symbolic qubit state and the `"+"` state in terms of fidelity:

```wl
FullSimplify[#, {\[Theta], \[Phi]} \[Element] Reals] &[
 QuantumDistance[
  QuantumState[{Cos[\[Theta]/2], E^(I \[Phi]) Sin[\[Theta]/2]}], 
  QuantumState["+"]]]
```
<!-- => 1 - Sqrt[2 + 2 Cos[\[Phi]] Sin[\[Theta]]]/2 -->

## Applications

Change of the trace distance under quantum channels

Initial maximum distance for orthogonal states:

```wl
QuantumDistance[QuantumState["+"], QuantumState["-"], "Trace"]
```
<!-- => 1 -->

Apply the "AmplitudeDamping" channel to both states:

```wl
{qs1, qs2} = 
  QuantumChannel["AmplitudeDamping"[\[Gamma]]] /@ {QuantumState["+"], 
    QuantumState["-"]};
```

Compute the trace distance as function of the channel parameter:

```wl
Simplify[QuantumDistance[qs1, qs2, "Trace"], 0 < \[Gamma] < 1]
```
<!-- => Sqrt[1 - \[Gamma]] -->

## Properties and Relations

The relative entropy of a mixed state with respect to a pure state is not finite:

```wl
QuantumDistance[QuantumState["RandomMixed"], 
 QuantumState["0"], "RelativeEntropy"]
```
<!-- => Quantity[\[Infinity], "Bits"] -->

---

Relative entropy is asymmetric, reversing the order of the quantum states changes the value:

```wl
QuantumDistance[##, "RelativeEntropy"] & @@@ 
 Permutations@Table[QuantumState["RandomMixed"], 2]
```
<!-- => a pair of different Quantity[_, "Bits"] values (random mixed states) -->

## Possible Issues

These two Hermitian matrices have unit trace but each has a negative eigenvalue, so neither is a density matrix; a message is issued, and their trace distance exceeds 1, the largest trace distance two density matrices can have:

```wl
QuantumDistance[QuantumState[{{1/4, 1}, {1, 3/4}}], 
 QuantumState[{{1/2, 2}, {2, 1/2}}], "Trace"]
```
<!-- => QuantumDistance::notphysical message, then Sqrt[17]/4 -->

---

The identity matrix has trace 2, so it is divided by its trace and read as the maximally mixed state, and a message is issued:

```wl
QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{{1, 0}, {0, 1}}]]
```
<!-- => QuantumDistance::notnormalized message, then 1 - Sqrt[3/2]/2 - 1/(2 Sqrt[2]) -->

---

The `"Bloch"` distance is defined only for a single qubit, so for a qutrit it returns a failure:

```wl
QuantumDistance[QuantumState[{1, 0, 0}, 3], 
 QuantumState[{0, 1, 0}, 3], "Bloch"]
```
<!-- => Failure["NonQubitBloch", <|"MessageTemplate" :> QuantumDistance::bloch, "MessageParameters" -> {3}, "Dimension" -> 3|>] -->
