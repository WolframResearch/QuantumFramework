---
Template: Default
---

# Quantum in Finite Dimensions: Answers in WL and in QuantumFramework (Part 4)

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

## Part 4. Composite systems: tensor product and partial trace

### 4.1 [BSc] How do I form the tensor product of states and of operators?

Two subsystems combine by the tensor product: states as $|\psi\rangle\otimes|\phi\rangle$, operators as
$A\otimes B$, and the two are compatible, $(A\otimes B)(|\psi\rangle\otimes|\phi\rangle) =
(A|\psi\rangle)\otimes(B|\phi\rangle)$. In the computational basis the tensor product is the Kronecker
product, and the composite lives in the $d_1 d_2$-dimensional product space.

**WL** : the tensor product of two amplitude vectors is their (flattened) Kronecker product.

```wl
Flatten[KroneckerProduct[{a, b}, {c, d}]]
```

The tensor product of two operators is the Kronecker product of their matrices, here $X\otimes Z$.

```wl
KroneckerProduct[PauliMatrix[1], PauliMatrix[3]]
```

**QF** : `QuantumTensorProduct` combines two states into one two-qubit state object.

```wl
QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]
```

Its `"StateVector"` is that same Kronecker product.

```wl
QuantumTensorProduct[QuantumState[{a, b}], QuantumState[{c, d}]]["StateVector"] // Normal
```

It joins operators the same way; here $X\otimes Z$ kept as a `QuantumOperator`.

```wl
QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]
```

Both give the amplitude vector $\{ac,\,ad,\,bc,\,bd\}$ and the operator $X\otimes Z$: the composite state
factorizes into its parts, and the QF object additionally knows it is a two-qubit system (its factors,
order, and basis), a matrix being just one representation read from it.

### 4.2 [BSc] How do I take a partial trace of a two-party state and obtain a reduced density matrix?

The state of subsystem $A$ alone is the reduced density matrix $\rho_A = \mathrm{Tr}_B[\rho_{AB}]$, formed
by tracing out $B$. For a product state $\rho_A$ stays pure; for an *entangled* state it comes out mixed,
which is the operational signature of entanglement.

**WL** : trace qubit $B$ out of the Bell state $\rho = |\Phi^+\rangle\langle\Phi^+|$ by reshaping the
$4\times4$ matrix to a $2\times2\times2\times2$ tensor and contracting the two $B$ indices.

```wl
With[{bell = {1, 0, 0, 1}/Sqrt[2]},
 TensorContract[ArrayReshape[KroneckerProduct[bell, Conjugate[bell]], {2, 2, 2, 2}], {{2, 4}}]]
```

**QF** : `QuantumPartialTrace` traces out the listed qubit and returns the reduced state as an object.

```wl
QuantumPartialTrace[QuantumState["PhiPlus"], {2}]
```

Its `"Purity"` is below one, so the reduction of the entangled pair is mixed.

```wl
QuantumPartialTrace[QuantumState["PhiPlus"], {2}]["Purity"]
```

A *product* input, by contrast, reduces to a pure state (purity $1$).

```wl
QuantumPartialTrace[QuantumTensorProduct[QuantumState["0"], QuantumState["1"]], {2}]["Purity"]
```

Both give $\rho_A = I/2$ with purity $\tfrac12$: tracing a maximally entangled pair down to one qubit
yields the maximally mixed state, whereas a product state stays pure. The mixedness of the reduced state
is exactly the entanglement discarded along with $B$.
