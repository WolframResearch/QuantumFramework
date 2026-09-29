---
Template: TechNote
Name: Tutorial
Title: Quantum Computation
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/Tutorial
RelatedGuides: [WolframQuantumComputationFramework]
---

Quantum computation is the use of quantum mechanical systems to perform computations. Wolfram quantum framework aims to simulate a wide range of quantum computations in the Wolfram Mathematica.

## Basic quantum objects

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [QuantumBasis]() | basis vectors encoding quantum states and operators |
| [QuantumState]() | quantum discrete state, defined by an association, a complex vector, or a density matrix |
| [QuantumOperator]() | quantum discrete operator, defined by an association or a matrix representation |
| [QuantumMeasurementOperator]() | quantum measurement operator, defined by a matrix, a collection of positive matrices (POVM) or an eigenbasis (projective) |
| [QuantumMeasurement]() | description of possible measurement results, containing a probability distribution as well as possible quantum states after the measurement |

Basic objects of discrete quantum mechanics.

A fundamental object in our framework is [QuantumBasis](). Quantum states and operators are defined with respect to a basis. [QuantumBasis]() can be constructed for any number of qudits, with any dimensionality. There are also a set of named basis (e.g., PauliX or Bell) built-in into the quantum framework.

<!-- #| style: MathCaption -->
Define a quantum basis given dimensions (3x5):

```wl
QuantumBasis[{3, 5}]
```

<!-- #| style: MathCaption -->
Define a quantum basis by an association:

```wl
QuantumBasis[<|0 -> {1, 1}, 1 -> {1, -I}|>]
```

<!-- #| style: MathCaption -->
There are many named-basis built into the framework:

```wl
{QuantumBasis["I"], QuantumBasis["Y"], QuantumBasis["Bell"], 
 QuantumBasis["Dirac"], QuantumBasis["Schwinger"], 
 QuantumBasis["Wigner"]}
```

<!-- #| style: MathCaption -->
There are many properties that can be extracted from [QuantumBasis]() object. For example:

```wl
basis = QuantumBasis[3];
Normal /@ basis["ElementAssociation"]
```

<!-- #| style: MathCaption -->
[QuantumBasis]() supports cases with input and output:

```wl
QuantumBasis[{2}, {3}]
```

```wl
{#["Input"], #["Output"]} &@QuantumBasis[{2}, {3}]
```

A quantum state is represented by [QuantumState]() object and a quantum operator is represented by [QuantumOperator]().

<!-- #| style: MathCaption -->
Define a pure 2-dimensional quantum state (qubit) in Pauli-X basis:

```wl
QuantumState[{1, -I}, "X"]
```

```wl
%["Amplitudes"]
```

<!-- #| style: MathCaption -->
If the basis is not specified, the default is the computational basis:

```wl
state = QuantumState["RandomPure"[3]];
state["Formula"]
```

<!-- #| style: MathCaption -->
Many named states are available for easy access:

```wl
{QuantumState["UniformSuperposition"], QuantumState["PsiPlus"], 
 QuantumState["GHZ"]}
```

<!-- #| style: MathCaption -->
We can also define a qubit state by specifying a Bloch vector:

```wl
QuantumState["BlochVector"[{.1, .2, .3}]]
```

<!-- #| style: MathCaption -->
Define a quantum operator by a matrix, basis, and order. The order is the information about which subsystems the operator would act on. For example order {1,2} means it would act on subsystem 1 and 2. If the basis is not specified, the default is the computational basis.

```wl
QuantumOperator["CNOT", {2, 3}, "XI"]
```

<!-- #| style: MathCaption -->
Define operator by names:

```wl
{QuantumOperator["H", {2}], QuantumOperator["CNOT"], 
 QuantumOperator["Toffoli"]}
```

<!-- #| style: MathCaption -->
Extract properties of operator:

```wl
operator = QuantumOperator["CNOT", {1, 2}];
operator["Table"]
```

```wl
operator["HermitianQ"]
```

<!-- #| style: MathCaption -->
[QuantumOperator]() can operate on a [QuantumState]():

```wl
QuantumOperator["H"][QuantumState["0"]]
```

In Wolfram quantum framework, a measurement is represented by [QuantumMeasurementOperator]() and [QuantumMeasurement]().

<!-- #| style: MathCaption -->
Define a measurement operator on the 2nd qubit (in computational basis):

```wl
QuantumMeasurementOperator[{2}]
```

<!-- #| style: MathCaption -->
Define measurement operator by names:

```wl
QuantumMeasurementOperator["X"]
QuantumMeasurementOperator["QBismSICPOVM"]
```

<!-- #| style: MathCaption -->
Define measurement operator by a [QuantumBasis]() object (as an eigenbasis) and a list of eigenvalues:

```wl
qmo = QuantumMeasurementOperator["Bell" -> {1, 3, 4, 5}]
```

```wl
qmo[QuantumState[{1, 3, 0, I}, "Bell"]]["ProbabilityPlot", 
 AspectRatio -> 1/4]
```

<!-- #| style: MathCaption -->
Define a POVM measurement:

```wl
povm = QuantumMeasurementOperator["TetrahedronSICPOVM"]
```

```wl
povm[QuantumState["RandomPure"]]["Probabilities"]
```

<!-- #| style: MathCaption -->
When [QuantumMeasurementOperator]() acts on [QuantumState](), [QuantumMeasurement]() is generated:

```wl
result = QuantumMeasurementOperator["X"][QuantumState[{0.2, 0.6}]]
```

```wl
result["ProbabilityPlot", AspectRatio -> 1/4]
```

<!-- #| style: MathCaption -->
There are two ways we can perform multiqubit measurements. The first method is sequential measurement:

```wl
result1 = QuantumMeasurementOperator[][QuantumState["RandomPure"[2]]]
```

```wl
result2 = QuantumMeasurementOperator[{2}][result1]
```

```wl
result2["ProbabilityPlot", AspectRatio -> 1/4]
```

To define the state of non-interacting subsystem or to define local operations on composite system, [QuantumTensorProduct]() is used. It can performs tensor product on bases, states, or operators.

<!-- #| style: MathCaption -->
Tensor product of states and bases:

```wl
QuantumTensorProduct[{QuantumState["0"], QuantumState["1"]}]
```

```wl
QuantumTensorProduct[{QuantumBasis["X"], QuantumBasis["Y"]}]
```

<!-- #| style: MathCaption -->
Tensor product of operators ([QuantumOperator]() or [QuantumMeasurementOperator]()):

```wl
QuantumTensorProduct[QuantumOperator["Z"], QuantumOperator["Z"]]
```

```wl
QuantumTensorProduct[
 QuantumMeasurementOperator["X"], QuantumMeasurementOperator["Y"]]
```

The representation of quantum state and operator can be changed according to any basis transformation. Basis transformation is applicable to [QuantumState](), [QuantumOperator](), and [QuantumMeasurementOperator]() as demonstrated below:

<!-- #| style: MathCaption -->
Basis change for a quantum state:

```wl
QuantumState[QuantumState["1"], "Y"]
```

<!-- #| style: MathCaption -->
Basis change for a quantum operator:

```wl
QuantumOperator[QuantumOperator["X"], "Y"]
```

<!-- #| style: MathCaption -->
One can do more operations on quantum operators and the result will be a QuantumOperator:

```wl
op = Exp[-I \[Phi]/2 QuantumOperator["X"]]
```

```wl
op == QuantumOperator["RX"[\[Phi]]] // FullSimplify
```

## Advanced quantum objects

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [QuantumCircuitOperator]() | quantum circuit as a list of quantum operations |
| [QuantumPartialTrace]() | reduced state or reduced basis by tracing out subsystems |
| [QuantumPartialTranspose]() | partial transposition of a quantum state or operator |
| [QuantumDistance]() | distance between two quantum states |
| [QuantumEntangledQ]() | check whether a quantum state is entangled |
| [QuantumEntanglementMonotone]() | entanglement monotone that measures the entanglement of a quantum state |
| [QuantumShortcut]() | label list of a given quantum circuit |

Basic objects and operations in discrete quantum mechanics naturally lead to an application in quantum information and computation. Quantum Information is the study of information encoded in quantum systems. Meanwhile, quantum computation is the manipulation of quantum information to perform a computation.

In quantum information, it is often convenient to have a local state description of a multiqubit state. For instance, Alice and Bob shared a (possibly entangled) quantum state. But, since they are really far away from each other, they do not have access to each other state. In such scenario, Alice's description of her part of quantum state is described by a "reduced state" that can be obtained by [QuantumPartialTrace]() operation.

<!-- #| style: MathCaption -->
Trace out the second subsystem in a two-qubit state:

```wl
QuantumPartialTrace[QuantumState["PsiPlus"], {2}]
```

<!-- #| style: MathCaption -->
Partial trace can also be applied to [QuantumBasis]():

```wl
QuantumPartialTrace[QuantumBasis[2, 2], {2}]
```

Partial transpose is an operation defined for multiple subsystems state in which only a specific subsystem is transposed. For instance, if there are two subsystems, partial transpose with respect to the first subsystem is defined by a map $\rho_{1}\otimes \rho_{2}\to \, \rho_{1}^{T}\otimes \rho_{2}$. Partial transpose operation is represented by [QuantumPartialTranspose]().

<!-- #| style: MathCaption -->
Example of partial transpose:

```wl
QuantumPartialTranspose[QuantumState["PsiPlus"], {2}]
```

Although partial transpose is not a physical operation (it does not map a density matrix to another density matrix), this operation is quite useful for detecting entanglement. In particular, it is used to compute logarithmic negativity which is commonly used metric to measure entanglement between two states. There are several other metrics to measure entanglement such as concurrence and entanglement entropy. The calculation of entanglement measure is represented by [QuantumEntanglementMonotone]() function.

<!-- #| style: MathCaption -->
Plotting various entanglement measure for the state $\alpha |00\rangle +\, \sqrt{1-\alpha^{2}}|11\rangle $

```wl
Plot[{
  QuantumEntanglementMonotone[
   QuantumState[{\[Alpha], 0, 0, Sqrt[1 - \[Alpha]^2]}], "Concurrence"],
  	QuantumEntanglementMonotone[
   QuantumState[{\[Alpha], 0, 0, Sqrt[1 - \[Alpha]^2]}], 
   "LogNegativity"],
  	QuantityMagnitude@
   QuantumEntanglementMonotone[
    QuantumState[{\[Alpha], 0, 0, Sqrt[1 - \[Alpha]^2]}], 
    "EntanglementEntropy"]}, {\[Alpha], 0, 1},
 PlotLabel -> 
  "Entaglement monotone of \[Alpha]\!\(\*TemplateBox[{\"00\"},\n\"Ket\
\"]\)+ \!\(\*SqrtBox[\(1 - \*SuperscriptBox[\(\[Alpha]\), \(2\)]\)]\)\
\!\(\*TemplateBox[{\"11\"},\n\"Ket\"]\)",
 AxesLabel -> {"\[Alpha]", None}, 
 PlotLegends -> {"Concurrence", "LogNegativity", 
   "Entanglement Entropy"}]
```

To know whether a state is entangled or separable without computing its measure, [QuantumEntangledQ]() is used.

<!-- #| style: MathCaption -->
Checking whether a subsystem 1 and 3 is entangled in "W" state:

```wl
QuantumEntangledQ[QuantumState["W"], {{1}, {3}}]
```

In quantum information, there exist notions of distance between quantum state such as fidelity, trace distance, Bures angle, et cetera. One may compute distance between two quantum state using various metrics by [QuantumDistance]():

<!-- #| style: MathCaption -->
Measuring trace distance between a pure state and a mixed state:

```wl
QuantumDistance[QuantumState[{{1/4, 0}, {0, 3/4}}], 
 QuantumState[{1, 0}], "BuresAngle"]
```

Wolfram quantum framework also supports Schmidt decomposition and spectral decomposition.

<!-- #| style: MathCaption -->
Example of Schmidt decomposition:

```wl
state = QuantumState["PhiPlus"]["SchmidtBasis"];
state["Basis"]["Association"] // Map@Normal
state["Amplitudes"]
```

<!-- #| style: MathCaption -->
Example of spectral decomposition:

```wl
diagState = QuantumState[{{1/2, 1/2}, {1/2, 1/2}}]["SpectralBasis"];
diagState["Amplitudes"]
```

Lastly, one may create a list of many quantum operators. A quantum circuit is represented by the symbol [QuantumCircuitOperator]().

<!-- #| style: MathCaption -->
Example for the construction of quantum circuit without measurement:

```wl
qc = QuantumCircuitOperator[{"X", 1, 
    "CNOT" -> {3, 2}, {"R", \[Theta], "YY" -> {2, 3}}, "SWAP", 
    "SX" -> 3, "P"[\[Phi]], "T" -> 2, {"C", "NOT" -> 3, {1, 2}}, 
    "H" -> 2, "BitFlip"[p] -> {1}, "Braid" -> {1, 3}}, 
   "Parameters" -> {\[Theta], \[Phi], p}];
qc["Diagram"]
```

<!-- #| style: MathCaption -->
When there is no measurement in the circuit, its action on [QuantumState]() would return another [QuantumState]():

```wl
qc[\[Pi]/3, \[Pi], \[Pi]/5][]
```

<!-- #| style: MathCaption -->
Add measurement into the circuit:

```wl
qc2 = qc/*{{1, 3}};
qc2["Diagram"]
```

<!-- #| style: MathCaption -->
When there is measurement in the circuit, its action on [QuantumState]() would return [QuantumMeasurement]():

```wl
N[qc2[\[Pi]/3, \[Pi], \[Pi]/5]][]
```

<!-- #| style: MathCaption -->
Decomposition of CNOT:

```wl
qc = QuantumCircuitOperator[{"RY"[-\[Pi]/2] -> 2, "CZ", 
    "RY"[\[Pi]/2] -> 2}];
qc["Diagram"]
```

```wl
QuantumOperator[qc] == QuantumOperator["CNOT"]
```

```wl
qc = QuantumCircuitOperator[{"RY"[\[Pi]/2] -> 2, "RootSWAP", 
    "RZ"[\[Pi]], "RootSWAP", "RZ"[-(\[Pi]/2)] -> {1, 2}, 
    "RY"[-(\[Pi]/2)] -> 2}];
qc["Diagram"]
```

```wl
QuantumOperator[qc] == QuantumOperator["CNOT"]
```

<!-- #| style: MathCaption -->
A quantum gate for the magic basis transformation (transforming 2 qubit computational basis to the Bell basis):

```wl
qc = QuantumCircuitOperator["Magic"];
qc["Diagram"]
```

```wl
qc[QuantumState["00"]] == QuantumState["PhiPlus"]
```

```wl
qc[QuantumState["10"]] == I QuantumState["PsiPlus"]
```

```wl
qc[QuantumState["01"]] == I QuantumState["PhiMinus"]
```

```wl
qc[QuantumState["11"]] == QuantumState["PsiMinus"]
```
