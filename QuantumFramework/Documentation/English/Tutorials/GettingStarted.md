---
Template: TechNote
Name: GettingStarted
Title: Getting Started
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/GettingStarted
RelatedGuides: [WolframQuantumComputationFramework]
---

How to install and load the paclet

<!-- #| style: MathCaption -->
Install the paclet and load it:

```wl
#| eval: false
PacletInstall["Wolfram/QuantumFramework"]
Needs["Wolfram`QuantumFramework`"]
```

<!-- #| style: MathCaption -->
Check whether definitions are now available:

```wl
Names["Quantum*"]
```

<!-- #| style: MathCaption -->
A quantum gate for the magic basis transformation (transforming 2 qubit computational basis to the Bell basis):

```wl
qc = QuantumCircuitOperator["Magic"];
qc["Diagram"]
```

Test how above circuit transforms computational basis of 2-qubit into Bell states:

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

Generate corresponding tensor network of the circuit

```wl
qc["TensorNetwork", EdgeLabels -> Automatic]
```

Add measurements into above circuit

```wl
qc2 = {{1}, {2}}@*qc;
qc2["Diagram"]
```

Calculate the result of circuit on registered state:

```wl
mea = qc2[]
```

Represents the corresponding probabilities:

```wl
mea["ProbabilityPlot"]
```
