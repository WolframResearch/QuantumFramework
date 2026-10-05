---
Template: TechNote
Name: ExploringFundamentalsOfQuantumTheory
Title: Exploring Fundamentals of Quantum Theory
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/ExploringFundamentalsOfQuantumTheory
RelatedGuides: [WolframQuantumComputationFramework]
Typeset: _SuperDagger -> StandardForm
---

We will explore some important experiment in the foundation of quantum theory, using the Wolfram Quantum Framework. In this document, we will discuss quantum eraser, Elitzur-Vaidman bomb experiment, Hardy’s paradox, and quantum SWITCH.

## Quantum Eraser

If there is path information (stored somewhere, regardless of whether it is measured or not), there won’t be any interference. The mere presence of path information kills quantum interference. But by erasing the path information, one can recover quantum interference. In the following setup, clicking all detectors corresponds to a no-interference case, and clicking of only D1 and D4 corresponds to the interference case.

### With path information: no interference

Prepare a series of quantum gates representing the above setup:

```wl
ops = {"H", "CX", "H", {1}, {2}};
```

Construct the quantum circuit:

```wl
qc = QuantumCircuitOperator[ops];
qc["Diagram"]
```

Calculate the system evolution through each step:

```wl
steps = ComposeList[qc["Operators"], QuantumState["Register"[2]]];
```

Show the states prior to measurement:

```wl
Grid[Transpose[{Style[#, Bold] & /@ {"Initial state", 
     "After 1st Hadamard", "After CNOT", 
     "After 2nd Hadamard"}, #["Formula"] & /@ steps[[;; -3]]}], 
 Frame -> All, Alignment -> Left]
```

Get the probability of each outcome:

```wl
qc[QuantumState["Register"[2]]]["ProbabilityPlot"]
```

### Erasing path information: interference

Prepare a series of quantum gates representing above setup:

```wl
ops = {"H", "CX", "H", "H" -> 2, {1}, {2}};
```

Construct the quantum circuit:

```wl
qc = QuantumCircuitOperator[ops];
qc["Diagram"]
```

Note that the Hadamard gate acting on 2nd qubit serves as quantum eraser.

Calculate the system evolution through each step:

```wl
steps = ComposeList[qc["Operators"], QuantumState["Register"[2]]];
```

Show the states prior to measurement:

```wl
Grid[Transpose[{Style[#, Bold] & /@ {"Initial state", 
     "After Hadamard on 1st qubit", "After CNOT", 
     "After Hadamard on 1st qubit", 
     "After Hadamard on 2nd qubit"}, #["Formula"] & /@ steps[[;; -3]]}],
  Frame -> All, Alignment -> Left]
```

Get the probability of each outcome:

```wl
qc[QuantumState["Register"[2]]]["ProbabilityPlot"]
```

## Elitzur-Vaidman bomb

Given quantum interference in a Mach-Zehnder setup with only one beam-splitter, a detector click can be inferred as the presence of an object along one arm of interferometer. However, due to the interference, there is a 50% chance that the photon does not move through that arm, as if one can infer the presence of that object without any interaction. To dramatize the case, think about that object as a bomb! This thought experiment is also called the interaction-free measurement.

### Simple case: 50:50 beam splitter

Prepare a series of quantum gates representing the above setup:

```wl
ops = {"H", "CX", "H", {1}, {2}};
```

Construct the quantum circuit:

```wl
qc = QuantumCircuitOperator[ops];
qc["Diagram"]
```

Get the probability of each outcome:

```wl
prob = qc[]["Probabilities"]
```

Note that detecting $|01\rangle $ or $|11\rangle $ means that there is interaction with bomb (or a photon passed through the arm where the bomb is), while $|00\rangle $ means that, with the help of quantum interference, the bomb is detected with no interaction. Therefore, the efficiency rate $\eta$ can be expressed as $\eta = P_{00}/(P_{00}+P_{01}+P_{11})=P_{00}/(1-P_{10})$.

Calculate the probability:

```wl
prob[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]]/(1 - prob[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]])
```

### Beyond 50:50 beam splitter

The Hadamard gate acts like a 50:50 beam splitter. One can replace it with a $R_{y}(\theta)$ gate, which is rotates along the $y$-axis with an arbitrary angle θ:

```wl
QuantumOperator["RY"[\[Theta]]]
```

Compare the Hadamard gate with the $R_{y}(\theta)$ gate on the quantum state $|0\rangle $:

```wl
QuantumOperator["H"][QuantumState[{1, 0}]]["Amplitudes"]
```

```wl
FullSimplify /@ 
 QuantumOperator["RY"[\[Theta]]][QuantumState[{1, 0}]]["Amplitudes"]
```

Note that the reflection will turn the qubit to $|1\rangle $, meaning it will be in the arm where the bomb is located.

Define a new quantum circuit using $R_{y}(\theta)\, $:

```wl
qc = QuantumCircuitOperator[{"RY"[\[Theta]], 
    "CX", "RY"[\[Theta]], {1}, {2}}];
qc["Diagram"]
```

Get the probability of each outcome:

```wl
prob = FullSimplify[#, \[Theta] \[Element] Reals] & /@ 
  qc[QuantumState["Register"[2]]]["Probabilities"]
```

Return the efficiency rate:

```wl
prob[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]]/(1 - prob[
Wolfram`QuantumFramework`QuditName[{1, 0}, "Dual" -> False]])
```

Plot the efficiency rate:

```wl
Plot[prob[
Wolfram`QuantumFramework`QuditName[{0, 0}, "Dual" -> False]]/(1 - prob[
Wolfram`QuantumFramework`QuditName[{1, 0}, 
      "Dual" -> False]]), {\[Theta], 0, 2 \[Pi]}, 
 Ticks -> {{0, \[Pi]/2, \[Pi], (3 \[Pi])/2, 2 \[Pi]}, Automatic}, 
 GridLines -> {{0, \[Pi]/2, \[Pi], (3 \[Pi])/2, 2 \[Pi]}, Automatic}]
```

### Increasing efficiency by using a series of beam splitters

Define a new quantum circuit using $R_{y}(\theta)\, $:

```wl
qcN[\[Theta]_, m_Integer] := 
  With[{n = m - 1}, 
   QuantumCircuitOperator@
    Append[Table[{"RY"[\[Theta]], "CX" -> {1, i + 1}}, {i, 
       n}], "RY"[\[Theta]]]
   ];
```

```wl
qcN[\[Theta], 5]["Diagram"]
```

Return the efficiency rate $\eta = P_{00,...,0}/(1-P_{10,...,0})$:

```wl
\[Eta][\[Theta]_, n_] := 
 Module[{\[Psi]r, \[Psi]t, \[Psi]f, pDet, pNul},
  (*|000,...,00\[RightAngleBracket]*)
  \[Psi]r = QuantumState["Register"[n]];
  (*|100,...,00\[RightAngleBracket]*)
  \[Psi]t = QuantumOperator["X"]@QuantumState["Register"[n]];
  (*final state at the end of the cicruit*)
  \[Psi]f = qcN[\[Theta], n][\[Psi]r];
  (*inner products wrt final state*)
  {pDet, pNul} = 
   Abs[First@SuperDagger[#][\[Psi]f][
         "StateVector"]]^2 & /@ {\[Psi]r, \[Psi]t};
  (*\[Eta]*)
  pDet/(1 - pNul)
  ]
```

Plot the efficiency rate:

```wl
ListPlot[Table[{n, N@\[Eta][\[Pi]/n, n]}, {n, 2, 12}], 
 GridLines -> All, AxesLabel -> {"#BS", "\[Eta]"}, 
 PlotLabel -> "\[Eta] vs #Beam-Splitters, with \[Theta]=\[Pi]/#BS"]
```

Note that the dimension of matrices is $2^{\# \mathrm{BS}}$, so increasing the number of qubits can lengthen the computation time.

## Hardy’s paradox

Define a new quantum circuit using $R_{y}(\theta)$:

```wl
qcH[\[Theta]1_, \[Theta]2_] := 
  QuantumCircuitOperator[{QuantumOperator["YRotation"[\[Theta]1]], 
    QuantumOperator["YRotation"[\[Theta]2], {2}], 
    QuantumOperator["Toffoli"], 
    QuantumOperator["YRotation"[\[Pi] - \[Theta]1]], 
    QuantumOperator["YRotation"[\[Pi] - \[Theta]2], {2}]}];
```

```wl
qcH[\[Theta]1, \[Theta]2]["Diagram"]
```

Get the final state:

```wl
final = qcH[\[Theta]1, \[Theta]2][QuantumState["Register"[3]]];
```

Calculate the nonlocal probability for $|000\rangle $ (where one considers only cases where 3rd qubit is nonzero) as $P_{000}/(P_{000}+P_{100}+P_{010}+P_{110})$:

```wl
nonlocalProb = With[{prob = final["Amplitudes"]^2}, prob[
Wolfram`QuantumFramework`QuditName[{0, 0, 0}, 
      "Dual" -> False]]/(prob[
Wolfram`QuantumFramework`QuditName[{0, 0, 0}, "Dual" -> False]] + prob[
Wolfram`QuantumFramework`QuditName[{1, 0, 0}, "Dual" -> False]] + prob[
Wolfram`QuantumFramework`QuditName[{0, 1, 0}, "Dual" -> False]] + prob[
Wolfram`QuantumFramework`QuditName[{1, 1, 0}, "Dual" -> False]])] // 
  FullSimplify[#, \[Theta]1 \[Element] Reals && \[Theta]2 \[Element] 
      Reals] &
```

Plot the nonlocal probability in 3D:

```wl
Plot3D[100. nonlocalProb, {\[Theta]1, 0, 2 \[Pi]}, {\[Theta]2, 0, 
  2 \[Pi]}, 
 AxesLabel -> {"\[Theta]1", "\[Theta]2", "% nonlocal probability"}, 
 Ticks -> {{0, \[Pi]/2, \[Pi], 3 \[Pi]/2, 
    2 \[Pi]}, {0, \[Pi]/2, \[Pi], 3 \[Pi]/2, 2 \[Pi]}, Automatic}]
```

Show the nonlocal probability as a density plot:

```wl
DensityPlot[
 100. nonlocalProb, {\[Theta]1, 0, 2 \[Pi]}, {\[Theta]2, 0, 2 \[Pi]}, 
 AxesLabel -> {"\[Theta]1", "\[Theta]2", "% nonlocal probability"}, 
 FrameTicks -> {{0, \[Pi]/2, \[Pi], 3 \[Pi]/2, 
    2 \[Pi]}, {0, \[Pi]/2, \[Pi], 3 \[Pi]/2, 2 \[Pi]}}, 
 PlotLegends -> Automatic, PlotLabel -> {"% nonlocal probability"}]
```

## Quantum Switch

The following example is based on a paper published in [Nature Comm](https://www.nature.com/articles/ncomms8913.pdf).

### Commuting example

Create a random 2D unitary operator:

```wl
u = QuantumOperator["RandomUnitary"];
```

Create two random commuting operators:

```wl
A = QuantumOperator[
   u@QuantumOperator["Phase"[RandomReal[{0, 2 \[Pi]}]]]@SuperDagger[u],
    "Label" -> "A"];
B = QuantumOperator[
   u@QuantumOperator["Phase"[RandomReal[{0, 2 \[Pi]}]]]@SuperDagger[u],
    "Label" -> "B"];
```

Consider the case where [A,B]=0:

```wl
(A@B - B@A)["Matrix"] // Chop // MatrixForm
```

The action of following circuit on an initial state $|+\rangle $$|\psi\rangle $ is $\frac{1}{2}$$|0\rangle $$\{\mathrm{A}, \mathrm{B} \}$$|\psi\rangle $+$|1\rangle $$[\mathrm{A}, \mathrm{B}]$$|\psi\rangle $, which means if the qubit 1 is observed in the state $|0\rangle $ then it confirms that $[\mathrm{A}, \mathrm{B}]$$=0$. If the qubit-1 is observed in the state $|1\rangle $ then it confirms that $\{\mathrm{A}, \mathrm{B} \}$$=0$.

Construct the circuit based on a quantum switch (defined in the initialization):

```wl
qc = QuantumCircuitOperator[{QuantumOperator[{"Switch", A, B}], 
    QuantumOperator["H"], QuantumMeasurementOperator[]}];
qc["Diagram"]
```

Get the tensor product for $|+\rangle |RandomPure\rangle $:

```wl
\[Psi]0 = 
 QuantumTensorProduct[QuantumState["Plus"], QuantumState["RandomPure"]]
```

Perform the measurement:

```wl
qc[\[Psi]0]["ProbabilityPlot"]
```

The result confirms the commuting case.

### Anti-commuting example

Create a random 2D unitary operator:

```wl
u = QuantumOperator["RandomUnitary"];
```

Create two commuting operators:

```wl
A = QuantumOperator[u@QuantumOperator["X"]@SuperDagger[u], "Label" -> "A"];
B = QuantumOperator[u@QuantumOperator["Z"]@SuperDagger[u], "Label" -> "B"];
```

Note that we consider the case where {A,B}=0:

```wl
(A@B + B@A)["Matrix"] // Chop // MatrixForm
```

The action of following circuit on an initial state $|+\rangle $$|\psi\rangle $ is as follows: $\frac{1}{2}$$|0\rangle $$\{\mathrm{A}, \mathrm{B} \}$$|\psi\rangle $+$|1\rangle $$[\mathrm{A}, \mathrm{B}]$$|\psi\rangle $, which means if qubit 1 is observed in the state $|0\rangle $ then it confirms $[\mathrm{A}, \mathrm{B}]$$=0$; and if qubit 1 is observed in the state $|1\rangle $ then it confirms $\{\mathrm{A}, \mathrm{B} \}$$=0$.

Construct the circuit based on a quantum switch (defined in the initialization):

```wl
qc = QuantumCircuitOperator[{QuantumOperator[{"Switch", A, B}], 
    QuantumOperator["H"], QuantumMeasurementOperator[]}];
qc["Diagram"]
```

Get the tensor product for $|+\rangle |RandomPure\rangle $:

```wl
\[Psi]0 = 
 QuantumTensorProduct[QuantumState["Plus"], QuantumState["RandomPure"]]
```

Perform the measurement:

```wl
qc[\[Psi]0]["ProbabilityPlot"]
```

The result confirms the anti-commuting case.
