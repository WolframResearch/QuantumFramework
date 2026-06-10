---
Template: Symbol
Name: QuantumCircuitOperator
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumCircuitOperator
Keywords: [quantum circuit, gate sequence, circuit operator, quantum gate, circuit composition, Qiskit, QASM]
SeeAlso: [QuantumState, QuantumOperator, QuantumMeasurementOperator]
RelatedGuides: []
---

## Usage

<code>[QuantumCircuitOperator]()[{$\mathit{obj}_{1}$,$\mathit{obj}_{2}$,$\mathit{obj}_{3}$,...}]</code> represents a quantum circuit with a list of quantum objects $\mathit{obj}_{i}$, e.g., quantum operator, quantum channel, quantum state or quantum measurement operators.

<code>[QuantumCircuitOperator]()[*obj*]</code> wraps a single quantum object (operator, channel, state or measurement operator) as a one-element circuit.

<code>[QuantumCircuitOperator]()[*obj* -> *order*]</code> places *obj* on the qudit positions given by *order*.

<code>[QuantumCircuitOperator]()["*name*"]</code> or <code>[QuantumCircuitOperator]()["*name*"[*args*]]</code> gives a named circuit, such as `"Bell"`, `"Toffoli"` or `"Grover"[f]`.

<code>[QuantumCircuitOperator]()[*qs*]</code> gives a circuit that prepares the quantum state *qs*.

<code>[QuantumCircuitOperator]()["*string*"]</code> gives a tensor-product gate from a string of single-qubit gate letters, such as `"IXYZ"`.

<code>[QuantumCircuitOperator]()[*qasm*]</code> imports a circuit from an OpenQASM string or file.

<code>[QuantumCircuitOperator]()[*qc*, *opts*]</code> rebuilds an existing circuit *qc* with new options.

## Details & Options

- QuantumCircuitOperator[…]["*prop*"] will give the property *prop* of the circuit.

- The following properties *prop* are supported:

|   |   |
|---|---|
| `"Arity"` | the arity of the circuit: the number of input qudits it acts on |
| `"Width"` | the number of wires (qudits) the circuit spans |
| `"Depth"` | the depth of the circuit: the number of gate layers |
| `"OperatorCount"` or `"GateCount"` | the total number of operators (gates) in the circuit |
| `"InputDimensions"` / `"OutputDimensions"` | the list of input / output qudit dimensions of the circuit |
| `"Operators"` | the list of operators (gates) in the circuit |
| `"CircuitOperator"` | the circuit operator: the single quantum operator representing the overall transformation done by the circuit |
| `"Association"` | the underlying data of the circuit, as an Association |
| `"Orders"` | the list of orders (input/output qudits) of the operators in the circuit |
| `"Label"` | the label of the circuit |
| `"Picture"` | the picture of the circuit (e.g. Schrodinger, Heisenberg or PhaseSpace) |
| `"Parameters"` | the list of parameters of a parametrized circuit |
| `"ParameterArity"` | the number of parameters of the circuit |
| `"Properties"` | the list of supported properties |
| `"Diagram"` | the diagram of the quantum circuit |
| `"Topology"` | the circuit topology graph, with qudits/wires as vertices and gates as edges |
| `"TensorNetwork"` | the corresponding tensor network representing the circuit |
| `"TensorNetworkGraph"` | the graph of the circuit's tensor network |
| `"QASM"` | the circuit in OpenQASM 3.0 |
| `"QiskitCircuit"` or `"Qiskit"` | the corresponding circuit in Qiskit (note you need Python as external evaluator for this property) |
| `"Shift"`,*n* | shift the order of the circuit (the qudits it acts on) by *n* |

- There are many named-circuits, <code>[QuantumCircuitOperator]()[{"*circuitname*",...}]</code>, for example, Bell, Toffoli, Grover (and its many oracles), quantum phase estimation, quantum fourier transform, Bernstein-Vazirani algorithm (and its many oracles), Graph and etc.

|   |   |
|---|---|
| `"Bell"` | Bell-state preparation circuit ({H, CNOT}); prepares a Bell pair from the all-zeros state |
| `"GHZ"[n]` | preparation circuit for the n-qubit GHZ state (H, then a CNOT chain) |
| `"Graph"[g]` | graph-state preparation circuit for a graph g (Hadamard on each vertex, then an entangling gate along each edge) |
| `"QuantumState"[qs]` or `"StatePreparation"[qs]` | a circuit that prepares the quantum state qs (qubit states via the multiplexer / RZY decomposition, general qudits via a block-diagonal construction) |
| `"Number"[n,qubits]` | encodes the integer n as a computational basis state (X gates on the 1-bits of n) |
| `"PhaseNumber"[n,qubits,h]` | encodes the number n in qubits qubits as relative phases; with a leading Hadamard layer (one H per qubit) when h is True |
| `"Toffoli"` | Toffoli (CCNOT) gate as a Clifford+T decomposition (15 gates) |
| `"Fredkin"` | Fredkin (controlled-SWAP) gate as a decomposition (18 gates) |
| `"Magic"` | the magic-basis circuit mapping the computational basis to the Bell basis ({S, S, H, CNOT}) |
| `"Fourier"[n]` | n-qubit quantum Fourier transform (QFT) |
| `"InverseFourier"[n]` | n-qubit inverse quantum Fourier transform |
| `"BooleanOracle"[f]` | Boolean oracle for a Boolean function f (ESOP-decomposed into multi-controlled NOT gates) |
| `"BooleanOracleR"[f]` | Boolean oracle for a Boolean function f built from rotation gates (controlled rotations and CNOTs; default rotation `"RZ"[Pi]`) |
| `"PhaseOracle"[f]` | phase oracle for a Boolean function f (encodes f in the phase via kickback) |
| `"GrayOracle"[angles,prec]` | Gray-code-optimized rotation oracle encoding a list or function of angles to precision prec (default `"RY"` rotations) |
| `"Grover"[f]` | Grover operator (oracle and diffusion) for a Boolean function f |
| `"Grover0"[f]` | Grover operator for f, with the diffusion reflecting about the all-zeros state (controlled-on-0 gate) |
| `"GroverPhase"[f]` | Grover operator for f using a phase oracle |
| `"GroverPhase0"[f]` | Grover operator for f using a phase oracle, with controlled-on-0 diffusion |
| `"GroverDiffusion"[n]` | the Grover diffusion (amplitude-amplification) operator on n qubits |
| `"GroverDiffusion0"[n]` | Grover diffusion on n qubits, with a controlled-on-0 gate |
| `"GroverPhaseDiffusion"[n]` | Grover diffusion on n qubits, for a phase oracle |
| `"GroverPhaseDiffusion0"[n]` | Grover diffusion on n qubits, for a phase oracle and with a controlled-on-0 gate |
| `"Deutsch"[f]` | Deutsch circuit for a Boolean function f (single-qubit Deutsch-Jozsa) |
| `"DeutschPhase"[f]` | Deutsch circuit for f using a phase oracle |
| `"DeutschJozsa"[k,n]` or `"DeutschJozsa"[BooleanFunction[k-1,n]]` | Deutsch-Jozsa circuit for the *(k-1)*th Boolean function in *n* variables |
| `"DeutschJozsaPhase"[k,n]` or `"DeutschJozsaPhase"[BooleanFunction[k-1,n]]` | Deutsch-Jozsa circuit using the phase oracle of the *(k-1)*th Boolean function in *n* variables |
| `"DeutschJozsaBooleanOracle"[k,n]` or `"DeutschJozsaBooleanOracle"[BooleanFunction[k-1,n]]` | Deutsch-Jozsa Boolean oracle of the *(k-1)*th Boolean function in *n* variables |
| `"DeutschJozsaPhaseOracle"[k,n]` or `"DeutschJozsaPhaseOracle"[BooleanFunction[k-1,n]]` | Deutsch-Jozsa phase oracle of the *(k-1)*th Boolean function in *n* variables |
| `"BernsteinVazirani"[s]` | Bernstein-Vazirani circuit for a secret bit-string s |
| `"BernsteinVaziraniOracle"[s]` | Bernstein-Vazirani oracle for a secret binary vector or bit-string s |
| `"Simon"[s]` | Simon circuit for a secret binary vector or bit-string s |
| `"SimonOracle"[s]` | Simon oracle for a secret binary vector or bit-string s |
| `"PhaseEstimation"[u,n]` | quantum phase-estimation circuit for a unitary operator u using n counting qubits (qubit operators only) |
| `"Trotterization"[ops,order,reps,t]` | Trotter-Suzuki decomposition of a list of operators ops, for a given order (1, 2, 4, 6, …), number of repetitions reps, and time t (based on Suzuki's recursion relation) |
| `"Controlled"[qc,c1,c0]` or `"C"[qc,c1,c0]` | the controlled version of a circuit qc, controlled on the c1 qubits (value 1) and the c0 qubits (value 0) |
| `"Multiplexer"[op1,op2,...]` | a multiplexer (uniformly-controlled gate): applies the operator selected by the binary value of the control qubits |
| `"ControlledMultiplexer"[m,v]` | the multiplexer-based linear-solver circuit produced by QuantumLinearSolve for a matrix m and vector v |
| `"Switch"[a,b]` | quantum-switch circuit of two operators a and b (their order of application coherently controlled) |
| `"CHSH"` | CHSH Bell-test circuit with measurement angle $\theta$ (default $\pi/4$) |
| `"LeggettGarg"` | Leggett-Garg circuit for the Leggett-Garg inequality, with measurement angle $\theta$ (default $\pi/4$) |
| `"WignerCHSH"` | CHSH circuit in the Wigner (WignerMIC) basis, with measurement angle $\theta$ (default $\pi/4$) |

Tutorials

Related Workflows

XXXX

## Basic Examples

---

Create a quantum circuit composed of only single qubit gates:

```wl
qc = QuantumCircuitOperator[{QuantumOperator["X"], 
    QuantumOperator["Y", {2}]}];
```

Draw the associated circuit diagram:

```wl
qc["Diagram"]
```

Shorthand representation of the circuit:

```wl
QuantumCircuitOperator["XY"]["Diagram"]
```

Another shorter version of input, to create above circuit:

```wl
QuantumShortcut@qc
```

```wl
QuantumCircuitOperator[%]["Diagram"]
```

---

Create a circuit with multi-qubit gates:

```wl
qc = QuantumCircuitOperator[{QuantumOperator["CNOT"], 
    QuantumOperator["H" -> {1, 2}]}];
```

Associated circuit diagram:

```wl
qc["Diagram"]
```

The above circuit has a shorthand representation, too:

```wl
QuantumCircuitOperator[{"CNOT", "H" -> {1, 2}}] == qc
```

The dashed line between two Hadamards imply that they are defined as tensor product, although they are separable. Compared it with this case:

```wl
QuantumCircuitOperator[{"CNOT", "H", "H" -> {2}}]["Diagram"]
```

---

A quantum state is transformed by a quantum circuit

```wl
QuantumCircuitOperator[{"CNOT", "H", "S" -> {2}}][QuantumState["01"]]
```

If no input is given, the initial state is set as the corresponding register state

```wl
QuantumCircuitOperator[{"CNOT", "H", "S" -> {2}}][]
```

---

Quantum circuit can include many quantum objects:

```wl
qc = QuantumCircuitOperator[{QuantumState["Bell", "Label" -> "Bell"], 
    "XY", {"C", "RY"[\[Pi]/4]}, "BitFlip"[.3], "H" -> 2, 
    "P"[\[Pi]/3] -> 2, {1}, "M"["X"] -> 2}];
qc["Diagram"]
```

```wl
qc[]["ProbabilityPlot"]
```

## Scope

### Common quantum gates

`QuantumCircuitOperator` supports a large family of common quantum gates.

Some single-qubit quantum gates:

```wl
QuantumCircuitOperator[{"X", "Y", "Z", "P"[\[Theta]], 
   "RX"[Subscript[\[Phi], x]], "RY"[Subscript[\[Phi], y]], 
   "RZ"[Subscript[\[Phi], z]], "H" -> 2, "S" -> 2, "T" -> 2, 
   "V" -> 2, "NOT" -> 2, "RootNOT" -> 2, "GlobalPhase"[\[CapitalPhi]],
    "PhaseShift"[2] -> 2}]["Diagram"]
```

Some multi-qubit quantum gates:

```wl
QuantumCircuitOperator[{"CH", "CS", "CNOT", "CT", "CX", "CY", "CZ", 
   "CP"[\[Theta]], {"C", "RX"[Subscript[\[Theta], x]] -> 4, {3}}, {"C",
     "RY"[Subscript[\[Theta], y]] -> 4, {3}}, {"C", 
    "RZ"[Subscript[\[Theta], z]] -> 4, {3}}, {"R", Subscript[\[Theta],
      xx], "XX"} -> {3, 4}, {"R", Subscript[\[Theta], yy], 
     "YY"} -> {3, 4}, {"R", Subscript[\[Theta], zz], "ZZ"} -> {3, 4}, 
   "SWAP" -> {3, 4}, "RootSWAP" -> {3, 4}, 
   "Permutation" -> {3, 4} -> {4, 3}, 
   "CSWAP", {"C", "RandomUnitary", {3}, {4}}}]["Diagram"]
```

### Named built-in quantum circuits

[QuantumCircuitOperator]() supports a large family of common quantum circuit.

#### Bell, Cup and Cap circuits

Bell circuit:

```wl
QuantumCircuitOperator["Bell"]["Diagram"]
```

The action of Bell circuit on the register state generates the Bell state:

```wl
QuantumCircuitOperator["Bell"][] == QuantumState["Bell"]
```

Bell circuit is similar to "Cup", with $\sqrt{2}$ factor

```wl
QuantumCircuitOperator[{"Cup"}]["Diagram"]
```

With a $\sqrt{2}$ factor difference, Cup generates a Bell state:

```wl
QuantumCircuitOperator[{"Cup"}][]
```

The cap circuit is the inverse of Cup:

```wl
QuantumCircuitOperator[{"Cap"}]["Diagram"]
```

#### Toffoli circuit

The Toffoli circuit as a combination of one-qubit and two-qubit gates:

```wl
QuantumCircuitOperator["Toffoli"]["Diagram"]
```

```wl
QuantumCircuitOperator["Toffoli"]["CircuitOperator"] == 
 QuantumOperator["Toffoli"]
```

#### Fourier and InverseFourier circuit

For the quantum Fourier circuit, we have <code>[QuantumCircuitOperator]()["Fourier"[n]]</code> with *n* the number of qubits. For example, the Fourier circuit of 2-qubits:

```wl
QuantumCircuitOperator["Fourier"]["Diagram"]
```

Fourier circuit transforms computational basis into momentum (Pauli-X) basis. Check $|00\rangle$ is transformed into $|x_{+}x_{+}\rangle$:

```wl
QuantumState[QuantumCircuitOperator["Fourier"][], "XX"]
```

Fourier circuit of 5-qubits:

```wl
QuantumCircuitOperator["Fourier"[5]]["Diagram"]
```

InverseFourier circuit of 5-qubits:

```wl
QuantumCircuitOperator["InverseFourier"[5]]["Diagram"]
```

InverseFourier circuit transforms from the momentum (Pauli-X) basis into the computational basis into . Check $|x_{+}x_{+}\rangle$ is transformed into $|00\rangle$:

```wl
QuantumCircuitOperator["InverseFourier"][QuantumState["++"]]
```

InverseFourier is the same as conjugate transpose of Fourier:

```wl
QuantumCircuitOperator["InverseFourier"] == 
 QuantumCircuitOperator["Fourier"]["Dagger"]
```

QFT is the same as basis transformation into the Fourier basis:

```wl
n = 3;
SparseArray[
  QuantumCircuitOperator["Fourier"[n]][#]["StateVector"] & /@ 
   QuantumBasis[2, n]["BasisStates"]] == 
 QuantumBasis["Fourier"[2^n]]["Elements"]
```

#### Phase estimation circuit

The quantum phase estimation <code>[QuantumCircuitOperator]()["PhaseEstimation"[u,n]]</code> takes two input arguments: a unitary operator u, and an integer n. The integer n specifies the number of qubits and $\mathit{controlled}-u^{j}$ operators in the circuit, with `j=0,1,…,n-1`. The accuracy of phase estimation and the success probability depends on n.

Generate a random unitary operator:

```wl
u = QuantumOperator["RandomUnitary"];
```

Create the corresponding QPE circuit to find one of its eigenvalues:

```wl
qc = QuantumCircuitOperator["PhaseEstimation"[u]];
qc["Diagram"]
```

Calculate corresponding measurement probabilities:

```wl
N[qc][]["ProbabilityPlot", "LabelsAngle" -> \[Pi]/2, 
 AspectRatio -> 1/4]
```

Given the result with max probability, estimate the eigenvalue:

```wl
Keys[TakeLargest[N[qc][]["Probabilities"], 1]][[1]]["Name"]
est1 = FromDigits[%, 2]/2^4 // N
```

Generate another QPE circuit, to obtain the other eigenvalue:

```wl
qc2 = qc[[2 ;;]];
qc2["Flatten"]["Diagram"]
```

Note the only difference is dropping the Pauli-X gate.

Calculate corresponding measurement probabilities:

```wl
N[qc2][]["ProbabilityPlot", "LabelsAngle" -> \[Pi]/2, 
 AspectRatio -> 1/4]
```

Given the result with max probability, estimate the eigenvalue:

```wl
Keys[TakeLargest[N[qc2][]["Probabilities"], 1]][[1]]["Name"]
est2 = FromDigits[%, 2]/2^4 // N
```

Find corresponding eigenvalues (note they should be in the form $e^{i\, 2\pi \, \lambda }$):

```wl
\[Lambda] = Mod[Arg[u["Eigenvalues"]]/(2 \[Pi]), 1]
```

Calculate the estimated eigenvalues, compared to the expected ones and find the error%:

```wl
100. (1 - Sort[{est1, est2}]/Sort[\[Lambda]])
```

#### Graph circuit

Given a graph, the corresponding circuit (which generates quantum graph/cluster states) is a series of controlled-Z operators on edges of graph:

```wl
g = Graph[RandomGraph[{4, 5}], VertexLabels -> "Name"]
```

```wl
QuantumCircuitOperator["Graph"[g]]["Diagram"]
```

Graph circuit generates graph state (also called as cluster state):

```wl
QuantumCircuitOperator["Graph"[g]][] == QuantumState["Graph"[g]]
```

Check qubits on the same edge are entangled:

```wl
With[{st = QuantumState["Graph"[g]]}, 
 QuantumEntangledQ[st, #] &@*List @@@ EdgeList[g]]
```

#### Bernstein-Vazirani circuit

For the Bernstein-Vazirani algorithm, one needs to specify a secret string bit st as st as follows: for the oracle we have <code>[QuantumCircuitOperator]()[[BernsteinVaziraniOracle]()[st]]</code>, and for the whole circuit e <code>[QuantumCircuitOperator]()[[BernsteinVazirani]()[st]]</code>. For example, the Bernstein-Vazirani oracle for the secret bit of 101:

```wl
QuantumCircuitOperator["BernsteinVaziraniOracle"["101"]]["Diagram"]
```

The Bernstein-Vazirani circuit for the secret bit of 101:

```wl
QuantumCircuitOperator["BernsteinVazirani"["101"]]["Diagram"]
```

Measurement results of the Bernstein-Vazirani circuit for the secret bit of 101:

```wl
QuantumCircuitOperator[
   "BernsteinVazirani"["101"]][]["ProbabilitiesPlot"]
```

BernsteinVaziraniOracle circuit can take a binary vector, or bit-strings

```wl
QuantumCircuitOperator["BernsteinVaziraniOracle"[{1, 0, 1}]] == 
 QuantumCircuitOperator["BernsteinVaziraniOracle"["101"]]
```

#### Simon circuit

Simon’s algorithm is a quantum algorithm designed to solve the Simon’s problem efficiently. It aims to find a hidden bit-string by querying a black-box function.

Create corresponding Simon oracle:

```wl
so = QuantumCircuitOperator["SimonOracle"[{1, 0, 1, 1}]];
so["Diagram"]
```

A function that return the Boolean table of a Simon oracle:

```wl
boolTable[oracle_] := 
 With[{n = oracle["Width"]/2}, 
  KeyValueMap[TraditionalForm /@ {#1, #2} &]@
   AssociationMap[
    QuantumOperator["Discard"[2] -> Range[n]][oracle[#]] &, 
    QuantumBasis[2, n]["BasisStates"]]]
```

Return Boolean table of above Simon oracle:

```wl
TableForm[boolTable[so], TableDepth -> 2, 
 TableHeadings -> {None, {"bit", "oracle[bit]"}}]
```

Simon circuit for a given hidden sequence:

```wl
SeedRandom[0];
simon = QuantumCircuitOperator["Simon"[seq = {1, 0, 0, 1}]];
simon["Diagram"]
```

The measurement probabilities given the Simon circuit for the binary string of 1001:

```wl
simonMeas = N[simon][]
```

Plot probabilities:

```wl
simonMeas["ProbabilitiesPlot", AspectRatio -> 1/4, 
 "LabelsAngle" -> \[Pi]/2]
```

List of all potential results:

```wl
m = Normal /@ Keys@Select[simonMeas["Probabilities"], Positive]
```

Find the null vector (`matrix·m=0 mod 2`):

```wl
NullSpace[m, Modulus -> 2]
```

We found the secret bit:

```wl
seq == %[[1]]
```

Note, the Simon circuit and its corresponding Boolean function is not uniquely determined. For each evaluation, one gets a different oracle that can solve the Simon problem:

```wl
simon2 = QuantumCircuitOperator["Simon"[seq]];
simon == simon2
```

With the new circuit, one can find the hidden sequence:

```wl
NullSpace[Normal /@ Keys@Select[simon2[]["Probabilities"], Positive], 
   Modulus -> 2][[1]] == seq
```

#### GrayOracle

It takes a list of angles or a function, number of ancillary qubits and target qubit, to generate the corresponding oracle based on [Gray code](https://arxiv.org/pdf/quant-ph/0404089.pdf).

Given 4-qubits as ancillary, generate angles for all possible bit-strings, and prepend zero angle, too:

```wl
prec = 4;
angles = 
  1/2^prec 1/FromDigits[#, 2] & /@ Rest@Tuples[{0, 1}, prec] // N // 
   Prepend[0];
```

Generate the corresponding GrayOracle:

```wl
QuantumCircuitOperator["GrayOracle"[angles, prec, 1]]["Diagram", 
 "ShowGateLabels" -> False]
```

It is the same as the following oracle, using a function as input:

```wl
QuantumCircuitOperator[
  "GrayOracle"[If[# == 0, 0, 1./#] &, prec, 1]]["Diagram", 
 "ShowGateLabels" -> False]
```

GrayOracle is the same as Multiplexer with angles:

```wl
inv = QuantumCircuitOperator["Multiplexer" @@ "RY" /@ angles, 
   RotateLeft[Range[prec + 1]]];
inv["Diagram", "ShowGateLabels" -> False]
```

Check the Multiplexer is the same as GrayOracle:

```wl
(inv["Matrix"] - 
    QuantumCircuitOperator[
      "GrayOracle"[If[# == 0, 0, 1./#] &, prec, 1]]["Matrix"])[
  "ExplicitValues"] // Chop
```

#### Grover circuit (and family of Boolean oracles)

There are a collection of named-circuits for the Grover search algorithm, and its oracles (e.g., Boolean, or Phase oracles) or amplification. The class of Grover-related named-circuits are <code>[QuantumCircuitOperator]()[name[boo]]</code> with *boo* a Boolean function.

Boolean oracle for a given Boolean function:

```wl
bo = QuantumCircuitOperator["BooleanOracle"[a && ! b]];
bo["Diagram"]
```

Boolean oracle will flip the ancillary qubit (in above example, 3rd qubit) only if the quantum state is a solution of the Boolean function. For example, `{True,False}` (10) is a solution of `a&&!b`, thus the ancillary qubit will be flipped:

```wl
bo[QuantumState["100"]]["Formula"]
```

One can also get the decomposition of a Boolean oracle, in terms of CNOTs and Z-rotation:

```wl
QuantumCircuitOperator["BooleanOracleR"[a && ! b]]["Diagram"]
```

For the diffusion part of the Grover circuit (also called, the amplification part), one only needs to specify the number of qubits:

```wl
QuantumCircuitOperator["GroverDiffusion"[3]]["Diagram"]
```

```wl
QuantumCircuitOperator["GroverDiffusion0"[3]]["Diagram"]
```

Another way of storing the solution of a Boolean function is to save it as a phase, rather than on an ancillary qubit. This can be archived by using a phase oracle, that add an overall phase of π if the quantum state is a solution of the Boolean function.

```wl
po = QuantumCircuitOperator["PhaseOracle"[a && ! b]];
po["Diagram"]
```

```wl
po[QuantumState["10"]]["Formula"]
```

The corresponding diffusion part, given a phase oracle:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion"[2]]["Diagram"]
```

One can also create the overall Grover circuit, using controlled-1 gates and a Boolean oracle:

```wl
QuantumCircuitOperator["Grover"[a && ! b]]["Diagram"]
```

One can also create the overall Grover circuit, using controlled-0 gates and a Boolean oracle:

```wl
QuantumCircuitOperator["Grover0"[a && ! b]]["Diagram"]
```

One can also create the overall Grover circuit, using controlled-1 gates and a phase oracle:

```wl
QuantumCircuitOperator["GroverPhase"[a && ! b]]["Diagram"]
```

One can also create the overall Grover circuit, using controlled-0 gates and a phase oracle:

```wl
QuantumCircuitOperator["GroverPhase0"[a && ! b]]["Diagram"]
```

#### CHSH circuit

CHSH quantum circuit:

```wl
QuantumCircuitOperator["CHSH"]["Diagram", Sequence[
 "ShowMeasurementWire" -> False, "WireLabels" -> {
Placed["a", Right], 
Placed["x", Right], 
Placed["y", Right], 
Placed["b", Right]}, Epilog -> {
Opacity[0.3], 
Text[
Style["a\[Xor]b\!\(\*OverscriptBox[\(=\), \(?\)]\)x\[And]y", Black, 
      20, Bold], {9., -2.35}]}, PlotRange -> {{0, 10}, {-5, 0}}, 
  ImageSize -> Large]]
```

Calculate the probability of winning:

```wl
NProbability[
 BitAnd[x, y] == BitXor[a, b], {a, x, y, b} \[Distributed] 
  QuantumCircuitOperator["CHSH"][]["MultivariateDistribution"]]
```

#### Trotter-Suzuki

One can find Trotter-Suzuki decomposition, to approximate a time evolution.

First order, only:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator[
      "X", {1, 2}], QuantumOperator["X", {2, 3}]}, 1, 1, 1}];
qc["Diagram"]
```

2nd order:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator[
      "X", {1, 2}], QuantumOperator["X", {2, 3}]}, 2, 1, 1}];
qc["Diagram"]
```

2nd order, with 4 steps:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator[
      "X", {1, 2}], QuantumOperator["X", {2, 3}]}, 2, 4, 1}];
qc["Diagram", FontSize -> 8]
```

4th order, only 1 step:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator[
      "X", {1, 2}], QuantumOperator["X", {2, 3}]}, 4, 1, 1}];
qc["Diagram", "ShowGateLabels" -> False]
```

#### PhaseNumber circuit

PhaseNumber circuit encodes a number in n-qubit.

Encode 5 in 4-qubits:

```wl
num = 5;
QuantumCircuitOperator["PhaseNumber"[num, 4]]["Diagram"]
```

Remove Hadamards:

```wl
QuantumCircuitOperator["PhaseNumber"[num, 4, False]]["Diagram"]
```

PhaseNumber followed by InverseFourier generates a circuit that can be used to calculate the encoded number:

```wl
num = 5;
qc = QuantumCircuitOperator[{"PhaseNumber"[num, 6], 
    "InverseFourier"[6], Range[6]}];
qc["Diagram"]
```

Corresponding measurement probabilities:

```wl
m = N[qc][]
```

Given the result with max probability, find the corresponding integer from bit-string:

```wl
FromDigits[Keys[TakeLargest[m["Probabilities"], 1]][[1]]["Name"], 2]
```

Check it is the same as initial encoded number:

```wl
% == num
```

The PhaseNumber circuits follow the usual addition:

```wl
qc1 = QuantumCircuitOperator[{"PhaseNumber", 5, 6}];
qc2 = QuantumCircuitOperator[{"PhaseNumber", 3, 6, False}];
(qc1/*qc2)["Diagram"]
```

Test that 5+3=8, in terms of their circuits:

```wl
(qc1/*qc2)["CircuitOperator"] == 
 QuantumCircuitOperator[{"PhaseNumber", 5 + 3, 6}]["CircuitOperator"]
```

#### Deutsch-Jozsa circuit

For the Deutsch-Jozsa algorithm, the goal is to determine whether a Boolean function is constant (producing the same output regardless of the input, either 0 or 1) or balanced (yielding an equal number of 0s and 1s).

Define a balanced Boolean function:

```wl
bf = BooleanFunction[2, 2];
BooleanConvert[bf][\[FormalX], \[FormalY]]
```

Show its Boolean table (note same number of True and False)

```wl
BooleanTable[bf]
```

```wl
QuantumCircuitOperator["DeutschJozsa"[bf]]["Diagram"]
```

One can also get the Deutsch-Jozsa with a Phase oracle:

```wl
QuantumCircuitOperator["DeutschJozsaPhase"[bf]]["Diagram"]
```

Also, “DeutschJozsaPhase” can take two integers, with the first one $k+1$ and second one as $n$, where the Boolean function is defined as BooleanFunction[k,n]

```wl
QuantumCircuitOperator["DeutschJozsaPhase"[3, 2]]["Diagram"]
```

Given all Boolean functions with two variables, show number of True and False in their solutions, and also corresponding Deutsch-Jozsa probabilities:

```wl
Grid[Prepend[
  Table[With[{bf = 
      BooleanFunction[f, 2]}, {BooleanConvert[
       bf][\[FormalX], \[FormalY]], Count[True][BooleanTable[bf]], 
     Count[False][BooleanTable[bf]], 
     With[{rs = 
        QuantumCircuitOperator["DeutschJozsa"[bf]][]["Probabilities"]},
       BarChart[rs, ChartLabels -> Keys[rs], AspectRatio -> 1/3, 
       ImageSize -> Small]]}], {f, 0, 15}], 
  Style[#, Bold] & /@ {"Boolean function", "True", "False", 
    "Deutsch-Jozsa probabilities"}], Frame -> All, 
 Alignment -> {{Left, Center, Center, Center}}]
```

#### Leggett-Garg:

```wl
QuantumCircuitOperator["LeggettGarg"]["Diagram"]
```

#### State-preparation

Generate a random pure state of 3-qubits:

```wl
state = QuantumState["RandomPure"[3]]
```

Generate a circuit whose outcome is above state:

```wl
qc = QuantumCircuitOperator["QuantumState"[state]];
qc["Diagram", "ShowGateLabels" -> False]
```

Check the execution of the circuit gives the expected state:

```wl
qc[] == state
```

One can also use "StatePreparation":

```wl
QuantumCircuitOperator["StatePreparation"[state]]["Diagram", 
 "ShowGateLabels" -> False]
```

### Quantum circuit diagram

One can customize many features of a circuit diagram. See [Circuit Diagram]() documentation page, for more details.

```wl
QuantumCircuitOperator[{"H", "CY", "H", {1}}]["Diagram", 
 "WireLabels" -> {{Placed["1st qubit", Left], 
    Placed["Eigenvector of Y", Right]}, {Placed["2nd qubit", Left]}}, 
 "MeasurementWireLabel" -> "Measurement Wire"]
```

Also, the overall size of the circuit can be adjusted:

```wl
superdense = 
  QuantumCircuitOperator[{"H", "CNOT", "CX" -> {4, 1}, "CZ" -> {3, 1},
     "CNOT", "H", {1}, {2}}];
superdense["Diagram", "Size" -> .5]
```

Additionally, one can customize the horizontal and vertical gaps:

```wl
superdense["Diagram", "HorizontalGapSize" -> 3]
```

Customize the vertical gap:

```wl
superdense["Diagram", "VerticalGapSize" -> 2]
```

Another interesting feature for the Diagram is showing the connectors as dot, specifying explicitly on what wires a gate acts:

```wl
QuantumCircuitOperator[{QuantumOperator[{"R", \[Theta], "XX"}, {1, 
     3}]}]["Diagram", "ShowConnectors" -> True]
```

But a portion of a circuit in a gray box, to emphasizing one portion:

```wl
QuantumCircuitOperator[{"CNOT", 
   QuantumCircuitOperator[{QuantumOperator["SWAP", {2, 4}], 
     QuantumOperator[{"Phase", \[Gamma]}]}, "MyCir"]}]["Diagram"]
```

Drop the gate labels (if you hover you mouse over the diagram, the labels will be shown):

```wl
QuantumCircuitOperator[{"H", "CY", "H", {1}}]["Diagram", 
 "ShowGateLabels" -> False]
```

### Tensor network of a quantum circuit

Corresponding tensor network of a quantum circuit can be generated as a property. See [Tensor Network]() documentation page, for more details.

```wl
tn = QuantumCircuitOperator["Grover"[a && ! b || c]]["TensorNetwork"]
```

Corresponding index graph:

```wl
TensorNetworkIndexGraph[tn, 
 GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}]
```

### Parametrized quantum circuit

Define a parametrized quantum circuit, with two parameters:

```wl
qc = QuantumCircuitOperator[{"RY"[\[Theta]], "RZ"[\[Phi]]}, 
  "Parameters" -> {\[Theta], \[Phi]}]
```

Assign numerical values to parameters:

```wl
qc[<|\[Theta] -> \[Pi]/3, \[Phi] -> \[Pi]/4|>]["Diagram"]
```

### Properties and downvalues

A `QuantumCircuitOperator` exposes a rich set of properties, accessed as `qc["property"]` (its downvalues). Define a sample circuit:

```wl
qc = QuantumCircuitOperator[{"H", "CNOT", "RZ"[Pi/4] -> 2}]
```

The full list of supported properties is returned by `"Properties"`:

```wl
qc["Properties"]
```

Sizes and counts of the circuit (`"GateCount"` is an alias of `"OperatorCount"`):

```wl
Dataset @ AssociationMap[qc, {"Arity", "Width", "Depth", "OperatorCount", 
   "GateCount"}]
```

The input and output qudit dimensions:

```wl
{qc["InputDimensions"], qc["OutputDimensions"]}
```

The operators (gates) that make up the circuit:

```wl
qc["Operators"]
```

The order (input and output qudits) of each operator:

```wl
qc["Orders"]
```

The circuit label and its picture:

```wl
{qc["Label"], qc["Picture"]}
```

The underlying data, as an association:

```wl
qc["Association"]
```

The whole circuit collapsed into a single quantum operator:

```wl
qc["CircuitOperator"]
```

The circuit diagram:

```wl
qc["Diagram"]
```

The circuit topology, with qudits/wires as vertices and gates as edges:

```wl
qc["Topology"]
```

The tensor network representing the circuit:

```wl
qc["TensorNetwork"]
```

The graph of that tensor network:

```wl
qc["TensorNetworkGraph"]
```

Export the circuit to OpenQASM 3.0:

```wl
qc["QASM"]
```

The corresponding Qiskit circuit (requires Python as an external evaluator):

```wl
qc["QiskitCircuit"]
```

Shift the qudits the circuit acts on by an offset (here, by 2):

```wl
qc["Shift", 2]["Orders"]
```

For a parametrized circuit, `"Parameters"` lists the declared parameters and `"ParameterArity"` counts them:

```wl
pqc = QuantumCircuitOperator[{"RY"[\[Theta]], "RZ"[\[Phi]] -> 2}, 
   "Parameters" -> {\[Theta], \[Phi]}]
```

```wl
{pqc["Parameters"], pqc["ParameterArity"]}
```

## Generalizations & Extensions

---

`QuantumCircuitOperator` also supports controlled operators with more than one control or targets.

```wl
QuantumCircuitOperator[
  QuantumOperator[{"C", "X", {1, 4}, {2}}, {3}]]["Diagram"]
```

---

Decompose Toffoli gate:

```wl
v = QuantumOperator[MatrixPower[PauliMatrix[1], 1/2], "Label" -> "V"];
```

```wl
qc = QuantumCircuitOperator[{QuantumOperator[{"Controlled", v}, {2, 3}],
     QuantumOperator["CX"], 
    QuantumOperator[{"Controlled", v["Dagger"]}, {2, 3}], 
    QuantumOperator["CX"], QuantumOperator[{"Controlled", v}, {1, 3}]}];
qc["Diagram"]
```

Show above implementation is the same as built-in Toffoli gate:

```wl
QuantumOperator["Toffoli"] == qc["CircuitOperator"]
```

---

Quantum circuits consisting of qudit gates are supported. Define a qudit gate of dimensionality 3:

```wl
QuantumCircuitOperator[{"RandomUnitary"[3], 
   "Spider"[QuantumBasis[{3, 2}, {3, 2}]] -> {1, 2} -> {2, 
      1}}]["Diagram", "ShowWireDimensions" -> True]
```

---

One can shift a circuit, which means changing its overall order

```wl
qc = QuantumCircuitOperator[{"Fourier", 3}];
qc["Diagram"]
```

Shift above circuit by 2 qudits:

```wl
qc["Shift", 2]["Diagram"]
```

## Options

### Method

The default method creates the tensor network and then contracts it to get the results:

```wl
QuantumCircuitOperator[{"Bell", {1, 2}}][]
```

The tensor-network engine can also be requested explicitly with `Method -> "TensorNetwork"`; tensor-network options are passed as `Method -> {"TensorNetwork", opts}`:

```wl
QuantumCircuitOperator[{"H", "CNOT"}][QuantumState["00"], Method -> "TensorNetwork"]
```

With `Method -> "Schrodinger"`, the circuit is applied gate by gate (a fold over the operators, simplifying after each step). It is slower than the tensor-network engine, but keeps the result in closed symbolic form, which is convenient for parametrized circuits:

```wl
QuantumCircuitOperator[{"H", "RX"[\[Theta]] -> 2, "CNOT"}][
  QuantumState["00"], Method -> "Schrodinger"]
```

For Clifford (stabilizer) circuits, `Method -> "Stabilizer"` simulates the circuit in the stabilizer formalism and returns a PauliStabilizer, scaling polynomially in the number of qubits:

```wl
QuantumCircuitOperator[{"H", "CNOT"}][Method -> "Stabilizer"]
```

With `Method -> "QuEST"`, the circuit is simulated through the external QuEST high-performance simulator (qubit circuits only, requires the QuESTLink paclet):

```wl
QuantumCircuitOperator[{"H", "CNOT"}][Method -> "QuEST"]
```

Convert the circuit into Qiskit, and simulate it using the Qiskit default backend such as Aer (the results will be quantum objects in Wolfram quantum framework):

```wl
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], 
   Range[5]}][Method -> "Qiskit"]
```

IBMProvider as Qiskit provider:

```wl
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], 
   Range[5]}][Method -> {"Qiskit", "Provider" -> "IBMProvider"}]
```

IBMProvider as Qiskit provider on a specific backend:

```wl
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], 
   Range[5]}][Method -> {"Qiskit", "Provider" -> "IBMProvider", 
   "Backend" -> backend}]
```

A delay can be added into Qiskit circuit too:

```wl
QuantumCircuitOperator[{"I" -> "Delay"[2], "M"}]["Qiskit"]
```

Return the list of IBM-quantum backends and their queue:

```wl
ServiceConnect["IBMQ", "New"]["BackendQueue"]
```

AWSBraket as Qiskit provider on a specific backend:

```wl
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], 
   Range[5]}][Method -> {"Qiskit", "Provider" -> "AWSBraket", 
   "Backend" -> backend}]
```

Return the list of AWSBraket backends:

```wl
aws = ServiceConnect["AWS"];
braket = ServiceExecute[aws, "GetService", {"Name" -> "Braket"}];
devices = braket["SearchDevices", "Filters" -> {}]["Devices", All]
```

For some circuit creation, users can connect to Classiq:

```wl
qc = QuantumCircuitOperator[
   QuantumState[RandomReal[{0, 1}, 2^3]]["Normalized"], 
   Method -> "Classiq"];
qc["Diagram", "ShowGateLabels" -> False]
```

## Applications

### POVM in a quantum circuit

Collective measurement as well as POVM measurement are supported.

Define a POVM measurement in a minimal information-ally complete basis:

```wl
povm = QuantumMeasurementOperator["WignerMICPOVM"]
```

Define a circuit consisting of the above measurements:

```wl
qc = QuantumCircuitOperator[{QuantumState["RandomPure", 
    "Label" -> "RandomState"], povm}]; qc["Diagram"]
```

Show

```wl
qc[]["Probabilities"]
```

### Quantum teleportation

We shall implement a quantum circuit for teleporting a qubit. The 1st and 2nd qubits represent Alice’s system, while the last one is Bob’s. The double lines carry classical bits. The state to be teleported is $|\psi _{1}\rangle=\alpha |0\rangle+\beta |1\rangle$, where α and β are unknown amplitudes.

```wl
circuit = 
  QuantumCircuitOperator[{QuantumCircuitOperator[{QuantumState[\
{\[Alpha], \[Beta]}], "Cup" -> {2, 3}}, "State Prep."], "CX", "H", 
    "CX" -> {2, 3}, "CZ" -> {1, 3}, {1}, {2}}];
```

```wl
circuit["Diagram"]
```

Trace out 1st and 2nd qubit on possible post-measurement states, and compare the reduced state of only qubit 3:

```wl
Thread[QuantumState[{\[Alpha], \[Beta]}] == 
   QuantumOperator["Trace" -> {1, 2}] /@ 
    circuit[]["States"]] // FullSimplify
```

### Quantum key distribution as in BB84 protocol

As an example of a quantum key distribution (QKD), consider the BB84 protocol. Alice generates random bits and for each, randomly chooses either a Z-basis or X-basis of measurement. Consider a random sequence of 9 bits and corresponding basis:

```wl
aliceBits = RandomChoice[{0, 1}, 9]
```

```wl
aliceBases = RandomChoice[{"Z", "X"}, 9]
```

Depending on her choice of bit and base, she follows the following protocol for sending a qubit to Bob:

```wl
aliceProtocol = <|{0, "Z"} -> QuantumState[{1, 0}, "Z"], {0, "X"} -> 
    QuantumState[{1, 0}, "X"], {1, "Z"} -> 
    QuantumState[{0, 1}, "Z"], {1, "X"} -> QuantumState[{0, 1}, "X"]|>;
#["Formula"] & /@ aliceProtocol // Normal // Column
```

With Alice's bits and bases having been generated randomly, she will send Bob these qubits:

```wl
aliceQubits = aliceProtocol /@ Thread[List[aliceBits, aliceBases]];
#["Formula"] & /@ aliceQubits
```

Bob randomly generates a sequence of either Z-bases or X-bases of measurement:

```wl
bobBases = RandomChoice[{"Z", "X"}, 9]
```

Thus, Bob applies the following measurement operators on the qubits he got from Alice:

```wl
ops = QuantumMeasurementOperator /@ bobBases;
```

Post-measurement, it is the case that, whenever Bob's measurement basis is not the same as Alice's, he gets two random results:

```wl
mea = MapThread[#1[#2] &, {ops, aliceQubits}];
```

Return one possible realization of Bob' s series of measurements (noting that some results will be consistent and other random, depending on Alice' s qubit and Bob's measurement):

```wl
bobRealization = #["SimulatedMeasurement"] & /@ mea
```

Bob assigns the bit 0 whenever he gets $|z_{+}\rangle$ or $|x_{+}\rangle$, and the bit 1 for $|z_{-}\rangle$ or $|x_{-}\rangle$:

```wl
bobBits = 
 bobRealization /. 
  Flatten[Alternatives @@ #[[;; 2]] -> #[[-1]] & /@ 
    Transpose[{QuditBasis["Z"]["Names"], 
      QuditBasis["X"]["Names"], {0, 1}}]]
```

Finally, Alice and Bob publicly share with each other the measurement basis that they used (in other words, they conduct a classical communication). If the basis is same, they each keep the bit; if not, they each disregard the corresponding bit. (Note that each does so separately, and each one obtains a secret key, which will be the same as the other's):

```wl
aliceSecretKey = 
 Select[Transpose[{aliceBits, aliceBases, 
     bobBases}], #[[2]] == #[[3]] &][[All, 1]]
```

```wl
bobSecretKey = 
 Select[Transpose[{bobBits, aliceBases, 
     bobBases}], #[[2]] == #[[3]] &][[All, 1]]
```

```wl
bobSecretKey == aliceSecretKey
```

Alice and Bob now share a secret key. Note that the key length is less than the number of qubits Alice sent to Bob (because Bob picks the measurement basis randomly, and he might not get the same outcome as Alice). Alice can keep sending Bob qubits until they reach a desired length. This can be implemented as follows:

```wl
add[{x_, y_}] := 
 Module[{aliceBit, aliceBase, aliceQubit, bobBase, bobBit},
  aliceBit = RandomChoice[{0, 1}];
  aliceBase = RandomChoice[{"Z", "X"}];
  aliceQubit = aliceProtocol[{aliceBit, aliceBase}];
  bobBase = RandomChoice[{"Z", "X"}];
  bobBit = 
   QuantumMeasurementOperator[bobBase][aliceQubit][
     "SimulatedMeasurement"] /. 
    Flatten[Alternatives @@ #[[;; 2]] -> #[[-1]] & /@ 
      Transpose[{QuditBasis["Z"]["Names"], 
        QuditBasis["X"]["Names"], {0, 1}}]];
  If[aliceBase == bobBase, {Append[x, aliceBit], 
    Append[y, bobBit]}, {x, y}]]
```

For example, the key length is set to 5 bits below:

```wl
NestWhile[add, {{}, {}}, Length[#[[1]]] < 5 &]
```

### Parity test

Quantum circuit for parity test. If the number of 1 in a bit string (e.g. 0001101..) is even, the parity is assigned to 0 (even parity), if it is odd, the parity is assigned to 1 (odd parity). A quantum circuit to measure parity of an n-qubit state is given below:

```wl
parityMeasurementCircuit[n_] := 
 QuantumCircuitOperator[
  Append["CNOT" -> {#, n + 1} & /@ Range[n], {n + 1}]]
```

Example of a parity measurement circuit for 3-qubit state:

```wl
parityMeasurementCircuit[3]["Diagram"]
```

Apply this circuit to $|000\rangle\otimes |0\rangle$ (even parity) and $|001\rangle\otimes |0\rangle$ (odd parity):

```wl
{parityMeasurementCircuit[3][QuantumState["0000"]]["ProbabilityPlot"],
  parityMeasurementCircuit[3][QuantumState["0010"]]["ProbabilityPlot"]}
```

### Quantum Approximate Optimization Algorithm  (QAOA)

One can find a binary vector $x$ that minimizes (or maximizes) the quadratic objective $x.q.x+c.x$, with no constrains. This is usually called as Quadratic Unconstrained Binary Optimization (QUBO). The linear term $c$ can be dropped, without losing generality.

Define a $4\times 4$ real-valued upper triangular matrix, as the objective matrix:

```wl
q = UpperTriangularize[RandomReal[{-5, 5}, {4, 4}]];
q // MatrixForm
```

Generate the corresponding objective function:

```wl
Expand@Block[{x = Array[Indexed[\[FormalX], #] &, Length[q]]}, 
  x . q . x]
```

Maximize and minimize the objective function (classical approach):

```wl
solMax = 
 NArgMax[{x . q . x, 0 \[VectorLess] x^2 \[VectorLessEqual] 1}, 
  x \[Element] Vectors[Length[q], Integers]]
solMin = 
 NArgMin[{x . q . x, 0 \[VectorLess] x^2 \[VectorLessEqual] 1}, 
  x \[Element] Vectors[Length[q], Integers]]
```

Define a quantum circuit for any objective function:

```wl
ClearAll[\[Gamma], \[Beta], qaoa]
qaoa[obj_?SquareMatrixQ, p_ : 1] := With[{l = Length[obj]},
  QuantumCircuitOperator[Join[
     Table["+" -> i, {i, l}],
     Catenate@Table[
       {"Barrier",
        Catenate@Table[
          
          If[obj[[i, j]] == 0, Nothing, 
           "R"[2 \[Gamma][m] obj[[i, j]], "ZZ"] -> {i, j}],
          {i, l}, {j, i + 1, l}
          ],
        "Barrier",
        Table["R"[ 2 \[Beta][m], "X"] -> i, {i, l}]}
       , {m, p}]]
    ]["Flatten"]]
```

Generate the circuit for the objective matrix, with only two layers:

```wl
p = 2;
qaoa[q, p]["Diagram", FontSize -> 9, "ShowGateLabels" -> False]
```

Define variable (angles) that should be optimized

```wl
var = Join[Table[\[Gamma][i], {i, p}], Table[\[Beta][i], {i, p}]]
```

Find the corresponding probabilities in the computational basis:

```wl
weight = Abs[qaoa[q, p]["Amplitudes"]]^2;
```

Find $\langle \gamma ,\beta |H|\gamma ,\beta \rangle$ with $H=\sum _{i>j}q_{\mathit{ij}}\sigma _{z,i}\sigma _{z,j}$ and $|\gamma ,\beta \rangle$ the output of above circuit

```wl
weightedEnergy = 
  Association@
   KeyValueMap[#1 -> #2 Total[
        UpperTriangularize[q, 1] UpperTriangularize[
          KroneckerProduct[(-1)^#1["Name"], (-1)^#1["Name"]], 1], 2] &,
     weight];
```

Maximize and minimize the total weighted-mean value

```wl
optMax = 
 Thread[var -> 
   NArgMax[{Total[weightedEnergy], 
     0 \[VectorLessEqual] var \[VectorLessEqual] 2 \[Pi]}, var]]
optMin = 
 Thread[var -> 
   NArgMin[{Total[weightedEnergy], 
     0 \[VectorLessEqual] var \[VectorLessEqual] 2 \[Pi]}, var]]
```

Plot the results given the optimal values:

```wl
Row[{With[{data = N[weightedEnergy /. optMax]}, 
   BarChart[data, 
    ChartLabels -> (Rotate[
         If[(solMax /. -1 -> 0) == #["Name"], 
          Style[#["Name"], Red], #["Name"]], Pi/2] & /@ Keys[data]), 
    PlotLabel -> "Maximizing objective function", ImageSize -> 200]], 
  With[{data = N[weightedEnergy /. optMin]}, 
   BarChart[data, 
    ChartLabels -> (Rotate[
         If[(solMin /. -1 -> 0) == #["Name"], 
          Style[#["Name"], Red], #["Name"]], Pi/2] & /@ Keys[data]), 
    PlotLabel -> "Minimizing objective function", ImageSize -> 200]]}]
```

Two solutions is due to the fact that NOT of solution is also a max/min ($(-x).q.(-x)=x.q.x$)

Check the classical solution is the same as quantum one:

```wl
MemberQ[#["Name"] & /@ 
  Keys[TakeLargestBy[Normal@N[weightedEnergy /. optMax], Last, 2]], 
 solMax /. -1 -> 0]
```

```wl
MemberQ[#["Name"] & /@ 
  Keys[TakeSmallestBy[Normal@N[weightedEnergy /. optMin], Last, 2]], 
 solMin /. -1 -> 0]
```

### Harrow-Hassidim-Lloyd (HHL) algorithm

Given a matrix A and a vector b, we want to solve $A\cdot x=b$

Generate HHL circuit that takes *A*, *b*, and number of ancillary qubits:

```wl
ClearAll[HHL]
HHL[A_?SquareMatrixQ, b_?VectorQ, prec : _Integer?Positive : 4] :=
 Block[{d, qpe, t0, 
   yTilde, \[Lambda]Tilde, \[Theta]\[Lambda], \[Gamma]},
  d = Log2@Length[b];
  t0 = Floor@\[Pi]/Norm[A];
  yTilde = 
   If[0 <= # <= 2^(prec - 1), #, # - 2^prec] & /@ Range[0, 2^prec - 1];
  \[Lambda]Tilde = (2 \[Pi])/t0 yTilde/2^prec;
  \[Gamma] = Min[DeleteCases[Abs[\[Lambda]Tilde], 0. | 0]];
  \[Theta]\[Lambda] = 
   If[# == 0 || # == 0., 0., 2 ArcSin[\[Gamma]/#]] & /@ \[Lambda]Tilde;
  QuantumCircuitOperator[{
     QuantumState[Normalize[b]] -> 1 + prec + Range[d],
     "Barrier",
     qpe = 
      QuantumCircuitOperator[
        "PhaseEstimation"[QuantumOperator[MatrixExp[t0 I A]], prec], 
        1 + Range[prec + d]][[d + 1 ;; -prec - 1]],
     "Barrier",
     QuantumCircuitOperator["GrayOracle"[\[Theta]\[Lambda], prec, 1], 
      "Eigenvalue Inversion"],
     "Barrier",
     qpe["Dagger"],
     "Barrier",
     QuantumOperator[
      QuantumState["1" <> StringJoin[ConstantArray["0", prec]]][
       "Dagger"]]
     }
    ] /; IntegerQ[d]
  ]
```

Generate a complex-valued random Hermitian of *n*-qubits (a matrix of *n×n*)

```wl
n = 4;
a = QuantumOperator["RandomHermitian", n]["Matrix"];
HermitianMatrixQ[a]
```

Generate a complex-valued random vector *b*:

```wl
b = QuantumOperator["RandomUnitary", n]["Matrix"][[All, 1]];
Norm[b]
```

Create HHL circuit:

```wl
hhl = N@HHL[a, b, 6];
hhl["Diagram", "ShowGateLabels" -> False]
```

Given the result of above circuit, find corresponding amplitudes, and then normalize them to get solution of HHL solver:

```wl
xHHL = hhl[]["Normalize"]["AmplitudesList"] // Chop
```

Find the corresponding solution using classical approach:

```wl
xClassical = Normalize@LinearSolve[a, b]
```

Compare classical result with the quantum one (using the inner product as a measure; the closer to one the better)"

```wl
Abs[Conjugate[xClassical] . xHHL]^2
```

As mentioned, the number of phase-estimation qubits (3rd argument) determines the precision of final outcome. Execute HHL with less phase-estimation qubits::

```wl
xHHLLowPrecision = 
 N[HHL[a, b, 3]][]["Normalize"]["AmplitudesList"] // Chop
```

Quantify the good ness of HHL result in compared with the classical result:

```wl
Abs[Conjugate[xClassical] . xHHLLowPrecision]^2
```

### Quantum error correction with the three-qubit code

A logical qubit can be shielded from bit-flip noise by spreading it across three physical qubits. The encoder $|\psi\rangle|00\rangle \mapsto \alpha|000\rangle+\beta|111\rangle$ copies the computational value with two CNOTs, and a later majority vote, implemented by a single Toffoli, restores the logical qubit as long as at most one of the three qubits flips.

We build the encoder and the decoder. The decoder undoes the two CNOTs and then applies a Toffoli whose controls are the two redundancy qubits, flipping the data qubit back exactly when the majority disagrees with it:

```wl
encode = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}}];
decode = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}, 
    "Toffoli" -> {2, 3, 1}}];
```

Noise enters as a quantum channel. A bit-flip channel `"BitFlip"[p]` applies an X with probability $p$; dropped into a circuit as an ordinary element it turns the pure state into a mixed density matrix. We place an independent copy on every physical qubit:

```wl
noise[p_] := 
 QuantumCircuitOperator[QuantumChannel["BitFlip"[p]] -> # & /@ {1, 2, 3}]
```

We assemble the full pipeline, encode then noise then decode, for an arbitrary input qubit, and trace out the two redundancy qubits to read off the recovered logical state:

```wl
psi = QuantumState[{2, 1}/Sqrt[5]];
protected[p_] := 
 QuantumPartialTrace[
  QuantumCircuitOperator[{psi -> 1, encode, noise[p], decode}][], {2, 3}]
```

Compare the fidelity of the recovered logical qubit with that of the same qubit sent unprotected through one bit-flip channel:

```wl
{QuantumSimilarity[psi, protected[0.1], "Fidelity"], 
 QuantumSimilarity[psi, QuantumChannel["BitFlip"[0.1]][psi], "Fidelity"]}
```

The encoded qubit stays closer to the original. Because the code removes every single-qubit error, the logical infidelity grows as $p^2$ rather than $p$, so encoding wins for every $p<1/2$.

Tabulate both fidelities across the noise range and plot them:

```wl
data = Table[{p, QuantumSimilarity[psi, protected[p], "Fidelity"], 
    QuantumSimilarity[psi, QuantumChannel["BitFlip"[p]][psi], 
     "Fidelity"]}, {p, 0, 0.5, 0.05}];
ListLinePlot[{data[[All, {1, 2}]], data[[All, {1, 3}]]}, 
 PlotLegends -> {"encoded", "unencoded"}, 
 AxesLabel -> {"p", "fidelity"}, PlotMarkers -> Automatic]
```

The encoded curve sits above the unencoded one throughout, the gap widening as the redundancy suppresses the leading-order error.

### Indefinite causal order: the quantum switch

Every circuit so far fixes the order in which its gates act. The quantum switch removes even that assumption: a control qubit coherently decides whether operation $A$ precedes $B$ or the reverse, a resource with indefinite causal order that no fixed circuit can reproduce. `QuantumCircuitOperator` provides it as the named circuit `"Switch"[A, B]`.

Build the switch of two single-qubit operations and inspect its structure:

```wl
QuantumCircuitOperator[
  "Switch"[QuantumOperator["X"], QuantumOperator["Z"]]]["Diagram"]
```

With the control prepared in $|+\rangle$, the switch outputs $\frac{1}{\sqrt{2}}(|0\rangle BA + |1\rangle AB)|\psi\rangle$. Rotating the control back with a Hadamard and measuring it interferes the two orders: the control collapses to $|0\rangle$ when the operations commute ($AB=BA$) and to $|1\rangle$ when they anticommute ($AB=-BA$).

We assemble the full single-query test, preparing the control in $|+\rangle$, applying the switch, rotating back, and measuring the control:

```wl
switchTest[A_, B_] := 
 QuantumCircuitOperator[{"H" -> 1, 
   QuantumCircuitOperator["Switch"[A, B]], "H" -> 1, {1}}]
```

Two operators that commute (here X with itself) send the control to 0 with certainty:

```wl
switchTest[QuantumOperator["X"], QuantumOperator["X"]][]["Probabilities"]
```

Two that anticommute (X and Z) send it to 1 with certainty:

```wl
switchTest[QuantumOperator["X"], QuantumOperator["Z"]][]["Probabilities"]
```

Each run uses $A$ and $B$ exactly once, yet decides commutation versus anticommutation deterministically. No fixed-order circuit can do this with a single use of each operation: the single-query discrimination is a genuine signature of indefinite causal order.

When the operations neither commute nor anticommute (X and H), the control emerges maximally mixed, reflecting the intermediate algebra:

```wl
switchTest[QuantumOperator["X"], QuantumOperator["H"]][]["Probabilities"]
```

## Possible Issues

Using Qiskit as Method, if too many qubits, the result will be an association and not a quantum measurement object:

```wl
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[10]], 
   Range[10]}][Method -> "Qiskit"]
```

## Neat Examples

Design Discussion

Application Notes

XXXX

Test Cases

Function Essay

XXXX
