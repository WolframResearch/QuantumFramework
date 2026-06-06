---
Template: Symbol
Name: QuantumCircuitOperator
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/QuantumCircuitOperator
Keywords: [quantum circuit, gate sequence, circuit operator, quantum gate, circuit composition, Qiskit, QASM]
SeeAlso: [QuantumState, QuantumOperator, QuantumMeasurementOperator]
RelatedGuides: [WolframQuantumComputationFramework]
---

## Usage

<code>[QuantumCircuitOperator]()[{*obj*$_1$, *obj*$_2$, *obj*$_3$, …}]</code> represents a quantum circuit built from a list of quantum objects *obj*$_i$, such as a [QuantumOperator](), a [QuantumChannel](), a [QuantumState]() or a [QuantumMeasurementOperator]().

<code>[QuantumCircuitOperator]()["*name*"]</code> gives the named built-in circuit "*name*", for example `"Bell"` or `"Toffoli"`.

<code>[QuantumCircuitOperator]()[*name*[*args*]]</code> gives a parametrized named circuit, such as `"Fourier"[*n*]` or `"Grover"[*f*]`.

<code>*circ*[*state*]</code> applies the circuit *circ* to a [QuantumState]() *state*.

<code>*circ*[]</code> applies the circuit to the corresponding register (computational-basis) state.

<code>*circ*["*prop*"]</code> gives the property *prop* of the circuit.

## Details & Options

- A `QuantumCircuitOperator` represents a quantum circuit as an ordered composition of quantum objects. Each object can be a [QuantumOperator]() (a gate), a [QuantumChannel]() (a noisy operation), a [QuantumState]() (a state preparation) or a [QuantumMeasurementOperator]() (a measurement).
- Many gates have a shorthand representation, so that `QuantumCircuitOperator[{"CNOT", "H" -> {1, 2}}]` is equivalent to building the same circuit from explicit [QuantumOperator]() objects. Use <code>[QuantumShortcut]()</code> to recover the shorthand of an existing circuit.
- Applying a circuit to a [QuantumState]() evolves the state through the circuit. With no argument, `circ[]` uses the corresponding register state as the input.
- <code>*circ*["*prop*"]</code> gives the property *prop* of the circuit. The following properties are supported:

| Property | Result |
|---|---|
| `"Arity"` | the arity of the circuit |
| `"Dimensions"` | the dimensionality of the circuit |
| `"Depth"` | the depth of the circuit |
| `"Qiskit"` or `"QiskitCircuit"` | the corresponding circuit in Qiskit (needs Python as an external evaluator) |
| `"Operators"` | the list of operators whose composition is the circuit |
| `"GateCount"` or `"OperatorCount"` | the total number of operators |
| `"Diagram"` | the circuit diagram |
| `"Label"` | the label of the circuit |
| `"QASM"` | the circuit as an OpenQASM 3.0 string |
| `"TensorNetwork"` | the tensor network of tensors representing the circuit |
| `"Topology"` | the circuit topology, with qudits/wires as vertices and gates as edges |
| `"CircuitOperator"` | the single operator giving the overall transformation of the circuit |
| `"Width"` | the number of wires in the circuit |
| `"Shift", *n*` | shift the order of the circuit (the qudits it acts on) by *n* |

- `QuantumCircuitOperator` supports a large family of named built-in circuits, given as `QuantumCircuitOperator["*name*"]` or `QuantumCircuitOperator[*name*[*args*]]`, including Bell, Toffoli, Grover (and its many oracles), quantum phase estimation, the quantum Fourier transform, the Bernstein-Vazirani algorithm (and its oracles), and graph circuits, among others:

| Named circuit | Result |
|---|---|
| `"Bell"` | Bell circuit |
| `"Toffoli"` | Toffoli circuit |
| `"Graph"[*g*]` | graph circuit of a given graph *g* |
| `"Grover"[*f*]` | Grover circuit of a Boolean function *f* |
| `"Grover0"[*f*]` | Grover circuit of a Boolean function *f*, with the diffusion part as a controlled-0 gate |
| `"BooleanOracle"[*f*]` | Boolean oracle for a Boolean function *f* |
| `"BooleanOracleR"[*f*]` | Boolean oracle for a Boolean function *f*, decomposed as CNOTs and Z-rotation operators |
| `"PhaseOracle"[*f*]` | phase oracle for a Boolean function *f* |
| `"GroverPhase"[*f*]` | Grover circuit with a phase oracle |
| `"GroverPhase0"[*f*]` | Grover circuit with a phase oracle, and a diffusion part with a controlled-0 gate |
| `"GroverDiffusion"[*n*]` | diffusion part of a Grover oracle acting on *n* qubits |
| `"GroverDiffusion0"[*n*]` | diffusion part of a Grover oracle acting on *n* qubits, with a controlled-0 gate |
| `"GroverPhaseDiffusion"[*n*]` | diffusion part of a Grover oracle acting on *n* qubits, for a phase oracle |
| `"GroverPhaseDiffusion0"[*n*]` | diffusion part of a Grover oracle acting on *n* qubits, for a phase oracle and with a controlled-0 gate |
| `"Fourier"[*n*]` | quantum Fourier transform of *n* qubits |
| `"InverseFourier"[*n*]` | quantum inverse Fourier transform of *n* qubits |
| `"BernsteinVazirani"[*s*]` | Bernstein-Vazirani circuit for a given secret string *s* |
| `"BernsteinVaziraniOracle"[*s*]` | Bernstein-Vazirani oracle for a given secret binary vector or bit-string *s* |
| `"SimonOracle"[*s*]` | Simon oracle for a given secret binary vector or bit-string *s* |
| `"Simon"[*s*]` | Simon circuit for a given secret binary vector or bit-string *s* |
| `"PhaseEstimation"[*u*, *n*]` | quantum phase estimation circuit for a given unitary operator *u* and *n* qubits |
| `"DeutschJozsaBooleanOracle"[*k*, *n*]` | Deutsch-Jozsa Boolean oracle of the $(k-1)^{\text{th}}$ Boolean function in *n* variables |
| `"DeutschJozsa"[*k*, *n*]` | Deutsch-Jozsa circuit with the Boolean oracle of the $(k-1)^{\text{th}}$ Boolean function in *n* variables |
| `"DeutschJozsaPhaseOracle"[*k*, *n*]` | Deutsch-Jozsa phase oracle of the $(k-1)^{\text{th}}$ Boolean function in *n* variables |
| `"DeutschJozsaPhase"[*k*, *n*]` | Deutsch-Jozsa circuit with the phase oracle of the $(k-1)^{\text{th}}$ Boolean function in *n* variables |
| `"Trotterization"[*ops*, *order*, *steps*, *time*]` | Trotter-Suzuki decomposition for a list of operators *ops*, for a given *order* (1, 2, 4, 6, …), number of *steps* and *time* (based on Suzuki's recursion relation) |
| `"PhaseNumber"[*n*, *m*, *boo*]` | a circuit encoding the number *n* within *m* qubits; if *boo* is True, with the leading Hadamards, and if False, without them |
| `"LeggettGarg"` | Leggett-Garg quantum circuit for the Leggett-Garg inequality |

- The `"DeutschJozsa"` family also accepts a [BooleanFunction]() directly, as in `"DeutschJozsa"[BooleanFunction[*k*-1, *n*]]`.
- Several further named circuits are demonstrated below, including `"Cup"` and `"Cap"`, `"GrayOracle"`, `"Multiplexer"`, `"CHSH"`, `"StatePreparation"` and `"QuantumState"` (state-preparation circuits).
- With the `Method` option, a circuit can be evaluated by contracting its tensor network (the default), or by converting it to Qiskit (`Method -> "Qiskit"`), optionally using a provider and backend, or by connecting to Classiq (`Method -> "Classiq"`).

## Basic Examples

Create a quantum circuit composed of only single-qubit gates:

```wl
qc = QuantumCircuitOperator[{QuantumOperator["X"], QuantumOperator["Y", {2}]}];
```

Draw the associated circuit diagram:

```wl
qc["Diagram"]
```

---

A shorthand representation of the circuit:

```wl
QuantumCircuitOperator["XY"]["Diagram"]
```

---

Another shorter form of input that creates the above circuit:

```wl
QuantumShortcut[qc]
```

<!-- => {"X" -> {1}, "Y" -> {2}} -->

```wl
QuantumCircuitOperator[%]["Diagram"]
```

---

Create a circuit with multi-qubit gates:

```wl
qc = QuantumCircuitOperator[{QuantumOperator["CNOT"], QuantumOperator["H" -> {1, 2}]}];
```

The associated circuit diagram:

```wl
qc["Diagram"]
```

---

The above circuit has a shorthand representation, too:

```wl
QuantumCircuitOperator[{"CNOT", "H" -> {1, 2}}] == qc
```

<!-- => True -->

---

The dashed line between the two Hadamards implies that they are defined as a tensor product, although they are separable. Compare it with this case:

```wl
QuantumCircuitOperator[{"CNOT", "H", "H" -> {2}}]["Diagram"]
```

---

A quantum state is transformed by a quantum circuit:

```wl
QuantumCircuitOperator[{"CNOT", "H", "S" -> {2}}][QuantumState["01"]]
```

---

If no input is given, the initial state is set to the corresponding register state:

```wl
QuantumCircuitOperator[{"CNOT", "H", "S" -> {2}}][]
```

---

A quantum circuit can include many quantum objects:

```wl
qc = QuantumCircuitOperator[{QuantumState["Bell", "Label" -> "Bell"], "XY", {"C", "RY"[π/4]}, "BitFlip"[.3], "H" -> 2, "P"[π/3] -> 2, {1}, "M"["X"] -> 2}];
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
QuantumCircuitOperator[{"X", "Y", "Z", "P"[θ], "RX"[Subscript[ϕ, x]], "RY"[Subscript[ϕ, y]], "RZ"[Subscript[ϕ, z]], "H" -> 2, "S" -> 2, "T" -> 2, "SX" -> 2, "NOT" -> 2, "RootNOT" -> 2, "GlobalPhase"[Φ], "PhaseShift"[2] -> 2}]["Diagram"]
```

---

Some multi-qubit quantum gates:

```wl
QuantumCircuitOperator[{"CH", "CS", "CNOT", "CT", "CX", "CY", "CZ", "CP"[θ], {"C", "RX"[Subscript[θ, x]] -> 4, {3}}, {"C", "RY"[Subscript[θ, y]] -> 4, {3}}, {"C", "RZ"[Subscript[θ, z]] -> 4, {3}}, {"R", Subscript[θ, xx], "XX"} -> {3, 4}, {"R", Subscript[θ, yy], "YY"} -> {3, 4}, {"R", Subscript[θ, zz], "ZZ"} -> {3, 4}, "SWAP" -> {3, 4}, "RootSWAP" -> {3, 4}, "Permutation" -> {3, 4} -> {4, 3}, "CSWAP", {"C", "RandomUnitary", {3}, {4}}}]["Diagram"]
```

### Named built-in quantum circuits

`QuantumCircuitOperator` supports a large family of common quantum circuits.

#### Bell, Cup and Cap circuits

The Bell circuit:

```wl
QuantumCircuitOperator["Bell"]["Diagram"]
```

---

The action of the Bell circuit on the register state generates the Bell state:

```wl
QuantumCircuitOperator["Bell"][] == QuantumState["Bell"]
```

<!-- => True -->

---

The Bell circuit is similar to `"Cup"`, up to a $\sqrt{2}$ factor:

```wl
QuantumCircuitOperator["Cup"]["Diagram"]
```

---

Up to a $\sqrt{2}$ factor, `"Cup"` generates a Bell state:

```wl
QuantumCircuitOperator["Cup"][]
```

<!-- => |00⟩ + |11⟩ -->

---

The Cap circuit is the inverse of Cup:

```wl
QuantumCircuitOperator["Cap"]["Diagram"]
```

#### Toffoli circuit

The Toffoli circuit as a combination of one-qubit and two-qubit gates:

```wl
QuantumCircuitOperator["Toffoli"]["Diagram"]
```

---

Its overall operation is the built-in Toffoli operator:

```wl
QuantumCircuitOperator["Toffoli"]["CircuitOperator"] == QuantumOperator["Toffoli"]
```

<!-- => True -->

#### Fourier and inverse Fourier circuit

For the quantum Fourier circuit, use `QuantumCircuitOperator["Fourier"[*n*]]` with *n* the number of qubits. For example, the Fourier circuit of 2 qubits:

```wl
QuantumCircuitOperator["Fourier"]["Diagram"]
```

---

The Fourier circuit transforms the computational basis into the momentum (Pauli-X) basis. Check that $|00\rangle$ is transformed into $|x_+ x_+\rangle$:

```wl
QuantumState[QuantumCircuitOperator["Fourier"][], "XX"]
```

<!-- => |x_+ x_+⟩ -->

---

The Fourier circuit of 5 qubits:

```wl
QuantumCircuitOperator["Fourier"[5]]["Diagram"]
```

---

The inverse Fourier circuit of 5 qubits:

```wl
QuantumCircuitOperator["InverseFourier"[5]]["Diagram"]
```

---

The inverse Fourier circuit transforms from the momentum (Pauli-X) basis into the computational basis. Check that $|x_+ x_+\rangle$ is transformed into $|00\rangle$:

```wl
QuantumCircuitOperator["InverseFourier"][QuantumState["++"]]
```

<!-- => |00⟩ -->

---

The inverse Fourier circuit is the same as the conjugate transpose of the Fourier circuit:

```wl
QuantumCircuitOperator["InverseFourier"] == QuantumCircuitOperator["Fourier"]["Dagger"]
```

<!-- => True -->

---

The QFT is the same as a basis transformation into the Fourier basis:

```wl
n = 3;
SparseArray[QuantumCircuitOperator["Fourier"[n]][#]["StateVector"] & /@ QuantumBasis[2, n]["BasisStates"]] == QuantumBasis["Fourier"[2^n]]["Elements"]
```

<!-- => True -->

#### Phase estimation circuit

The quantum phase estimation circuit `QuantumCircuitOperator["PhaseEstimation"[*u*, *n*]]` takes two arguments: a unitary operator *u* and an integer *n*. The integer *n* specifies the number of qubits and `controlled-`*u*$^j$ operators in the circuit, with $j = 0, 1, \dots, n-1$. The accuracy of phase estimation and the success probability depend on *n*.

Generate a random unitary operator:

```wl
u = QuantumOperator["RandomUnitary"];
```

---

Create the corresponding QPE circuit to find one of its eigenvalues:

```wl
qc = QuantumCircuitOperator["PhaseEstimation"[u]];
qc["Diagram"]
```

---

Calculate the corresponding measurement probabilities:

```wl
N[qc][]["ProbabilityPlot", "LabelsAngle" -> π/2, AspectRatio -> 1/4]
```

---

Given the result with maximum probability, estimate the eigenvalue:

```wl
Keys[TakeLargest[N[qc][]["Probabilities"], 1]][[1]]["Name"]
est1 = N[FromDigits[%, 2]/2^4]
```

<!-- => {0, 1, 1, 1} -->
<!-- => 0.4375 -->

---

Generate another QPE circuit to obtain the other eigenvalue (the only difference is dropping the Pauli-X gate):

```wl
qc2 = qc[[2 ;;]];
qc2["Flatten"]["Diagram"]
```

---

Calculate the corresponding measurement probabilities:

```wl
N[qc2][]["ProbabilityPlot", "LabelsAngle" -> π/2, AspectRatio -> 1/4]
```

---

Given the result with maximum probability, estimate the eigenvalue:

```wl
Keys[TakeLargest[N[qc2][]["Probabilities"], 1]][[1]]["Name"]
est2 = N[FromDigits[%, 2]/2^4]
```

<!-- => {0, 0, 1, 1} -->
<!-- => 0.1875 -->

---

Find the corresponding eigenvalues (note they should be of the form $e^{i 2\pi\lambda}$):

```wl
λ = Mod[Arg[u["Eigenvalues"]]/(2 π), 1]
```

---

Calculate the percentage error of the estimated eigenvalues against the expected ones:

```wl
100. (1 - Sort[{est1, est2}]/Sort[λ])
```

#### Graph circuit

Given a graph, the corresponding circuit (which generates quantum graph/cluster states) is a series of controlled-Z operators on the edges of the graph:

```wl
g = Graph[RandomGraph[{4, 5}], VertexLabels -> "Name"]
```

```wl
QuantumCircuitOperator["Graph"[g]]["Diagram"]
```

---

The graph circuit generates the graph state (also called the cluster state):

```wl
QuantumCircuitOperator["Graph"[g]][] == QuantumState["Graph"[g]]
```

<!-- => True -->

---

Check that qubits on the same edge are entangled:

```wl
With[{st = QuantumState["Graph"[g]]}, QuantumEntangledQ[st, #] & @* List @@@ EdgeList[g]]
```

<!-- => {True, True, True, True, True} -->

#### Bernstein-Vazirani circuit

For the Bernstein-Vazirani algorithm, specify a secret bit-string *s*: the oracle is `QuantumCircuitOperator["BernsteinVaziraniOracle"[*s*]]` and the whole circuit is `QuantumCircuitOperator["BernsteinVazirani"[*s*]]`. For example, the Bernstein-Vazirani oracle for the secret bit-string 101:

```wl
QuantumCircuitOperator["BernsteinVaziraniOracle"["101"]]["Diagram"]
```

---

The Bernstein-Vazirani circuit for the secret bit-string 101:

```wl
QuantumCircuitOperator["BernsteinVazirani"["101"]]["Diagram"]
```

---

Measurement results of the Bernstein-Vazirani circuit for the secret bit-string 101:

```wl
QuantumCircuitOperator["BernsteinVazirani"["101"]][]["ProbabilitiesPlot"]
```

---

The Bernstein-Vazirani oracle can take a binary vector or a bit-string:

```wl
QuantumCircuitOperator["BernsteinVaziraniOracle"[{1, 0, 1}]] == QuantumCircuitOperator["BernsteinVaziraniOracle"["101"]]
```

<!-- => True -->

#### Simon circuit

Simon's algorithm is a quantum algorithm designed to solve Simon's problem efficiently. It aims to find a hidden bit-string by querying a black-box function.

Create the corresponding Simon oracle:

```wl
so = QuantumCircuitOperator["SimonOracle"[{1, 0, 1, 1}]];
so["Diagram"]
```

---

A function that returns the Boolean table of a Simon oracle:

```wl
boolTable[oracle_] := With[{n = oracle["Width"]/2}, KeyValueMap[TraditionalForm /@ {#1, #2} & ][AssociationMap[QuantumOperator["Discard"[2] -> Range[n]][oracle[#1]] & , QuantumBasis[2, n]["BasisStates"]]]]
```

---

Return the Boolean table of the above Simon oracle:

```wl
TableForm[boolTable[so], TableDepth -> 2, TableHeadings -> {None, {"bit", "oracle[bit]"}}]
```

---

The Simon circuit for a given hidden sequence:

```wl
SeedRandom[0];
simon = QuantumCircuitOperator["Simon"[seq = {1, 0, 0, 1}]];
simon["Diagram"]
```

---

The measurement probabilities given the Simon circuit for the binary string 1001:

```wl
simonMeas = N[simon][]
```

---

Plot the probabilities:

```wl
simonMeas["ProbabilitiesPlot", AspectRatio -> 1/4, "LabelsAngle" -> π/2]
```

---

The list of all potential results:

```wl
m = Normal /@ Keys @ Select[simonMeas["Probabilities"], Positive]
```

---

Find the null vector ($\text{matrix} \cdot m = 0 \bmod 2$):

```wl
NullSpace[m, Modulus -> 2]
```

<!-- => {{1, 0, 0, 1}} -->

---

The secret bit-string is recovered:

```wl
seq == %[[1]]
```

<!-- => True -->

---

The Simon circuit and its corresponding Boolean function are not uniquely determined. For each evaluation, one gets a different oracle that can solve Simon's problem:

```wl
simon2 = QuantumCircuitOperator["Simon"[seq]];
simon == simon2
```

<!-- => False -->

---

With the new circuit, one can still find the hidden sequence:

```wl
NullSpace[Normal /@ Keys @ Select[simon2[]["Probabilities"], Positive], Modulus -> 2][[1]] == seq
```

<!-- => True -->

#### Gray oracle

`"GrayOracle"` takes a list of angles or a function, a number of ancillary qubits and a target qubit, and generates the corresponding oracle based on [Gray code](https://arxiv.org/pdf/quant-ph/0404089.pdf).

Given 4 qubits as ancillary, generate angles for all possible bit-strings, prepending a zero angle:

```wl
prec = 4;
angles = Prepend[0][N[((1/2^prec)*(1/FromDigits[#1, 2]) & ) /@ Rest[Tuples[{0, 1}, prec]]]];
```

---

Generate the corresponding Gray oracle:

```wl
QuantumCircuitOperator["GrayOracle"[angles, prec, 1]]["Diagram", "ShowGateLabels" -> False]
```

---

It is the same as the following oracle, using a function as input:

```wl
QuantumCircuitOperator["GrayOracle"[If[#1 == 0, 0, 1./#1] & , prec, 1]]["Diagram", "ShowGateLabels" -> False]
```

---

The Gray oracle is the same as a Multiplexer with angles:

```wl
inv = QuantumCircuitOperator["Multiplexer" @@ "RY" /@ angles, RotateLeft[Range[prec + 1]]];
inv["Diagram", "ShowGateLabels" -> False]
```

---

Check that the Multiplexer is the same as the Gray oracle:

```wl
(inv["Matrix"] - QuantumCircuitOperator["GrayOracle"[If[#1 == 0, 0, 1./#1] & , prec, 1]]["Matrix"])["ExplicitValues"] // Chop
```

<!-- => {0, 0, …, 0} (64 zeros) -->

#### Grover circuit (and family of Boolean oracles)

There is a collection of named circuits for the Grover search algorithm and its oracles (for example Boolean or phase oracles) and amplification. The class of Grover-related named circuits is `QuantumCircuitOperator[*name*[*boo*]]`, with *boo* a Boolean function.

The Boolean oracle for a given Boolean function:

```wl
bo = QuantumCircuitOperator["BooleanOracle"[a && ! b]];
bo["Diagram"]
```

---

A Boolean oracle flips the ancillary qubit (here the 3rd qubit) only if the quantum state is a solution of the Boolean function. For example, `{True, False}` (10) is a solution of `a && ! b`, so the ancillary qubit is flipped:

```wl
bo[QuantumState["100"]]["Formula"]
```

<!-- => |101⟩ -->

---

One can also get the decomposition of a Boolean oracle in terms of CNOTs and Z-rotations:

```wl
QuantumCircuitOperator["BooleanOracleR"[a && ! b]]["Diagram"]
```

---

For the diffusion part of the Grover circuit (also called the amplification part), one only needs to specify the number of qubits:

```wl
QuantumCircuitOperator["GroverDiffusion"[3]]["Diagram"]
```

```wl
QuantumCircuitOperator["GroverDiffusion0"[3]]["Diagram"]
```

---

Another way to store the solution of a Boolean function is to save it as a phase rather than on an ancillary qubit. This is achieved with a phase oracle, which adds an overall phase of $\pi$ if the quantum state is a solution of the Boolean function:

```wl
po = QuantumCircuitOperator["PhaseOracle"[a && ! b]];
po["Diagram"]
```

```wl
po[QuantumState["10"]]["Formula"]
```

<!-- => -|10⟩ -->

---

The corresponding diffusion part, given a phase oracle:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion"[2]]["Diagram"]
```

---

The overall Grover circuit, using controlled-1 gates and a Boolean oracle:

```wl
QuantumCircuitOperator["Grover"[a && ! b]]["Diagram"]
```

---

The overall Grover circuit, using controlled-0 gates and a Boolean oracle:

```wl
QuantumCircuitOperator["Grover0"[a && ! b]]["Diagram"]
```

---

The overall Grover circuit, using controlled-1 gates and a phase oracle:

```wl
QuantumCircuitOperator["GroverPhase"[a && ! b]]["Diagram"]
```

---

The overall Grover circuit, using controlled-0 gates and a phase oracle:

```wl
QuantumCircuitOperator["GroverPhase0"[a && ! b]]["Diagram"]
```

#### CHSH circuit

The CHSH quantum circuit:

```wl
QuantumCircuitOperator["CHSH"]["Diagram", "ShowMeasurementWire" -> False, "WireLabels" -> {Placed["a", Right], Placed["x", Right], Placed["y", Right], Placed["b", Right]}, Epilog -> {Opacity[0.3], Text[Style["a⊻b\!\(\*OverscriptBox[\(=\), \(?\)]\)x∧y", GrayLevel[0], 20, Bold], {9., -2.35}]}, PlotRange -> {{0, 10}, {-5, 0}}, ImageSize -> Large]
```

---

Calculate the probability of winning:

```wl
NProbability[BitAnd[x, y] == BitXor[a, b], Distributed[{a, x, y, b}, QuantumCircuitOperator["CHSH"][]["MultivariateDistribution"]]]
```

<!-- => 0.853553390593274 -->

#### Trotter-Suzuki

One can find a Trotter-Suzuki decomposition to approximate a time evolution.

First order only:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator["X", {1, 2}], QuantumOperator["X", {2, 3}]}, 1, 1, 1}];
qc["Diagram"]
```

---

Second order:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator["X", {1, 2}], QuantumOperator["X", {2, 3}]}, 2, 1, 1}];
qc["Diagram"]
```

---

Second order, with 4 steps:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator["X", {1, 2}], QuantumOperator["X", {2, 3}]}, 2, 4, 1}];
qc["Diagram", FontSize -> 8]
```

---

Fourth order, only 1 step:

```wl
qc = QuantumCircuitOperator[{"Trotterization", {QuantumOperator["X", {1, 2}], QuantumOperator["X", {2, 3}]}, 4, 1, 1}];
qc["Diagram", "ShowGateLabels" -> False]
```

#### PhaseNumber circuit

The PhaseNumber circuit encodes a number in *n* qubits.

Encode 5 in 4 qubits:

```wl
num = 5;
QuantumCircuitOperator["PhaseNumber"[num, 4]]["Diagram"]
```

---

Remove the Hadamards:

```wl
QuantumCircuitOperator["PhaseNumber"[num, 4, False]]["Diagram"]
```

---

PhaseNumber followed by InverseFourier gives a circuit that can be used to recover the encoded number:

```wl
num = 5;
qc = QuantumCircuitOperator[{"PhaseNumber"[num, 6], "InverseFourier"[6], Range[6]}];
qc["Diagram"]
```

---

The corresponding measurement probabilities:

```wl
m = N[qc][]
```

---

Given the result with maximum probability, find the corresponding integer from the bit-string:

```wl
FromDigits[Keys[TakeLargest[m["Probabilities"], 1]][[1]]["Name"], 2]
```

<!-- => 5 -->

---

Check that it is the same as the initially encoded number:

```wl
% == num
```

<!-- => True -->

---

The PhaseNumber circuits follow the usual addition:

```wl
qc1 = QuantumCircuitOperator[{"PhaseNumber", 5, 6}];
qc2 = QuantumCircuitOperator[{"PhaseNumber", 3, 6, False}];
(qc1 /* qc2)["Diagram"]
```

---

Test that $5 + 3 = 8$ in terms of their circuits:

```wl
(qc1 /* qc2)["CircuitOperator"] == QuantumCircuitOperator[{"PhaseNumber", 5 + 3, 6}]["CircuitOperator"]
```

<!-- => True -->

#### Deutsch-Jozsa circuit

For the Deutsch-Jozsa algorithm, the goal is to determine whether a Boolean function is constant (producing the same output regardless of the input, either 0 or 1) or balanced (yielding an equal number of 0s and 1s).

Define a balanced Boolean function:

```wl
bf = BooleanFunction[2, 2];
BooleanConvert[bf][\[FormalX], \[FormalY]]
```

<!-- => ! x && y -->

---

Show its Boolean table (note the same number of True and False):

```wl
BooleanTable[bf]
```

<!-- => {False, False, True, False} -->

```wl
QuantumCircuitOperator["DeutschJozsa"[bf]]["Diagram"]
```

---

One can also get the Deutsch-Jozsa circuit with a phase oracle:

```wl
QuantumCircuitOperator["DeutschJozsaPhase"[bf]]["Diagram"]
```

---

`"DeutschJozsaPhase"` can also take two integers, the first being $k+1$ and the second $n$, where the Boolean function is defined as `BooleanFunction[k, n]`:

```wl
QuantumCircuitOperator["DeutschJozsaPhase"[3, 2]]["Diagram"]
```

---

Given all Boolean functions with two variables, show the number of True and False in their solutions and the corresponding Deutsch-Jozsa probabilities:

```wl
Grid[Prepend[Table[With[{bf = BooleanFunction[f, 2]}, {BooleanConvert[bf][\[FormalX], \[FormalY]], Count[True][BooleanTable[bf]], Count[False][BooleanTable[bf]], With[{rs = QuantumCircuitOperator["DeutschJozsa"[bf]][]["Probabilities"]}, BarChart[rs, ChartLabels -> Keys[rs], AspectRatio -> 1/3, ImageSize -> Small]]}], {f, 0, 15}], Style[#1, Bold] & /@ {"Boolean function", "True", "False", "Deutsch-Jozsa probabilities"}], Frame -> All, Alignment -> {{Left, Center, Center, Center}}]
```

#### Leggett-Garg

The Leggett-Garg circuit:

```wl
QuantumCircuitOperator["LeggettGarg"]["Diagram"]
```

#### State preparation

Generate a random pure state of 3 qubits:

```wl
state = QuantumState["RandomPure"[3]]
```

---

Generate a circuit whose outcome is the above state:

```wl
qc = QuantumCircuitOperator["QuantumState"[state]];
qc["Diagram", "ShowGateLabels" -> False]
```

---

Check that executing the circuit gives the expected state:

```wl
qc[] == state
```

<!-- => True -->

---

One can also use `"StatePreparation"`:

```wl
QuantumCircuitOperator["StatePreparation"[state]]["Diagram", "ShowGateLabels" -> False]
```

### Quantum circuit diagram

One can customize many features of a circuit diagram. See the [Circuit Diagram](paclet:Wolfram/QuantumFramework/tutorial/CircuitDiagram) documentation page for more details.

```wl
QuantumCircuitOperator[{"H", "CY", "H", {1}}]["Diagram", "WireLabels" -> {{Placed["1st qubit", Left], Placed["Eigenvector of Y", Right]}, {Placed["2nd qubit", Left]}}, "MeasurementWireLabel" -> "Measurement Wire"]
```

---

The overall size of the circuit can be adjusted:

```wl
superdense = QuantumCircuitOperator[{"H", "CNOT", "CX" -> {4, 1}, "CZ" -> {3, 1}, "CNOT", "H", {1}, {2}}];
superdense["Diagram", "Size" -> .5]
```

---

One can customize the horizontal gap:

```wl
superdense["Diagram", "HorizontalGapSize" -> 3]
```

---

Customize the vertical gap:

```wl
superdense["Diagram", "VerticalGapSize" -> 2]
```

---

Show the connectors as dots, specifying explicitly which wires a gate acts on:

```wl
QuantumCircuitOperator[{QuantumOperator[{"R", θ, "XX"}, {1, 3}]}]["Diagram", "ShowConnectors" -> True]
```

---

Put a portion of a circuit in a gray box to emphasize it:

```wl
QuantumCircuitOperator[{"CNOT", QuantumCircuitOperator[{QuantumOperator["SWAP", {2, 4}], QuantumOperator[{"Phase", γ}]}, "MyCir"]}]["Diagram"]
```

---

Drop the gate labels (hovering over the diagram still shows the labels):

```wl
QuantumCircuitOperator[{"H", "CY", "H", {1}}]["Diagram", "ShowGateLabels" -> False]
```

### Tensor network of a quantum circuit

The corresponding tensor network of a quantum circuit can be generated as a property. See the [Tensor Network](paclet:Wolfram/QuantumFramework/tutorial/TensorNetwork) documentation page for more details.

```wl
tn = QuantumCircuitOperator["Grover"[a && ! b || c]]["TensorNetwork"]
```

---

The corresponding index graph:

```wl
TensorNetworkIndexGraph[tn, GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}]
```

### Parametrized quantum circuit

Define a parametrized quantum circuit with two parameters:

```wl
qc = QuantumCircuitOperator[{"RY"[θ], "RZ"[\[Phi]]}, "Parameters" -> {θ, \[Phi]}]
```

---

Assign numerical values to the parameters:

```wl
qc[<|θ -> π/3, \[Phi] -> π/4|>]["Diagram"]
```

### Generalizations and extensions

`QuantumCircuitOperator` also supports controlled operators with more than one control or target:

```wl
QuantumCircuitOperator[QuantumOperator[{"C", "X", {1, 4}, {2}}, {3}]]["Diagram"]
```

---

Decompose the Toffoli gate:

```wl
v = QuantumOperator[MatrixPower[PauliMatrix[1], 1/2], "Label" -> "V"];
```

```wl
qc = QuantumCircuitOperator[{QuantumOperator[{"Controlled", v}, {2, 3}], QuantumOperator["CX"], QuantumOperator[{"Controlled", v["Dagger"]}, {2, 3}], QuantumOperator["CX"], QuantumOperator[{"Controlled", v}, {1, 3}]}];
qc["Diagram"]
```

---

Show that the above implementation is the same as the built-in Toffoli gate:

```wl
QuantumOperator["Toffoli"] == qc["CircuitOperator"]
```

<!-- => True -->

---

Quantum circuits consisting of qudit gates are supported. Define a qudit gate of dimensionality 3:

```wl
QuantumCircuitOperator[{"RandomUnitary"[3], "Spider"[QuantumBasis[{3, 2}, {3, 2}]] -> {1, 2} -> {2, 1}}]["Diagram", "ShowWireDimensions" -> True]
```

---

One can shift a circuit, which means changing its overall order:

```wl
qc = QuantumCircuitOperator[{"Fourier", 3}];
qc["Diagram"]
```

---

Shift the above circuit by 2 qudits:

```wl
qc["Shift", 2]["Diagram"]
```

## Options

### Method

The default method creates the tensor network and then contracts it to get the results:

```wl
QuantumCircuitOperator[{"Bell", {1, 2}}][]
```

---

Convert the circuit into Qiskit and simulate it using the Qiskit default backend, such as Aer (the results are quantum objects in the Wolfram quantum framework):

```wl
#| eval: false
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], Range[5]}][Method -> "Qiskit"]
```

---

IBMProvider as the Qiskit provider:

```wl
#| eval: false
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], Range[5]}][Method -> {"Qiskit", "Provider" -> "IBMProvider"}]
```

---

IBMProvider as the Qiskit provider on a specific *backend*:

```wl
#| eval: false
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], Range[5]}][Method -> {"Qiskit", "Provider" -> "IBMProvider", "Backend" -> backend}]
```

---

A delay can be added into a Qiskit circuit, too:

```wl
#| eval: false
QuantumCircuitOperator[{"I" -> "Delay"[2], "M"}]["Qiskit"]
```

---

Return the list of IBM-quantum backends and their queue:

```wl
#| eval: false
ServiceConnect["IBMQ", "New"]["BackendQueue"]
```

---

AWSBraket as the Qiskit provider on a specific *backend*:

```wl
#| eval: false
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[5]], Range[5]}][Method -> {"Qiskit", "Provider" -> "AWSBraket", "Backend" -> backend}]
```

---

Return the list of AWSBraket backends:

```wl
#| eval: false
aws = ServiceConnect["AWS"];
braket = ServiceExecute[aws, "GetService", {"Name" -> "Braket"}];
devices = braket["SearchDevices", "Filters" -> {}]["Devices", All]
```

---

For some circuit creation, users can connect to Classiq:

```wl
#| eval: false
qc = QuantumCircuitOperator[QuantumState[RandomReal[{0, 1}, 2^3]]["Normalized"], Method -> "Classiq"];
qc["Diagram", "ShowGateLabels" -> False]
```

## Applications

### POVM in a quantum circuit

Collective measurements as well as POVM measurements are supported. Define a POVM measurement in a minimal informationally complete basis:

```wl
povm = QuantumMeasurementOperator["WignerMICPOVM"]
```

---

Define a circuit consisting of the above measurement:

```wl
qc = QuantumCircuitOperator[{QuantumState["RandomPure", "Label" -> "RandomState"], povm}];
qc["Diagram"]
```

---

Show the outcome probabilities:

```wl
qc[]["Probabilities"]
```

### Quantum teleportation

Implement a quantum circuit for teleporting a qubit. The 1st and 2nd qubits represent Alice's system, while the last is Bob's. The double lines carry classical bits. The state to be teleported is $|\psi_1\rangle = \alpha|0\rangle + \beta|1\rangle$, where $\alpha$ and $\beta$ are unknown amplitudes:

```wl
circuit = QuantumCircuitOperator[{QuantumCircuitOperator[{QuantumState[{α, β}], "Cup" -> {2, 3}}, "State Prep."], "CX", "H", "CX" -> {2, 3}, "CZ" -> {1, 3}, {1}, {2}}];
```

```wl
circuit["Diagram"]
```

---

Trace out the 1st and 2nd qubits over the possible post-measurement states, and compare the reduced state of qubit 3:

```wl
FullSimplify[Thread[QuantumState[{α, β}] == QuantumOperator["Trace" -> {1, 2}] /@ circuit[]["States"]]]
```

<!-- => {{True, True, True, True}, {True, True, True, True}, {True, True, True, True}, {True, True, True, True}} -->

### Quantum key distribution as in the BB84 protocol

As an example of quantum key distribution (QKD), consider the BB84 protocol. Alice generates random bits and, for each, randomly chooses either a Z-basis or X-basis of measurement. Consider a random sequence of 9 bits and corresponding bases:

```wl
aliceBits = RandomChoice[{0, 1}, 9]
```

```wl
aliceBases = RandomChoice[{"Z", "X"}, 9]
```

---

Depending on her choice of bit and basis, she follows this protocol for sending a qubit to Bob:

```wl
aliceProtocol = Association[{0, "Z"} -> QuantumState[{1, 0}, "Z"], {0, "X"} -> QuantumState[{1, 0}, "X"], {1, "Z"} -> QuantumState[{0, 1}, "Z"], {1, "X"} -> QuantumState[{0, 1}, "X"]];
Column[Normal[(#1["Formula"] & ) /@ aliceProtocol]]
```

---

With Alice's bits and bases generated randomly, she will send Bob these qubits:

```wl
aliceQubits = aliceProtocol /@ Thread[{aliceBits, aliceBases}];
(#1["Formula"] & ) /@ aliceQubits
```

---

Bob randomly generates a sequence of either Z-bases or X-bases of measurement:

```wl
bobBases = RandomChoice[{"Z", "X"}, 9]
```

---

Thus, Bob applies the following measurement operators on the qubits he got from Alice:

```wl
ops = QuantumMeasurementOperator /@ bobBases;
```

---

Post-measurement, whenever Bob's measurement basis is not the same as Alice's, he gets a random result:

```wl
mea = MapThread[#1[#2] & , {ops, aliceQubits}];
```

---

Return one possible realization of Bob's series of measurements (some results are consistent and others random, depending on Alice's qubit and Bob's measurement):

```wl
bobRealization = (#1["SimulatedMeasurement"] & ) /@ mea
```

---

Bob assigns the bit 0 whenever he gets $|z_+\rangle$ or $|x_+\rangle$, and the bit 1 for $|z_-\rangle$ or $|x_-\rangle$:

```wl
bobBits = bobRealization /. Flatten[(Alternatives @@ #1[[1 ;; 2]] -> #1[[-1]] & ) /@ Transpose[{QuditBasis["Z"]["Names"], QuditBasis["X"]["Names"], {0, 1}}]]
```

---

Finally, Alice and Bob publicly share the measurement bases they used (a classical communication). If a basis is the same, they each keep the bit; if not, they each discard the corresponding bit. Each does so separately, and each obtains a secret key that matches the other's:

```wl
aliceSecretKey = Select[Transpose[{aliceBits, aliceBases, bobBases}], #1[[2]] == #1[[3]] & ][[All, 1]]
```

```wl
bobSecretKey = Select[Transpose[{bobBits, aliceBases, bobBases}], #1[[2]] == #1[[3]] & ][[All, 1]]
```

```wl
bobSecretKey == aliceSecretKey
```

<!-- => True -->

---

The key length is less than the number of qubits Alice sent to Bob, because Bob picks his measurement basis randomly. Alice can keep sending qubits until they reach a desired length:

```wl
add[{x_, y_}] := Module[{aliceBit, aliceBase, aliceQubit, bobBase, bobBit}, aliceBit = RandomChoice[{0, 1}]; aliceBase = RandomChoice[{"Z", "X"}]; aliceQubit = aliceProtocol[{aliceBit, aliceBase}]; bobBase = RandomChoice[{"Z", "X"}]; bobBit = QuantumMeasurementOperator[bobBase][aliceQubit]["SimulatedMeasurement"] /. Flatten[(Alternatives @@ #1[[1 ;; 2]] -> #1[[-1]] & ) /@ Transpose[{QuditBasis["Z"]["Names"], QuditBasis["X"]["Names"], {0, 1}}]]; If[aliceBase == bobBase, {Append[x, aliceBit], Append[y, bobBit]}, {x, y}]]
```

---

For example, set the key length to 5 bits:

```wl
NestWhile[add, {{}, {}}, Length[#1[[1]]] < 5 & ]
```

### Parity test

A quantum circuit for a parity test. If the number of 1s in a bit-string (for example 0001101…) is even, the parity is 0 (even parity); if odd, the parity is 1 (odd parity). A quantum circuit to measure the parity of an *n*-qubit state:

```wl
parityMeasurementCircuit[n_] := QuantumCircuitOperator[Append[("CNOT" -> {#1, n + 1} & ) /@ Range[n], {n + 1}]]
```

---

Example of a parity measurement circuit for a 3-qubit state:

```wl
parityMeasurementCircuit[3]["Diagram"]
```

---

Apply this circuit to $|000\rangle \otimes |0\rangle$ (even parity) and $|001\rangle \otimes |0\rangle$ (odd parity):

```wl
{parityMeasurementCircuit[3][QuantumState["0000"]]["ProbabilityPlot"], parityMeasurementCircuit[3][QuantumState["0010"]]["ProbabilityPlot"]}
```

### Quantum Approximate Optimization Algorithm (QAOA)

One can find a binary vector $x$ that minimizes (or maximizes) the quadratic objective $x \cdot q \cdot x + c \cdot x$, with no constraints. This is usually called Quadratic Unconstrained Binary Optimization (QUBO). The linear term $c$ can be dropped without loss of generality.

Define a $4 \times 4$ real-valued upper-triangular matrix as the objective matrix:

```wl
q = UpperTriangularize[RandomReal[{-5, 5}, {4, 4}]];
MatrixForm[q]
```

---

Generate the corresponding objective function:

```wl
Expand[Block[{x = Array[Indexed[\[FormalX], #1] & , Length[q]]}, x . q . x]]
```

---

Maximize and minimize the objective function (classical approach):

```wl
solMax = NArgMax[{x . q . x, 0 \[VectorLess] x^2 \[VectorLessEqual] 1}, x \[Element] Vectors[Length[q], Integers]]
solMin = NArgMin[{x . q . x, 0 \[VectorLess] x^2 \[VectorLessEqual] 1}, x \[Element] Vectors[Length[q], Integers]]
```

---

Define a quantum circuit for any objective function:

```wl
ClearAll[γ, β, qaoa]
qaoa[(obj_)?SquareMatrixQ, p_:1] := With[{l = Length[obj]}, QuantumCircuitOperator[Join[Table["+" -> i, {i, l}], Catenate[Table[{"Barrier", Catenate[Table[If[obj[[i, j]] == 0, Nothing, "R"[2*γ[m]*obj[[i, j]], "ZZ"] -> {i, j}], {i, l}, {j, i + 1, l}]], "Barrier", Table["R"[2*β[m], "X"] -> i, {i, l}]}, {m, p}]]]]["Flatten"]]
```

---

Generate the circuit for the objective matrix, with only two layers:

```wl
p = 2;
qaoa[q, p]["Diagram", FontSize -> 9, "ShowGateLabels" -> False]
```

---

Define the variables (angles) to be optimized:

```wl
var = Join[Table[γ[i], {i, p}], Table[β[i], {i, p}]]
```

---

Find the corresponding probabilities in the computational basis:

```wl
weight = Abs[qaoa[q, p]["Amplitudes"]]^2;
```

---

Find $\langle\gamma,\beta|H|\gamma,\beta\rangle$ with $H = \sum_{i>j} q_{ij}\,\sigma_{z,i}\,\sigma_{z,j}$ and $|\gamma,\beta\rangle$ the output of the above circuit:

```wl
weightedEnergy = Association[KeyValueMap[#1 -> #2*Total[UpperTriangularize[q, 1]*UpperTriangularize[KroneckerProduct[(-1)^#1["Name"], (-1)^#1["Name"]], 1], 2] & , weight]];
```

---

Maximize and minimize the total weighted-mean value:

```wl
optMax = Thread[var -> NArgMax[{Total[weightedEnergy], 0 \[VectorLessEqual] var \[VectorLessEqual] 2 π}, var]]
optMin = Thread[var -> NArgMin[{Total[weightedEnergy], 0 \[VectorLessEqual] var \[VectorLessEqual] 2 π}, var]]
```

---

Plot the results given the optimal values:

```wl
Row[{With[{data = N[weightedEnergy /. optMax]}, BarChart[data, ChartLabels -> (Rotate[If[(solMax /. -1 -> 0) == #1["Name"], Style[#1["Name"], Red], #1["Name"]], Pi/2] & ) /@ Keys[data], PlotLabel -> "Maximizing objective function", ImageSize -> 200]], With[{data = N[weightedEnergy /. optMin]}, BarChart[data, ChartLabels -> (Rotate[If[(solMin /. -1 -> 0) == #1["Name"], Style[#1["Name"], Red], #1["Name"]], Pi/2] & ) /@ Keys[data], PlotLabel -> "Minimizing objective function", ImageSize -> 200]]}]
```

Two solutions arise because the negation of a solution is also a max/min ($(-x) \cdot q \cdot (-x) = x \cdot q \cdot x$).

---

Check that the classical solution matches the quantum one:

```wl
MemberQ[(#1["Name"] & ) /@ Keys[TakeLargestBy[Normal[N[weightedEnergy /. optMax]], Last, 2]], solMax /. -1 -> 0]
```

<!-- => True -->

```wl
MemberQ[(#1["Name"] & ) /@ Keys[TakeSmallestBy[Normal[N[weightedEnergy /. optMin]], Last, 2]], solMin /. -1 -> 0]
```

<!-- => True -->

### Harrow-Hassidim-Lloyd (HHL) algorithm

Given a matrix $A$ and a vector $b$, solve $A \cdot x = b$. Generate an HHL circuit that takes *A*, *b*, and a number of ancillary qubits:

```wl
ClearAll[HHL]
HHL[(A_)?SquareMatrixQ, (b_)?VectorQ, prec : _Integer?Positive : 4] := Block[{d, qpe, t0, yTilde, λTilde, θλ, γ}, d = Log2[Length[b]]; t0 = Floor[Pi]/Norm[A]; yTilde = (If[0 <= #1 <= 2^(prec - 1), #1, #1 - 2^prec] & ) /@ Range[0, 2^prec - 1]; λTilde = ((2*Pi)/t0)*(yTilde/2^prec); γ = Min[DeleteCases[Abs[λTilde], 0. | 0]]; θλ = (If[#1 == 0 || #1 == 0., 0., 2*ArcSin[γ/#1]] & ) /@ λTilde; QuantumCircuitOperator[{QuantumState[Normalize[b]] -> 1 + prec + Range[d], "Barrier", qpe = QuantumCircuitOperator["PhaseEstimation"[QuantumOperator[MatrixExp[t0*I*A]], prec], 1 + Range[prec + d]][[d + 1 ;; -prec - 1]], "Barrier", QuantumCircuitOperator["GrayOracle"[θλ, prec, 1], "Eigenvalue Inversion"], "Barrier", qpe["Dagger"], "Barrier", QuantumOperator[QuantumState[StringJoin["1", StringJoin[ConstantArray["0", prec]]]]["Dagger"]]}] /; IntegerQ[d]]
```

---

Generate a complex-valued random Hermitian matrix of *n* qubits (an $n \times n$ matrix):

```wl
n = 4;
a = QuantumOperator["RandomHermitian", n]["Matrix"];
HermitianMatrixQ[a]
```

<!-- => True -->

---

Generate a complex-valued random vector *b*:

```wl
b = QuantumOperator["RandomUnitary", n]["Matrix"][[All, 1]];
Norm[b]
```

<!-- => 1. -->

---

Create the HHL circuit:

```wl
hhl = N[HHL[a, b, 6]];
hhl["Diagram", "ShowGateLabels" -> False]
```

---

Given the result of the above circuit, find the corresponding amplitudes, then normalize them to get the solution of the HHL solver:

```wl
xHHL = Chop[hhl[]["Normalize"]["AmplitudesList"]]
```

---

Find the corresponding solution using the classical approach:

```wl
xClassical = Normalize[LinearSolve[a, b]]
```

---

Compare the classical result with the quantum one (using the inner product as a measure; the closer to one, the better):

```wl
Abs[Conjugate[xClassical] . xHHL]^2
```

<!-- => 0.998… -->

---

The number of phase-estimation qubits (3rd argument) determines the precision of the final outcome. Execute HHL with fewer phase-estimation qubits:

```wl
xHHLLowPrecision = Chop[N[HHL[a, b, 3]][]["Normalize"]["AmplitudesList"]]
```

---

Quantify the goodness of the HHL result compared with the classical result:

```wl
Abs[Conjugate[xClassical] . xHHLLowPrecision]^2
```

<!-- => 0.716… -->

## Possible Issues

Using Qiskit as the method, if there are too many qubits, the result is an association and not a quantum measurement object:

```wl
#| eval: false
QuantumCircuitOperator[{QuantumOperator["RandomUnitary", Range[10]], Range[10]}][Method -> "Qiskit"]
```

<!-- => <|"0000111101" -> 1, "0110010010" -> 1, …|> (an Association of bit-string counts) -->
