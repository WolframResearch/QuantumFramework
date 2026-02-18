# QuantumFramework Comprehensive Function Report

A detailed reference guide for all functions in the Wolfram QuantumFramework package.

---

## Table of Contents

1. [Core Quantum Objects](#part-i-core-quantum-objects)
   - [QuantumState](#quantumstate)
   - [QuantumOperator](#quantumoperator)
   - [QuantumBasis](#quantumbasis)
   - [QuditBasis](#quditbasis)
   - [QuantumCircuitOperator](#quantumcircuitoperator)
   - [QuantumChannel](#quantumchannel)
   - [QuantumMeasurement](#quantummeasurement)
   - [QuantumMeasurementOperator](#quantummeasurementoperator)

2. [Transformations & Operations](#part-ii-transformations--operations)
   - [QuantumTensorProduct](#quantumtensorproduct)
   - [QuantumPartialTrace](#quantumpartialtrace)
   - [QuantumWignerTransform](#quantumwignertransform)
   - [QuantumWeylTransform](#quantumweyltransform)
   - [QuantumPhaseSpaceTransform](#quantumphasespace-transform)
   - [QuantumWignerMICTransform](#quantumwignermictransform)

3. [Measures & Entanglement](#part-iii-measures--entanglement)
   - [QuantumDistance](#quantumdistance)
   - [QuantumEntangledQ](#quantumentangledq)
   - [QuantumEntanglementMonotone](#quantumentanglementmonotone)

4. [Evolution & Dynamics](#part-iv-evolution--dynamics)
   - [QuantumEvolve](#quantumevolve)

5. [State Estimation](#part-v-state-estimation)
   - [QuantumStateEstimate](#quantumstateestimate)
   - [QuantumMeasurementSimulation](#quantummeasurementsimulation)

6. [Optimization Functions](#part-vi-optimization-functions)
   - [GradientDescent](#gradientdescent)
   - [QuantumNaturalGradientDescent](#quantumnaturalgradientdescent)
   - [SPSRGradientValues](#spsrgradientvalues)
   - [FubiniStudyMetricTensor](#fubinistudymetrictensor)
   - [QuantumLinearSolve](#quantumlinearsolve)
   - [ParametrizedLayer](#parametrizedlayer)
   - [EntanglementLayer](#entanglementlayer)

7. [Second Quantization (Quantum Optics)](#part-vii-second-quantization)
   - [FockState](#fockstate)
   - [CoherentState](#coherentstate)
   - [AnnihilationOperator](#annihilationoperator)
   - [DisplacementOperator](#displacementoperator)
   - [SqueezeOperator](#squeezeoperator)
   - [BeamSplitterOperator](#beamsplitteroperator)
   - [WignerRepresentation](#wignerrepresentation)
   - [HusimiQRepresentation](#husimiqrepresentation)

8. [Advanced Features](#part-viii-advanced-features)
   - [QuantumMPS](#quantummps)
   - [PauliStabilizer](#paulistabilizer)

---

## Part I: Core Quantum Objects

---

### QuantumState

A quantum state represents the state of a quantum system, either as a pure state (state vector) or mixed state (density matrix).

#### Summary

`QuantumState` is the fundamental object for representing quantum states. It supports both pure states (represented as vectors) and mixed states (represented as density matrices).

#### Constructor Patterns

```mathematica
(* From amplitude vector *)
QuantumState[{a1, a2, ...}]

(* From density matrix *)
QuantumState[{{rho11, rho12}, {rho21, rho22}}]

(* Named states *)
QuantumState["StateName"]
QuantumState["StateName"[parameters...]]

(* With explicit basis *)
QuantumState[amplitudes, basis]
QuantumState[amplitudes, QuantumBasis[...]]

(* With options *)
QuantumState[..., "Label" -> label]
```

#### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| amplitudes | List | State vector amplitudes or density matrix |
| basis | QuantumBasis/QuditBasis | The basis for the state representation |

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"Label"` | None | Label for display purposes |
| `"Picture"` | "Schrodinger" | Quantum picture ("Schrodinger", "Heisenberg", "Interaction", "PhaseSpace") |

#### Properties (50+)

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"StateType"` | String | "Vector" or "Matrix" |
| `"State"` | List | Raw state data |
| `"Basis"` | QuantumBasis | The quantum basis |
| `"Amplitudes"` | Association | State amplitudes with basis labels |
| `"Amplitude"[label]` | Complex | Amplitude for specific basis element |
| `"Weights"` | Association | Squared amplitudes (weights) |
| `"Probabilities"` | Association | Measurement probabilities |
| `"Probability"[label]` | Real | Probability for specific outcome |
| `"StateVector"` | List | State vector representation |
| `"DensityMatrix"` | Matrix | Density matrix representation |
| `"AmplitudesList"` | List | List of amplitudes |
| `"ProbabilitiesList"` | List | List of probabilities |
| `"NormalizedState"` | QuantumState | Normalized quantum state |
| `"NormalizedAmplitudes"` | Association | Normalized amplitudes |
| `"NormalizedStateVector"` | List | Normalized state vector |
| `"NormalizedDensityMatrix"` | Matrix | Normalized density matrix |
| `"Entropy"` | Real | Von Neumann entropy |
| `"VonNeumannEntropy"` | Real | Von Neumann entropy |
| `"LogicalEntropy"` | Real | Logical entropy |
| `"Purity"` | Real | State purity (Tr[rho^2]) |
| `"Type"` | String | "Pure" or "Mixed" |
| `"PureStateQ"` | Boolean | True if pure state |
| `"MixedStateQ"` | Boolean | True if mixed state |
| `"MatrixQ"` | Boolean | True if density matrix representation |
| `"VectorQ"` | Boolean | True if vector representation |
| `"PhysicalQ"` | Boolean | True if valid physical state |
| `"NumericQ"` | Boolean | True if all numeric values |
| `"Kind"` | String | Kind of quantum object |
| `"Scalar"` | Complex | Scalar factor |
| `"Norm"` | Real | State norm |
| `"TraceNorm"` | Real | Trace norm |
| `"NormalizedQ"` | Boolean | True if normalized |
| `"BlochSphericalCoordinates"` | List | {theta, phi} Bloch sphere angles |
| `"BlochCartesianCoordinates"` | List | {x, y, z} Bloch vector |
| `"Projector"` | Matrix | Projection operator |
| `"Eigenvalues"` | List | Eigenvalues of density matrix |
| `"Eigenvectors"` | List | Eigenvectors of density matrix |
| `"Eigenstates"` | List | Eigenstates as QuantumState objects |
| `"Computational"` | QuantumState | State in computational basis |
| `"SchmidtBasis"` | QuantumState | State in Schmidt basis |
| `"SpectralBasis"` | QuantumState | State in spectral basis |
| `"PrimeBasis"` | QuantumState | State in prime basis |
| `"StateTensor"` | Array | Tensor representation |
| `"StateMatrix"` | Matrix | Matrix representation |
| `"Tensor"` | Array | Tensor representation |
| `"Matrix"` | Matrix | Matrix form |
| `"Purify"` | QuantumState | Purification of mixed state |
| `"Unpurify"` | QuantumState | Partial trace to unpurify |
| `"Bend"` | QuantumState | Channel-state duality map |
| `"Double"` | QuantumState | Doubled Hilbert space representation |
| `"Pure"` | QuantumState | Pure state projection |
| `"Mixed"` | QuantumState | Maximally mixed state |
| `"Trace"` | Complex | Trace of state |
| `"Transpose"` | QuantumState | Transposed state |
| `"Conjugate"` | QuantumState | Complex conjugate |
| `"ConjugateTranspose"` | QuantumState | Conjugate transpose |
| `"Inverse"` | QuantumState | Inverse state |
| `"ReverseOutput"` | QuantumState | Reverse output qudits |
| `"ReverseInput"` | QuantumState | Reverse input qudits |
| `"Reverse"` | QuantumState | Reverse all qudits |
| `"Split"` | List | Split into subsystems |
| `"Bipartition"` | List | Bipartite split |
| `"Disentangle"` | QuantumState | Disentangled state |
| `"Decompose"` | List | Decomposition |
| `"SchmidtDecompose"` | Association | Schmidt decomposition |
| `"Formula"` | Expression | Symbolic formula |
| `"Simplify"` | QuantumState | Simplified state |
| `"FullSimplify"` | QuantumState | Fully simplified state |
| `"Diagram"` | Graphics | Visual diagram |
| `"CircuitDiagram"` | Graphics | Circuit diagram |
| `"BlochPlot"` | Graphics | Bloch sphere visualization |
| `"AmplitudePlot"` | Graphics | Amplitude bar chart |
| `"ProbabilityPlot"` | Graphics | Probability bar chart |
| `"PieChart"` | Graphics | Pie chart of probabilities |
| `"SectorChart"` | Graphics | Sector chart |

#### Named States ($QuantumStateNames)

| Name | Description |
|------|-------------|
| `"0"` / `"Zero"` / `"Up"` | Computational basis state \|0> |
| `"1"` / `"One"` / `"Down"` | Computational basis state \|1> |
| `"Plus"` | Superposition state \|+> = (\|0> + \|1>)/sqrt(2) |
| `"Minus"` | Superposition state \|-> = (\|0> - \|1>)/sqrt(2) |
| `"Left"` | Circular polarization \|L> = (\|0> + i\|1>)/sqrt(2) |
| `"Right"` | Circular polarization \|R> = (\|0> - i\|1>)/sqrt(2) |
| `"PhiPlus"` | Bell state \|Phi+> = (\|00> + \|11>)/sqrt(2) |
| `"PhiMinus"` | Bell state \|Phi-> = (\|00> - \|11>)/sqrt(2) |
| `"PsiPlus"` | Bell state \|Psi+> = (\|01> + \|10>)/sqrt(2) |
| `"PsiMinus"` | Bell state \|Psi-> = (\|01> - \|10>)/sqrt(2) |
| `"BasisState"[n, dim]` | Computational basis state \|n> in dimension dim |
| `"Register"[n, val]` | n-qubit register with value val |
| `"UniformSuperposition"[n]` | Equal superposition of all n-qubit states |
| `"UniformMixture"[n]` | Maximally mixed state on n qubits |
| `"RandomPure"[n]` | Random pure state on n qubits |
| `"RandomMixed"[n]` | Random mixed state on n qubits |
| `"GHZ"[n]` | n-qubit GHZ state |
| `"Bell"` | Bell state (same as PhiPlus) |
| `"Dicke"[n, k]` | Dicke state with n qubits and k excitations |
| `"W"[n]` | n-qubit W state |
| `"Werner"[p]` | Werner state with parameter p |
| `"Graph"[g]` | Graph state from graph g |
| `"BlochVector"[{x,y,z}]` | State from Bloch vector coordinates |

#### Sample Code

```mathematica
(* Basic state creation *)
state1 = QuantumState[{1, 0}]                (* |0> state *)
state2 = QuantumState[{1, 1}/Sqrt[2]]        (* |+> state *)

(* Named states *)
bell = QuantumState["PhiPlus"]               (* Bell state *)
ghz = QuantumState["GHZ"[3]]                 (* 3-qubit GHZ *)

(* Density matrix *)
mixed = QuantumState[{{1/2, 0}, {0, 1/2}}]   (* Maximally mixed *)

(* Access properties *)
bell["Probabilities"]                         (* Measurement probabilities *)
bell["Entropy"]                              (* Von Neumann entropy *)
bell["BlochCartesianCoordinates"]            (* Bloch vector *)

(* Visualizations *)
state1["BlochPlot"]                          (* Bloch sphere plot *)
bell["ProbabilityPlot"]                      (* Probability bar chart *)

(* Operations *)
bell["Conjugate"]                            (* Complex conjugate *)
bell["Purify"]                               (* Purification *)
```

---

### QuantumOperator

A quantum operator represents a linear operator acting on quantum states, including unitary gates, Hamiltonians, and general operators.

#### Summary

`QuantumOperator` represents linear operators on quantum Hilbert spaces. It is used for quantum gates, Hamiltonians, and other linear transformations.

#### Constructor Patterns

```mathematica
(* From matrix *)
QuantumOperator[matrix]
QuantumOperator[matrix, {outputOrder, inputOrder}]

(* Named operators *)
QuantumOperator["OperatorName"]
QuantumOperator["OperatorName"[parameters...]]
QuantumOperator[{"OperatorName", parameters...}]

(* With specific order *)
QuantumOperator[matrix, order]
QuantumOperator["Gate", {qubit1, qubit2}]

(* Controlled gates *)
QuantumOperator[{"C", "X"}, {control, target}]

(* With options *)
QuantumOperator[..., "Label" -> label]
```

#### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| matrix | Matrix | Operator matrix representation |
| order | List | `{outputOrder, inputOrder}` or just `inputOrder` |
| name | String | Named operator identifier |

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"Label"` | Automatic | Display label for the operator |
| `"Picture"` | "Schrodinger" | Quantum picture |
| `"ParameterSpec"` | {} | Parameter specifications for parametric operators |

#### Properties (34+)

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"Order"` | List | `{outputOrder, inputOrder}` |
| `"InputOrder"` | List | Input qudit indices |
| `"OutputOrder"` | List | Output qudit indices |
| `"ControlOrder"` | List | Control qudit indices (for controlled gates) |
| `"TargetOrder"` | List | Target qudit indices |
| `"MatrixRepresentation"` | Matrix | Matrix in computational basis |
| `"Matrix"` | Matrix | Matrix representation |
| `"TensorRepresentation"` | Array | Tensor representation |
| `"Tensor"` | Array | Tensor form |
| `"Ordered"` | QuantumOperator | Ordered operator |
| `"OrderedInput"` | QuantumOperator | Input-ordered operator |
| `"OrderedOutput"` | QuantumOperator | Output-ordered operator |
| `"SortInput"` | QuantumOperator | Sort input indices |
| `"SortOutput"` | QuantumOperator | Sort output indices |
| `"Sort"` | QuantumOperator | Sorted operator |
| `"SortedQ"` | Boolean | True if sorted |
| `"ReverseOutput"` | QuantumOperator | Reverse output order |
| `"ReverseInput"` | QuantumOperator | Reverse input order |
| `"Reverse"` | QuantumOperator | Reverse all orders |
| `"Shift"[n]` | QuantumOperator | Shift qudit indices by n |
| `"Arity"` | Integer | Number of input qudits |
| `"MaxArity"` | Integer | Maximum arity |
| `"FullArity"` | Integer | Full arity including gaps |
| `"TargetArity"` | Integer | Target arity |
| `"Range"` | Integer | Range of qudit indices |
| `"HermitianQ"` | Boolean | True if Hermitian |
| `"UnitaryQ"` | Boolean | True if unitary |
| `"Eigenvalues"` | List | Eigenvalues |
| `"Eigenvectors"` | List | Eigenvectors |
| `"Eigensystem"` | List | `{eigenvalues, eigenvectors}` |
| `"Projectors"` | List | Projection operators |
| `"Transpose"` | QuantumOperator | Transposed operator |
| `"ConjugateTranspose"` | QuantumOperator | Conjugate transpose (dagger) |
| `"Dagger"` | QuantumOperator | Conjugate transpose |
| `"Dual"` | QuantumOperator | Dual operator |
| `"TraceNorm"` | Real | Trace norm |
| `"PauliDecompose"` | Association | Pauli decomposition |
| `"CircuitDiagram"` | Graphics | Circuit diagram |
| `"Simplify"` | QuantumOperator | Simplified operator |
| `"FullSimplify"` | QuantumOperator | Fully simplified operator |

#### Named Operators ($QuantumOperatorNames) - 80+ Gates

**Single-Qubit Gates:**

| Name | Description | Matrix |
|------|-------------|--------|
| `"Identity"` / `"I"` | Identity gate | I |
| `"X"` / `"PauliX"` / `"NOT"` | Pauli X (NOT gate) | [[0,1],[1,0]] |
| `"Y"` / `"PauliY"` | Pauli Y | [[0,-i],[i,0]] |
| `"Z"` / `"PauliZ"` | Pauli Z | [[1,0],[0,-1]] |
| `"H"` / `"Hadamard"` | Hadamard gate | (1/sqrt(2))[[1,1],[1,-1]] |
| `"S"` | S gate (Phase gate) | [[1,0],[0,i]] |
| `"T"` | T gate | [[1,0],[0,e^(i*pi/4)]] |
| `"V"` | V gate (sqrt(X)) | (1/2)[[1+i,1-i],[1-i,1+i]] |

**Rotation Gates:**

| Name | Description | Parameters |
|------|-------------|------------|
| `"RX"[theta]` | X-rotation | Angle theta |
| `"RY"[theta]` | Y-rotation | Angle theta |
| `"RZ"[theta]` | Z-rotation | Angle theta |
| `"Phase"[phi]` / `"P"[phi]` | Phase gate | Angle phi |
| `"U"[theta, phi, lambda]` | Universal single-qubit gate | Three Euler angles |
| `"U2"[phi, lambda]` | Two-parameter universal gate | Two angles |
| `"U3"[theta, phi, lambda]` | Three-parameter universal gate | Three angles |
| `"R"[axis, theta]` | Rotation about axis | Axis vector, angle |
| `"GlobalPhase"[phi]` | Global phase | Phase angle |

**Two-Qubit Gates:**

| Name | Description |
|------|-------------|
| `"CNOT"` / `"CX"` | Controlled-NOT |
| `"CY"` | Controlled-Y |
| `"CZ"` | Controlled-Z |
| `"SWAP"` | SWAP gate |
| `"RootSWAP"` | Square root of SWAP |
| `"CSWAP"` / `"Fredkin"` | Controlled-SWAP (Fredkin gate) |
| `"Braid"` | Braiding operator |

**Multi-Qubit Gates:**

| Name | Description |
|------|-------------|
| `"Toffoli"` | Toffoli gate (CCNOT) |
| `"Deutsch"[theta]` | Deutsch gate |
| `"C"["Gate", controls, targets]` | General controlled gate |
| `"C0"["Gate", controls, targets]` | Control on 0 state |
| `"Controlled"[gate]` | Controlled version of gate |
| `"Controlled0"[gate]` | Controlled on 0 |

**Fourier & Special Gates:**

| Name | Description |
|------|-------------|
| `"Fourier"[n]` | n-qubit Quantum Fourier Transform |
| `"InverseFourier"[n]` | Inverse QFT |
| `"Permutation"[perm]` | Permutation operator |

**ZX Calculus & Spiders:**

| Name | Description |
|------|-------------|
| `"Spider"[...]` | Generic spider |
| `"ZSpider"[phase, m, n]` | Z spider with phase |
| `"XSpider"[phase, m, n]` | X spider with phase |
| `"WSpider"[m, n]` | W spider |

**Measurement & Encoding:**

| Name | Description |
|------|-------------|
| `"Measure"` | Measurement operator |
| `"Encode"[state]` | State encoding |
| `"Copy"` | Copy operator |
| `"Decohere"` | Decoherence operator |
| `"Marginal"` | Marginal operator |
| `"Discard"` | Discard (partial trace) |
| `"Cup"` | Cup operator (bend) |
| `"Cap"` | Cap operator |
| `"Trace"` | Trace operator |
| `"Reset"` | Reset to |0> |

**Hamiltonians & Liouvillians:**

| Name | Description |
|------|-------------|
| `"Hamiltonian"[H]` | Hamiltonian operator |
| `"Liouvillian"[H, L]` | Lindbladian superoperator |

**Spin Operators:**

| Name | Description |
|------|-------------|
| `"JX"[j]` | Spin-j X operator |
| `"JY"[j]` | Spin-j Y operator |
| `"JZ"[j]` | Spin-j Z operator |
| `"J+"` / `"J-"` | Raising/lowering operators |

**Random & Miscellaneous:**

| Name | Description |
|------|-------------|
| `"RandomUnitary"[n]` | Random unitary on n qubits |
| `"RandomHermitian"[n]` | Random Hermitian operator |
| `"WignerD"[j, alpha, beta, gamma]` | Wigner D-matrix |
| `"Zero"` | Zero operator |
| `"Multiplexer"[...]` | Multiplexer gate |
| `"Switch"[...]` | Quantum switch |

#### Sample Code

```mathematica
(* Basic operators *)
X = QuantumOperator["X"]                        (* Pauli X *)
H = QuantumOperator["H"]                        (* Hadamard *)

(* Operators with order *)
cnot = QuantumOperator["CNOT", {1, 2}]          (* CNOT on qubits 1,2 *)
cz = QuantumOperator["CZ", {2, 3}]              (* CZ on qubits 2,3 *)

(* Parametric gates *)
rx = QuantumOperator[{"RX", Pi/4}]              (* X rotation by Pi/4 *)
phase = QuantumOperator[{"Phase", theta}]       (* Parametric phase gate *)

(* Controlled gates *)
ccnot = QuantumOperator["Toffoli", {1, 2, 3}]   (* Toffoli gate *)
cx = QuantumOperator[{"C", "X"}, {{1}, {2}}]    (* Controlled-X *)

(* Apply to state *)
state = QuantumState[{1, 0}]
result = H[state]                               (* Apply Hadamard *)

(* Chain operators *)
circuit = cnot @ (H \[CircleTimes] QuantumOperator["I"])
circuit["MatrixRepresentation"]

(* Properties *)
H["UnitaryQ"]                                   (* True *)
H["Eigenvalues"]                                (* {1, -1} *)
H["PauliDecompose"]                             (* Pauli decomposition *)

(* Circuit diagram *)
cnot["CircuitDiagram"]
```

---

### QuantumBasis

A quantum basis represents the basis for a Hilbert space, including input and output spaces.

#### Summary

`QuantumBasis` defines the complete basis structure for quantum states and operators, combining input and output `QuditBasis` objects.

#### Constructor Patterns

```mathematica
(* Default computational basis *)
QuantumBasis[]
QuantumBasis[dimension]
QuantumBasis[{dim1, dim2, ...}]

(* Named basis *)
QuantumBasis["BasisName"]
QuantumBasis[{"BasisName", parameters...}]

(* Custom basis *)
QuantumBasis[<|label1 -> vector1, ...|>]

(* With input/output specification *)
QuantumBasis["Output" -> outBasis, "Input" -> inBasis]
```

#### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| dimension | Integer | Dimension of the Hilbert space |
| dimensions | List | List of qudit dimensions |
| elements | Association | Named basis vectors |

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"Label"` | None | Basis label |
| `"Picture"` | "Schrodinger" | Quantum picture |
| `"Output"` | QuditBasis[] | Output space basis |
| `"Input"` | QuditBasis[] | Input space basis |

#### Properties

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"Dimension"` | Integer | Total dimension |
| `"Dimensions"` | List | Qudit dimensions |
| `"Qudits"` | Integer | Number of qudits |
| `"Output"` | QuditBasis | Output basis |
| `"Input"` | QuditBasis | Input basis |
| `"Elements"` | List | Basis elements |
| `"Names"` | List | Element names |
| `"Picture"` | String | Quantum picture |

#### Sample Code

```mathematica
(* Computational basis *)
basis2 = QuantumBasis[2]                        (* 2D computational basis *)
basis4 = QuantumBasis[{2, 2}]                   (* 2-qubit basis *)

(* Named basis *)
bellBasis = QuantumBasis["Bell"]                (* Bell basis *)

(* Custom basis *)
customBasis = QuantumBasis[<|
    "+" -> {1, 1}/Sqrt[2],
    "-" -> {1, -1}/Sqrt[2]
|>]

(* Access properties *)
basis4["Dimension"]                             (* 4 *)
basis4["Qudits"]                                (* 2 *)
```

---

### QuditBasis

A qudit basis represents the basis for a single qudit (d-dimensional quantum system).

#### Summary

`QuditBasis` is the building block for quantum bases, representing the basis states for individual qudits.

#### Constructor Patterns

```mathematica
(* Computational basis *)
QuditBasis[]
QuditBasis[dimension]
QuditBasis[{"Computational", dimension}]

(* Named basis *)
QuditBasis["BasisName"]
QuditBasis[{"BasisName", parameters...}]

(* Custom basis *)
QuditBasis[<|label1 -> vector1, ...|>]
```

#### Named Bases ($QuditBasisNames)

| Name | Description |
|------|-------------|
| `"Computational"` / `"Identity"` / `"I"` | Standard computational basis |
| `"PauliX"` / `"X"` | Pauli X eigenbasis (\|+>, \|->) |
| `"PauliY"` / `"Y"` | Pauli Y eigenbasis (\|L>, \|R>) |
| `"PauliZ"` / `"Z"` | Pauli Z eigenbasis (\|0>, \|1>) |
| `"JX"[j]` / `"JY"[j]` / `"JZ"[j]` | Spin-j angular momentum bases |
| `"J"[j]` / `"JI"[j]` | Spin-j identity basis |
| `"Bell"` | Bell basis |
| `"Fourier"[d]` | Fourier basis of dimension d |
| `"Schwinger"[d]` | Schwinger basis (d^2 elements) |
| `"Dirac"` | Dirac gamma matrix basis |
| `"Ivanovic"[d, n]` | Mutually unbiased basis n for dimension d |
| `"Wigner"[d]` | Wigner basis for phase space |
| `"WignerMIC"[d]` | Wigner MIC (minimal informationally complete) basis |
| `"Pauli"` | Pauli matrix basis |
| `"GellMann"[d]` | Generalized Gell-Mann basis |
| `"GellMannMIC"[d]` | Gell-Mann MIC basis |
| `"Bloch"[d]` | Bloch basis |
| `"BlochSphere"[d]` | Bloch sphere basis |
| `"Wootters"[d]` | Wootters basis |
| `"Feynman"[d]` | Feynman basis |
| `"Tetrahedron"` | SIC-POVM tetrahedron basis |
| `"RandomMIC"[d]` | Random MIC basis |
| `"RandomHaarMIC"[d]` | Random Haar MIC basis |
| `"RandomBlochMIC"[d]` | Random Bloch MIC basis |
| `"QBismSIC"[d]` | QBism SIC-POVM basis |
| `"HesseSIC"` | Hesse SIC (dimension 3) |
| `"HoggarSIC"` | Hoggar SIC (dimension 8) |

#### Sample Code

```mathematica
(* Standard bases *)
comp = QuditBasis[2]                            (* 2D computational *)
fourier = QuditBasis["Fourier"[4]]              (* 4D Fourier basis *)

(* Pauli eigenbases *)
xBasis = QuditBasis["PauliX"]                   (* |+>, |-> *)
yBasis = QuditBasis["PauliY"]                   (* |L>, |R> *)

(* Spin bases *)
spinHalf = QuditBasis["JZ"[1/2]]                (* Spin-1/2 *)
spinOne = QuditBasis["JZ"[1]]                   (* Spin-1 *)

(* Phase space bases *)
wigner = QuditBasis["Wigner"[3]]                (* 3D Wigner basis *)

(* MIC bases for tomography *)
sicPovm = QuditBasis["QBismSIC"[2]]             (* Qubit SIC-POVM *)
```

---

### QuantumCircuitOperator

A quantum circuit operator represents a sequence of quantum gates applied to qubits.

#### Summary

`QuantumCircuitOperator` represents a quantum circuit as a composition of quantum operators. It provides circuit visualization, compilation, and various analysis tools.

#### Constructor Patterns

```mathematica
(* From list of operators *)
QuantumCircuitOperator[{op1, op2, op3, ...}]

(* Named circuits *)
QuantumCircuitOperator["CircuitName"]
QuantumCircuitOperator["CircuitName"[parameters...]]

(* With options *)
QuantumCircuitOperator[ops, "Label" -> label]
```

#### Arguments

| Argument | Type | Description |
|----------|------|-------------|
| operators | List | List of QuantumOperator objects |
| name | String | Named circuit identifier |

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"Label"` | None | Circuit label |
| `"Expand"` | True | Whether to expand subcircuits in diagrams |

#### Properties

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"Association"` | Association | Internal data association |
| `"Operators"` | List | List of operators (excluding barriers) |
| `"Elements"` | List | All circuit elements including barriers |
| `"Diagram"` | Graphics | Circuit diagram visualization |
| `"OperatorCount"` | Integer | Number of operators |
| `"Orders"` | List | List of operator orders |
| `"CircuitOperator"` | QuantumOperator | Compiled circuit as single operator |
| `"QuantumOperator"` | QuantumOperator | Same as CircuitOperator |
| `"Compile"` | QuantumOperator | Compile to single operator |
| `"QiskitCircuit"` | String | Qiskit circuit representation |
| `"Label"` | Any | Circuit label |
| `"Depth"` | Integer | Circuit depth |
| `"Arity"` | Integer | Number of qubits |
| `"Width"` | Integer | Circuit width |
| `"TensorNetwork"` | TensorNetwork | Tensor network representation |
| `"TensorNetworkGraph"` | Graph | Tensor network graph |
| `"Topology"` | Graph | Circuit topology graph |
| `"Picture"` | String | Quantum picture |
| `"Parameters"` | List | Circuit parameters |
| `"ParameterArity"` | Integer | Number of parameters |
| `"InputOrder"` | List | Input qudit order |
| `"OutputOrder"` | List | Output qudit order |
| `"InputDimensions"` | List | Input dimensions |
| `"OutputDimensions"` | List | Output dimensions |
| `"Flatten"[n]` | QuantumCircuitOperator | Flatten nested circuits to depth n |
| `"Sort"` | QuantumCircuitOperator | Sorted circuit |
| `"Dagger"` | QuantumCircuitOperator | Adjoint circuit |
| `"Conjugate"` | QuantumCircuitOperator | Conjugated circuit |
| `"Double"` | QuantumCircuitOperator | Doubled circuit |
| `"Bend"` | QuantumCircuitOperator | Bent circuit |
| `"Layers"` | List | Circuit layers |
| `"Measurements"` | Integer | Number of measurements |
| `"Channels"` | Integer | Number of channels |
| `"QASM"` | String | OpenQASM representation |

#### Named Circuits ($QuantumCircuitOperatorNames)

| Name | Description |
|------|-------------|
| `"Bell"` | Bell state preparation circuit |
| `"GHZ"[n]` | n-qubit GHZ state preparation |
| `"Graph"[g]` | Graph state preparation |
| `"Fourier"[n]` | n-qubit Quantum Fourier Transform |
| `"InverseFourier"[n]` | Inverse QFT |
| `"PhaseEstimation"[U, n]` | Quantum phase estimation |
| `"GroverDiffusion"[n]` | Grover diffusion operator |
| `"GroverDiffusion0"[n]` | Grover diffusion about |0> |
| `"GroverPhaseDiffusion"[n]` | Phase version of diffusion |
| `"BooleanOracle"[f]` | Boolean function oracle |
| `"PhaseOracle"[f]` | Phase oracle |
| `"GrayOracle"[...]` | Gray code oracle |
| `"DeutschJozsaPhaseOracle"[f]` | Deutsch-Jozsa oracle |
| `"DeutschJozsa"[f]` | Deutsch-Jozsa algorithm |
| `"Deutsch"[f]` | Deutsch algorithm |
| `"SimonOracle"[s]` | Simon's algorithm oracle |
| `"Simon"[s]` | Simon's algorithm |
| `"Grover"[oracle, n]` | Grover's search algorithm |
| `"GroverPhase"[oracle, n]` | Phase version of Grover |
| `"Toffoli"` | Toffoli gate decomposition |
| `"Fredkin"` | Fredkin gate decomposition |
| `"BernsteinVazirani"[s]` | Bernstein-Vazirani algorithm |
| `"Trotterization"[H, t, n]` | Trotter decomposition of exp(-iHt) |
| `"Magic"` | Magic state preparation |
| `"Multiplexer"[...]` | Multiplexer circuit |
| `"ControlledMultiplexer"[...]` | Controlled multiplexer |
| `"CHSH"` | CHSH inequality test circuit |
| `"LeggettGarg"` | Leggett-Garg inequality test |
| `"WignerCHSH"` | Wigner CHSH circuit |

#### Sample Code

```mathematica
(* Create circuit from operators *)
circuit = QuantumCircuitOperator[{
    QuantumOperator["H", {1}],
    QuantumOperator["CNOT", {1, 2}]
}]

(* Named circuits *)
bellCircuit = QuantumCircuitOperator["Bell"]
qft = QuantumCircuitOperator["Fourier"[3]]
grover = QuantumCircuitOperator["Grover"[oracle, 3]]

(* View diagram *)
circuit["Diagram"]

(* Compile to single operator *)
U = circuit["CircuitOperator"]
U["MatrixRepresentation"]

(* Apply to state *)
input = QuantumState["00"]
output = circuit[input]

(* Circuit properties *)
circuit["Depth"]                                (* Circuit depth *)
circuit["OperatorCount"]                        (* Gate count *)
circuit["QASM"]                                 (* OpenQASM output *)

(* Decompose *)
circuit["Flatten"]["Operators"]                 (* All gates *)
circuit["Layers"]                               (* Layer by layer *)
```

---

### QuantumChannel

A quantum channel represents a completely positive trace-preserving (CPTP) map.

#### Summary

`QuantumChannel` represents quantum operations that include noise, decoherence, and other non-unitary evolution. Channels are specified via their Kraus operators, Choi matrix, or other representations.

#### Constructor Patterns

```mathematica
(* From Kraus operators *)
QuantumChannel[{K1, K2, ...}]

(* Named channels *)
QuantumChannel["ChannelName"]
QuantumChannel["ChannelName"[parameters...]]

(* From operator *)
QuantumChannel[quantumOperator]
```

#### Named Channels ($QuantumChannelNames)

| Name | Parameters | Description |
|------|------------|-------------|
| `"BitFlip"[p]` | p: probability | Bit flip with probability p |
| `"PhaseFlip"[p]` | p: probability | Phase flip with probability p |
| `"BitPhaseFlip"[p]` | p: probability | Bit-phase flip with probability p |
| `"Depolarizing"[p]` | p: probability | Depolarizing channel |
| `"AmplitudeDamping"[gamma]` | gamma: damping rate | Amplitude damping (T1 decay) |
| `"GeneralizedAmplitudeDamping"[p, gamma]` | p, gamma | Generalized amplitude damping |
| `"PhaseDamping"[lambda]` | lambda: damping rate | Phase damping (T2 decay) |
| `"ResetError"[p0, p1]` | p0, p1: reset probabilities | Reset error channel |

#### Properties

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"QuantumOperator"` | QuantumOperator | Underlying operator representation |
| `"Operator"` | QuantumOperator | Same as QuantumOperator |
| `"TraceOrder"` | List | Trace qudit indices |
| `"TraceQudits"` | Integer | Number of trace qudits |
| `"TraceDimensions"` | List | Dimensions of traced qudits |
| `"Trace"` | QuantumOperator | Trace map |
| `"TracePreservingQ"` | Boolean | True if trace-preserving |
| `"Dagger"` | QuantumChannel | Adjoint channel |
| `"Conjugate"` | QuantumChannel | Conjugate channel |
| `"Dual"` | QuantumChannel | Dual channel |
| `"Double"` | QuantumChannel | Doubled channel |
| `"Computational"` | QuantumChannel | In computational basis |
| `"Sort"` | QuantumChannel | Sorted channel |
| `"CircuitDiagram"` | Graphics | Circuit diagram |
| `"EvolutionChannel"` | QuantumChannel | Evolution channel |

#### Sample Code

```mathematica
(* Named channels *)
depol = QuantumChannel["Depolarizing"[0.1]]     (* 10% depolarizing *)
ampDamp = QuantumChannel["AmplitudeDamping"[0.05]]
phaseDamp = QuantumChannel["PhaseDamping"[0.1]]

(* Apply channel to state *)
state = QuantumState["Plus"]
noisyState = depol[state]

(* Check properties *)
depol["TracePreservingQ"]                       (* True *)

(* Compose channels *)
combined = depol @ ampDamp                      (* Channel composition *)

(* View representation *)
depol["QuantumOperator"]["MatrixRepresentation"]
```

---

### QuantumMeasurement

A quantum measurement represents a measurement outcome with associated probabilities.

#### Summary

`QuantumMeasurement` represents the result of measuring a quantum state, including the measurement outcomes and their probabilities.

#### Constructor Patterns

```mathematica
(* From measurement operator application *)
measurement = measurementOp[state]

(* From probabilities *)
QuantumMeasurement[<|outcome1 -> prob1, ...|>]
```

#### Properties

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"Probabilities"` | Association | Outcome probabilities |
| `"States"` | Association | Post-measurement states |
| `"Outcomes"` | List | Possible outcomes |

#### Sample Code

```mathematica
(* Create measurement operator *)
meas = QuantumMeasurementOperator["ComputationalBasis"]

(* Measure a state *)
state = QuantumState["Plus"]
result = meas[state]

(* Get probabilities *)
result["Probabilities"]                         (* <|0 -> 1/2, 1 -> 1/2|> *)
```

---

### QuantumMeasurementOperator

A quantum measurement operator represents a POVM (positive operator-valued measure) for quantum measurements.

#### Summary

`QuantumMeasurementOperator` defines the measurement basis or POVM elements for performing quantum measurements.

#### Constructor Patterns

```mathematica
(* Computational basis measurement *)
QuantumMeasurementOperator[]
QuantumMeasurementOperator["ComputationalBasis"]

(* Named POVMs *)
QuantumMeasurementOperator["POVMName"]
QuantumMeasurementOperator["POVMName"[parameters...]]

(* Custom POVM *)
QuantumMeasurementOperator[{E1, E2, ...}]
```

#### Named Measurement Operators ($QuantumMeasurementOperatorNames)

| Name | Description |
|------|-------------|
| `"M"` | Standard computational basis measurement |
| `"RandomHermitian"` | Random Hermitian POVM |
| `"WignerMICPOVM"` | Wigner MIC POVM |
| `"GellMannMICPOVM"[d]` | Gell-Mann MIC POVM |
| `"TetrahedronSICPOVM"` | Tetrahedron SIC-POVM |
| `"QBismSICPOVM"[d]` | QBism SIC-POVM |
| `"HesseSICPOVM"` | Hesse SIC-POVM (d=3) |
| `"HoggarSICPOVM"` | Hoggar SIC-POVM (d=8) |
| `"RandomPOVM"` | Random POVM elements |

#### Sample Code

```mathematica
(* Standard measurement *)
meas = QuantumMeasurementOperator[]

(* SIC-POVM measurement *)
sicMeas = QuantumMeasurementOperator["QBismSICPOVM"[2]]

(* Apply measurement *)
state = QuantumState["Plus"]
result = meas[state]
result["Probabilities"]

(* POVM elements *)
sicMeas["POVMElements"]
```

---

## Part II: Transformations & Operations

---

### QuantumTensorProduct

Computes the tensor product of quantum objects.

#### Summary

`QuantumTensorProduct` computes the tensor (Kronecker) product of quantum states, operators, bases, or circuits.

#### Constructor Patterns

```mathematica
QuantumTensorProduct[obj1, obj2, ...]
QuantumTensorProduct[{obj1, obj2, ...}]
obj1 \[CircleTimes] obj2                        (* Infix notation *)
```

#### Sample Code

```mathematica
(* Tensor product of states *)
state1 = QuantumState[{1, 0}]
state2 = QuantumState[{0, 1}]
combined = QuantumTensorProduct[state1, state2]  (* |01> *)

(* Using infix notation *)
twoQubit = state1 \[CircleTimes] state2

(* Tensor product of operators *)
hh = QuantumTensorProduct[
    QuantumOperator["H"],
    QuantumOperator["H"]
]                                               (* H tensor H *)

(* Tensor product of bases *)
basis = QuantumTensorProduct[QuditBasis[2], QuditBasis[3]]
```

---

### QuantumPartialTrace

Computes the partial trace over specified subsystems.

#### Summary

`QuantumPartialTrace` traces out (removes) specified qudits from a quantum state or operator, producing a reduced state.

#### Constructor Patterns

```mathematica
QuantumPartialTrace[state, {qudit1, qudit2, ...}]
QuantumPartialTrace[operator, {qudit1, qudit2, ...}]
```

#### Sample Code

```mathematica
(* Create Bell state *)
bell = QuantumState["PhiPlus"]

(* Trace out second qubit *)
reduced = QuantumPartialTrace[bell, {2}]
reduced["DensityMatrix"]                        (* Maximally mixed state *)

(* Partial trace of operator *)
cnot = QuantumOperator["CNOT", {1, 2}]
traced = QuantumPartialTrace[cnot, {2}]
```

---

### QuantumWignerTransform

Transforms quantum states to/from Wigner phase space representation.

#### Summary

`QuantumWignerTransform` converts quantum states between Hilbert space and discrete Wigner function representations.

#### Constructor Patterns

```mathematica
QuantumWignerTransform[state]
QuantumWignerTransform[state, dimension]
```

#### Sample Code

```mathematica
(* Transform to Wigner representation *)
state = QuantumState["Plus"]
wigner = QuantumWignerTransform[state]

(* Properties of Wigner function *)
wigner["Probabilities"]
```

---

### QuantumWeylTransform

Performs the Weyl transform between operators and phase space functions.

#### Summary

`QuantumWeylTransform` implements the Weyl transform, mapping between operators and their phase space representations.

#### Constructor Patterns

```mathematica
QuantumWeylTransform[operator]
```

---

### QuantumPhaseSpaceTransform

General phase space transformation for quantum states.

#### Summary

`QuantumPhaseSpaceTransform` provides a general interface for various phase space representations including Wigner, Husimi Q, and others.

#### Constructor Patterns

```mathematica
QuantumPhaseSpaceTransform[state, "Method" -> method]
```

---

### QuantumWignerMICTransform

Transforms to Wigner MIC (Minimal Informationally Complete) representation.

#### Summary

`QuantumWignerMICTransform` converts quantum states to the Wigner MIC basis, useful for state tomography.

#### Constructor Patterns

```mathematica
QuantumWignerMICTransform[state]
```

---

## Part III: Measures & Entanglement

---

### QuantumDistance

Computes various distance measures between quantum states.

#### Summary

`QuantumDistance` calculates distance or distinguishability measures between two quantum states using various metrics.

#### Constructor Patterns

```mathematica
QuantumDistance[state1, state2]
QuantumDistance[state1, state2, "DistanceType"]
```

#### Distance Types ($QuantumDistances)

| Name | Description | Range |
|------|-------------|-------|
| `"Fidelity"` | Fidelity F(rho, sigma) | [0, 1] |
| `"RelativeEntropy"` | Quantum relative entropy S(rho \|\| sigma) | [0, inf) |
| `"RelativePurity"` | Relative purity | [0, 1] |
| `"Trace"` | Trace distance D(rho, sigma) = (1/2)\|rho - sigma\|_1 | [0, 1] |
| `"Bures"` | Bures distance | [0, sqrt(2)] |
| `"BuresAngle"` | Bures angle (arccos of fidelity) | [0, pi/2] |
| `"HilbertSchmidt"` | Hilbert-Schmidt distance | [0, sqrt(2)] |
| `"Bloch"` | Bloch vector distance (qubits) | [0, 2] |

#### Sample Code

```mathematica
(* Create two states *)
state1 = QuantumState["0"]
state2 = QuantumState["Plus"]

(* Compute fidelity *)
QuantumDistance[state1, state2, "Fidelity"]     (* 1/Sqrt[2] *)

(* Trace distance *)
QuantumDistance[state1, state2, "Trace"]

(* Compare all distances *)
Table[
    dist -> QuantumDistance[state1, state2, dist],
    {dist, {"Fidelity", "Trace", "Bures", "HilbertSchmidt"}}
]
```

---

### QuantumEntangledQ

Tests whether a quantum state is entangled.

#### Summary

`QuantumEntangledQ` determines if a bipartite quantum state is entangled (non-separable).

#### Constructor Patterns

```mathematica
QuantumEntangledQ[state]
QuantumEntangledQ[state, {subsystem1, subsystem2}]
```

#### Sample Code

```mathematica
(* Test Bell state *)
bell = QuantumState["PhiPlus"]
QuantumEntangledQ[bell]                         (* True *)

(* Test product state *)
product = QuantumTensorProduct[
    QuantumState["0"],
    QuantumState["1"]
]
QuantumEntangledQ[product]                      (* False *)
```

---

### QuantumEntanglementMonotone

Computes various entanglement measures for quantum states.

#### Summary

`QuantumEntanglementMonotone` calculates quantitative measures of entanglement using various monotones.

#### Constructor Patterns

```mathematica
QuantumEntanglementMonotone[state]
QuantumEntanglementMonotone[state, "MonotoneName"]
QuantumEntanglementMonotone[state, "MonotoneName", {subsystems}]
```

#### Entanglement Monotones ($QuantumEntanglementMonotones)

| Name | Description | Applicable States |
|------|-------------|-------------------|
| `"Concurrence"` | Concurrence | 2-qubit states |
| `"Negativity"` | Negativity | Bipartite states |
| `"LogNegativity"` | Logarithmic negativity | Bipartite states |
| `"EntanglementEntropy"` | Entanglement entropy (von Neumann) | Pure bipartite states |
| `"RenyiEntropy"[alpha]` | Renyi entanglement entropy | Pure bipartite states |
| `"Realignment"` | Realignment criterion | Bipartite states |
| `"MutualInformationI"` | Mutual information (type I) | Bipartite states |
| `"MutualInformationJ"` | Mutual information (type J) | Bipartite states |
| `"Discord"` | Quantum discord | Bipartite states |

#### Sample Code

```mathematica
(* Bell state entanglement *)
bell = QuantumState["PhiPlus"]

(* Concurrence (maximal for Bell states) *)
QuantumEntanglementMonotone[bell, "Concurrence"]    (* 1 *)

(* Entanglement entropy *)
QuantumEntanglementMonotone[bell, "EntanglementEntropy"]  (* Log[2] *)

(* Negativity *)
QuantumEntanglementMonotone[bell, "Negativity"]

(* Werner state with varying entanglement *)
werner = QuantumState["Werner"[0.5]]
QuantumEntanglementMonotone[werner, "Negativity"]
```

---

## Part IV: Evolution & Dynamics

---

### QuantumEvolve

Evolves quantum states under Hamiltonian or Lindbladian dynamics.

#### Summary

`QuantumEvolve` simulates the time evolution of quantum states under Schrodinger equation (unitary) or Lindblad master equation (open system) dynamics.

#### Constructor Patterns

```mathematica
(* Unitary evolution *)
QuantumEvolve[hamiltonian, initialState, time]
QuantumEvolve[hamiltonian, initialState, {t, tmin, tmax}]

(* Lindblad evolution *)
QuantumEvolve[{hamiltonian, {L1, L2, ...}}, initialState, time]

(* Evolution operator only *)
QuantumEvolve[hamiltonian, None, time]
```

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"AdditionalEquations"` | {} | Extra differential equations |
| `"ReturnEquations"` | False | Return the differential equations |
| `"ReturnSolution"` | False | Return the full NDSolve solution |
| `"MergeInterpolatingFunctions"` | True | Merge interpolating functions |

(Also accepts all options from `NDSolveValue` and `DSolveValue`)

#### Sample Code

```mathematica
(* Define Hamiltonian *)
H = QuantumOperator[{{0, 1}, {1, 0}}]           (* sigma_x *)

(* Initial state *)
psi0 = QuantumState["0"]

(* Evolve for time t *)
psiT = QuantumEvolve[H, psi0, t]
psiT /. t -> Pi/2                               (* Evolved state *)

(* Time-dependent evolution *)
psiT = QuantumEvolve[H, psi0, {t, 0, 2*Pi}]
psiT[Pi/4]["Probabilities"]                     (* State at t=Pi/4 *)

(* Lindblad evolution with decay *)
gamma = 0.1
L = Sqrt[gamma] * QuantumOperator[{{0, 1}, {0, 0}}]  (* Lowering operator *)
rhoT = QuantumEvolve[{H, {L}}, psi0, {t, 0, 10}]

(* Get evolution operator *)
U = QuantumEvolve[H, None, t]
U["MatrixRepresentation"]
```

---

## Part V: State Estimation

---

### QuantumStateEstimate

Reconstructs a quantum state from measurement data.

#### Summary

`QuantumStateEstimate` performs quantum state tomography, reconstructing the density matrix from measurement outcomes.

#### Constructor Patterns

```mathematica
QuantumStateEstimate[measurementData, "Method" -> method]
QuantumStateEstimate[frequencies, basis]
```

#### Methods

- `"MaximumLikelihood"` - Maximum likelihood estimation
- `"LinearInversion"` - Linear inversion tomography
- `"Bayesian"` - Bayesian estimation

#### Sample Code

```mathematica
(* Simulate measurements *)
state = QuantumState["Plus"]
meas = QuantumMeasurementOperator["QBismSICPOVM"[2]]
data = QuantumMeasurementSimulation[meas, state, 1000]

(* Reconstruct state *)
estimated = QuantumStateEstimate[data]
QuantumDistance[estimated, state, "Fidelity"]   (* Should be close to 1 *)
```

---

### QuantumMeasurementSimulation

Simulates quantum measurements with statistical sampling.

#### Summary

`QuantumMeasurementSimulation` generates simulated measurement data by sampling from the theoretical probability distribution.

#### Constructor Patterns

```mathematica
QuantumMeasurementSimulation[measurementOp, state, numSamples]
```

#### Sample Code

```mathematica
(* Simulate 1000 measurements of Bell state *)
bell = QuantumState["PhiPlus"]
meas = QuantumMeasurementOperator[]
samples = QuantumMeasurementSimulation[meas, bell, 1000]

(* Check frequencies *)
Counts[samples]                                 (* Should show ~500 each for 00, 11 *)
```

---

## Part VI: Optimization Functions

---

### GradientDescent

Performs gradient descent optimization for variational quantum algorithms.

#### Summary

`GradientDescent` implements gradient descent optimization for parameterized quantum circuits, supporting various gradient estimation methods.

#### Constructor Patterns

```mathematica
GradientDescent[costFunction, initialParams, options]
```

#### Options

| Option | Default | Description |
|--------|---------|-------------|
| `"LearningRate"` | 0.1 | Step size for updates |
| `"MaxIterations"` | 100 | Maximum number of iterations |
| `"Tolerance"` | 10^-6 | Convergence tolerance |
| `"GradientMethod"` | "ParameterShift" | Method for gradient computation |

---

### QuantumNaturalGradientDescent

Natural gradient descent using the quantum Fisher information metric.

#### Summary

`QuantumNaturalGradientDescent` performs optimization using the natural gradient, which accounts for the geometry of the parameter space.

#### Constructor Patterns

```mathematica
QuantumNaturalGradientDescent[costFunction, initialParams, options]
```

---

### SPSRGradientValues

Computes gradients using the Stochastic Parameter Shift Rule.

#### Summary

`SPSRGradientValues` estimates gradients of parameterized quantum circuits using the stochastic parameter shift rule.

#### Constructor Patterns

```mathematica
SPSRGradientValues[circuit, parameters, observable]
```

---

### ASPSRGradientValues

Computes gradients using the Analytical Stochastic Parameter Shift Rule.

#### Summary

`ASPSRGradientValues` provides analytical gradient computation for certain circuit structures.

---

### FubiniStudyMetricTensor

Computes the Fubini-Study metric tensor for a parameterized state.

#### Summary

`FubiniStudyMetricTensor` calculates the quantum geometric tensor (Fubini-Study metric), which quantifies distinguishability in parameter space.

#### Constructor Patterns

```mathematica
FubiniStudyMetricTensor[circuit, parameters]
FubiniStudyMetricTensor[parameterizedState]
```

---

### QuantumLinearSolve

Solves linear systems using quantum algorithms.

#### Summary

`QuantumLinearSolve` implements quantum algorithms for solving linear systems of equations (HHL algorithm variants).

#### Constructor Patterns

```mathematica
QuantumLinearSolve[matrix, vector, options]
```

---

### ParametrizedLayer

Creates a parametrized layer of single-qubit rotations.

#### Summary

`ParametrizedLayer` constructs a layer of parameterized single-qubit gates for variational circuits.

#### Constructor Patterns

```mathematica
ParametrizedLayer[numQubits, "GateType"]
ParametrizedLayer[numQubits, {angles...}]
```

#### Sample Code

```mathematica
(* Create parametrized rotation layer *)
layer = ParametrizedLayer[3, "RY"]

(* With specific angles *)
layer = ParametrizedLayer[3, {theta1, theta2, theta3}]
```

---

### EntanglementLayer

Creates an entangling layer for variational circuits.

#### Summary

`EntanglementLayer` constructs a layer of entangling gates (typically CNOTs) for variational quantum circuits.

#### Constructor Patterns

```mathematica
EntanglementLayer[numQubits]
EntanglementLayer[numQubits, "Pattern" -> pattern]
```

#### Sample Code

```mathematica
(* Linear entanglement pattern *)
entLayer = EntanglementLayer[4, "Pattern" -> "Linear"]

(* Circular pattern *)
entLayer = EntanglementLayer[4, "Pattern" -> "Circular"]
```

---

### GenerateParameters

Generates initial parameters for variational circuits.

#### Summary

`GenerateParameters` creates initial parameter values for variational quantum algorithms with various initialization strategies.

#### Constructor Patterns

```mathematica
GenerateParameters[circuit]
GenerateParameters[numParameters, "Method" -> method]
```

---

## Part VII: Second Quantization

Functions for quantum optics and bosonic systems using Fock space representations.

---

### SetFockSpaceSize

Sets the truncation dimension for Fock space calculations.

#### Summary

`SetFockSpaceSize` sets the maximum photon number (dimension) for Fock space truncation.

#### Constructor Patterns

```mathematica
SetFockSpaceSize[n]
```

#### Sample Code

```mathematica
SetFockSpaceSize[20]                            (* Truncate at 20 photons *)
```

---

### FockState

Creates a Fock (number) state.

#### Summary

`FockState` creates a quantum state with a definite number of photons/excitations.

#### Constructor Patterns

```mathematica
FockState[n]                                    (* |n> state *)
FockState[n, dimension]                         (* With explicit dimension *)
```

#### Sample Code

```mathematica
vacuum = FockState[0]                           (* Vacuum state |0> *)
single = FockState[1]                           (* Single photon |1> *)
fock5 = FockState[5]                            (* 5 photon state *)
```

---

### CoherentState

Creates a coherent state.

#### Summary

`CoherentState` creates a coherent state, the closest quantum analog to a classical electromagnetic wave.

#### Constructor Patterns

```mathematica
CoherentState[alpha]                            (* Coherent state |alpha> *)
CoherentState[alpha, dimension]
```

#### Sample Code

```mathematica
(* Create coherent state with amplitude 2 *)
coh = CoherentState[2]

(* Properties *)
coh["Probabilities"]                            (* Poisson distribution *)
```

---

### AnnihilationOperator

Creates the bosonic annihilation (lowering) operator.

#### Summary

`AnnihilationOperator` returns the photon annihilation operator 'a' that removes one quantum of excitation.

#### Constructor Patterns

```mathematica
AnnihilationOperator[]
AnnihilationOperator[dimension]
```

#### Sample Code

```mathematica
a = AnnihilationOperator[]
adag = a["Dagger"]                              (* Creation operator *)

(* Number operator *)
n = adag @ a

(* Apply to Fock state *)
fock3 = FockState[3]
a[fock3]                                        (* Sqrt[3] |2> *)
```

---

### DisplacementOperator

Creates the displacement operator D(alpha).

#### Summary

`DisplacementOperator` creates the unitary operator that displaces a state in phase space, generating coherent states from vacuum.

#### Constructor Patterns

```mathematica
DisplacementOperator[alpha]
```

#### Sample Code

```mathematica
D = DisplacementOperator[2]

(* Coherent state from vacuum *)
vacuum = FockState[0]
coherent = D[vacuum]                            (* Same as CoherentState[2] *)
```

---

### SqueezeOperator

Creates the squeezing operator S(r, phi).

#### Summary

`SqueezeOperator` creates a unitary operator that produces squeezed states with reduced uncertainty in one quadrature.

#### Constructor Patterns

```mathematica
SqueezeOperator[r]                              (* Real squeezing parameter *)
SqueezeOperator[r, phi]                         (* With phase *)
```

#### Sample Code

```mathematica
S = SqueezeOperator[0.5]

(* Squeezed vacuum *)
squeezed = S[FockState[0]]
```

---

### ThermalState

Creates a thermal (Bose-Einstein) state.

#### Summary

`ThermalState` creates a mixed state representing thermal equilibrium at a given mean photon number.

#### Constructor Patterns

```mathematica
ThermalState[nbar]                              (* Mean photon number *)
```

---

### CatState

Creates a Schrodinger cat state.

#### Summary

`CatState` creates a superposition of two coherent states, a macroscopic quantum superposition.

#### Constructor Patterns

```mathematica
CatState[alpha]                                 (* Even cat state *)
CatState[alpha, "Parity" -> parity]             (* Even or odd *)
```

---

### QuadratureOperators

Returns the position and momentum quadrature operators.

#### Summary

`QuadratureOperators` returns the X and P quadrature operators for a harmonic oscillator.

#### Constructor Patterns

```mathematica
{X, P} = QuadratureOperators[]
```

---

### BeamSplitterOperator

Creates a beam splitter operator.

#### Summary

`BeamSplitterOperator` creates a two-mode unitary operator representing a beam splitter.

#### Constructor Patterns

```mathematica
BeamSplitterOperator[theta]                     (* Mixing angle *)
BeamSplitterOperator[theta, phi]                (* With phase *)
```

#### Sample Code

```mathematica
BS = BeamSplitterOperator[Pi/4]                 (* 50:50 beam splitter *)

(* Two-mode state *)
input = QuantumTensorProduct[FockState[1], FockState[0]]
output = BS[input]                              (* Creates superposition *)
```

---

### PhaseShiftOperator

Creates a phase shift operator.

#### Summary

`PhaseShiftOperator` creates an operator that applies a phase shift proportional to photon number.

#### Constructor Patterns

```mathematica
PhaseShiftOperator[phi]
```

---

### WignerRepresentation

Computes the Wigner function for a quantum state.

#### Summary

`WignerRepresentation` calculates the Wigner quasi-probability distribution for visualizing quantum states in phase space.

#### Constructor Patterns

```mathematica
WignerRepresentation[state, {x, y}]
WignerRepresentation[state, {xmin, xmax}, {ymin, ymax}]
```

#### Sample Code

```mathematica
(* Wigner function of coherent state *)
coh = CoherentState[2]
WignerRepresentation[coh, {x, -4, 4}, {p, -4, 4}]

(* Plot Wigner function *)
DensityPlot[
    WignerRepresentation[coh, {x, p}],
    {x, -4, 4}, {p, -4, 4}
]
```

---

### HusimiQRepresentation

Computes the Husimi Q function for a quantum state.

#### Summary

`HusimiQRepresentation` calculates the Husimi Q quasi-probability distribution, a smoothed version of the Wigner function.

#### Constructor Patterns

```mathematica
HusimiQRepresentation[state, alpha]
HusimiQRepresentation[state, {x, y}]
```

---

### G2Coherence

Computes the second-order coherence function g^(2)(0).

#### Summary

`G2Coherence` calculates the degree of second-order coherence, which characterizes photon statistics (bunching, antibunching, coherent).

#### Constructor Patterns

```mathematica
G2Coherence[state]
```

#### Sample Code

```mathematica
(* Coherent state: g^(2) = 1 *)
G2Coherence[CoherentState[2]]                   (* 1 *)

(* Fock state: g^(2) < 1 (antibunching) *)
G2Coherence[FockState[1]]                       (* 0 *)

(* Thermal state: g^(2) > 1 (bunching) *)
G2Coherence[ThermalState[1]]                    (* 2 *)
```

---

### OperatorVariance

Computes the variance of an operator in a given state.

#### Summary

`OperatorVariance` calculates Var(O) = <O^2> - <O>^2 for a quantum state.

#### Constructor Patterns

```mathematica
OperatorVariance[operator, state]
```

---

## Part VIII: Advanced Features

---

### QuantumMPS

Matrix Product State representation for quantum states.

#### Summary

`QuantumMPS` represents quantum states using Matrix Product States (MPS), an efficient representation for many-body systems with limited entanglement.

#### Constructor Patterns

```mathematica
QuantumMPS[state]
QuantumMPS[tensors]
```

#### Properties

| Property | Return Type | Description |
|----------|-------------|-------------|
| `"Tensors"` | List | List of MPS tensors |
| `"BondDimensions"` | List | Bond dimensions |
| `"QuantumState"` | QuantumState | Reconstructed quantum state |

---

### QuantumMPO

Matrix Product Operator representation.

#### Summary

`QuantumMPO` represents quantum operators using Matrix Product Operators (MPO), enabling efficient representation of many-body operators.

#### Constructor Patterns

```mathematica
QuantumMPO[operator]
QuantumMPO[tensors]
```

---

### PauliStabilizer

Pauli stabilizer formalism for Clifford circuits.

#### Summary

`PauliStabilizer` implements the stabilizer formalism for efficient simulation of Clifford circuits.

#### Constructor Patterns

```mathematica
PauliStabilizer[generators]
PauliStabilizer["NamedState"]
```

#### Sample Code

```mathematica
(* Bell state stabilizer *)
stab = PauliStabilizer["Bell"]

(* GHZ stabilizer *)
ghzStab = PauliStabilizer["GHZ"[4]]
ghzStab["Generators"]
```

---

## Appendix A: Complete Property Reference

### All QuantumState Properties

```
"StateType", "State", "Basis",
"Amplitudes", "Amplitude", "Weights", "Probabilities", "Probability",
"StateVector", "DensityMatrix",
"AmplitudesList", "ProbabilitiesList",
"NormalizedState", "NormalizedAmplitudes", "NormalizedStateVector", "NormalizedDensityMatrix",
"Entropy", "VonNeumannEntropy", "LogicalEntropy",
"Purity", "Type", "PureStateQ", "MixedStateQ", "MatrixQ", "VectorQ",
"PhysicalQ", "NumericQ", "Kind", "Scalar", "Norm", "TraceNorm", "NormalizedQ",
"BlochSphericalCoordinates", "BlochCartesianCoordinates",
"Projector", "Eigenvalues", "Eigenvectors", "Eigenstates",
"Computational", "SchmidtBasis", "SpectralBasis", "PrimeBasis", "UniformBasis",
"StateTensor", "StateMatrix", "Tensor", "Matrix", "Table",
"Purify", "Unpurify", "Bend", "BendDual", "Unbend", "Double",
"Pure", "Mixed", "Trace", "Transpose", "Conjugate", "ConjugateTranspose", "Inverse",
"Physical", "ReverseOutput", "ReverseInput", "Reverse",
"Split", "SplitDual", "Bipartition", "Disentangle", "Decompose", "SchmidtDecompose",
"Formula", "Simplify", "FullSimplify",
"Diagram", "CircuitDiagram", "BlochPlot", "AmplitudePlot", "ProbabilityPlot",
"PieChart", "SectorChart", "PauliTree"
```

### All QuantumOperator Properties

```
"Order", "InputOrder", "OutputOrder", "ControlOrder", "TargetOrder",
"MatrixRepresentation", "Matrix", "TensorRepresentation", "Tensor",
"Ordered", "OrderedInput", "OrderedOutput",
"SortInput", "SortOutput", "Sort", "SortedQ",
"ReverseOutput", "ReverseInput", "Reverse", "Shift",
"OrderedMatrixRepresentation", "OrderedMatrix",
"OrderedTensorRepresentation", "OrderedTensor",
"Arity", "MaxArity", "FullArity", "TargetArity",
"Range", "FullInputOrder", "FullOutputOrder", "FullOrder",
"InputQuditOrder", "OutputQuditOrder", "QuditOrder",
"HermitianQ", "UnitaryQ", "Eigenvalues", "Eigenvectors", "Eigensystem", "Projectors",
"Transpose", "ConjugateTranspose", "UnstackOutput", "UnstackInput",
"QuantumOperator", "Operator", "Computational", "Diagonalize",
"Dagger", "Dual", "TraceNorm", "PauliDecompose",
"CircuitDiagram", "OrderedFormula",
"Simplify", "FullSimplify", "Chop", "ComplexExpand"
```

### All QuantumCircuitOperator Properties

```
"Association", "Operators", "Diagram", "OperatorCount", "Orders",
"CircuitOperator", "QiskitCircuit", "Label",
"Depth", "Arity", "Width", "TensorNetwork", "TensorNetworkGraph", "Topology",
"Properties", "Picture", "Parameters", "ParameterArity",
"InputOrder", "OutputOrder", "InputDimensions", "OutputDimensions",
"Flatten", "Sort", "Dagger", "Conjugate", "Double", "Bend",
"Layers", "Measurements", "Channels", "QASM"
```

---

## Appendix B: Quick Reference Tables

### Distance Measures

| Measure | Formula | Range | Notes |
|---------|---------|-------|-------|
| Fidelity | F = (Tr[sqrt(sqrt(rho) sigma sqrt(rho))])^2 | [0,1] | 1 = identical |
| Trace | D = (1/2) Tr[\|rho - sigma\|] | [0,1] | 1 = orthogonal |
| Bures | D_B = sqrt(2(1-sqrt(F))) | [0,sqrt(2)] | Metric |
| Hilbert-Schmidt | D_HS = sqrt(Tr[(rho-sigma)^2]) | [0,sqrt(2)] | Easy to compute |

### Entanglement Monotones

| Monotone | States | Range | Notes |
|----------|--------|-------|-------|
| Concurrence | 2-qubit | [0,1] | Analytic formula |
| Negativity | Bipartite | [0,1/2] | Computable |
| Log Negativity | Bipartite | [0,log(d)] | Additive |
| Entanglement Entropy | Pure bipartite | [0,log(d)] | Von Neumann |

### Common Gate Matrices

| Gate | Matrix |
|------|--------|
| X | [[0,1],[1,0]] |
| Y | [[0,-i],[i,0]] |
| Z | [[1,0],[0,-1]] |
| H | (1/sqrt(2))[[1,1],[1,-1]] |
| S | [[1,0],[0,i]] |
| T | [[1,0],[0,e^(i*pi/4)]] |
| CNOT | [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]] |

---

## Appendix C: Common Workflows

### Workflow 1: State Preparation and Measurement

```mathematica
(* 1. Create initial state *)
initial = QuantumState["0"]

(* 2. Build circuit *)
circuit = QuantumCircuitOperator[{
    QuantumOperator["H", {1}],
    QuantumOperator["CNOT", {1, 2}]
}]

(* 3. Apply circuit *)
final = circuit[QuantumState["00"]]

(* 4. Measure *)
meas = QuantumMeasurementOperator[]
result = meas[final]
result["Probabilities"]
```

### Workflow 2: Variational Quantum Eigensolver (VQE)

```mathematica
(* 1. Define Hamiltonian *)
H = QuantumOperator[{{1, 0}, {0, -1}}]

(* 2. Create parameterized ansatz *)
ansatz[theta_] := QuantumCircuitOperator[{
    QuantumOperator[{"RY", theta}, {1}]
}]

(* 3. Cost function *)
cost[theta_] := Module[{state},
    state = ansatz[theta][QuantumState["0"]];
    Re[Tr[H["Matrix"] . state["DensityMatrix"]]]
]

(* 4. Optimize *)
result = FindMinimum[cost[theta], {theta, 0}]
```

### Workflow 3: Quantum State Tomography

```mathematica
(* 1. Prepare unknown state *)
unknown = QuantumState["Plus"]

(* 2. Define POVM *)
povm = QuantumMeasurementOperator["QBismSICPOVM"[2]]

(* 3. Simulate measurements *)
data = QuantumMeasurementSimulation[povm, unknown, 10000]

(* 4. Reconstruct *)
reconstructed = QuantumStateEstimate[data]

(* 5. Verify *)
QuantumDistance[unknown, reconstructed, "Fidelity"]
```

### Workflow 4: Open System Dynamics

```mathematica
(* 1. System Hamiltonian *)
H = QuantumOperator["Z"]

(* 2. Jump operators (decay) *)
gamma = 0.1
L = Sqrt[gamma] * QuantumOperator[{{0, 1}, {0, 0}}]

(* 3. Initial state *)
psi0 = QuantumState["Plus"]

(* 4. Evolve *)
rho[t_] = QuantumEvolve[{H, {L}}, psi0, {t, 0, 50}]

(* 5. Plot evolution *)
Plot[rho[t]["Purity"], {t, 0, 50}]
```

---

*Report generated for QuantumFramework*
*Total Functions Documented: 50+*
*Total Named Instances: 150+*
*Total Properties: 100+*
