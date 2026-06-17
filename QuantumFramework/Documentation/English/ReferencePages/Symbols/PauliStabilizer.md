---
Template: Symbol
Name: PauliStabilizer
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/ref/PauliStabilizer
Keywords: [stabilizer, Clifford, tableau, Aaronson-Gottesman, Pauli, quantum error correction, graph state, random Clifford, Mallows sampler, SymPhase]
SeeAlso: [StabilizerFrame, StabilizerStateQ, GraphState, LocalComplement, CliffordChannel, QuantumState, QuantumOperator, QuantumCircuitOperator, QuantumMeasurementOperator, QuantumChannel, PauliMatrix, KroneckerProduct]
RelatedGuides: [WolframQuantumComputationFramework]
---

## Usage

<code>[PauliStabilizer]()[]</code> represents an empty stabilizer state, equivalent to the single-qubit register $|0\rangle$.

<code>[PauliStabilizer]()[*n*]</code> represents the *n*-qubit $|0\cdots 0\rangle$ register state, with stabilizers $Z_{i}$.

<code>[PauliStabilizer]()[{$stab_{1}$, $stab_{2}$, …}]</code> constructs a stabilizer state from a list of Pauli strings.

<code>[PauliStabilizer]()[{*stabs*}, {*destabs*}]</code> constructs a stabilizer state with explicit stabilizers and destabilizers.

<code>[PauliStabilizer]()[*name*]</code> returns a named stabilizer state. Names: `"5QubitCode"`, `"SteaneCode"`, `"9QubitCode"`, their `"·1"` siblings, and `"Random"`.

<code>[PauliStabilizer]()["Random", *n*]</code> returns a uniformly random *n*-qubit Clifford state via the Mallows sampler.

<code>[PauliStabilizer]()[*qs*]</code> converts a [QuantumState](), [QuantumOperator](), or [QuantumCircuitOperator]() to a stabilizer state (Clifford only).

<code>*ps*[*method*,…]</code> applies a dispatched gate, measurement, or property accessor.

<code>*ps*[$ps_{2}$]</code> composes two stabilizer states (apply $\mathit{ps}_{2}$ first, then *ps*).

## Details & Options

- A [PauliStabilizer]() encodes an *n*-qubit pure stabilizer state by a tableau of Pauli generators plus signs, requiring $O(n^{2})$ memory.

- The internal representation is an Association `<|"Tableau"->t,"Signs"->s|>` where *t* is a rank-3 binary array of shape $\{2,\, n,\, 2n\}$ (X-block, Z-block; *n* qubits; $2n$ rows: destabilizers then stabilizers) and *s* is a length-$2n$ list of $\pm 1$ signs.

- Equivalent serialization keys are also accepted: `"Phase"` (= $1-Signs/2$) and `"Matrix"` (the flat $2n\times 2n$ binary tableau).

- Stabilizer rows are the *n* parity-check generators of the codespace. Destabilizer rows complete the symplectic basis and are required for measurement and `"Dagger"`. String-list constructors auto-generate destabilizers via `Reverse`; supply them explicitly with the two-list form to avoid singular fixtures.

- Named codes follow the convention that `"5QubitCode"` encodes $|0_{L}\rangle$ (the +1 eigenstate of logical $Z$), and the `"·1"` sibling encodes $|1_{L}\rangle$.

- <code>[PauliStabilizer]()["Random", *n*]</code> uses the Bravyi-Maslov / Koenig-Smolin (arXiv:1406.2170) Mallows sampler for uniform sampling over the $n$-qubit Clifford group modulo Paulis.

- A stabilizer state *ps* supports two access patterns: properties via <code>*ps*["Property"]</code> (e.g. `"Stabilizers"`, `"Tableau"`, `"State"`) and methods via <code>*ps*["Method",…]</code> (e.g. Clifford gates, measurement, `"InnerProduct"`). The full list of accepted strings is available via <code>*ps*["Properties"]</code>.

- Clifford gate updates run in $O(n^{2})$ time using the Aaronson-Gottesman tableau (arXiv:quant-ph/0406196). Single-qubit gates: `"H"`, `"S"`, $"S"^{\dagger }$, `"X"`, `"Y"`, `"Z"`, `"V"`. Two-qubit gates: `"CNOT"` (or `"CX"`), `"CZ"`, `"SWAP"`.

- Non-Clifford gates $"P"[\theta ]$, `"T"`, and $"T"^{\dagger }$ return a [StabilizerFrame]() (a coherent superposition of stabilizer states) rather than a [PauliStabilizer](). The frame closes under further Clifford gates.

- Symbolic measurement (Fang-Ying 2023, arXiv:2311.03906) via `"SymbolicMeasure"` allocates a fresh $s[k]$ outcome symbol per non-deterministic measurement, yielding a single PauliStabilizer with symbolic phase. Resolve via `"SubstituteOutcomes"` or `"SampleOutcomes"`.

- The option `Method` on <code>*ps*["InnerProduct", *other*]</code> selects between exact state-vector materialization (`"Direct"`, default) and the closed-form Garcia-Markov-Cross algorithm (`"ClosedForm"`, arXiv:1210.6646).

| Option | Default value | Description |
|---|---|---|
| `"Method"` | `"Direct"` | the algorithm to use when computing the inner product |

- Possible `"Method"` values:

| Value | Default | Description |
|---|---|---|
| `"Direct"` | `default` | materialize state vectors and compute exactly ($O(2^{n})$, recovers complex phase) |
| `"ClosedForm"` |  | closed-form $O(n^{3})$ magnitude only, no complex phase (concrete signs only) |

## Basic Examples

The empty form constructs the single-qubit "|0〉" register:

```wl
PauliStabilizer[]
```

---

For a multi-qubit register, pass the qubit count:

```wl
PauliStabilizer[3]
```

---

Construct a Bell state from two Pauli-string stabilizers:

```wl
PauliStabilizer[{"XX", "ZZ"}]
```

---

Construct a named stabilizer code:

```wl
PauliStabilizer["5QubitCode"]
```

---

Sample a random Clifford state via the Mallows sampler:

```wl
Module[{}, SeedRandom[42]; PauliStabilizer["Random", 3]]
```

---

Build a 3-qubit GHZ state by applying Clifford gates:

```wl
PauliStabilizer[3]["H", 1]["CNOT", 1, 2]["CNOT", 2, 3]
```

## Scope

### Constructors

The empty form returns the single-qubit register state:

```wl
PauliStabilizer[]
```

---

Specify a number of qubits to get an n-qubit register:

```wl
PauliStabilizer[5]
```

---

Build a state from a list of Pauli strings:

```wl
PauliStabilizer[{"XX", "ZZ"}]
```

---

Inspect the stabilizers via the property accessor:

```wl
PauliStabilizer[{"XX", "ZZ"}]["Stabilizers"]
```

---

Provide explicit destabilizers as a second list:

```wl
PauliStabilizer[{"XX", "ZZ"}, {"IX", "ZI"}]
```

---

Construct from a QuantumState via 4^n Pauli tomography:

```wl
PauliStabilizer[QuantumState[{1, 0, 0, 1}/Sqrt[2]]]
```

---

Construct from a QuantumOperator that is Clifford:

```wl
PauliStabilizer[QuantumOperator["H", {1}]]
```

---

Fold a QuantumCircuitOperator over the n-qubit register:

```wl
PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]
```

### Named stabilizer codes

The 5-qubit perfect code:

```wl
PauliStabilizer["5QubitCode"]
```

---

The Steane CSS code on 7 qubits:

```wl
PauliStabilizer["SteaneCode"]
```

---

The 9-qubit Shor code:

```wl
PauliStabilizer["9QubitCode"]
```

---

The "·1" siblings encode the |1〉 logical (sign-flipped Z\[Macron]):

```wl
PauliStabilizer["5QubitCode1"]
```

---

Read off the sign-flipped logical stabilizer:

```wl
Last[PauliStabilizer["5QubitCode1"]["Stabilizers"]]
```

### Random Clifford states

Sample a random 4-qubit Clifford state via the Mallows sampler:

```wl
Module[{}, SeedRandom[42]; PauliStabilizer["Random", 4]]
```

---

Read off the random stabilizers:

```wl
Module[{}, SeedRandom[42]; PauliStabilizer["Random", 4]["Stabilizers"]]
```

---

The discrete spectrum of squared overlaps for 3-qubit pairs:

```wl
Module[{products}, SeedRandom[42]; products = Table[Abs[PauliStabilizer["Random", 3]["InnerProduct", PauliStabilizer["Random", 3]]]^2, {200}]; Tally[Round[products, 1/8]]]
```

### Properties and tableau views

List the available property names:

```wl
Length[PauliStabilizer[2]["Properties"]]
```

---

The first eight property names:

```wl
Take[PauliStabilizer[2]["Properties"], 8]
```

---

The rank-3 binary tableau dimensions:

```wl
Dimensions[PauliStabilizer["5QubitCode"]["Tableau"]]
```

---

Stabilizers and destabilizers separately:

```wl
Module[{ps}, ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]; {ps["Stabilizers"], ps["Destabilizers"]}]
```

---

Materialize a dense state vector:

```wl
Normal[PauliStabilizer[{"XX", "ZZ"}]["State"]["StateVector"]]
```

### Single-qubit Clifford gates

Apply Hadamard:

```wl
PauliStabilizer[2]["H", 1]["Stabilizers"]
```

---

Apply S, X, Y, Z, V to qubit 1:

```wl
Module[{ps}, ps = PauliStabilizer[2]; AssociationThread[{"S", "X", "Y", "Z", "V"} -> {ps["S", 1]["Stabilizers"], ps["X", 1]["Stabilizers"], ps["Y", 1]["Stabilizers"], ps["Z", 1]["Stabilizers"], ps["V", 1]["Stabilizers"]}]]
```

### Two-qubit Clifford gates

Build a Bell state with H + CNOT:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
```

---

CZ acts symmetrically on both qubits:

```wl
PauliStabilizer[2]["H", 1]["H", 2]["CZ", 1, 2]["Stabilizers"]
```

---

SWAP via column permutation:

```wl
PauliStabilizer[3]["X", 1]["SWAP", 1, 3]["Stabilizers"]
```

### Permutation and padding

Permute qubit labels with a Cycles spec:

```wl
PauliStabilizer[3]["X", 1]["PermuteQudits", Cycles[{{1, 3}}]]["Stabilizers"]
```

---

Pad to more qubits via tensor product with identity:

```wl
PauliStabilizer[{"X"}]["PadRight", 3]["Stabilizers"]
```

### Dagger / Inverse

Compute the inverse Clifford via the symplectic-matrix dagger:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Dagger"]["Stabilizers"]
```

### Measurement

Single-qubit Z-basis measurement returns an Association keyed by outcome:

```wl
Keys[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", 1]]
```

---

Multi-qubit measurement via a list of qubits:

```wl
Keys[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", {1, 2}]]
```

---

Measure an arbitrary Pauli-string observable:

```wl
Keys[PauliStabilizer["5QubitCode"]["M", "XZZXI"]]
```

---

Symbolic measurement embeds a fresh outcome symbol in the phase:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["Stabilizers"]
```

---

Substitute the symbol back to recover concrete branches:

```wl
Module[{psSym}, psSym = PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]; psSym["SubstituteOutcomes", s[_] -> 0]["Stabilizers"]]
```

### Inner product and expectation

Self-inner-product is 1:

```wl
Module[{bell}, bell = PauliStabilizer[{"XX", "ZZ"}]; bell["InnerProduct", bell]]
```

---

Sign-flipped Bell states are orthogonal:

```wl
PauliStabilizer[{"XX", "ZZ"}]["InnerProduct", PauliStabilizer[{"-XX", "ZZ"}]]
```

---

Expectation of a stabilizer-group element returns +1:

```wl
PauliStabilizer[{"XX", "ZZ"}]["Expectation", "XX"]
```

---

An anticommuting Pauli gives expectation 0:

```wl
PauliStabilizer[{"XX", "ZZ"}]["Expectation", "XI"]
```

### Composition

Compose two stabilizer states (apply second, then first):

```wl
Module[{ps1, ps2}, ps1 = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1}]]; ps2 = PauliStabilizer[QuantumCircuitOperator[{"S" -> 1}]]; ps1[ps2]["Stabilizers"]]
```

---

Apply a QuantumOperator via UpValue:

```wl
Module[{qo}, qo = QuantumOperator["H", {1}]; qo[PauliStabilizer[1]]["Stabilizers"]]
```

### Round-trips and interop

Convert a stabilizer state to a QuantumState:

```wl
Normal[QuantumState[PauliStabilizer[{"XX", "ZZ"}]]["StateVector"]]
```

---

Apply a Pauli QuantumMeasurementOperator natively (fast path):

```wl
Module[{bell}, bell = PauliStabilizer[{"XX", "ZZ"}]; Keys[QuantumMeasurementOperator["XX", {1, 2}][bell]]]
```

## Generalizations & Extensions

### Non-Clifford promotion to StabilizerFrame

Applying T promotes a stabilizer state to a StabilizerFrame:

```wl
Module[{psT}, psT = PauliStabilizer[1]["T", 1]; Head[psT]]
```

---

Successive non-Clifford gates double the frame size:

```wl
PauliStabilizer[1]["T", 1]["T", 1]["Length"]
```

## Options

### Method

The default Method is "Direct" (state-vector materialization):

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["InnerProduct", PauliStabilizer[2]]
```

---

Explicit "Direct" recovers the full complex amplitude:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["InnerProduct", PauliStabilizer[2], Method -> "Direct"]
```

---

"ClosedForm" returns the magnitude in O(n³) (Garcia-Markov-Cross 2012):

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["InnerProduct", PauliStabilizer[2], Method -> "ClosedForm"]
```

### "H"

Hadamard maps X $\to$ Z:

```wl
PauliStabilizer[{"X"}]["H", 1]["Stabilizers"]
```

### "S"

S maps X $\to$ Y:

```wl
PauliStabilizer[{"X"}]["S", 1]["Stabilizers"]
```

### SuperDagger["S"]

$S^\dagger$ maps Y $\to$ X:

```wl
PauliStabilizer[{"Y"}][SuperDagger["S"], 1]["Stabilizers"]
```

### "X"

X anticommutes with Z:

```wl
PauliStabilizer[{"Z"}]["X", 1]["Stabilizers"]
```

### "Y"

Y anticommutes with Z:

```wl
PauliStabilizer[{"Z"}]["Y", 1]["Stabilizers"]
```

### "Z"

Z anticommutes with X:

```wl
PauliStabilizer[{"X"}]["Z", 1]["Stabilizers"]
```

### "V"

V = $\sqrt{X}$ maps Z $\to$ -Y:

```wl
PauliStabilizer[{"Z"}]["V", 1]["Stabilizers"]
```

### SuperDagger["V"]

$V^\dagger$ maps Z $\to$ Y:

```wl
PauliStabilizer[{"Z"}][SuperDagger["V"], 1]["Stabilizers"]
```

### "CNOT"

Build a Bell state with H + CNOT:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Stabilizers"]
```

### "CX"

CX is an alias for CNOT:

```wl
PauliStabilizer[2]["H", 1]["CX", 1, 2]["Stabilizers"]
```

### "CZ"

CZ propagates X-on-target into Z-on-control:

```wl
PauliStabilizer[{"IX", "XI"}]["CZ", 1, 2]["Stabilizers"]
```

### "SWAP"

SWAP exchanges qubit labels:

```wl
PauliStabilizer[{"IX", "ZI"}]["SWAP", 1, 2]["Stabilizers"]
```

### "Permute"

Permute swaps tableau rows:

```wl
PauliStabilizer[{"XI", "IZ"}]["Permute", Cycles[{{1, 2}}]]["Stabilizers"]
```

### "PermuteQudits"

PermuteQudits swaps qubit columns:

```wl
PauliStabilizer[{"XI", "IZ"}]["PermuteQudits", Cycles[{{1, 2}}]]["Stabilizers"]
```

### "PadRight"

Tensor identity on the right:

```wl
PauliStabilizer[{"X"}]["PadRight", 3]["Stabilizers"]
```

### "PadLeft"

Tensor identity on the left:

```wl
PauliStabilizer[{"X"}]["PadLeft", 3]["Stabilizers"]
```

### "Dagger"

F⪡2⪢ symplectic inverse:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Dagger"]["Stabilizers"]
```

### "Inverse"

Inverse is an alias for Dagger:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Inverse"]["Stabilizers"]
```

### "P"[$\theta$]

Phase rotation P[$\theta$] returns a StabilizerFrame:

```wl
Head[PauliStabilizer[1]["P"[Pi/4], 1]]
```

### "T"

T = P[$\pi$/2] returns a StabilizerFrame:

```wl
Head[PauliStabilizer[1]["T", 1]]
```

### SuperDagger["T"]

$T^\dagger$ returns a StabilizerFrame:

```wl
Head[PauliStabilizer[1][SuperDagger["T"], 1]]
```

### "M" (single qubit)

Z-basis measurement on qubit q:

```wl
Keys[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", 1]]
```

### "M" (qubit list)

Joint Z measurement on a list of qubits:

```wl
Keys[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", {1, 2}]]
```

### "M" (Pauli string)

Measurement of an arbitrary Pauli observable:

```wl
Keys[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["M", "ZI"]]
```

### "SymbolicMeasure"

Allocates a fresh s[k] outcome symbol:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["Phase"]
```

### "SubstituteOutcomes"

Resolve outcome symbols to concrete values:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["SubstituteOutcomes", s[_] -> 1]["Phase"]
```

### "SampleOutcomes"

Draw n random samples by independent substitution:

```wl
Length[PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["SymbolicMeasure", 1]["SampleOutcomes", 4]]
```

### "RowSum"

AG row-multiplication primitive:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["RowSum", 1, 2]["Stabilizers"]
```

### "InnerProduct"

Inner product against another stabilizer state:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["InnerProduct", PauliStabilizer[2]]
```

### "Expectation"

Expectation value of a Pauli observable:

```wl
PauliStabilizer[2]["H", 1]["CNOT", 1, 2]["Expectation", "XX"]
```

### Property accessors

PauliStabilizer exposes ~30 string-keyed property accessors. The full list is available via <code>*ps*["Properties"]</code>. A representative subset:

```wl
Module[{ps}, ps = PauliStabilizer[{"XX", "ZZ"}]; {ps["Qubits"], ps["GeneratorCount"], ps["Stabilizers"], ps["Destabilizers"], Dimensions[ps["Tableau"]], Dimensions[ps["Matrix"]], Head[ps["State"]]}]
```

## Applications

### Stabilizer-state tomography

Recover a stabilizer state from a dense state vector via 4^n tomography and round-trip back:

```wl
Module[{qs, ps}, qs = QuantumState[{1, 0, 0, 1}/Sqrt[2]]; ps = PauliStabilizer[qs]; ps["Stabilizers"]]
```

### Random Clifford circuits and statistical fingerprint

Sample 200 random 5-qubit Cliffords and tabulate their stabilizer-weight distribution:

```wl
Module[{weights}, SeedRandom[42]; weights = Table[Total[(StringCount[#1, "X" | "Y" | "Z"] & ) /@ PauliStabilizer["Random", 5]["Stabilizers"]], {200}]; Sort[Tally[weights]]]
```

### Quantum error correction with the 5-qubit code

Encode |0_L〉 in the 5-qubit code, inject an X error on qubit 3, and observe a non-trivial syndrome:

```wl
Module[{ps0, psErr, syndrome}, ps0 = PauliStabilizer["5QubitCode"]; psErr = ps0["X", 3]; syndrome = Sign /@ Table[psErr["Expectation", g], {g, ps0["Stabilizers"]}]; syndrome]
```

## Properties & Relations

### State conversion preserves the stabilizer up to global phase

The QuantumState round-trip is exact for circuit-built fixtures:

```wl
Module[{ps, qs, ps2}, ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]; qs = ps["State"]; ps2 = PauliStabilizer[qs]; ps["Stabilizers"] === ps2["Stabilizers"]]
```

### Inner-product squared on a discrete spectrum

For random pairs the squared overlap is in {0, 1/2^k}:

```wl
Module[{vals}, SeedRandom[1]; vals = Table[Abs[PauliStabilizer["Random", 3]["InnerProduct", PauliStabilizer["Random", 3]]]^2, {50}]; Union[Round[vals, 1/8]]]
```

### StabilizerStateQ on a PauliStabilizer is True

The predicate accepts PauliStabilizer values:

```wl
StabilizerStateQ[PauliStabilizer[{"XX", "ZZ"}]]
```

### Bitflip channel as a tableau-mixture list

A BitFlip channel returns a probability-weighted mixture:

```wl
Length[QuantumChannel["BitFlip"[1/3], {1}][PauliStabilizer[1]]]
```

### Inverse of a Clifford state

Composing ps with its Dagger yields a tableau equal to the identity register up to row reordering:

```wl
Module[{ps, dag}, ps = PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]; dag = ps["Dagger"]; Sort[ps[dag]["Stabilizers"]]]
```

## Possible Issues

### PauliStabilizer::nonclifford

The PauliStabilizer tableau formalism tracks only Clifford gates. Applying a non-Clifford operator such as a continuous rotation "RX"[$\theta$] (via the public QuantumOperator-on-PauliStabilizer dispatch) cannot be folded into the tableau, emits PauliStabilizer::nonclifford, and returns $Failed. The fold short-circuits, so the message fires only once even in a multi-gate circuit:

```wl
QuantumOperator["RX"[Pi/3], {1}][PauliStabilizer[1]]
```

Resolution: evaluate non-Clifford operators via the dense state-vector path: apply the QuantumOperator to a QuantumState instead of a PauliStabilizer. This sidesteps the tableau and supports arbitrary unitaries:

```wl
QuantumOperator["RX"[Pi/3], {1}][QuantumState["0"]]
```

### PauliStabilizer::singular

PauliStabilizer[{p1, p2, ...}] from stabilizer rows alone fills the destabilizer rows with a Reverse-padded default, which may yield a symplectic matrix that is singular mod 2. Dagger requires the full tableau to be invertible, so this case emits PauliStabilizer::singular and returns $Failed:

```wl
PauliStabilizer[{"XX", "ZZ"}]["Dagger"]
```

Resolution: build the stabilizer from a Clifford circuit. The constructor PauliStabilizer[qco] runs the circuit on the all-zero state and sets up consistent destabilizers, so Dagger is well-defined:

```wl
PauliStabilizer[QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}]]["Dagger"]
```

### PauliStabilizer::expectationdim

ps["Expectation", P] requires the Pauli string P to have one letter per qubit. A length mismatch emits PauliStabilizer::expectationdim and returns $Failed:

```wl
PauliStabilizer[2]["Expectation", "X"]
```

Resolution: pass a Pauli string of length n where n is the qubit count. For a 2-qubit stabilizer, use a 2-letter string from the alphabet {I, X, Y, Z}:

```wl
PauliStabilizer[2]["Expectation", "XX"]
```

### PauliStabilizer::partition

Measurement and entropy queries require qubit indices in {1, ..., n}. Out-of-range indices emit PauliStabilizer::partition and return $Failed:

```wl
PauliStabilizer[2]["M", {1, 5}]
```

Resolution: use indices within {1, ..., n}. For a 2-qubit stabilizer, valid indices are 1 and 2:

```wl
PauliStabilizer[2]["M", {1, 2}]
```

### PauliStabilizer::nonpaulibasis

The PauliStabilizer fast-path for QuantumMeasurementOperator and QuantumChannel requires a Pauli-basis label. A non-Pauli basis (here a diagonal operator with eigenvalues {1, 2}) emits PauliStabilizer::nonpaulibasis and falls back to dense state-vector evaluation: the head of the result becomes QuantumMeasurement:

```wl
QuantumMeasurementOperator[QuantumOperator[DiagonalMatrix[{1, 2}]], {1}][PauliStabilizer[1]]
```

Resolution: use a Pauli-basis QuantumMeasurementOperator (a Pauli string label or one of "X", "Y", "Z"). The fast path applies the corresponding Pauli measurement directly to the tableau, keeping evaluation on the O(n^2) path and preserving the PauliStabilizer head:

```wl
QuantumMeasurementOperator["X", {1}][PauliStabilizer[1]]
```

## Neat Examples

### 50-qubit random Clifford tableau as MatrixPlot

Visualize the symplectic structure of a large random Clifford:

```wl
Module[{ps}, SeedRandom[1]; ps = PauliStabilizer["Random", 50]; MatrixPlot[ps["Matrix"], FrameTicks -> None, ImageSize -> 250]]
```
