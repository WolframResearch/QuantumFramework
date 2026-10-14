---
Template: Default
Title: Named Quantum Circuits and Operators
Context: Wolfram`QuantumFramework`
---

# Named Quantum Circuits and Operators

A single catalogue of every named construction `QuantumCircuitOperator` accepts, each shown as `QuantumCircuitOperator[{"Name"[...]}]["Diagram"]` with a one-line caption stating what that signature does, and an action cell on the canonical examples. The families are ordered by impact for a quantum audience, from complete algorithms down to the primitive gate set. Every cell is verified against the kernel.

**Tier 1 - Complete algorithms.**

## Grover

`QuantumCircuitOperator` accepts a library of named constructions. A bare string
builds the default instance; `"Name"[...]` passes arguments. This page walks the
**Grover** family through its main signatures, one line each, showing the
circuit diagram.

The names split into two layers: the diffusion block alone, and the full search
step $G = D\,O$ that wraps a diffusion block around an oracle.

### The diffusion block

The diffusion reflection on *n* wires, $I - 2|{+}\cdots{+}{-}\rangle\langle{+}\cdots{+}{-}|$: every input wire in $|+\rangle$ and the flag wire in $|-\rangle$:

```wl
QuantumCircuitOperator["GroverDiffusion"[4]]["Diagram"]
```

---

The block placed on a chosen (possibly non-contiguous) set of wires, marking onto the highest:

```wl
QuantumCircuitOperator["GroverDiffusion"[{1, 3}]]["Diagram"]
```

---

A supplied operator sets the marking gate: its wire (2) becomes the target, the other wires its controls:

```wl
QuantumCircuitOperator["GroverDiffusion"[3, QuantumOperator["Z", {2}]]]["Diagram"]
```

---

The same reflection built from zero-controls instead of an X-sandwich: an identical operator in a shallower circuit (depth 3 rather than 5):

```wl
QuantumCircuitOperator["GroverDiffusion0"[3]]["Diagram"]
```

---

The phase-oracle diffusion, $I - 2|{+}\cdots{+}\rangle\langle{+}\cdots{+}|$: a reflection about the all-plus state with no flag wire, distinct from `GroverDiffusion`:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion"[3]]["Diagram"]
```

---

The same phase reflection from zero-controls: an identical operator to `GroverPhaseDiffusion` in a shallower circuit:

```wl
QuantumCircuitOperator["GroverPhaseDiffusion0"[3]]["Diagram"]
```

---

The same block under its amplitude-named synonym, which builds only in call form (even `[]`), not as the bare string:

```wl
QuantumCircuitOperator["GroverAmplification"[3]]["Diagram"]
```

### The full search operator

The search step marking the states that satisfy a boolean formula:

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2] && ! \[FormalX][3]]]["Diagram"]
```

---

A variable-to-wire list fixing which wire each variable occupies (both this and the bare two-variable form use three wires, with wire 3 the target):

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2], {1, 2}]]["Diagram"]
```

---

A second argument that is an integer sets the oracle's target wire; choosing 1 collides with a variable, so the variables shift to 2,3 and the target lands on 1:

```wl
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2], 1]]["Diagram"]
```

---

A `BooleanFunction` object accepted directly in place of a symbolic formula:

```wl
QuantumCircuitOperator["Grover"[BooleanFunction[7, 2]]]["Diagram"]
```

---

A prebuilt phase oracle wrapped with the matching phase diffusion, no formula compilation (here `CZ`, which sign-flips the marked $|11\rangle$; pairing a phase oracle with the plain `Grover` instead would mismatch the diffusion):

```wl
QuantumCircuitOperator["GroverPhase"[QuantumOperator["CZ"]]]["Diagram"]
```

---

The phase-oracle search step, which marks by sign flip on the data wires and so needs no ancilla:

```wl
QuantumCircuitOperator["GroverPhase"[\[FormalX][1] && \[FormalX][2]]]["Diagram"]
```

---

The search step built with the shallower zero-control diffusion in place of the X-sandwich form (the same reflection):

```wl
QuantumCircuitOperator["Grover0"[\[FormalX][1] && \[FormalX][2]]]["Diagram"]
```

### Does it search?

Two phase-Grover iterations on the uniform superposition concentrate the probability on the marked $|111\rangle$:

```wl
Nest[QuantumCircuitOperator["GroverPhase"[\[FormalX][1] && \[FormalX][2] && \[FormalX][3]]], QuantumState["+++"], 2]["Probabilities"]
```

## Deutsch and Deutsch-Jozsa

Deutsch-Jozsa decides, for a function *promised* to be either constant or
balanced, which of the two it is in a single query. The named circuit assembles
the oracle, the surrounding Hadamard layers, and the final measurement. The
function argument does not enforce the promise: a formula or an integer indexes
a general truth table, and `Automatic` draws from all Boolean functions (for
three inputs only 72 of 256 are constant or balanced), so a random draw can be a
valid circuit on an invalid problem instance. The explicit examples below pick
balanced or constant functions; the `Automatic` example does not enforce the
promise. Each algorithm comes in a *boolean* form (bit-flip oracle plus a
$|-\rangle$ ancilla) and a *phase* form (phase oracle, no ancilla); both read the
input register out in the X basis.

### Deutsch: one qubit

The balanced case, with the oracle `f(x) = x`:

```wl
QuantumCircuitOperator[{"Deutsch"[\[FormalX][1]]}]["Diagram"]
```

---

The constant case, selected by the integer index (here `f ≡ 0`):

```wl
QuantumCircuitOperator[{"Deutsch"[1]}]["Diagram"]
```

---

The phase form, which drops the ancilla and measures in the X basis (width 2):

```wl
QuantumCircuitOperator[{"DeutschPhase"[\[FormalX][1]]}]["Diagram"]
```

### Deutsch-Jozsa: n qubits

A balanced function on three inputs, here the parity `x1 ⊻ x2 ⊻ x3`:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Xor[\[FormalX][1], \[FormalX][2], \[FormalX][3]]]}]["Diagram"]
```

---

A constant function on three inputs, selected by the integer index:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[1, 3]}]["Diagram"]
```

---

The decision, read out: the parity function sends all the weight onto the all-minus outcome, so the all-plus string that would certify a constant function has probability zero (the constant case puts probability 1 there instead):

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Xor[\[FormalX][1], \[FormalX][2], \[FormalX][3]]]}][]["Probabilities"]
```

---

A second argument fixes the number of input qubits when the function is left to a random draw:

```wl
QuantumCircuitOperator[{"DeutschJozsa"[Automatic, 3]}]["Diagram"]
```

---

The phase form on two inputs, here the balanced equivalence `¬(x1 ⊻ x2)`: phase oracle, no ancilla:

```wl
QuantumCircuitOperator[{"DeutschJozsaPhase"[! Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

### The oracles alone

The bit-flip oracle for a balanced function, flipping the ancilla where the function is true:

```wl
QuantumCircuitOperator[{"DeutschJozsaBooleanOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

The phase oracle for the same function, marking by sign with no ancilla:

```wl
QuantumCircuitOperator[{"DeutschJozsaPhaseOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

The constant oracle on three inputs, selected by the integer index:

```wl
QuantumCircuitOperator[{"DeutschJozsaBooleanOracle"[1, 3]}]["Diagram"]
```

## Simon's Algorithm

Simon's algorithm finds a hidden nonzero mask *s* of a function promised to be
two-to-one with that period: `f(x) = f(x')` exactly when `x' = x` or `x' = x ⊕ s`.
(The bare period `f(x) = f(x ⊕ s)` alone is weaker: a constant function satisfies
it for every *s*.) It uses exponentially fewer queries than any classical method.
The named circuit is the *sampling* step: oracle, Hadamard layers on the input
register, and a measurement that returns one string `y` with `y · s = 0`;
recovering *s* needs repeated samples and classical linear algebra on top. The
mask is given as a bit string, a bit list, or left to a random draw. One caveat
verified in the kernel: the generated function respects the period `f(x) = f(x ⊕ s)`
but is not guaranteed strictly two-to-one (the outputs are drawn at random), so
it realizes the period rather than the full Simon promise.

### Sampling circuit: `Simon`

The mask given as a bit string, on a 3-wire input register feeding a 3-wire output register:

```wl
QuantumCircuitOperator[{"Simon"["101"]}]["Diagram"]
```

---

The mask given as a bit list:

```wl
QuantumCircuitOperator[{"Simon"[{1, 0}]}]["Diagram"]
```

---

A random mask, defaulting to three bits (so the register is fixed at three input and three output wires):

```wl
QuantumCircuitOperator[{"Simon"[]}]["Diagram"]
```

### Oracle alone: `SimonOracle`

The two-register oracle for a chosen mask, mapping a 3-wire input register into a 3-wire output register:

```wl
QuantumCircuitOperator[{"SimonOracle"["110"]}]["Diagram"]
```

---

The same mask given as a bit list:

```wl
QuantumCircuitOperator[{"SimonOracle"[{1, 1, 0}]}]["Diagram"]
```

## Bernstein-Vazirani

Bernstein-Vazirani recovers a hidden bit string *a* in a single query, where the
oracle computes the dot product `f(x) = a · x (mod 2)`. The named circuit builds
the oracle (a `CNOT` from each wire whose secret bit is 1 onto the ancilla), the
surrounding Hadamard layers, and the measurement, whose outcome reads the secret
back directly.

### The algorithm

The secret given as a bit string:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["101"]}]["Diagram"]
```

---

Running the full circuit, a single measurement returns the secret `101` with certainty:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["101"]}][]["Probabilities"]
```

---

The secret given as a bit list:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"[{1, 1, 0, 1}]}]["Diagram"]
```

---

A second argument shifts every wire up by one, leaving wire 1 idle with the register on wires 2-5 and the ancilla on 6:

```wl
QuantumCircuitOperator[{"BernsteinVazirani"["1101", 1]}]["Diagram"]
```

### The oracle alone

The oracle for a chosen secret: one `CNOT` onto the ancilla from each wire whose bit is 1:

```wl
QuantumCircuitOperator[{"BernsteinVaziraniOracle"["1011"]}]["Diagram"]
```

---

The same oracle from a bit list:

```wl
QuantumCircuitOperator[{"BernsteinVaziraniOracle"[{1, 0, 1, 1}]}]["Diagram"]
```

## Fourier Transform and Phase Estimation

The quantum Fourier transform moves a register between the computational and the
Fourier basis; phase estimation uses its inverse to read the eigenphase of a
unitary into a counting register. The Fourier transform and its inverse also
exist as single `QuantumOperator` matrices, so under the circuit form here the
diagram is the gate *decomposition*, Hadamards and controlled phase rotations,
not an atomic box.

### Quantum Fourier transform

The transform on three qubits: a Hadamard and a ladder of controlled phase shifts on each wire, then the wire-reversing swap:

```wl
QuantumCircuitOperator[{"Fourier"[3]}]["Diagram"]
```

---

A second argument offsets the transform up by that many leading wires:

```wl
QuantumCircuitOperator[{"Fourier"[3, 1]}]["Diagram"]
```

---

The inverse transform, the dagger of the forward one:

```wl
QuantumCircuitOperator[{"InverseFourier"[3]}]["Diagram"]
```

### Phase estimation

The eigenphase of the `T` gate read into a three-qubit counting register; the constructor prepares the target in $|1\rangle$, selecting its eigenphase 1/8 (both `T` and `S` have eigenphase 0 on $|0\rangle$), and ends in an inverse Fourier transform and measurement:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}]["Diagram"]
```

---

Running it, the counting register reads `001` = 1/8, the eigenphase of `T`, with certainty:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3]}][]["Probabilities"]
```

---

The same for the `S` gate (phase 1/4), a different target unitary:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["S"], 3]}]["Diagram"]
```

---

A larger counting register, five qubits: more resolution for a general phase, though the dyadic 1/8 here is already exact with three:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 5]}]["Diagram"]
```

---

The `PowerExpand` option replaces each controlled `U^(2^k)` by `2^k` repeated controlled-`U` gates:

```wl
QuantumCircuitOperator[{"PhaseEstimation"[QuantumOperator["T"], 3, "PowerExpand" -> True]}]["Diagram"]
```

**Tier 2 - Oracles.**

## Oracles

An oracle is the black box the search and query algorithms compile to: it marks
the basis states that satisfy a Boolean condition. The named oracles differ in
*how* they mark: flip a target qubit, flip a sign, rotate a target, or load a
whole function of the input. The first three take a Boolean formula and
variable-placement arguments; `GrayOracle` instead takes a function or angle
list, a precision, a target, and a rotation axis.

### Flag oracle: `BooleanOracle`

Flipping the target on the states that satisfy a formula, here a single multi-controlled gate for the one product term:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && \[FormalX][2] && \[FormalX][3]]}]["Diagram"]
```

---

A sum of exclusive-or terms compiles to one controlled gate per term (the ESOP form), so the gate count tracks the formula:

```wl
QuantumCircuitOperator[{"BooleanOracle"[Xor[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

A variable-to-wire list placing the variables on wires 1 and 3, leaving wire 2 idle and the target on wire 4:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] || \[FormalX][2], {1, 3}]}]["Diagram"]
```

---

An integer in the second slot sets the target wire, here wire 5, leaving the variables on 1,2 and wires 3,4 idle:

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && ! \[FormalX][2], 5]}]["Diagram"]
```

---

An operator in the second slot becomes the marking gate, turning the oracle into a controlled-*U* that applies an arbitrary unitary to the target on the satisfying states:

```wl
QuantumCircuitOperator[{"BooleanOracle"[! \[FormalX][1] && \[FormalX][2], QuantumOperator["RandomUnitary", "Label" -> "U"]]}]["Diagram"]
```

### Phase oracle: `PhaseOracle`

Flipping the sign of the satisfying states with no ancilla (width 2, not 3):

```wl
QuantumCircuitOperator[{"PhaseOracle"[Nand[\[FormalX][1], \[FormalX][2]]]}]["Diagram"]
```

---

On the uniform state $|{+}{+}\rangle$, the oracle for $x_1 \wedge x_2$ flips the sign of the marked $|11\rangle$ amplitude alone:

```wl
Normal[QuantumCircuitOperator[{"PhaseOracle"[\[FormalX][1] && \[FormalX][2]]}][QuantumState["++"]]["StateVector"]]
```

---

A variable-to-wire list fixing the wire of each variable:

```wl
QuantumCircuitOperator[{"PhaseOracle"[Implies[\[FormalX][1], \[FormalX][2]], {1, 2}]}]["Diagram"]
```

---

An integer in the second slot offsets every wire up by *m*, here placing the variables on wires 3,4 with 1,2 idle:

```wl
QuantumCircuitOperator[{"PhaseOracle"[\[FormalX][1] && \[FormalX][2], 2]}]["Diagram"]
```

### Rotation oracle: `BooleanOracleR`

Marking the states that satisfy a formula by a target rotation (`RZ` by π) rather than a flip:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[\[FormalX][1] && \[FormalX][2]]}]["Diagram"]
```

---

A variable-to-wire list fixing each variable's wire:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[Xor[\[FormalX][1], \[FormalX][2]], {1, 2}]}]["Diagram"]
```

---

A fourth argument sets the axis and angle of the marking rotation:

```wl
QuantumCircuitOperator[{"BooleanOracleR"[! \[FormalX][1] && ! \[FormalX][2], Automatic, Automatic, "RY"[Pi/2]]}]["Diagram"]
```

### Gray-code angle oracle: `GrayOracle`

A pure function loaded as target rotations $\theta_x = f(x)/2^{\text{prec}}$, sampled over the input values (default precision 4, so $R_Y(x^2/16)$ here):

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &)]}]["Diagram"]
```

---

An explicit list of angles in place of a function, padded with zeros to the $2^{\text{prec}}$ the default four input wires expect (16 angles):

```wl
QuantumCircuitOperator[{"GrayOracle"[{0, 1, 2, 3}]}]["Diagram"]
```

---

A precision argument sets the number of input wires, here 3 instead of 4, so the angles become $f(x)/8$ and the ladder halves:

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &), 3]}]["Diagram"]
```

---

A rotation-axis argument sets which rotation the angles drive, here `RZ`:

```wl
QuantumCircuitOperator[{"GrayOracle"[(#^2 &), 3, Automatic, "RZ"]}]["Diagram"]
```

**Tier 3 - Composition and simulation.**

## Composition and Simulation

These constructions take operators or whole circuits as input and assemble a
larger circuit: controlling a circuit on extra wires, multiplexing a choice of
gate on a control register, routing two unitaries in indefinite order, solving a
linear system, and compiling a Hamiltonian into a product-formula circuit.

### Controlling a circuit: `Controlled`

A whole circuit controlled on an extra wire, here a two-qubit Fourier transform on wires 2,3 fired by wire 1 (the transform is offset to 2,3 so it does not collide with the control):

```wl
QuantumCircuitOperator[{"Controlled"[QuantumCircuitOperator["Fourier"[2, 1]], {1}]}]["Diagram"]
```

---

Two control wires (3 and 4) gating a GHZ preparation on wires 1,2, so the pair is prepared only when both controls fire:

```wl
QuantumCircuitOperator[{"Controlled"[QuantumCircuitOperator["GHZ"[2]], {3, 4}]}]["Diagram"]
```

### Multiplexing a gate: `Multiplexer`

A gate chosen by the control register, one target gate per control value, here three gates on a two-wire control (shown as its gate decomposition):

```wl
QuantumCircuitOperator[{"Multiplexer"["X", "Y", "Z"]}]["Diagram"]
```

---

Two gates on a single control wire:

```wl
QuantumCircuitOperator[{"Multiplexer"["X", "H"]}]["Diagram"]
```

---

Its matrix is block-diagonal: $X$ on the control-0 block, $H$ on the control-1 block:

```wl
Normal[QuantumCircuitOperator[{"Multiplexer"["X", "H"]}]["Matrix"]]
```

### Indefinite causal order: `Switch`

The quantum switch of two unitaries, applying them in an order set by an ancilla; preparing that ancilla in superposition puts the order in superposition (shown as its gate decomposition):

```wl
QuantumCircuitOperator[{"Switch"[QuantumOperator["X"], QuantumOperator["Z"]]}]["Diagram"]
```

### Linear systems: `ControlledMultiplexer`

The variational linear-solver construction for `A x = b` (coefficient-state preparation, an ansatz, multiplexing, and an ancilla contraction, returned before the fidelity cost and optimization), here a diagonal `A` and a two-component `b`:

```wl
QuantumCircuitOperator[{"ControlledMultiplexer"[{{1, 0}, {0, 2}}, {1, 1}]}]["Diagram"]
```

### Hamiltonian simulation: `Trotterization`

A second-order Trotter step for the sum `X + Z`, symmetrically ordered:

```wl
QuantumCircuitOperator[{"Trotterization"[{"X", "Z"}, 2]}]["Diagram"]
```

---

A first-order step for `X + Y + Z` repeated three times, which shrinks the step error:

```wl
QuantumCircuitOperator[{"Trotterization"[{"X", "Y", "Z"}, 1, 3]}]["Diagram"]
```

**Tier 4 - Beyond unitary gates.**

## Nonclassicality Tests

These circuits stage the standard tests that separate quantum correlations from
classical ones: a spatial Bell test, its temporal analogue, and a phase-space
version. Each is parametrized by a single measurement angle, with the Alice,
Bob, and Charlie roles laid out as labelled sub-circuits, barrier-separated in
the spatial and temporal tests. The page builds the measurement circuits; the
correlation value and its classical bound (for CHSH, the $2\sqrt2$ Tsirelson
violation of the bound 2) are a classical post-processing step on the outcome
statistics, not computed here.

### Spatial Bell test: `CHSH`

The CHSH Bell test: a shared entangled pair, Charlie's preparation, and Alice and Bob's measurements, with Bob's rotation carrying the angle, ending in four readouts:

```wl
QuantumCircuitOperator[{"CHSH"[]}]["Diagram"]
```

---

The measurement angle is the single knob; setting it changes Bob's rotation:

```wl
QuantumCircuitOperator[{"CHSH"[Pi/3]}]["Diagram"]
```

### Temporal Bell test: `LeggettGarg`

The temporal analogue of CHSH, correlating measurements of one system at different times rather than two systems in space:

```wl
QuantumCircuitOperator[{"LeggettGarg"[]}]["Diagram"]
```

### Phase-space Bell test: `WignerCHSH`

The CHSH test built on the Wigner minimal informationally complete basis, a discrete phase-space *representation* of the qubit, rather than a physical measurement of phase space:

```wl
QuantumCircuitOperator[{"WignerCHSH"[]}]["Diagram"]
```

## Non-Unitary Elements

Not every circuit element is a unitary gate. A circuit can also measure a wire,
reset it, trace it out, or discard it. Each of these is a named element whose
action on a state is described below, not just drawn as a diagram.

### Measurement: `Measurement`

An unrecorded projective measurement on wire 1; applied to $|+\rangle$ it returns the maximally mixed state $I/2$:

```wl
QuantumCircuitOperator[{"Measurement"[1]}]["Diagram"]
```

---

Acting on $|+\rangle$, the density matrix is $I/2$: the coherence is gone:

```wl
Normal[QuantumCircuitOperator[{"Measurement"[1]}][QuantumState["+"]]["DensityMatrix"]]
```

### Reset: `Reset`

Resets a wire to $|0\rangle$ regardless of its input, so $|1\rangle$ becomes $|0\rangle\langle0|$:

```wl
QuantumCircuitOperator[{"Reset"[]}]["Diagram"]
```

### Partial trace: `Trace`

The partial trace over a wire, keeping the reduced state of the rest; tracing one half of a Bell pair leaves the other maximally mixed:

```wl
QuantumCircuitOperator[{"Trace"[2]}]["Diagram"]
```

### Discard: `Discard`

The effect $\langle0| + \langle1|$ that contracts a wire's amplitudes, not a trace-preserving removal: it sends a Bell pair's kept wire to $|+\rangle$ but maps $|-\rangle$ to zero, so it carries a normalization or postselection caveat:

```wl
QuantumCircuitOperator[{"Discard"[2]}]["Diagram"]
```

### Channel: `Channel`

A completely positive map applied directly in the circuit (the strength nests inside the name, `"Depolarizing"[p]`); the depolarizing channel at 3/10 partially mixes a pure input:

```wl
QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}]["Diagram"]
```

---

Acting on $|+\rangle$ it damps the coherence rather than erasing it, giving $\begin{pmatrix}1/2 & 7/20\\ 7/20 & 1/2\end{pmatrix}$:

```wl
Normal[QuantumCircuitOperator[{"Channel"["Depolarizing"[3/10]]}][QuantumState["+"]]["DensityMatrix"]]
```

### Generators and the diagonal map

Two of these names are the one open-system generator in its two conventions.
`Hamiltonian[H, {jumps}]` is the form `QuantumEvolve` consumes: the commutator
$[H, \cdot]$ for a closed system, and $i\mathcal L$, with $\mathcal L$ the
Lindblad generator, once jump operators are given. `Liouvillian[H, {jumps}]` is
$\mathcal L$ itself, the bare generator of $\dot\rho = \mathcal L[\rho]$, equal
to $-i\,$`Hamiltonian[H, {jumps}]`; reach for it when $\mathcal L$ is wanted (its
spectrum, say) rather than for `QuantumEvolve`, which exponentiates $-i$ times
its argument. Evolving with `Hamiltonian[H, {jumps}]` and the separate
Hamiltonian-and-jump-operator form `QuantumEvolve[H, {jumps}, ...]` give the same
result. A third, `Decohere`, is not a generator at all but a
dimension-reducing *representation* map: the $2\times4$ operator
$\{\{1,0,0,0\},\{0,0,0,1\}\}$ that pulls the diagonal from a doubled
representation (4 amplitudes in, 2 out). Because its input and output dimensions
differ it cannot serve as an evolution generator, and it is not a
trace-preserving dephasing channel.

## Tensor-Network Atoms

These are the wire-level building blocks of the ZX and tensor-network calculus:
bending a wire into a state or an effect, crossing two wires, copying a basis
state, and the phase spider. Each caption states what the atom's matrix does,
verified in the kernel.

### Cup: `Cup`

The cup bends a wire from two outputs, giving the unnormalized Bell state $|00\rangle + |11\rangle$:

```wl
QuantumCircuitOperator[{"Cup"[2]}]["Diagram"]
```

---

It has no inputs; its output state has amplitudes $(1, 0, 0, 1)$:

```wl
Normal[QuantumCircuitOperator[{"Cup"[2]}][]["StateVector"]]
```

### Cap: `Cap`

The cap is the matching effect $\langle00| + \langle11|$, a wire bent into two inputs that contracts a pair:

```wl
QuantumCircuitOperator[{"Cap"[2]}]["Diagram"]
```

### Braid: `Braid`

The braid crosses two wires: a swap carrying a $-1$ phase on $|11\rangle$:

```wl
QuantumCircuitOperator[{"Braid"[2]}]["Diagram"]
```

### Copy: `Copy`

The copy spider clones a basis state onto two wires, $|0\rangle \mapsto |00\rangle$ and $|1\rangle \mapsto |11\rangle$:

```wl
QuantumCircuitOperator[{"Copy"[2]}]["Diagram"]
```

### Phase spider: `ZSpider`

A one-in one-out Z-spider carries a phase on the $|1\rangle$ leg; at phase $\pi$ it is $\mathrm{diag}(1, -1) = Z$:

```wl
QuantumCircuitOperator[{"ZSpider"[Pi]}]["Diagram"]
```

---

At phase $\pi/2$ the same spider is $\mathrm{diag}(1, i) = S$:

```wl
QuantumCircuitOperator[{"ZSpider"[Pi/2]}]["Diagram"]
```

## Qudits and Angular Momentum

Beyond qubits, the library carries the clock-and-shift operators of a
$d$-level qudit, the qudit entangling gate, and angular-momentum operators whose
dimension is set by the spin. Each caption states the matrix or the action.

### Clock shift: `Shift`

The clock shift on a three-level qudit, cycling $|j\rangle \to |j+1 \bmod 3\rangle$:

```wl
QuantumCircuitOperator[{"Shift"[3]}]["Diagram"]
```

---

Acting on $|0\rangle$ it advances the qutrit to $|1\rangle$, amplitudes $(0, 1, 0)$:

```wl
Normal[QuantumCircuitOperator[{"Shift"[3]}][QuantumState[{1, 0, 0}, 3]]["StateVector"]]
```

### Clock phase: `ShiftPhase`

The companion clock phase, a diagonal of the cube roots of unity on a qutrit:

```wl
QuantumCircuitOperator[{"ShiftPhase"[3]}]["Diagram"]
```

### Qudit sum gate: `SUM`

The qudit analogue of CNOT, adding the control to the target modulo 3:

```wl
QuantumCircuitOperator[{"SUM"[3]}]["Diagram"]
```

### Sign flip: `FlipSign`

Flips the sign of the all-ones basis state on three qutrits, leaving every other state fixed:

```wl
QuantumCircuitOperator[{"FlipSign"[3]}]["Diagram"]
```

### Angular momentum: `JZ`

The spin-1 angular-momentum z-component; in its stored eigenbasis $\mathrm{diag}(-1, 0, 1)$ (the computational representation is $\mathrm{diag}(1, 0, -1)$):

```wl
QuantumCircuitOperator[{"JZ"[1]}]["Diagram"]
```

### Raising operator: `J+`

The spin-1 raising operator, $J_+|m\rangle = \sqrt{(j-m)(j+m+1)}\,|m+1\rangle$ with $J_+|j,j\rangle = 0$; for spin 1 the nonzero coefficients are $\sqrt2$:

```wl
QuantumCircuitOperator[{"J+"[1]}]["Diagram"]
```

### Spin rotation: `WignerD`

A spin rotation by Euler angles, here a $\pi/2$ z-rotation on a spin-1/2:

```wl
QuantumCircuitOperator[{"WignerD"[{Pi/2, 0, 0}]}]["Diagram"]
```

**Tier 5 - Gate gadgets and encoders.**

## Gate Gadgets and Encoders

These are the small fixed circuits that recur as building blocks: entangled-state
preparations, the three-qubit permutation gates, integer encoders, graph states,
and the general state-preparation compiler.

### Bell pair: `Bell`

A Hadamard then a CNOT, preparing the Bell state $(|00\rangle + |11\rangle)/\sqrt2$:

```wl
QuantumCircuitOperator[{"Bell"}]["Diagram"]
```

---

Acting on $|00\rangle$ it produces the Bell state, amplitudes $(1, 0, 0, 1)/\sqrt2$:

```wl
Normal[QuantumCircuitOperator[{"Bell"}][QuantumState["00"]]["StateVector"]]
```

### GHZ state: `GHZ`

The n-qubit GHZ state $(|0\cdots0\rangle + |1\cdots1\rangle)/\sqrt2$, here on four wires:

```wl
QuantumCircuitOperator[{"GHZ"[4]}]["Diagram"]
```

### Doubly-controlled NOT: `Toffoli`

The doubly-controlled NOT, shown as its 15-gate Clifford+T decomposition:

```wl
QuantumCircuitOperator[{"Toffoli"}]["Diagram"]
```

### Controlled-SWAP: `Fredkin`

The controlled-SWAP, shown as its elementary-gate decomposition:

```wl
QuantumCircuitOperator[{"Fredkin"}]["Diagram"]
```

### Magic basis: `Magic`

The two-qubit magic-basis circuit, `S` on each wire then `H` and `CNOT`:

```wl
QuantumCircuitOperator[{"Magic"}]["Diagram"]
```

### Integer encoding: `Number`

Encodes an integer in the computational basis, here 5 as the bit string `101` with an `X` on each 1:

```wl
QuantumCircuitOperator[{"Number"[5]}]["Diagram"]
```

---

A second argument pads the encoding to a fixed number of qubits:

```wl
QuantumCircuitOperator[{"Number"[5, 4]}]["Diagram"]
```

### Fourier-basis integer: `PhaseNumber`

Encodes an integer in the Fourier basis, Hadamards followed by phase shifts rather than bit flips:

```wl
QuantumCircuitOperator[{"PhaseNumber"[5]}]["Diagram"]
```

### Graph state: `Graph`

The graph state of a given graph: a Hadamard on every vertex and a controlled-Z on every edge, here a four-vertex path:

```wl
QuantumCircuitOperator[{"Graph"[PathGraph[Range[4]]]}]["Diagram"]
```

### State preparation: `QuantumState`

Compiles a circuit that prepares a given pure state (mixed inputs are rejected), here the Bell state:

```wl
QuantumCircuitOperator[{"QuantumState"[QuantumState["Bell"]]}]["Diagram"]
```

---

The same constructor under its other name, `StatePreparation`, here preparing a GHZ state:

```wl
QuantumCircuitOperator[{"StatePreparation"[QuantumState["GHZ"]]}]["Diagram"]
```

**Tier 6 - Primitive gate set.**

## Primitive Gates

The atomic gates every circuit is built from. These are commodity, so this page
is a compact reference: one cell per distinct gate, with trivial name aliases
noted rather than repeated.

### Pauli gates: `X`, `Y`, `Z`

The bit flip $X$, swapping $|0\rangle \leftrightarrow |1\rangle$:

```wl
QuantumCircuitOperator[{"X"}]["Diagram"]
```

---

The $Y$ gate, a combined bit-and-phase flip sending $|0\rangle \to i|1\rangle$:

```wl
QuantumCircuitOperator[{"Y"}]["Diagram"]
```

---

The sign flip $Z$, negating $|1\rangle$:

```wl
QuantumCircuitOperator[{"Z"}]["Diagram"]
```

The aliases `PauliX`, `PauliY`, `PauliZ` build the same gates, as does `NOT` for `X`.

### Hadamard: `H`

The Hadamard, mapping $|0\rangle \to |+\rangle$ and $|1\rangle \to |-\rangle$ (alias `Hadamard`):

```wl
QuantumCircuitOperator[{"H"}]["Diagram"]
```

---

Acting on $|0\rangle$ it gives the equal superposition $(|0\rangle + |1\rangle)/\sqrt2$:

```wl
Normal[QuantumCircuitOperator[{"H"}][QuantumState["0"]]["StateVector"]]
```

### Phase gates: `S`, `T`

The `S` gate, a $\pi/2$ phase on $|1\rangle$, $\mathrm{diag}(1, i)$:

```wl
QuantumCircuitOperator[{"S"}]["Diagram"]
```

---

The `T` gate, a $\pi/4$ phase on $|1\rangle$, $\mathrm{diag}(1, e^{i\pi/4})$:

```wl
QuantumCircuitOperator[{"T"}]["Diagram"]
```

### Square root of NOT: `V`

The square root of `NOT`, with $V^2 = X$ (note: bare `SX` is the two-qubit $S \otimes X$, while `SX[]` resolves to the single-qubit root):

```wl
QuantumCircuitOperator[{"V"}]["Diagram"]
```

### Rotations: `RX`, `RY`, `RZ`

Rotation about the X axis by an angle:

```wl
QuantumCircuitOperator[{"RX"[Pi/2]}]["Diagram"]
```

---

Rotation about the Y axis:

```wl
QuantumCircuitOperator[{"RY"[Pi/3]}]["Diagram"]
```

---

Rotation about the Z axis:

```wl
QuantumCircuitOperator[{"RZ"[Pi/4]}]["Diagram"]
```

### General single-qubit gate: `U`

The general single-qubit gate from three Euler angles; here $U[\pi/2, 0, \pi]$ is the Hadamard:

```wl
QuantumCircuitOperator[{"U"[Pi/2, 0, Pi]}]["Diagram"]
```

### Phase gate: `P`

A phase gate $\mathrm{diag}(1, e^{i\theta})$ by an explicit angle (the integer form `PhaseShift[k]` gives a $2\pi/2^k$ phase for positive $k$, and $2\pi\,\mathrm{Sign}(k)/2^{|k|}$ in general):

```wl
QuantumCircuitOperator[{"P"[Pi/4]}]["Diagram"]
```

### Two-qubit gates: `CNOT`, `CZ`, `SWAP`

Controlled NOT (alias `CX`):

```wl
QuantumCircuitOperator[{"CNOT"}]["Diagram"]
```

---

Controlled Z:

```wl
QuantumCircuitOperator[{"CZ"}]["Diagram"]
```

---

The swap of two wires:

```wl
QuantumCircuitOperator[{"SWAP"}]["Diagram"]
```
