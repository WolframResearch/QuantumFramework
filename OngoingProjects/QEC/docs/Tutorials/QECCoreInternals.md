---
Template: TechNote
Name: QECCoreInternals
Title: The Circuit-Level Layer, Function by Function
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/QECCoreInternals
Keywords: [quantum error correction, syndrome extraction circuit, Pauli frame, detector error model, memory experiment, circuit-level noise, hook error, Stim, logical error rate, fault mechanism, minimum weight decoder]
RelatedGuides: [WolframQuantumComputationFramework]
---

The companion note *Stabilizer Codes in the Quantum Framework* reads the code layer from the
outside: codes in, error rates out. This one opens it. It walks the four modules that turn a
code into a noisy experiment — `Circuit.wl`, `DetectorModel.wl`, `Memory.wl` and `Stim.wl` —
and calls **every** function in them on a real code, including the ones that are internal to
the package.

The reason to document internals at all is that they are where the physics is. A syndrome is a
matrix product only as long as the measurement is assumed perfect; once a check is a subcircuit,
a fault in it can spread, and the object that records how is the detector error model. Reading
the layer from the outside you see a polynomial in $p$; reading it from the inside you see why
that polynomial has a linear term.

The chain each section follows:

    QECCode  ->  instructions + Pauli frame       (Circuit.wl)
             ->  fault -> detectors, observables  (DetectorModel.wl)
             ->  decode, and a rate               (Memory.wl)
             ->  the same experiment in Stim      (Stim.wl)

## Loading the layer

The code layer is a development package under `OngoingProjects/QEC/`, so this page loads both the
framework and the layer itself, exactly as the companion note does. Symbols introduced here stay
bound for the whole page:

```wl
Needs["Wolfram`QuantumFramework`"];

qecCore = SelectFirst[
    {
        Quiet @ Check[FileNameJoin[{ParentDirectory[NotebookDirectory[], 2], "QECCore", "QECCore.wl"}], $Failed],
        FileNameJoin[{ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
            "OngoingProjects", "QEC", "QECCore", "QECCore.wl"}],
        FileNameJoin[{Directory[], "QECCore", "QECCore.wl"}]
    },
    StringQ[#] && FileExistsQ[#] &,
    $Failed
];

If[ qecCore === $Failed,
    Failure["QECCoreNotFound", <|"MessageTemplate" ->
        "Could not locate QECCore/QECCore.wl. Get it by hand, then run the rest of the page."|>],
    Get[qecCore];
    MemberQ[$ContextPath, "Wolfram`QuantumFramework`QEC`"]
]
```

<!-- => True -->

### Reaching the internals

Twenty-one symbols are `PackageExport` and land on the context path. Everything else this page
touches is `PackageScope`, which means it exists but is not exported, and is reached by its full
context name. This is the same idiom the test suite uses:

```wl
qec = "Wolfram`QuantumFramework`QEC`PackageScope`";
sym[n_String] := Symbol[qec <> n];

{codeData, codeStabCount, codeLogQubits, codeDistanceF, codeLabelMatrix} =
    sym /@ {"codeData", "codeStabilizerCount", "codeLogicalQubits", "codeDistance", "codeLabelMatrix"};

{circuitData, circuitInstr, circuitQubits, codeCircInstr, genInstr} =
    sym /@ {"circuitData", "circuitInstructions", "circuitQubitCount",
            "codeCircuitInstructions", "generatorInstructions"};

{framePropagate, frameZero, engineGates, oneQubitOps} =
    sym /@ {"framePropagate", "frameZero", "instructionEngineGates", "$circuitOneQubitOps"};

{circuitSchedule, circuitIdleQubits, circuitInstructionLayers, circuitIdleAnchor,
 circuitIdleSlots} =
    sym /@ {"circuitSchedule", "circuitIdleQubits", "circuitInstructionLayers",
            "circuitIdleAnchor", "circuitIdleSlots"};

{codeDetModel, faultMechanisms, faultEffect, recordDetectors, roundLength, defaultRounds} =
    sym /@ {"codeDetectorModel", "circuitFaultMechanisms", "faultEffect",
            "recordDetectors", "roundLength", "defaultRounds"};

{demRowKeys, demGroups, demDecoderTable, demExactFailure, demSampled, demRate} =
    sym /@ {"demRowKeys", "demGroups", "demDecoderTable", "demExactFailure",
            "demSampledFailure", "demRate"};

{noiseLevel, noiseCircRates, noiseMeasErr, noiseSymbolicQ, noiseLevels, circLocations} =
    sym /@ {"noiseLevel", "noiseCircuitRates", "noiseMeasurementError",
            "noiseSymbolicQ", "$QECNoiseLevels", "$QECCircuitLocations"};

{stimGateName, codeStimCircuit} = sym /@ {"stimGateName", "codeStimCircuit"};

Length[Names["Wolfram`QuantumFramework`QEC`*"]]
```

A third tier is invisible even here: `letterAt`, `basisGate`, `unbasisGate`, `rotation`,
`oddParityRate`, `bitsToInteger`, `encoderLines` and friends are private to their file and cannot
be called by name at all. They appear below through their effect — the basis rotations show up in
the emitted instruction list, `oddParityRate` through the `"DetectorRates"` property.

Two codes carry the page. The bit-flip code is small enough for every exact route; the five-qubit
code has generators carrying $Y$, so it exercises the $\sqrt{X}$ basis change:

```wl
bitFlip = QECCode["BitFlipCode"];
five = QECCode["5QubitCode"];
{bitFlip["Parameters"], five["Parameters"]}
```

---

The internal views of the same object. `codeData` unwraps the association; the label matrix has
$m + 2k$ rows, the syndrome half stacked on the logical-class half:

```wl
{codeStabCount[codeData[bitFlip]], codeLogQubits[codeData[bitFlip]],
 codeDistanceF[codeData[bitFlip]], Dimensions[codeLabelMatrix[codeData[bitFlip]]]}
```

The distance is worth a pause: it is **1**, not 3. The classical repetition code corrects three
bit flips, but as a quantum code a lone $Z$ commutes with both checks and is not in the
stabilizer group, so it is an undetectable weight-one logical. The layer reports the
mathematics, not the folklore.

## Circuit.wl — the syndrome-extraction circuit

Up to here a check was measured by fiat. That is the right object at code-capacity level, where
measurement is assumed perfect. It is also the thing that has to go if the noise is going to be
honest, because in a real device a check is not a matrix product — it is a subcircuit, and every
gate and every readout in it can fail.

The construction is one ancilla per generator, which is Gottesman's section 12.1.1, figure 12.1a
— a section titled *Non-Fault-Tolerant Measurement of Paulis*, which is the honest name for what
this is.

### One generator, by hand

An instruction is a plain list `{op, qubits...}`, deliberately without a head: the frame
propagator walks millions of them and a `Switch` on a string part is the cheapest dispatch
available. The one-qubit vocabulary:

```wl
oneQubitOps
```

---

Take the five-qubit code's first generator. To measure a Pauli $g$, rotate each data qubit in
$g$'s support so that $g$'s letter there becomes $Z$, CNOT each of those into the ancilla, rotate
back, and measure the ancilla in $Z$:

$$X \to H \quad (HXH = Z), \qquad Y \to V \quad (VYV^{-1} = Z), \qquad Z \to \text{nothing}$$

```wl
g1 = First[five["CheckMatrix"]];
QECPauliString[Append[g1, 0]]
```

---

`generatorInstructions` emits that for one generator, one ancilla, one round. Here `XZZXI` has
support $\{1,2,3,4\}$, so qubits 1 and 4 get an `H` on the way in and out and qubits 2 and 3 get
nothing:

```wl
genInstr[g1, 5, 6] // Column
```

---

A generator carrying a $Y$ instead brings in `V` and `Vdg`. The two act identically on the
symplectic part — $V^2$ is the Pauli $X$ — and are separate operations only because the noise
model counts instructions, and a rotation is one location whatever the hardware needs to realise
it:

```wl
gY = First[QECCode[{"YYI", "IYY"}]["CheckMatrix"]];
genInstr[gY, 3, 4] // Column
```

### The whole circuit

Data qubits are $1 \ldots n$ and the ancilla of generator $j$ is $n + j$. The ancillas are reset
and reused every round, so an $r$-round circuit still runs on $n + m$ qubits, not $n + rm$:

```wl
{Length[codeCircInstr[codeData[five], 1]], Length[codeCircInstr[codeData[five], 3]]}
```

---

The public constructor wraps that with its metadata:

```wl
circ = QECSyndromeCircuit[five]
```

---

```wl
circ3 = QECSyndromeCircuit[five, 3]
```

---

The properties:

```wl
circ["Properties"]
```

---

```wl
AssociationMap[circ, {"DataQubits", "Ancillas", "Qubits", "Rounds",
    "Depth", "InstructionCount", "MeasurementCount"}]
```

---

```wl
circ["GateCounts"]
```

---

The record is filled in emission order, so measurement index $(r-1)m + j$ is generator $j$ of
round $r$. `"MeasurementLabels"` names them:

```wl
circ3["MeasurementLabels"]
```

---

The `PackageScope` accessors reach the same association without the property dispatch:

```wl
{Keys[circuitData[circ]], Length[circuitInstr[circuitData[circ]]], circuitQubits[circuitData[circ]]}
```

### The Pauli frame

Here is the move that makes everything downstream viable. The circuit is Clifford and the noise
is Pauli, so a fault does not have to be simulated as a state — it can be carried as a Pauli
*frame* and pushed through the gates by conjugation, which over GF(2) is four XORs:

| gate | action on the frame |
|---|---|
| `H` | $x \leftrightarrow z$ |
| `S` | `z ^= x` |
| `V`, `Vdg` | `x ^= z` |
| `CNOT` | `x_t ^= x_c`, `z_c ^= z_t` |
| `R` | $x, z := 0$ |
| `M` | record $x$ |

Phases never enter: the frame's contribution to an outcome is whether it anticommutes with the
measured $Z$, which is the ancilla's $x$ bit.

An empty frame:

```wl
frameZero[9]
```

---

With no faults, every ancilla reads 0 and the data frame is empty. That is what makes a detector
"deterministically zero when nothing goes wrong":

```wl
framePropagate[circ["Instructions"], circ["Qubits"], {}]
```

---

A fault is `{instructionIndex, qubit, {x, z}}`, applied *after* the instruction at that index.
Index 0 means before the circuit starts, which is where a code-capacity error on the data lands.
Put an $X$ on data qubit 1 and read the measurement record:

```wl
framePropagate[circ["Instructions"], circ["Qubits"], {{0, 1, {1, 0}}}]["Record"]
```

---

The oracle: `codeSyndrome` computes the same thing by matrix product and knows nothing about
circuits. They agree, and that agreement over every weight-one and weight-two error is what
`Tests/Circuit.wlt` checks:

```wl
five["Syndrome", QECPauliString["XIIII"]]
```

---

Same for a $Z$, which the five-qubit code also sees:

```wl
{framePropagate[circ["Instructions"], circ["Qubits"], {{0, 3, {0, 1}}}]["Record"],
 five["Syndrome", QECPauliString["IIZII"]]}
```

### The hook error

The book's figure 12.1b puts the ancilla in $|+\rangle$ as the *control*. Here the ancilla is the
*target*, and the mechanism is mirrored. Getting this backwards is the natural mistake, so it is
worth watching happen.

Find the CNOT ladder of generator 1, which runs into ancilla 6:

```wl
ladder = Select[Range[circ["InstructionCount"]], MatchQ[circ["Instructions"][[#]], {"CNOT", _, 6}] &]
```

---

An $X$ on the ancilla mid-ladder reaches no data qubit at all: CNOT carries $X$ from control to
target, and the ancilla *is* the target. It sits there until the measurement and flips that one
outcome. A pure readout error:

```wl
framePropagate[circ["Instructions"], circ["Qubits"], {{ladder[[2]], 6, {1, 0}}}]
```

---

A $Z$ is the one that spreads. CNOT carries $Z$ from target back to control, and the ancilla
keeps its $Z$, so it lands on every data qubit the ladder has *still to touch*. After the
un-rotation the residual is exactly the un-extracted tail of the generator:

```wl
Grid[
    Prepend[
        Table[
            With[{run = framePropagate[circ["Instructions"], circ["Qubits"],
                    {{ladder[[k]], 6, {0, 1}}}]},
                {k, Take[run["Frame"][[1]], 5], Take[run["Frame"][[2]], 5], run["Record"]}],
            {k, Length[ladder]}],
        {"after CNOT #", "x on data", "z on data", "record"}],
    Frame -> All, Alignment -> Left]
```

`XZZXI` leaves `IZZXI`, `IIZXI` or `IIIXI` depending on where the fault struck — and being a
sub-product of the generator, it does **not** fire that generator's own check in that round. A
$Z$ before the first CNOT leaves the whole generator, which is a stabilizer and harmless; a $Z$
after the last leaves nothing. Only mid-ladder hurts.

This is what turns a distance-three code's logical error rate from order $p^2$ into order $p$,
and it is a defect of the *gadget*, not a property of circuit-level noise. A gadget satisfying
the book's gate and error-correction propagation properties keeps the $p^2$.

### The schedule, and idling

The instruction list is flat and sequential. The *circuit* is not: two instructions on
disjoint qubits happen in the same time step. `circuitSchedule` recovers those layers
from the flat list by ASAP list scheduling — each instruction goes into the earliest
layer in which none of its qubits is already busy:

```wl
{Length[circ["Instructions"]], Length[circuitSchedule[circ["Instructions"], circ["Qubits"]]]}
```

---

The schedule is *derived*, never stored. The instruction list stays canonical, so every
consumer that does not care about time — the frame propagator, the Stim writer, the tests
— is untouched. What it is for is the fourth kind of fault location: a qubit that no
instruction of a layer touches is **waiting**, and waiting is a location in its own right.
Section 10.1.1's Definition 10.1 lists preparation, gate, *wait* and measurement.

```wl
bfCirc = QECSyndromeCircuit[bitFlip];
bfLayers = circuitSchedule[bfCirc["Instructions"], bfCirc["Qubits"]];
Grid[
    Prepend[
        Table[{L, bfLayers[[L]], circuitIdleQubits[bfCirc["Instructions"], 5, bfLayers[[L]]]},
            {L, Length[bfLayers]}],
        {"layer", "instructions", "waiting"}],
    Frame -> All, Alignment -> Left]
```

Eight instructions in six layers: `{R,4}` and `{R,5}` share the first, `{M,4}` and
`{CNOT,2,5}` share the fourth. And 18 qubit-steps of waiting, against 28 if every
instruction had its own step — which is the point of section 15.5.1, where partial
parallelism multiplies the storage rate and *"to have a threshold, it is essential to do
parallel gates, at least when the storage error rate is non-zero"*.

---

Where an idle fault goes is the part that is easy to get wrong, and the wrong answer is
tempting. "The end of the layer" is not a point in time, because **the instruction list is
not ordered by layer**: layer 1 holds instructions 1 and 5 while layer 2 holds instruction
2, so in emission order instruction 5 comes *after* instruction 2.

What is well defined is per qubit. The frame propagator walks the flat list, so for a fault
on qubit $q$ the only thing that matters is which instructions touching $q$ run after it.
A qubit idle in layer $L$ has no instruction in that layer, so the fault belongs after the
last instruction touching $q$ in an *earlier* layer — and at index 0, before the circuit,
when none has run yet:

```wl
With[{layerOf = circuitInstructionLayers[bfCirc["Instructions"], 5]},
    circuitIdleAnchor[bfCirc["Instructions"], layerOf, 1, #] & /@ Range[6]]
```

Data qubit 1 is touched exactly once, by the CNOT at instruction 2 in layer 2. So it is
idle in layers 1, 3, 4, 5 and 6, anchored before the circuit for the first and after that
one instruction for the other four. Several layers sharing an anchor is correct rather than
a collapse: they are independent chances to decohere during a wait the propagator sees as a
single gap.

---

`circuitIdleSlots` is the single source of truth, and it is shared on purpose — the
detector model turns each slot into a mechanism and the Stim writer turns each slot into a
`DEPOLARIZE1` at the same anchor, so the exported circuit cannot drift from the model it
exists to cross-check:

```wl
Length[circuitIdleSlots[bfCirc["Instructions"], 5]]
```

---

One consequence worth seeing, because it is why the anchors are computed over the whole
circuit rather than per round: **rounds pipeline**. An ancilla's reset for round 2 schedules
while the other ancilla of round 1 is still being measured, so $r$ rounds are not $r$ copies
of one round's schedule:

```wl
Table[Length[circuitSchedule[codeCircInstr[codeData[bitFlip], r], 5]], {r, 1, 3}]
```

### Handing the circuit to the engine

For validation only: the engine simulates states, so this is how one checks that the emitted
circuit really measures the generators we think it does. `R` has no unitary form and is dropped,
`M` is not a gate either, and `Vdg` expands to three `V`:

```wl
Take[engineGates[circ["Instructions"]], UpTo[10]]
```

## Noise.wl — the three levels

A code means nothing until you say what it is protecting against. Three levels, each dropping an
assumption the one before it made:

```wl
noiseLevels
```

---

Circuit-level noise splits by location type:

```wl
circLocations
```

---

Code capacity: the data is hit once, the checks are perfect.

```wl
nCap = QECNoiseModel["BitFlip", p]
```

---

Phenomenological: the data is hit between rounds and the readout can lie, but the extraction
itself is a black box.

```wl
nPhen = QECNoiseModel["BitFlip", p, "MeasurementError" -> q]
```

---

Circuit level: every gate, reset and readout can fail.

```wl
nCirc = QECNoiseModel["Circuit", p]
```

---

Rates may be given per location type, and anything omitted is filled with zero:

```wl
nMix = QECNoiseModel["Circuit", <|"OneQubit" -> p, "TwoQubit" -> 2 p, "Measurement" -> q|>];
noiseCircRates[First[nMix]]
```

---

The accessors that the rest of the layer dispatches on:

```wl
{noiseLevel[First[nCap]], noiseLevel[First[nPhen]], noiseLevel[First[nCirc]],
 noiseMeasErr[First[nPhen]], noiseMeasErr[First[nCap]],
 noiseSymbolicQ[First[nCirc]], noiseSymbolicQ[First[QECNoiseModel["Circuit", 0.001]]]}
```

`Idle` is **modelled** — see the schedule above — but it still **defaults to zero**, and that
default is optimistic rather than neutral.
The generators are extracted one after another rather than interleaved, so there is *more*
waiting, not less. What excuses it is that a fixed-size code extracted serially is a constant
number of steps — not the claim that there is no time step to idle in, which section 15.5.2 of
the book explicitly forbids.

## DetectorModel.wl — from faults to detectors

Once the checks are measured by a circuit, a syndrome bit is no longer trustworthy on its own: it
can be wrong because the data was hit or because the readout lied, and one round cannot tell
those apart. The fix used here is to decode the *differences*:

$$D_1 = S_1, \qquad D_r = S_r \oplus S_{r-1}, \qquad D_{R+1} = S_{\text{final}} \oplus S_R$$

A detector is a parity that is deterministically 0 when nothing goes wrong, so a fired detector is
unambiguous evidence of a fault, and a readout error fires two detectors in a row rather than
looking like a data error.

Whose fix this is matters. The book's answer to an untrustworthy syndrome is repeat *and agree*:
trust a run of $t+1$ consecutive identical syndromes. Differencing consecutive syndromes and
decoding the whole spacetime history is a different procedure, from the topological-code
literature — Dennis, Kitaev, Landahl and Preskill — and as an explicit object with an error model
attached it is Gidney's Stim. Section 12.2.2 of the book raises the strategy and declines to
analyse it: *"it is difficult to analyze such strategies in the completely general case."*

### The default number of rounds

```wl
bfData = codeData[bitFlip];
{defaultRounds[bfData], defaultRounds[codeData[five]], roundLength[bfData, 1]}
```

The default is the code distance, which is the topological-memory convention: with repeated noisy
measurement the history is a lattice in one more dimension, and taking fewer rounds than the
distance protects the time direction less than the space ones. It is deliberately *not* the
book's fault-tolerant criterion, and the two should not be conflated.

### The mechanism list

A mechanism is `<|"Faults" -> {{index, qubit, {x, z}}, ...}, "Probability" -> p, "Location" -> label|>`.
The `"Location"` is the key to everything downstream: it groups mechanisms that are alternatives
of the same physical location.

```wl
mechC = faultMechanisms[bfData, First[QECNoiseModel["Circuit", p]], 1];
{Length[mechC], First[mechC]}
```

---

Three per one-qubit gate at $p/3$, fifteen per CNOT at $p/15$, one $X$ before each measurement,
one after each reset:

```wl
Counts[First /@ Lookup[mechC, "Location"]]
```

---

Waiting is the fifth kind, and it is the one that only appears if you ask for it. The `"Idle"`
rate defaults to zero, and at zero `idleMechanisms` returns nothing at all rather than a list of
zero-probability rows — so every number in this note is what it was before the schedule existed.
Turn it on and three Paulis appear per idle slot:

```wl
idleMech = sym["idleMechanisms"];
withIdle = QECNoiseModel["Circuit",
    <|"OneQubit" -> p, "TwoQubit" -> p, "Measurement" -> p, "Reset" -> p, "Idle" -> pIdle|>];
{
  Length[mechC],
  Length[faultMechanisms[bfData, First[withIdle], 1]],
  3 Length[circuitIdleSlots[codeCircInstr[bfData, 1], 5]]
}
```

The bit-flip extraction circuit has no one-qubit gates at all — both checks are pure $Z$, so
nothing needs rotating — which makes it a clean place to see the idle count on its own.

---

Phenomenological instead puts three Paulis on each data qubit at the start of each round, plus one
flip per readout:

```wl
mechP = faultMechanisms[bfData, First[QECNoiseModel["BitFlip", p, "MeasurementError" -> q]], 1];
{Length[mechP], Counts[First /@ Lookup[mechP, "Location"]]}
```

The book counts the identity as an admissible outcome at a faulty location; the $p/3$ and $p/15$
here exclude it, which makes this marginally *less* conservative. And listing the three outcomes
of a depolarizing channel as three independent mechanisms is an approximation — the same one Stim
makes. Sampling draws the true channel, so simulated numbers are exact; only the probabilities
attached to the rows are approximate, and the minimum-weight decoder does not read them.

### What one fault does

`faultEffect` runs the frame propagator and reads the residual through the code's label matrix.
One GF(2) matrix product returns both halves: the syndrome a perfect final round would see, and
the logical class.

```wl
bfInstr = codeCircInstr[bfData, 1];
faultEffect[bfData, bfInstr, 5, 2, 1, {{0, 1, {1, 0}}}]
```

An $X$ on data qubit 1 fires the checks it anticommutes with, and flips nothing logical.

---

Now a $Z$ on the same qubit. It is invisible to a bit-flip code's checks — but it *is* a logical:

```wl
faultEffect[bfData, bfInstr, 5, 2, 1, {{0, 1, {0, 1}}}]
```

No detector, one observable. A weight-one undetectable fault, which is the distance-1 result from
the top of the page seen from inside the circuit.

---

`recordDetectors` is the piece that turns a raw measurement record into that detector vector.
Given two rounds of two checks plus the final syndrome:

```wl
recordDetectors[{1, 0, 1, 0}, {1, 0}, 2, 2]
```

### The model object

```wl
dem = QECDetectorModel[bitFlip, QECNoiseModel["Circuit", p], 1]
```

---

```wl
dem["Properties"]
```

---

```wl
AssociationMap[dem, {"Faults", "Detectors", "Observables", "Rounds", "Checks"}]
```

---

```wl
dem["LocationCounts"]
```

---

Rows are faults; columns are detectors, then logical observables. Every noisy run of the
experiment is a GF(2) sum of rows:

```wl
{Dimensions[dem["DetectorMatrix"]], Dimensions[dem["ObservableMatrix"]]}
```

---

```wl
Grid[Prepend[
    Table[{dem["Locations"][[i]], dem["DetectorMatrix"][[i]],
           dem["ObservableMatrix"][[i]], dem["Probabilities"][[i]]}, {i, Min[8, dem["Faults"]]}],
    {"location", "detectors", "observables", "probability"}],
  Frame -> All, Alignment -> Left]
```

---

Faults no detector can see. Their existence is not a bug: they are the circuit-level analogue of a
logical operator, and how many it takes to build one is the circuit-level distance.

```wl
dem["UndetectableFaults"]
```

### Exact firing rates

A detector fires when an odd number of the mechanisms touching it fire, and odd-parity
probabilities multiply through the $(1 - 2p)$ form, one factor per independent location:

$$P(\text{odd}) = \frac{1}{2}\left(1 - \prod_{\text{locations}} \Bigl(1 - 2\!\!\sum_{i \in \text{loc},\, i \,\text{touches}\, d}\!\! p_i\Bigr)\right)$$

```wl
Simplify[dem["DetectorRates"]]
```

---

```wl
Simplify[dem["ObservableRates"]]
```

Exact, symbolic if the rates are, and the natural quantity to check an external simulator
against: it is basis-free and every detector is a separate test. This is what `StimCrossCheck/`
compares, detector by detector.

---

Numerically:

```wl
QECDetectorModel[bitFlip, QECNoiseModel["Circuit", 0.001], 1]["DetectorRates"]
```

### The rows that get dropped

`codeDetectorModel` is the raw builder underneath the constructor. It drops two kinds of row:

```wl
nP = First[QECNoiseModel["BitFlip", p, "MeasurementError" -> q]];
raw = codeDetModel[bfData, nP, 1];
{Length[faultMechanisms[bfData, nP, 1]], Length[raw["Probabilities"]]}
```

First, faults that fire no detector *and* flip no observable — the ones a reset wipes before they
reach anything. Keeping them only inflates the subset search a decoder has to do.

Second, faults the channel cannot actually produce, meaning probability literally zero: bit-flip
noise never makes a $Y$ or a $Z$. These matter more than they look, because the minimum-weight
decoder does not read probabilities. A zero-probability row left in the table can claim a detector
pattern that a real fault should have had, and quietly make the decoder worse than the noise
deserves. A symbolic rate is never dropped, since nothing can be concluded about it.

## Memory.wl — the memory experiment

Prepare a codeword, run $r$ rounds of noisy syndrome extraction, read the data out perfectly at
the end, decode the whole detector history, and ask whether the logical qubit survived.

Two limits bound the two routes:

```wl
{$QECExactDetectorLimit, $QECDecoderSubsetLimit}
```

### Packing and grouping

Detector and observable patterns are only ever XORed and compared, so they are carried as
integers: one `BitXor` replaces a vector `Mod`.

```wl
a = codeDetModel[bfData, First[QECNoiseModel["Circuit", p]], 1];
Take[#, UpTo[10]] & /@ demRowKeys[a]
```

---

Mechanisms belonging to the same physical location are alternatives, not independent coins: at
most one of the three Paulis of a depolarizing channel happens. The location label already carries
that grouping.

```wl
{Length[demGroups[a]], Take[demGroups[a], UpTo[4]]}
```

### The decoder

Detector pattern to the observable flip implied by the lightest fault set that explains it. Because
subsets are enumerated by increasing weight, the first set to produce a key claims it — which is
what makes the table minimum-weight.

```wl
t1 = demDecoderTable[a, 1];
{Length[t1], Take[Normal[t1], UpTo[6]]}
```

---

A larger reach explains more patterns. A pattern the table never reaches is undecodable and counts
as a failure:

```wl
{Length[demDecoderTable[a, 1]], Length[demDecoderTable[a, 2]]}
```

This is a different kind of decoder from the book's, not a variant of it. Gottesman's uses the
repetitions only as a certificate: find the last run of $t+1$ agreeing syndromes, conclude one of
them was taken with no faults, discard everything else, and look that single syndrome up. What
happens here is an optimisation over the whole spacetime history: nothing is discarded, all rounds
are used jointly, and the answer is the lightest fault configuration consistent with every
detector. It inherits none of that section's guarantees — what is exact here is the error rate of
a stated protocol, not a fault-tolerance proof.

### The exact route

Rather than enumerate $2^{\#\text{faults}}$ configurations, carry the *distribution over effects*
— a vector of length $2^{\text{detectors} + \text{observables}}$ — and fold the coins in one at a
time:

$$P'[s] = P[s]\Bigl(1 - \sum_i p_i\Bigr) + \sum_i p_i\, P[s \oplus \sigma_i]$$

the sum running over the mutually exclusive outcomes of one location, so a three-way depolarizing
channel is handled exactly rather than as three independent coins. The state space is what bounds
it:

```wl
2^(a["Detectors"] + a["Observables"])
```

---

Cost is (locations) $\times$ (state space), not $2^{\#\text{faults}}$, and with symbolic rates the
answer is a circuit-level logical error rate in closed form:

```wl
Simplify[demExactFailure[a, 2]]
```

This is the move a sampling simulator structurally cannot make.

### The sampled route

One draw per location, all shots at once: the group's own probabilities plus the leftover mass for
"nothing happened". Mechanisms dropped from the model because they had no effect are part of that
leftover, which is why the weights still add up.

```wl
aN = codeDetModel[bfData, First[QECNoiseModel["Circuit", 0.01]], 1];
{demSampled[aN, 200000, 2], N[Simplify[demExactFailure[aN, 2]]]}
```

The same experiment by two routes; the sampled one should sit within a few standard errors of the
exact one.

### The entry point

`count === None` picks the exact fold; an integer picks the sampler:

```wl
{Simplify[demRate[bfData, First[QECNoiseModel["Circuit", p]], 1, None, 2]],
 demRate[bfData, First[QECNoiseModel["Circuit", 0.01]], 1, 100000, 2]}
```

## The public entry point

`QECLogicalErrorRate` dispatches on the noise level: code capacity keeps the original route
(exact enumeration over cosets), and the two noisy levels are routed to `Memory.wl`.

```wl
Options[QECLogicalErrorRate]
```

---

```wl
Expand[QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", p]]]
```

### The consistency identity

A phenomenological model with perfect readout over one round *is* code-capacity noise. Two
completely different code paths — enumeration over cosets, and a fold over detectors — therefore
have to produce the same polynomial:

```wl
Simplify[
    QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", p]] ==
    QECLogicalErrorRate[bitFlip, QECNoiseModel["BitFlip", p, "MeasurementError" -> 0],
        "Rounds" -> 1]
]
```

---

Let the readout lie and a second parameter appears:

```wl
Simplify[QECLogicalErrorRate[bitFlip,
    QECNoiseModel["BitFlip", p, "MeasurementError" -> q], "Rounds" -> 1]]
```

### What circuit-level noise costs

```wl
Series[QECLogicalErrorRate[bitFlip, QECNoiseModel["Circuit", p], "Rounds" -> 1], {p, 0, 1}]
```

The leading term is $O(p)$, not $O(p^2)$. That is the hook error from the top of the page, priced.
And it is the *gadget* that costs the distance, not the noise level: a fault-tolerant extraction
under the same circuit-level noise would keep the $p^2$.

### A duality that survives one level and breaks at the next

The bit-flip and phase-flip codes are related by a Hadamard on every qubit, so at code capacity
with the dual channel they give the same polynomial. At circuit level they do not:

```wl
{Series[QECLogicalErrorRate[bitFlip, QECNoiseModel["Circuit", p], "Rounds" -> 1], {p, 0, 1}],
 Series[QECLogicalErrorRate[QECCode["PhaseFlipCode"], QECNoiseModel["Circuit", p],
    "Rounds" -> 1], {p, 0, 1}]}
```

---

The whole difference is that the $X$-type checks need eight extra Hadamards. Switch off
one-qubit-gate noise and the duality comes back:

```wl
nNoGates = QECNoiseModel["Circuit", <|"TwoQubit" -> p, "Measurement" -> p, "Reset" -> p|>];
{Series[QECLogicalErrorRate[bitFlip, nNoGates, "Rounds" -> 1], {p, 0, 1}],
 Series[QECLogicalErrorRate[QECCode["PhaseFlipCode"], nNoGates, "Rounds" -> 1], {p, 0, 1}]}
```

That is a statement about a circuit, not about a code, and it is awkward to make any other way.

---

Sampling, when the state space is too large or the rates are numeric anyway:

```wl
QECLogicalErrorRate[bitFlip, QECNoiseModel["Circuit", 0.01], 100000, "Rounds" -> 1]
```

## Stim.wl — the bridge

Stim is the community's reference stabilizer simulator, and the plan is explicit that we connect
to it rather than compete with it: it samples fast and large, this layer derives exactly and
symbolically. Writing the memory experiment out in its syntax buys an independent check on
everything above, and a path to the decoders built on top of it.

The gate map:

```wl
stimGateName /@ {"H", "S", "V", "Vdg", "CNOT", "CZ", "SWAP", "T"}
```

`V` is the engine's $\sqrt{X}$; Stim's `SQRT_X` is the same Clifford up to a Pauli, and a Pauli
difference is invisible to detectors, which are defined relative to the circuit's own noiseless
run. `T` is not Clifford and has no entry.

---

Idling crosses the bridge too, and by the shared primitive rather than by a second computation:
one `DEPOLARIZE1` per `circuitIdleSlots` entry, at the same anchor the detector model uses. So
the exported circuit carries exactly one idle channel per slot where the model carries three
mechanisms, and the two cannot drift:

```wl
noisyIdle = QECNoiseModel["Circuit", <|"OneQubit" -> 1/1000, "TwoQubit" -> 1/1000,
    "Measurement" -> 1/1000, "Reset" -> 1/1000, "Idle" -> 1/500|>];
{
  Length[circuitIdleSlots[codeCircInstr[bfData, 2], 5]],
  StringCount[QECStimCircuit[bitFlip, noisyIdle, 2], "DEPOLARIZE1(0.002)"]
}
```

No `TICK` is emitted. `TICK` would claim a time step, and these lines are in *emission* order
rather than layer order — the two disagree, and reordering to fix that would move measurements
relative to each other and silently invalidate every `rec[-k]` index the detectors are built
from. The final noiseless round carries no idling either, since the residual it exists to reveal
must not be corrupted by the act of revealing it.

---

A full memory experiment. The encoder comes first and is not optional: Stim starts in
$|0\ldots0\rangle$, which is not a codeword for a general code, and a memory experiment has to
start inside the code space or its first round of detectors means nothing.

```wl
Column[StringSplit[QECStimCircuit[bitFlip, QECNoiseModel["Circuit", 0.001], 2], "\n"]]
```

---

Phenomenological emits a single Pauli channel across the data at each round instead of per-gate
noise, because there the extraction is a black box:

```wl
Column[StringSplit[
    QECStimCircuit[bitFlip, QECNoiseModel["BitFlip", 0.01, "MeasurementError" -> 0.02], 1],
    "\n"]]
```

---

The option chooses which logical operator is kept. Our own logical error rate asks whether the
decoder got the residual's *class* right, which covers logical $X$ and $Z$ damage at once. A Stim
memory experiment is single-basis by construction: prepare $|0_L\rangle$, keep it, measure
$\bar{Z}$. The detector part is basis-independent, and is where the two models are compared.

```wl
StringTake[QECStimCircuit[QECCode["SteaneCode"], QECNoiseModel["Circuit", 0.001], 3,
    "Observable" -> "X"], -60]
```

---

The shorthand on the code object, and the internal builder:

```wl
{StringCount[bitFlip["StimCircuit", QECNoiseModel["Circuit", 0.001]], "\n"],
 StringLength[codeStimCircuit[bfData, First[QECNoiseModel["Circuit", 0.001]], 1, "Z"]]}
```

---

Code capacity is refused, because Stim has no notion of a model whose checks are perfect by
assumption:

```wl
QECStimCircuit[bitFlip, QECNoiseModel["BitFlip", 0.01]]
```

## Where this comes from

| claim | source |
|---|---|
| one ancilla per check, and that it is not fault tolerant | Gottesman, *Stabilizer Codes and Quantum Error Correction*, sec. 12.1.1, fig. 12.1a-b |
| using it anyway, as the surface-code community does | ibid., sec. 12.5.1 |
| the Pauli frame | ibid., sec. 12.5.2 |
| locations, faults and one rate per location type | ibid., sec. 10.1.1, Definitions 10.1-10.3 |
| repeat-and-agree as the book's answer to noisy syndromes | ibid., sec. 12.2.2, 12.2.3 |
| differencing set aside as hard to analyse in general | ibid., sec. 12.2.2 |
| ancilla reuse via reset, and the residual as a preparation error | ibid., sec. 15.4, 15.4.3 |
| parallelism is required for a threshold; ill-defined time steps are *defined*, not dropped | ibid., sec. 15.5.1, 15.5.2 |
| waiting is a fault location alongside preparation, gate and measurement | ibid., sec. 10.1.1, Definition 10.1 |
| detectors as differences of consecutive syndromes | Dennis, Kitaev, Landahl, Preskill, *Topological quantum memory*, arXiv:quant-ph/0110143 |
| the detector error model as an object with an error model attached | Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum 5, 497 (2021), arXiv:2103.02202 |
| the tableau engine underneath | Aaronson, Gottesman, Phys. Rev. A 70, 052328 (2004), arXiv:quant-ph/0406196 |
