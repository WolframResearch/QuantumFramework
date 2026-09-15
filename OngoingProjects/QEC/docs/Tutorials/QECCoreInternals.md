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

Twenty-five symbols are `PackageExport` and land on the context path. Everything else this page
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

A third tier is invisible even here: `oddParityRate`, `bitsToInteger`, `encoderLines` and friends
are private to their file and cannot be called by name at all. They appear below through their
effect — `oddParityRate` through the `"DetectorRates"` property. The basis rotations
(`letterAt`, `basisGate`, `unbasisGate`, `rotation`) used to be in that tier and were promoted to
`PackageScope` when the fault-tolerant measurement needed the same letter vocabulary: two
constructions deriving it separately is how they drift apart.

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

## Cat.wl — a probe a single fault cannot spread

Everything above measures the bare-ancilla circuit, and the last two sections priced what it
costs: order $p$ instead of $p^2$, and 288 faults collapsing onto 71 signatures. The diagnosis in
section 12.1.2 is precise — one ancilla shared by every data qubit of a check lets one fault reach
several of them — and so is the fix: make the controlled-$P$ transversal, one ancilla qubit per
letter, which needs those $m$ qubits in a coherent superposition of "do nothing" and "apply $P$ to
all of them". That is the cat state $|0\ldots0\rangle + |1\ldots1\rangle$.

```wl
{catPrep, catInstr, catXPatterns, catPairs, catCoverQ, catCanonical} =
    sym /@ {"catPrepInstructions", "catInstructions", "catXPatterns",
            "catCheckPairs", "catChecksCoverQ", "catCanonicalPattern"};

catPrep[{1, 2, 3, 4}]
```

---

A cat built this way is useless unless it is checked, and the interesting question is *which* pairs
to check. The book says "if we do this on enough pairs of qubits, any single fault in the circuit
originally constructing the cat state will be picked up by the checks" and leaves the set open.
Here it is derived: propagate every single fault of the preparation and collect the $X$ patterns it
can leave.

Two reductions make that set small. $X^{\otimes m}$ stabilises the cat, so a pattern and its
complement are the same error; and a pattern of canonical weight one is a single data error, which
the code corrects. Only canonical weight two and above is dangerous:

```wl
{catCanonical[{1, 1, 0, 0}], catCanonical[{0, 0, 1, 1}], catXPatterns[{1, 2, 3, 4}, 5]}
```

<!-- => {{0, 0, 1, 1}, {0, 0, 1, 1}, {{0, 0, 1, 1}}} -->

That single pattern is the book's own figure 12.3b, recovered rather than transcribed. And the
count follows a rule the book does not state: a fault reaching qubit $j$ leaves a suffix of ones of
length $L$, which canonicalises to $\min(L, m - L)$, so the dangerous $L$ run from 2 to $m-2$ and
there are $m - 3$ of them. An $m = 4$ cat needs one check, not three.

```wl
Table[{m, Length[catXPatterns[Range[m], m + 1]], catPairs[m], m - 3}, {m, 2, 7}] // TableForm
```

---

`catChecksCoverQ` is the proof obligation rather than a comment: the chosen pairs have to detect
every dangerous pattern, and a `False` here says the gadget is not fault tolerant.

```wl
Table[catCoverQ[catPairs[m], catXPatterns[Range[m], m + 1]], {m, 2, 9}]
```

<!-- => {True, True, True, True, True, True, True, True} -->

---

And no pair is redundant — drop one and coverage fails, which is what says the default set is
minimal rather than merely sufficient:

```wl
With[{m = 6, pats = catXPatterns[Range[6], 7]},
    Table[catCoverQ[Drop[catPairs[m], {k}], pats], {k, Length[catPairs[m]]}]]
```

<!-- => {False, False, False} -->

---

The check outcomes are heralds, not syndrome bits, which is why the instruction language grew
`"MH"` two sections ago: a nonzero one means discard the attempt, and the price of that
post-selection is an acceptance rate to carry alongside the answer.

```wl
cat = QECCatState[5];
AssociationMap[cat, {"Size", "Pairs", "Heralds", "DangerousPatterns", "ChecksCoverQ"}]
```

## Measurement.wl — spending the cat

The cat is the ancilla; this is the gadget that consumes it. One attempt applies a transversal
controlled-$P$ and reads the eigenvalue out with a Hadamard transform, since

$$H^{\otimes m}\left(|0\ldots0\rangle + (-1)^b |1\ldots1\rangle\right) = \sum_{\mathrm{wt}(x) \equiv b} |x\rangle$$

so the *parity of the weight* of the measured string is $b$ — $m$ bits of record, one bit of
answer.

The sandwich that applies the controlled-$P$ is the same one `generatorInstructions` uses to
measure, in the other direction. To apply a controlled-$L$, conjugate by the $U$ that sends $L$ to
$Z$: $\text{controlled-}L = (I \otimes U^{-1})\, CZ\, (I \otimes U)$, which in time order is $U$,
then $CZ$, then $U^{-1}$. For $L = X$ that is $H \cdot CZ \cdot H = CNOT$, so nothing is wasted.
The rotations are shared rather than re-derived, which is why `basisGate` and `rotation` are
`PackageScope`:

```wl
{catMeasure, pauliReps, pauliOutcome} =
    sym /@ {"catMeasureInstructions", "pauliMeasureRepetitions", "pauliMeasureOutcome"};

catMeasure[QECPauliVector["XZZXI"], 5, {6, 7, 8, 9}] // Column
```

---

One attempt still fails the measurement correctness property — a bit flip on a cat qubit before
its readout flips the reported $b$. The answer is $2t+1$ attempts with a *fresh* cat each time and
a majority vote (Theorem 12.1), since one fault can corrupt only one cat:

```wl
{pauliReps /@ {1, 3, 5, 7}, QECPauliMeasurement[five, "ZZZZZ"]["Repetitions"]}
```

<!-- => {{1, 3, 5, 7}, 3} -->

---

Now the number the whole construction exists for, and it needs two corrections before it means
anything.

The residual is read **modulo $P$**. $X^{\otimes m}$ stabilises the cat, so a cat pattern and its
complement are one physical error; through the transversal $CZ$ those two differ on the data by a
$Z$ on every qubit of the support, which after the un-rotation is exactly $P$. The frame
propagator tracks literal Paulis and cannot know the ancilla's own stabilizer, so without the
quotient it reports a weight-$m$ residual for a fault that does nothing at all.

And only **accepted** runs count: a fault early in the cat preparation does propagate along the
chain, which is what checking is for, and the check fires.

```wl
mz = QECPauliMeasurement[five, "ZZZZZ"];
{mz["TransversalQ"], mz["DataWeights"]}
```

<!-- => {True, {0, 1}} -->

---

Against the bare ancilla measuring the *same* operator, in the same units:

```wl
pauliWeights = sym["pauliMeasureDataWeights"];
{
    Max @ pauliWeights[genInstr[First[five["CheckMatrix"]], 5, 6], 6, 5, ConstantArray[0, 10], {}],
    QECPauliMeasurement[five, "XZZXI", 1]["MaxDataWeight"]
}
```

<!-- => {4, 1} -->

---

Two checks that the quotient is forced rather than convenient. The accepted weights come in
complementary pairs $w$ and $m - w$, which is the fingerprint of $X^{\otimes m}$; and dropping the
herald filter brings weight two back, so the cat checks are visibly load-bearing:

```wl
mg = QECPauliMeasurement[five, "XZZXI", 1];
Sort @ pauliWeights[mg["Instructions"], mg["Qubits"], 5, ConstantArray[0, 10], {}]
```

<!-- => {0, 1, 3, 4} -->

---

What repetition does **not** fix is the honest boundary. If the codeword carries an error $E$ with
$EP = -PE$, every repetition reads the flipped eigenvalue and the majority is confidently wrong —
"no matter how many times we repeat it". Compare a fault inside one repetition against a data
error that anticommutes with $P$:

```wl
parities[m_, f_] := Mod[Total /@ Partition[framePropagate[m["Instructions"], m["Qubits"], f]["Record"], m["Weight"]], 2];

With[{i = FirstPosition[mz["Instructions"], {"M", 6}][[1]]},
    {parities[mz, {{i - 1, 6, {1, 0}}}], parities[mz, {{0, 1, {1, 0}}}], parities[mz, {{0, 1, {0, 1}}}]}]
```

<!-- => {{1, 0, 0}, {1, 1, 1}, {0, 0, 0}} -->

One repetition lost, majority intact; every repetition lost, majority wrong; and a commuting error
invisible, as it should be. That middle column is why Theorem 12.1 needs error correction
interspersed and not just more repetitions.

## ErrorCorrection.wl — Steane EC, and why it is the one that fits

Section 12.3 opens on cost: Shor EC "involves a lot of locations. Lots of locations means lots of
opportunities for errors, which eventually will translate into a lousy threshold." Steane EC moves
that work into preparing the ancilla, and the asymmetry that licenses it is the idea of the whole
chapter — an ancilla can be checked and a data block cannot, "the reason is that we know the
precise ancilla state that we are trying to create, but the state of the logical qubits somewhere
in the middle of a long computation is unknown".

It is CSS-only, and for a reason rather than by restriction: the construction runs on transversal
CNOT being the logical CNOT.

```wl
{steaneRegions, steaneWeights, logicalWires, encZero, encPlus, applyGates} =
    sym /@ {"steaneRegions", "steaneResidualWeights", "codeLogicalWires",
            "codeEncodedZeroInstructions", "codeEncodedPlusInstructions", "applyGates"};

ec = QECErrorCorrection[QECCode["SteaneCode"]];
AssociationMap[ec, {"DataQubits", "AncillaQubits", "Qubits", "Measurements", "TransversalQ"}]
```

---

Two halves. The bit-flip half puts the ancilla in encoded $|+\rangle$ and runs transversal CNOT
from data to ancilla: on the encoded states $CNOT|\psi\rangle|+\rangle = |\psi\rangle|+\rangle$, so
nothing happens logically, but the gate still propagates errors and bit flips land on the same
positions of the ancilla. The phase half uses encoded $|0\rangle$, the CNOT the other way, and a
transversal Hadamard so the same $Z$-basis readout sees phase errors:

```wl
{Take[ec["BitFlipInstructions"], -14], Take[ec["PhaseInstructions"], -21]} // Column
```

---

The ancilla states are checked against the engine rather than asserted. Which wire carries the
logical qubit going into the encoder is derived the same way:

```wl
prepared[instrs_, n_] := applyGates[PauliStabilizer[n], engineGates[instrs]];

With[{a = codeData[QECCode["SteaneCode"]], st = QECCode["SteaneCode"]},
    {
        logicalWires[a],
        With[{s = prepared[encZero[a, 0], 7]},
            {AllTrue[st["Generators"], s["Expectation", #] === 1 &], s["Expectation", First[st["LogicalZ"]]]}],
        With[{s = prepared[encPlus[a, 0], 7]},
            {AllTrue[st["Generators"], s["Expectation", #] === 1 &], s["Expectation", First[st["LogicalX"]]]}]
    }]
```

<!-- => {{7}, {True, 1}, {True, 1}} -->

---

Now the property that makes this the right sub-gadget for the measurement's slot. Section 12.3.3
asks whether Steane EC needs repeating and answers no, because the errors are **additive**: the
measured classical error is $e + f + g$, for $e$ the true data error, $f$ the ancilla's and $g$ the
measurement's, so correcting by $e + f + g$ leaves $f + g$ — "a single-qubit error in the ancilla
can only produce a single-qubit error in the final state. This is in contrast to Shor EC, where a
single ancilla error changing one bit of the error syndrome could totally change the error we
deduce."

That is a claim to check, not to quote. The gadget splits in two regions and reports each:

```wl
{ec["RepetitionsNeeded"], ec["DataWeights"], ec["PreparationDataWeights"]}
```

<!-- => {1, {0, 1}, {0, 1, 2, 3, 4}} -->

---

One fault in the interaction, one data error. One fault in the ancilla **preparation**, up to four
— because that preparation is the code's own non-fault-tolerant encoder. Got26 is equally explicit
about both halves of that: "How do we make the ancilla states and ensure they do not have too many
errors? That is the complicated part of Steane EC", and it is chapter 13 work, needed for a full FT
protocol whatever EC gadget is chosen. So the gadget names the assumption rather than implying it:

```wl
ec["OpenAssumptions"]
```

---

Filling the slot completes Theorem 12.1's structure. Without it the measurement handles ancilla
faults only; with it, the composed gadget still leaves one data error per fault, and inherits the
sub-gadget's open assumption rather than hiding it behind that `True`:

```wl
steane = QECCode["SteaneCode"];
{
    QECPauliMeasurement[steane, First[steane["LogicalZ"]]]["MeasurementCorrectQ"],
    QECPauliMeasurement[steane, First[steane["LogicalZ"]], "ErrorCorrection" -> "Steane"]["MeasurementCorrectQ"],
    QECPauliMeasurement[steane, First[steane["LogicalZ"]], "ErrorCorrection" -> "Steane"]["DataWeights"]
}
```

<!-- => {False, True, {0, 1}} -->

---

And the restriction is enforced rather than documented. The five-qubit code is not CSS, so it can
have a cat measurement but not Steane EC:

```wl
{QECCode["5QubitCode"]["CSSQ"], steane["CSSQ"]}
```

<!-- => {False, True} -->

## The extraction, and the one thing that cannot cross over

The three gadgets above are built and checked, but the rate machinery routes through
`codeCircuitInstructions`, which is the bare-ancilla circuit and nothing else. `codeExtraction` is
what makes that a choice: it returns the instructions, the qubit count, and the number of record
bits each check contributes.

```wl
{codeExtraction, recordSyndromes} = sym /@ {"codeExtraction", "recordSyndromes"};

{
    KeyDrop[codeExtraction[codeData[bitFlip], 1, "BareAncilla"], "Instructions"],
    KeyDrop[codeExtraction[codeData[bitFlip], 1, "Transversal"], "Instructions"]
}
```

<!-- => {<|"Qubits" -> 5, "BlockSizes" -> {1, 1}|>, <|"Qubits" -> 6, "BlockSizes" -> {2, 2}|>} -->

---

Only half of Theorem 12.1's gadget can cross over, and the obstruction is structural. **A detector
error model is a matrix** — the effect of a set of faults is the XOR of their rows — and the exact
fold, the decoder table and the Stim writer all rest on that. A **majority** vote is not linear, so
two faults whose individual effects are known do not have an effect given by their XOR once a vote
is in the path.

This is the same fork the package already stands on one side of: §12.2.2 answers a noisy syndrome
by repeat-and-agree, and `DetectorModel.wl` takes differencing instead. `QECPauliMeasurement` is
the book's gadget and keeps the vote; the extraction takes **transversality**, which is the half
that cured the hook error in the first place. Its readout is a *parity*, which is linear, and
`recordSyndromes` is where that reduction happens — the identity when a check contributes one bit,
so the bare-ancilla path is untouched:

```wl
{recordSyndromes[{1, 0, 1, 1}, {1, 1}, 2],
 recordSyndromes[{1, 1, 0, 1, 0, 0, 1, 0}, {2, 2}, 2],
 recordSyndromes[{1, 1, 1, 0, 0, 1}, {3, 1, 2}, 1]}
```

<!-- => {{1, 0, 1, 1}, {0, 1, 0, 1}, {1, 0, 1}} -->

---

That linearity is a claim, so it is checked rather than argued: the effect of two faults against
the XOR of their rows, over two hundred random pairs.

```wl
Module[{a, ex, instr, nq, blocks, eff, faults, pairs},
    a = codeData[bitFlip];
    ex = codeExtraction[a, 1, "Transversal"];
    instr = ex["Instructions"]; nq = ex["Qubits"]; blocks = ex["BlockSizes"];
    eff[f_] := faultEffect[a, instr, nq, 2, 1, f, blocks];
    faults = Flatten[Table[{i, q, pauli}, {i, 0, Length[instr]}, {q, nq},
        {pauli, {{1, 0}, {1, 1}, {0, 1}}}], 2];
    pairs = BlockRandom[SeedRandom[20260914]; RandomSample[Subsets[faults, {2}], 200]];
    AllTrue[pairs, With[{c = eff[#], u = eff[{#[[1]]}], v = eff[{#[[2]]}]},
        c[[1]] === BitXor[u[[1]], v[[1]]] && c[[2]] === BitXor[u[[2]], v[[2]]] &&
        c[[3]] === BitXor[u[[3]], v[[3]]]] &]]
```

<!-- => True -->

---

One change was needed on the decoder side before any of this showed up in a rate, and it is the
kind of bug that hides behind a plausible number. `demDecoderTable` built its hypothesis space from
every row of the model. But the decoder only ever runs on a shot that was **accepted**, so a fault
that trips a verification check cannot have happened in it; offering that fault as an explanation
lets the lightest-set rule claim a detector pattern on behalf of something the experiment already
threw away, and since the rule is first-come, the claim displaces the real explanation.

Before the fix the transversal extraction measured an exponent of 1.73 and looked like it had a
residual linear term. It did not. The table now drops the herald-tripping rows, and the same
condition is visible in the model itself:

```wl
Module[{dem, d, o, h, keep, grouped},
    dem = QECDetectorModel[five, QECNoiseModel["Circuit", 1/1000], 1, "Extraction" -> "Transversal"];
    d = dem["DetectorMatrix"]; o = dem["ObservableMatrix"]; h = dem["HeraldMatrix"];
    keep = Select[Range[Length[d]], Total[h[[#]]] === 0 &];
    grouped = GroupBy[keep, d[[#]] &, DeleteDuplicates[o[[#]] & /@ #] &];
    <|"all rows" -> Length[d], "accepted" -> Length[keep],
      "ambiguous among accepted" -> Count[grouped, alt_ /; Length[alt] > 1]|>]
```

<!-- => <|"all rows" -> 624, "accepted" -> 484, "ambiguous among accepted" -> 0|> -->

Zero ambiguous signatures among accepted faults is the fault-tolerance property stated in the
detector model's own language, and it is the same conditional statement `QECPauliMeasurement` makes
about its residual data weight. Both are wrong if the condition is dropped on one side only, which
is exactly what the decoder was doing.

## Register.wl — blocks, and the layout that has to be right

A logical qubit is a block, so a logical two-qubit gate is a gate between two blocks. `QECRegister`
fixes the layout — block $b$ owns qubits $(b-1)n+1 \ldots bn$ — and `registerIndex` is the one line
that says so, stated once so no consumer reconstructs it.

```wl
{registerIndex, registerLift, registerLabelMatrix, registerLogicalVectors,
 registerTransversalCNOT, registerConjugate} =
    sym /@ {"registerIndex", "registerLift", "registerLabelMatrix", "registerLogicalVectors",
            "registerTransversalCNOT", "registerConjugate"};

{registerIndex[7, 1, 3], registerIndex[7, 2, 3], registerLift[{{"CNOT", 1, 2}}, 7, 2]}
```

<!-- => {3, 10, {{"CNOT", 8, 9}}} -->

---

The label matrix is where this gets delicate, and the delicacy is inherited rather than invented.
A code's label matrix has its $X$ and $Z$ halves **exchanged**, so that a plain matrix product
computes symplectic products. Lifting a block's rows into a register therefore scatters each row
into **two windows**, one in each half of the register's columns, not into one contiguous run.
Getting that wrong produces a matrix of exactly the right shape that quietly reports the wrong
syndromes — which is the worst kind of wrong.

```wl
Module[{steane = QECCode["SteaneCode"]},
    {Dimensions[codeLabelMatrix[codeData[steane]]],
     Dimensions[registerLabelMatrix[codeData[steane], 2]],
     registerLabelMatrix[codeData[steane], 1] === codeLabelMatrix[codeData[steane]]}]
```

<!-- => {{8, 14}, {16, 28}, True} -->

One block is the identity, which is the regression that keeps the rest of the package out of this
file's way.

---

`registerConjugate` is how a claim about a gate gets checked instead of cited: put a Pauli in as the
frame before the circuit starts, push it through, read the frame afterwards. Phases are dropped, as
everywhere in the frame propagator, so what comes back is the action on the Pauli group modulo
sign — which is exactly what "is this the logical CNOT" asks.

```wl
Module[{reg = QECRegister[QECCode["SteaneCode"], 2]},
    reg["LogicalAction", 1, 2]]
```

Steane EC rested on that answer being the CNOT action, and the header of `ErrorCorrection.wl` cites
§12.3.1 for it. This is where it is computed.

---

And the second half, which is the one that actually decides. A gate can act correctly on the
logical operators and still take the state out of the code space; both have to hold. On the
five-qubit code the $Z$ half of the action comes out right and the $X$ half does not, so a test
that looked only at the logical operators would have passed half the time:

```wl
Module[{leaves},
    leaves[r_] := Module[{instr = r["TransversalCNOT", 1, 2], nq = r["Qubits"], group},
        group = QECCode[r["Generators"]];
        Count[QECPauliString[registerConjugate[instr, nq, QECPauliVector[#]]] & /@ r["Generators"],
              g_ /; ! group["StabilizerMemberQ", g]]];
    <|"Steane" -> leaves[QECRegister[QECCode["SteaneCode"], 2]],
      "five-qubit" -> leaves[QECRegister[QECCode["5QubitCode"], 2]]|>]
```

<!-- => <|"Steane" -> 0, "five-qubit" -> 8|> -->

## Transversal.wl — tables that are computed, not typed

The criterion of sec. 11.2 is cheap to state and easy to get subtly wrong: $U$ applied to every
qubit separately is a gate gadget when it maps the stabilizer group to itself, and the logical gate
is what it does to $\bar{X}$ and $\bar{Z}$. The subtlety is that both halves are statements **about
signs**, and every conjugation in the layer up to this file threw signs away.

So the first decision here is a negative one: `framePropagate` is not used. Rows carry their
$\mathbb{Z}_4$ phase (`Pauli.wl`), products go through `QECPauliProduct`, and membership is asked of
`codeStabilizerElement`, which returns the group element *with its true phase* rather than a yes
about its symplectic part. An image equal to minus a generator preserves the stabilizer as a set of
Paulis and destroys the code space; only the phase-aware test can tell those apart.

The second decision is that the action tables are derived. A hand-written `H: X -> Z, S: Y -> -X`
table is exactly the kind of thing that acquires one wrong sign and then reports plausible answers
forever, so `transversalAction` computes each table once from the gate's own matrix — conjugate each
Pauli, read the coefficient off a trace, and insist that it is a fourth root of unity:

```wl
{transversalAction, transversalConjugate, transversalPairCode, transversalGateNames,
 transversalCliffordWords} =
    sym /@ {"transversalAction", "transversalConjugate", "transversalPairCode",
            "transversalGateNames", "transversalCliffordWords"};

QECPauliString /@ Lookup[transversalAction["S", 1], {{1, 0}, {0, 1}, {1, 1}}]
```

<!-- => {"Y", "Z", "-X"} -->

Adding a gate to this file is adding a matrix. The matrices are the engine's, and the tests compare
them against `QuantumOperator[name]["MatrixRepresentation"]` so that a drift in either convention
fails loudly instead of quietly flipping every sign below.

---

Conjugating a whole row is then qubit by qubit: the factors live on different qubits, so they
commute, and their phases simply add.

```wl
{QECPauliString[transversalConjugate["S", QECPauliVector["IIIXXXX"]]],
 QECPauliString[transversalConjugate["S", QECPauliVector["XXXXXXX"]]],
 QECPauliString[transversalConjugate["H", QECPauliVector["XZZXI"]]]}
```

<!-- => {"IIIYYYY", "YYYYYYY", "ZXXZI"} -->

Note where the minus of sec. 11.3 is *not*: the physical image of $X^{\otimes 7}$ is $Y^{\otimes 7}$
with phase zero. The sign appears only when that image is divided by the logical representative,
because $\bar{Y}$ is $-Y^{\otimes 7}$ and not $Y^{\otimes 7}$. `transversalDecompose` is what does
the dividing: the exponents of $\bar{X}_j$ and $\bar{Z}_j$ come from commutation — anticommuting
with $\bar{Z}_j$ is what an $\bar{X}_j$ factor does — and the sign is whatever is left once the
logical part and the stabilizer part have been divided out.

---

A gate gadget is checked generator by generator, and the check is reported rather than reduced to a
boolean, because the interesting failures are partial:

```wl
Normal[QECTransversalGate[QECCode["SteaneCode"], "S"]["StabilizerImages"]][[1]]
```

<!-- => <|"Generator" -> "IIIXXXX", "Image" -> "IIIYYYY", "StabilizerQ" -> True|> -->

$X^{\otimes 4}$ goes to $Y^{\otimes 4}$, which is in the group with the right sign because every
element of this stabilizer has weight 4 — the book's "the phases will take care of themselves",
checked instead of assumed.

---

Two blocks reuse the machinery rather than duplicating it. `transversalPairCode` builds the direct
sum of two copies of the stabilizer as an ordinary code association, which is what lets
`codeStabilizerElement` answer membership for the pair with no second implementation, and the
logical operators are taken block by block rather than from the pair's own standard form — the
question is what the gate does to *this* block's qubit and *that* one's.

```wl
With[{pc = transversalPairCode[codeData[QECCode["SteaneCode"]]]},
    {Length[pc["CheckMatrix"]], pc["Qubits"]}]
```

<!-- => {12, 14} -->

The qubit layout is `Register.wl`'s, so a two-block gate emitted here drops straight into a
register.

---

Naming the result is a lookup, not a classification. The twenty-four one-qubit Cliffords are
enumerated as words in $H$ and $S$, keyed by what they do to $X$ and $Z$; the named gates are added
last so that they win the key they share with a word, and the transversal $S$ of the seven-qubit
code is reported as `"Sdg"` rather than as `"SSS"`.

```wl
{Length[transversalCliffordWords[1]], Length[transversalGateNames[1]], Length[transversalGateNames[2]]}
```

<!-- => {24, 24, 3} -->

Two consequences worth naming. A code with $k > 1$ gets its logical action reported and its logical
gate left as `Missing["NotNamed"]` — on the $[[4,2,2]]$ code the transversal $H$ is a Hadamard on
each logical qubit *and* a swap of the two, and there is no one-qubit name for that. And one
vocabulary change leaked out of this file: `"Sdg"` had to join the circuit ops, where it moves the
frame exactly as `"S"` does and is emitted to the engine as three `"S"`, the same treatment `"Vdg"`
already had.

## FaultTolerant.wl — the substitution, and the bookkeeping that goes with it

Def 10.6 is short, and almost all of this file is the bookkeeping it implies. Three parts of that
bookkeeping are where a hand-assembled protocol goes wrong.

**Which gadget.** Not the one with the same name. `ftGateWord` searches the transversal gates of the
code for the one whose *logical action* is the gate requested, named gates first so that a single
instruction wins over a three-instruction word that does the same thing:

```wl
{ftGateWord, ftRelocate, ftLocations, ftReadoutSupport} =
    sym /@ {"ftGateWord", "ftRelocate", "ftLocations", "ftReadoutSupport"};

{ftGateWord[codeData[QECCode["SteaneCode"]], "S"],
 ftGateWord[codeData[QECCode["SteaneCode"]], "Sdg"],
 ftGateWord[codeData[QECCode["SteaneCode"]], "T"]}
```

<!-- => {"Sdg", "S", Missing["NoGadget", "T"]} -->

The gadget for a logical $S$ is the transversal $S^{\dagger}$ and vice versa (sec. 11.3), and a gate
outside the Clifford group has no gadget at all, which is reported as a `Missing` rather than
approximated. This is the whole payoff of carrying phases through `Transversal.wl`: without them the
search cannot tell the two apart and the assembler ships the conjugate circuit.

---

**Where the gadget goes.** The error correction gadget is written for data on qubits $1 \ldots n$
with its ancilla at an offset, because that is all it ever needed. Here the data is block $b$ and
the ancilla is that block's own workspace. Every instruction of that gadget touches either a data
qubit or an ancilla qubit, so relocating it is a two-branch map rather than a rewrite:

```wl
ftRelocate[{{"CNOT", 1, 8}}, 7, 7, 7, 21]
```

<!-- => {{"CNOT", 8, 22}} -->

Data qubit 1 becomes qubit 1 of block 2; ancilla qubit 8 becomes qubit 1 of that block's workspace.
One workspace per block rather than one shared one, so that the corrections of a layer can run in
parallel — a shared ancilla would serialise them and inflate the depth overhead for no physical
reason.

---

**When it goes there.** The assembly walks the *schedule* of $C$, not its instruction list. Def 10.6
puts a correction between every adjacent pair of locations, and in a circuit with parallel gates
"adjacent" means a layer. Within a layer: emit each location's gadget; then give every live block
nobody touched a **storage** gadget, which emits nothing and is still a location; then put an error
correction gadget on every block still live at the end of the layer. A block measured in that layer
is done, and gets none.

That is also why the size overhead is computed with `ftLocations` rather than `Length`: a location
is an instruction *or* a wait, and the waits come from `circuitIdleSlots`, the same source the
detector model and the Stim writer use.

```wl
ftLocations[{{"R", 1}, {"R", 2}, {"H", 1}, {"CNOT", 1, 2}, {"M", 1}, {"M", 2}}, 2]
```

<!-- => 7 -->

Six instructions and one wait — block 2 idling through the Hadamard.

---

**What comes back out.** The measurement record is tracked while emitting rather than recovered
afterwards, because a measurement gadget's $n$ outcomes sit in the middle of a record otherwise made
of error-correction outcomes. `"Readouts"` gives, per logical qubit, the absolute positions whose
parity is the logical bit; their offsets inside the block are the support of $\bar{Z}$, which a
transversal $Z$ measurement is only able to read because that operator is a product of $Z$s:

```wl
ftReadoutSupport[codeData[QECCode["SteaneCode"]]]
```

<!-- => {2, 4, 6} -->

And the gadgets can be taken apart again: the table is in emission order and each row carries its
own length, so `ft["GadgetInstructions", "Gate"]` slices the assembled circuit by cumulative sums
and never has to recognise a gadget by looking at its gates. That is what lets a check run the gate
gadgets alone — preparation and gates, no corrections — and ask the engine whether the encoded
blocks came out in the state the ideal circuit would have produced.

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
