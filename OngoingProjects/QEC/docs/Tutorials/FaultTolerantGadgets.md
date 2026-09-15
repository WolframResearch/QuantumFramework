---
Template: TechNote
Name: FaultTolerantGadgets
Title: Fault-Tolerant Gadgets, from the Cat State to FT(C)
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/FaultTolerantGadgets
Keywords: [quantum error correction, fault tolerance, cat state, verified ancilla, fault-tolerant measurement, Steane error correction, transversal gate, logical gate, code block, register, fault-tolerant simulation, gadget, overhead, Clifford group]
RelatedGuides: [WolframQuantumComputationFramework]
---

A quantum error-correcting code tells you which errors you could detect. It does not tell you how
to detect them without causing more, and that gap is the whole subject of fault tolerance: the
circuit that measures a check is built from the same failing gates as everything else, and one bad
gate inside it can spread into many bad qubits on the data.

This note is about the six gadgets that close that gap, in the order they build on each other.
Each one is an object you can hold, ask questions of, and hand to the next:

| | gadget | what it is |
|---|---|---|
| 1 | `QECCatState` | an ancilla worth trusting, because it has been checked |
| 2 | `QECPauliMeasurement` | spend the ancilla: measure one Pauli, fault tolerantly |
| 3 | `QECErrorCorrection` | do that for a whole code at once — Steane EC |
| 4 | `QECRegister` | several blocks, so there is something to compute *with* |
| 5 | `QECTransversalGate` | which logical gates the code performs for free |
| 6 | `QECFaultTolerant` | all of the above, spent on a whole circuit: $FT(C)$ |

The thread running through them is that every claim is **measured rather than cited**. "This gadget
is fault tolerant" becomes a number: how many data errors one fault leaves behind. "This is the
logical CNOT" becomes a conjugation you can read off. "This protocol costs more" becomes the three
overheads of Definition 10.6. And where a gadget does *not* deliver what its name suggests, it says
so — one of the objects below reports `False` for exactly the reason the book gives.

## Definition

The code layer is a development package under `OngoingProjects/QEC/`, not yet part of the paclet,
so this page loads both the framework and the layer itself. Symbols introduced in this section stay
bound for the whole note.

```wl
Needs["Wolfram`QuantumFramework`"];

qecCore = SelectFirst[
    {
        (* beside this notebook: docs/Tutorials/ up two levels *)
        Quiet @ Check[FileNameJoin[{ParentDirectory[NotebookDirectory[], 2], "QECCore", "QECCore.wl"}], $Failed],
        (* from a repo checkout of the framework, when that is the one loaded *)
        FileNameJoin[{ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
            "OngoingProjects", "QEC", "QECCore", "QECCore.wl"}],
        (* from the working directory *)
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

Two codes carry the whole note, and they are chosen to disagree with each other:

```wl
steane = QECCode["SteaneCode"];
five = QECCode["5QubitCode"];
{steane["Parameters"], five["Parameters"], {steane["CSSQ"], five["CSSQ"]}}
```

<!-- => {{7, 1, 3}, {5, 1, 3}, {True, False}} -->

Same distance, same single logical qubit, and one of them is CSS. Most of what separates the two
below comes from that last bit.

## 1. `QECCatState` — an ancilla worth trusting

The problem is one line long. To measure a weight-$m$ check you couple $m$ data qubits to an
ancilla; if that ancilla is a single qubit, one fault on it spreads to every data qubit it touches
afterwards. The fix (sec. 12.1.2) is to spread the ancilla out first: $m$ qubits in the cat state
$|0\ldots0\rangle + |1\ldots1\rangle$, one per data qubit, so each ancilla qubit touches exactly one
of them.

That only moves the problem, because preparing the cat state is itself a circuit, and one fault in
it can produce exactly the correlated error the construction was meant to avoid. So the cat is
**verified**: parity checks between pairs of its qubits, each one a herald that rejects the shot.

```wl
cat = QECCatState[4]
```

<!-- => a QECCatState summary box: 4 cat qubits, check qubit 5, one check pair -->

The object carries the circuit and the bookkeeping around it:

```wl
AssociationMap[cat, {"Size", "CatQubits", "CheckQubit", "Qubits", "Pairs", "Repetitions"}]
```

<!-- => <|"Size" -> 4, "CatQubits" -> {1, 2, 3, 4}, "CheckQubit" -> 5, "Qubits" -> 5, "Pairs" -> {{2, 3}}, "Repetitions" -> 1|> -->

And the circuit is short enough to read in full:

```wl
cat["Instructions"]
```

<!-- => {{"R",1},{"R",2},{"R",3},{"R",4},{"H",1},{"CNOT",1,2},{"CNOT",2,3},{"CNOT",3,4},{"R",5},{"CNOT",2,5},{"CNOT",3,5},{"MH",5}} -->

Four resets, a Hadamard and a chain of CNOTs — that is the cat. Then a fifth qubit, two CNOTs from
the pair being compared, and an `MH`: a measurement whose outcome is a *herald*, deciding whether
the shot is kept at all.

---

Now the part that is derived rather than assumed. A single fault during preparation leaves one of a
small set of patterns on the cat, and only some of them are dangerous: $X^{\otimes m}$ stabilises
the cat state, so a pattern and its complement are the same physical error, and only patterns of
canonical weight $\geq 2$ can hurt. The object computes that set, and whether the checks it was
given actually cover it:

```wl
<|"dangerous" -> cat["DangerousPatterns"], "checks" -> cat["Pairs"],
  "covered" -> cat["ChecksCoverQ"]|>
```

<!-- => <|"dangerous" -> {{0, 0, 1, 1}}, "checks" -> {{2, 3}}, "covered" -> True|> -->

For $m = 4$ the one dangerous pattern is $\{0,0,1,1\}$ — exactly the error drawn in the book's
figure 12.3b — and comparing qubits 2 and 3 catches it. Sweeping the width gives the rule:

```wl
Dataset @ Association @ Table[
    m -> <|"checks" -> Length[QECCatState[m]["Pairs"]],
           "dangerous" -> Length[QECCatState[m]["DangerousPatterns"]],
           "covered" -> QECCatState[m]["ChecksCoverQ"]|>,
    {m, 2, 6}]
```

<!-- => m=2: 0 checks, 0 dangerous; m=3: 0, 0; m=4: 1, 1; m=5: 2, 2; m=6: 3, 3 — covered True throughout -->

$m - 3$ checks suffice, and $m = 2$ and $m = 3$ need none at all. A small result, and the kind this
layer exists for: the book says a few checks are enough, and the object says how few.

---

The gadget can be placed anywhere, which is what the next section needs — the cat has to sit
somewhere above the data qubits:

```wl
QECCatState[{8, 9, 10, 11}, 12]["Instructions"]
```

<!-- => the same twelve instructions on qubits 8..11 with 12 as the check qubit -->

## 2. `QECPauliMeasurement` — spending the cat

With a verified cat, measuring a Pauli $P$ of weight $m$ goes: prepare the cat, apply a
controlled-$P$ from each cat qubit to its own data qubit, read the cat out in the $X$ basis, and take
the **parity** of the outcomes — that parity is the eigenvalue (eqs. 12.1–12.3). Because a single
round can still lie, the whole thing is repeated $2t+1$ times and a majority taken (Theorem 12.1).

```wl
meas = QECPauliMeasurement[five, "XZZXI"]
```

<!-- => a QECPauliMeasurement summary box: the Pauli XZZXI on the five-qubit code -->

```wl
AssociationMap[meas, {"Pauli", "Support", "Weight", "Repetitions", "CatRepetitions",
    "DataQubits", "CatQubits", "CheckQubit", "Qubits", "Measurements", "Depth"}]
```

<!-- => Pauli "XZZXI", Support {1,2,3,4}, Weight 4, Repetitions 3, DataQubits 5, CatQubits {6,7,8,9}, CheckQubit 10, Qubits 10, Measurements 12, Depth 23 -->

Three repetitions, because the five-qubit code has distance 3, so $t = 1$ and $2t+1 = 3$. That count
is derived from the code rather than typed in — and it can be overridden when you want to see what a
single unrepeated round does:

```wl
{QECPauliMeasurement[five, "XZZXI"]["Repetitions"],
 QECPauliMeasurement[five, "XZZXI", 1]["Repetitions"]}
```

<!-- => {3, 1} -->

---

Here is the number the gadget exists for. `"DataWeights"` is the set of residual data errors left by
**every single fault** anywhere in the circuit, conditional on the heralds accepting, and read
modulo $P$ itself:

```wl
<|"data weights" -> meas["DataWeights"], "worst" -> meas["MaxDataWeight"],
  "transversal" -> meas["TransversalQ"], "heralds" -> meas["Heralds"]|>
```

<!-- => <|"data weights" -> {0, 1}, "worst" -> 1, "transversal" -> True, "heralds" -> 3|> -->

**One fault, at most one data error.** The bare-ancilla circuit measuring the same generator leaves
**four** — that is the hook error, a defect of the gadget rather than a fact about circuit-level
noise, and removing it is what this object is for. The three heralds are the price: a shot can be
rejected, so the acceptance probability belongs in any rate computed from this circuit.

---

And now the thing worth pausing on, because the object refuses to overclaim:

```wl
<|"majority is right" -> meas["MeasurementCorrectQ"],
  "error correction between repetitions" -> meas["ErrorCorrection"]|>
```

<!-- => <|"majority is right" -> False, "error correction between repetitions" -> None|> -->

`False`. Repetition plus majority handles faults in the *ancilla*, and nothing else. A data error
that anticommutes with $P$ flips the true answer, so every repetition sees the same wrong value and
the majority is confidently wrong (sec. 12.1.4). The object says which assumption it is standing on:

```wl
Column[meas["OpenAssumptions"]]
```

<!-- => the sec. 12.1.4 assumption: no error correction is spliced between repetitions -->

That is also what Theorem 12.1 says to do about it: interleave fault-tolerant error correction
between the repetitions. Splice one in — which needs a CSS code, for reasons the next section gives
— and the answer changes:

```wl
withEC = QECPauliMeasurement[steane, "IIIXXXX", "ErrorCorrection" -> "Steane"];
<|"majority is right" -> withEC["MeasurementCorrectQ"], "repetitions" -> withEC["Repetitions"],
  "data weights" -> withEC["DataWeights"], "qubits" -> withEC["Qubits"],
  "depth" -> withEC["Depth"]|>
```

<!-- => <|"majority is right" -> True, "repetitions" -> 3, "data weights" -> {0, 1}, "qubits" -> 19, "depth" -> 77|> -->

`True`, at the price of nineteen qubits and a depth of 77 for one measurement of one generator. The
open assumption that remains is the one every gadget in this note inherits, and it is the honest
one:

```wl
Column[withEC["OpenAssumptions"]]
```

<!-- => the chapter-13 assumption: fault-tolerant preparation of the encoded ancillas -->

## 3. `QECErrorCorrection` — Steane EC

Measuring one Pauli is not error correction: for that you need every generator and then a decoder.
Section 12.3 does all of them at once for a CSS code, and with a different trick — instead of a cat,
a whole **encoded ancilla block**.

The bit-flip half prepares an encoded $|+\rangle$, applies a transversal CNOT from data to ancilla
and measures the ancilla, so bit flips in the data copy over and are read there. The phase half is
the mirror image: an encoded $|0\rangle$, the CNOT the other way, and a transversal Hadamard so the
same $Z$-basis readout sees phase errors.

```wl
ec = QECErrorCorrection[steane]
```

<!-- => a QECErrorCorrection summary box: Steane, 7 data qubits, ancilla 8..14 -->

```wl
AssociationMap[ec, {"DataQubits", "AncillaQubits", "Offset", "Order", "Qubits",
    "Measurements", "Depth", "InstructionCount", "RepetitionsNeeded"}]
```

<!-- => DataQubits 7, AncillaQubits {8,...,14}, Offset 7, Order "BitFlipFirst", Qubits 14, Measurements 14, Depth 27, InstructionCount 100, RepetitionsNeeded 1 -->

`"RepetitionsNeeded"` is 1, and that is why this gadget and not Shor's fills the slot left open in
section 2: the errors are **additive**. The measured classical error is $e + f + g$ for $e$ the true
data error, $f$ the ancilla's and $g$ the measurement's, so correcting by $e+f+g$ leaves $f+g$
behind — a single fault can only leave a single error, with no repetition and no majority vote at
all (sec. 12.3.3).

---

Measured in the same units as section 2, and split in two, because the gadget is honest about which
half is which:

```wl
<|"regions" -> Map[Length, ec["Regions"]],
  "interaction" -> ec["DataWeights"],
  "preparation" -> ec["PreparationDataWeights"]|>
```

<!-- => <|"regions" -> <|"Preparation" -> 65, "Interaction" -> 36|>, "interaction" -> {0, 1}, "preparation" -> {0, 1, 2, 3, 4}|> -->

The **interaction** — the transversal CNOTs and the readout — leaves at most one data error per
fault, which is the additivity above. The **preparation** of the ancilla leaves up to four, because
that preparation is still the code's ordinary, non-fault-tolerant encoder. That is chapter 13's work
appearing as a number instead of a caveat, and it is the open assumption every gadget in this note
carries.

---

Each half is readable on its own:

```wl
Take[ec["BitFlipInstructions"], 10]
```

<!-- => {{"R",8},...,{"R",14},{"H",14},{"H",13},{"H",12}} — resets, then the ancilla's encoder -->

and the order of the two halves is a choice with a consequence:

```wl
{QECErrorCorrection[steane]["Order"],
 QECErrorCorrection[steane, "Order" -> "PhaseFirst"]["Order"]}
```

<!-- => {"BitFlipFirst", "PhaseFirst"} -->

Either order is fault tolerant — Theorem 12.5 does not care — but the errors left at the end tend to
be of the type corrected *first*, since more locations follow it.

The construction rests on the transversal CNOT being the logical CNOT, which holds for CSS codes and
fails in general, so `QECErrorCorrection[five]` is refused with a message rather than assembled into
something that silently is not error correction. The next two sections are where that claim stops
being a citation.

## 4. `QECRegister` — a logical qubit is a block

Everything so far lives on one block. A memory experiment needs no more, which is why the sections
above never mentioned layout. Anything that *computes* does: a logical two-qubit gate is a gate
between two blocks, and Theorem 13.2 needs $2m$ blocks for a gate touching $m$ of them.

`QECRegister` fixes the layout once — block $b$ owns qubits $(b-1)n+1 \ldots bn$ — so nothing
downstream has to guess it:

```wl
reg = QECRegister[steane, 2];
AssociationMap[reg, {"Blocks", "BlockQubits", "Qubits", "StabilizerCount", "LogicalQubits", "CSSQ"}]
```

<!-- => <|"Blocks" -> 2, "BlockQubits" -> 7, "Qubits" -> 14, "StabilizerCount" -> 12, "LogicalQubits" -> 2, "CSSQ" -> True|> -->

```wl
{reg["BlockRange", 2], reg["Index", 2, 3],
 reg["Lift", {{"R", 1}, {"H", 2}, {"CNOT", 1, 2}}, 2]}
```

<!-- => {{8,9,10,11,12,13,14}, 10, {{"R",8},{"H",9},{"CNOT",8,9}}} -->

A register of one block is the identity on everything — same generators, same label matrix, same
circuits — which is what keeps the rest of the layer unaffected by this object existing.

---

The gate between blocks is the transversal CNOT: qubit by qubit, so each qubit of one block touches
exactly the corresponding qubit of the other and nothing else. That is the property that stops a
single fault from spreading inside either block.

```wl
reg["TransversalCNOT", 1, 2]
```

<!-- => {{"CNOT",1,8},{"CNOT",2,9},{"CNOT",3,10},{"CNOT",4,11},{"CNOT",5,12},{"CNOT",6,13},{"CNOT",7,14}} -->

Steane EC was already resting on this being the logical CNOT. Here it is computed, by conjugating
the logical operators through the gate and reading the images off:

```wl
reg["LogicalAction", 1, 2]
```

<!-- => <|{"X",1,1} -> "IIXIXXIIIXIXXI", {"X",2,1} -> "IIIIIIIIIXIXXI", {"Z",1,1} -> "IZIZIZIIIIIIII", {"Z",2,1} -> "IZIZIZIIZIZIZI"|> -->

Read as logical operators: $\bar{X}_1 \to \bar{X}_1\bar{X}_2$ (the image spans both blocks),
$\bar{X}_2 \to \bar{X}_2$, $\bar{Z}_1 \to \bar{Z}_1$, $\bar{Z}_2 \to \bar{Z}_1\bar{Z}_2$. That is
the action of a CNOT on the logical pair.

On the five-qubit code the same gate does something else, and *how* it fails is worth keeping in
mind: the $Z$ half of the action still comes out right and only the $X$ half breaks, so a check
that looked only at the logical operators would have passed half the time. What decides is whether
the gate maps the stabilizer group to itself — the next section's question.

## 5. `QECTransversalGate` — the gates a code performs for free

A gate applied to every qubit of a block separately cannot spread one fault into two, so transversal
gates are the cheapest fault-tolerant gates there are. Which ones a code admits is a question about
the *symmetry* of its stabilizer: $U^{\otimes n}$ is a valid gadget exactly when it maps the
stabilizer group to itself, and the logical gate you get is whatever it does to $\bar{X}$ and
$\bar{Z}$.

```wl
QECTransversalGate[steane, "H"]["LogicalAction"]
```

<!-- => <|{"X", 1} -> "Z", {"Z", 1} -> "X"|> -->

Transversal Hadamard, logical Hadamard (eq. 11.18–11.19). Now the same question for $S$:

```wl
Module[{g = QECTransversalGate[steane, "S"]},
    <|"transversal" -> g["TransversalQ"], "action" -> g["LogicalAction"],
      "logical gate" -> g["LogicalGate"]|>]
```

<!-- => <|"transversal" -> True, "action" -> <|{"X",1} -> "-Y", {"Z",1} -> "Z"|>, "logical gate" -> "Sdg"|> -->

**A minus sign, and it is a different gate.** $\bar{X} \to -\bar{Y}$ with $\bar{Z}$ fixed is the
action of $S^{\dagger}$: the transversal $S$ performs the *inverse* of the gate it is made of. The
sign is one line of Pauli algebra:

```wl
QECPauliString[QECPauliProduct["XXXXXXX", "ZZZZZZZ"]]
```

<!-- => "iYYYYYYY" -->

$\bar{X}\bar{Z} = i\,Y^{\otimes 7}$ and $\bar{Y} = i\bar{X}\bar{Z}$, so $Y^{\otimes 7} = -\bar{Y}$
(eq. 11.22). Everywhere else in this layer a conjugation modulo sign is enough; here the sign *is*
the answer, so this object carries the $\mathbb{Z}_4$ phase through every step and asks membership
of a routine that knows a generator's true phase.

---

Being a gate gadget at all is the other half, and it is reported generator by generator rather than
as one boolean, because the interesting failures are partial:

```wl
QECTransversalGate[steane, "S"]["StabilizerImages"]
```

<!-- => six rows: IIIXXXX -> IIIYYYY, XIXIXIX -> YIYIYIY, IXXIIXX -> IYYIIYY, the Z generators fixed; StabilizerQ True throughout -->

Each $X^{\otimes 4}$ goes to $Y^{\otimes 4}$, which is in the group *with the right sign* because
every element of this stabilizer has weight 4 — the book's "the phases will take care of
themselves", checked instead of assumed.

---

Asking a code which Cliffords it admits is one call:

```wl
Dataset @ <|"Steane" -> Length[QECTransversalGate[steane, All]],
            "five-qubit" -> Length[QECTransversalGate[five, All]]|>
```

<!-- => <|"Steane" -> 24, "five-qubit" -> 12|> -->

All twenty-four one-qubit Cliffords on the seven-qubit code — which is what chapter 11 singles it
out for — and exactly half of them on the five-qubit code. Half is not none, and the ones it has are
worth seeing: no $H$, no $S$, but the cyclic Clifford $X \to Y \to Z \to X$, which is $S$ followed by
$H$:

```wl
Module[{g = QECTransversalGate[five, {"S", "H"}]},
    <|"H alone" -> QECTransversalGate[five, "H"]["TransversalQ"],
      "S then H" -> g["TransversalQ"], "action" -> g["LogicalAction"]|>]
```

<!-- => <|"H alone" -> False, "S then H" -> True, "action" -> <|{"X",1} -> "-Y", {"Z",1} -> "-X"|>|> -->

---

Between two blocks the answer is the CNOT, with the signs now carried, and two more gates come free
with the same symmetry:

```wl
Dataset @ AssociationMap[QECTransversalGate[steane, #]["LogicalGate"] &, {"CNOT", "CZ", "SWAP"}]
```

<!-- => <|"CNOT" -> "CNOT", "CZ" -> "CZ", "SWAP" -> "SWAP"|> -->

On the five-qubit code the transversal CNOT is not a gadget at all, which is section 11.4's argument
run backwards: demanding it forces every generator to be $X$-only or $Z$-only, which is what being
CSS means — and therefore why Steane EC, and the spliced-in correction of section 2, are CSS-only.
$H$, $S$ and CNOT generate the Clifford group, so the seven-qubit code performs the whole logical
Clifford group transversally, each gate arriving conjugated per the sign above. Nothing outside that
group is transversal on any code of this kind; that is chapter 13's problem, not this note's.

## 6. `QECFaultTolerant` — spending all of it on a circuit

Definition 10.6 says how to put the gadgets together. Take an ideal circuit $C$, replace each of its
qubits with a **block**, replace each of its locations with the corresponding gadget, and after every
preparation, gate and storage gadget put an error correction gadget on each block involved — never
after a measurement gadget, whose output is classical.

```wl
circuit = {{"R", 1}, {"R", 2}, {"H", 1}, {"CNOT", 1, 2}, {"M", 1}, {"M", 2}};
ft = QECFaultTolerant[circuit, steane]
```

<!-- => a QECFaultTolerant summary box: 2 logical qubits, 28 qubits, 13 gadgets -->

```wl
ft["Gadgets"]
```

<!-- => a 13-row Dataset: gadget kind, the location it replaces, the blocks, the instruction count -->

```wl
ft["GadgetCounts"]
```

<!-- => <|"Preparation" -> 2, "ErrorCorrection" -> 6, "Gate" -> 2, "Storage" -> 1, "Measurement" -> 2|> -->

Two rows deserve a look. The **storage** gadget is block 2 waiting while block 1 gets its Hadamard:
a wait is a location, the book's storage gadget is "just putting a wait location for all physical
qubits in the code", and it earns a correction like any other. It emits no instructions and it still
costs. And there are **six** corrections rather than four — one per live block per layer, stopping
when a block is measured.

---

Definition 10.6 also defines what that costs, in three ratios:

```wl
ft["Overheads"]
```

<!-- => <|"Size" -> 263, "Qubits" -> 14, "Depth" -> 77/4|> -->

```wl
{ft["CircuitLocations"], ft["Locations"], ft["Depth"], ft["InstructionCount"]}
```

<!-- => {7, 1841, 77, 692} -->

Seven locations become 1841 — counted on both sides as instructions *plus waits*, which is what
keeps the first ratio honest. Fourteen physical qubits per logical one: seven for the block and
seven for its own correction workspace, one workspace per block so the corrections of a layer run in
parallel instead of queueing. And nineteen and a quarter times the depth. Those are the numbers a
fault-tolerance threshold has to beat, and they are why the subject is about *rates* rather than
counts.

---

Now the detail that makes this more than bookkeeping. Asked for a logical $S$, the assembler emits a
transversal $S^{\dagger}$:

```wl
Union @ Map[First, QECFaultTolerant[{{"R", 1}, {"S", 1}}, steane]["GadgetInstructions", "Gate"]]
```

<!-- => {"Sdg"} -->

It never emits "the transversal version of the gate it was asked for" — it asks which word has the
requested *logical* action, which is section 5 doing work. A protocol assembled by name would
compute the complex conjugate of the intended circuit and look perfectly healthy doing it: on the
engine, assembled properly the encoded qubit ends in the $+1$ eigenstate of $\bar{Y}$, and assembled
by name in the $-1$ eigenstate.

---

The classical side is tracked too. A measurement gadget's outcomes sit in a record otherwise full of
error-correction syndromes, so a decoder has to be told which ones make the logical bit:

```wl
<|"measurements" -> ft["Measurements"], "logical bit of qubit 1" -> ft["Readouts"][1],
  "of qubit 2" -> ft["Readouts"][2]|>
```

<!-- => <|"measurements" -> 98, "logical bit of qubit 1" -> {86, 88, 90}, "of qubit 2" -> {93, 95, 97}|> -->

Three outcomes out of ninety-eight, being the support of $\bar{Z} = IZIZIZI$ inside that block's own
measurements.

---

And what the protocol is not, said by the object rather than left to the reader:

```wl
Column[ft["OpenAssumptions"], Spacings -> 1]
```

<!-- => three assumptions: the encoder is not fault tolerant, the gate set is Clifford only, corrections live in the Pauli frame -->

The preparation gadget is the code's own encoder, and so are the ancillas of every correction — the
same debt section 3 measured as `{0,1,2,3,4}`. The gate gadgets are transversal, so the gate set is
Clifford: Definition 10.5 asks for a *universal* set, and a circuit with a $T$ in it is refused
rather than approximated. And corrections stay in the Pauli frame, so nothing in the emitted circuit
is classically conditioned — which is what real experiments do, and one fewer location that can
fail.

## The six, side by side

| function | takes | gives back | the one thing it is for |
|---|---|---|---|
| `QECCatState[m]` | the width $m$, or explicit qubits and a check qubit; options `"Pairs"`, `"Repetitions"` | the preparation and its checks, the dangerous patterns, whether the checks cover them | an ancilla a single fault cannot corrupt undetected |
| `QECPauliMeasurement[code, P]` | a code and a Pauli; an optional repetition count; options `"Pairs"`, `"CatRepetitions"`, `"ErrorCorrection"` | the circuit, the heralds, the residual data weights, whether the majority is right | measuring one Pauli without spreading errors |
| `QECErrorCorrection[code]` | a CSS code; an optional ancilla offset; option `"Order"` | the two halves, the two regions, and the data weights of each | correcting a whole block, with no repetition needed |
| `QECRegister[code, b]` | a code and a number of blocks | the layout, the lift, the transversal CNOT and what it does logically | making a logical qubit addressable |
| `QECTransversalGate[code, gate]` | a code and a gate, a sequence of gates, or `All` | whether it is a gadget, the images of the generators, the logical action and its name | knowing which logical gate you actually performed |
| `QECFaultTolerant[circuit, code]` | an ideal Clifford circuit on logical qubits, and a CSS code | $FT(C)$, its gadget table, its readouts, and the three overheads | turning a circuit into its fault-tolerant simulation |

Every one of them answers `obj["Properties"]` with its full list, and every one that stands on
something it has not built says so in `obj["OpenAssumptions"]`.

## Where this comes from

| claim | source |
|---|---|
| cat states as a spread-out ancilla, and their verification | Gottesman, *Surviving as a Quantum Computer in a Classical World*, sec. 12.1.2–12.1.3 |
| the parity of the cat readout is the eigenvalue | ibid., eqs. 12.1–12.3 |
| $2t+1$ repetitions, a majority vote, and FTEC interleaved between them | ibid., sec. 12.1.4–12.1.5, Theorem 12.1 |
| Steane error correction, and why not Shor's | ibid., sec. 12.3, fig. 12.8–12.9 |
| the $e+f+g$ additivity, and therefore no repetition | ibid., sec. 12.3.3 |
| either order of the two halves is fault tolerant | ibid., Theorem 12.5 |
| transversal gates as symmetries of the stabilizer | ibid., sec. 11.2 |
| the 7-qubit code's transversal Cliffords, and $Y^{\otimes 7} = -\bar{Y}$ | ibid., sec. 11.3, eq. 11.18–11.26 |
| transversal CNOT forces a CSS code | ibid., sec. 11.4 |
| gadgets, protocols and $FT(C)$, with the three overheads | ibid., Definitions 10.4–10.6, fig. 10.2 |
| blocks as logical qubits, and $2m$ of them for a gate on $m$ | ibid., Theorem 13.2 |
| corrections tracked rather than applied | ibid., sec. 12.5.2 |
| fault-tolerant preparation of encoded states — the open assumption | ibid., ch. 13 |
