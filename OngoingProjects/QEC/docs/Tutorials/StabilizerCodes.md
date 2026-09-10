---
Template: TechNote
Name: StabilizerCodes
Title: Stabilizer Codes in the Quantum Framework
Context: Wolfram`QuantumFramework`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/StabilizerCodes
Keywords: [quantum error correction, stabilizer code, syndrome, decoder, logical operator, code distance, encoding circuit, CSS code, concatenation, noise model, logical error rate, syndrome extraction, detector error model, circuit-level noise, memory experiment, Stim]
RelatedGuides: [WolframQuantumComputationFramework]
---

The framework's stabilizer engine holds a *state*: an $n$-qubit state stored as the Pauli
operators that fix it, pushed through Clifford gates in the Heisenberg picture. This layer
holds a *code*: a subspace fixed by a commuting group of Pauli checks, chosen so that the
errors we care about push a state out of that subspace in a way we can detect and undo.

The difference matters. A stabilizer state is fixed by $n$ independent checks and is a single
point. A stabilizer code is fixed by $m < n$ of them, so it is a subspace of dimension $2^{n-m}$,
and the $k = n - m$ qubits' worth of room inside it is where the protected information lives.
Everything below follows from that one shift.

What this note covers, in order: the code object and how it is stored; the Pauli layer beneath
it; the structural quantities (logical operators, distance); syndromes and decoding; encoding
circuits; building codes out of other codes; noise; and the figure of merit the whole subject is
about — whether a logical qubit survives. Then the same question asked honestly: the circuit that
measures the checks, the faults inside it, detectors, and a memory experiment run over several
rounds — ending with the whole thing written out in Stim's syntax and checked against it.

One thing here is not available anywhere else, and it is worth naming up front. Because the noise
rate may be left symbolic, the logical error rate comes back as an **exact polynomial** rather
than a sampled estimate. A fast simulator answers "about 0.0022 at p = 0.01"; this answers
"$10p^2 - 200p^3/9 + 160p^4/9 - 128p^5/27$, for every p". That is what exactness and symbolic
generality buy, and it is the one move a Monte Carlo simulator structurally cannot make. It holds
all the way up: even with every gate, reset and readout of the syndrome-extraction circuit
failing, the answer is still a closed-form polynomial, and the last sections use that to say two
things about fault tolerance that are awkward to state any other way.

## Definition

The code layer is a development package under `OngoingProjects/QEC/`, not yet part of the
paclet, so this page loads both the framework and the layer itself. Symbols introduced in this
section stay bound for the whole note.

Two notes on why this cell exists at all, since documentation pages normally have no setup code.
A reference page gets an `ExampleInitialization` cell from its template, and the converter fills
it from the frontmatter `Context:`; the tech-note template has no such cell, so the converter has
nothing to fill and a reader opening the notebook cold would find no context loaded. And the code
layer is not in any paclet yet, so no context declaration could load it. Both lines go away when
the layer moves into the paclet.

Load the framework and locate the code layer. The cell reports whether it worked, because a
setup cell that fails quietly leaves every later cell showing an unevaluated symbol — which looks
like a broken package rather than a path that did not resolve:

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

A `True` there means the layer is loaded and its symbols are on the context path, so every cell
below runs. A `Failure` means only that the file was not where this page guessed: `Get` the
`QECCore/QECCore.wl` beside it by hand and continue from the next cell.

## The code object

A code is built from its stabilizer generators, given as Pauli strings. The constructor checks
what a code needs: equal length, pairwise commuting, independent.

The five-qubit code, from its four checks:

```wl
QECCode[{"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}]
```

---

The textbook codes are available by name:

```wl
QECCode["SteaneCode"]["Generators"]
```

<!-- => {"IIIXXXX", "XIXIXIX", "IXXIIXX", "IIIZZZZ", "ZIZIZIZ", "IZZIIZZ"} -->

---

Its parameters are the $[[n, k, d]]$ of the literature — physical qubits, logical qubits, distance:

```wl
QECCode["SteaneCode"]["Parameters"]
```

<!-- => {7, 1, 3} -->

A code is stored the way the engine stores a tableau, and for the same reason: not as Pauli
strings but as a check matrix of binary symplectic rows plus a separate list of phases. Strings
are an input and output convention, never the thing computations run on.

The check matrix of the Steane code, X half on the left, Z half on the right:

```wl
ArrayPlot[QECCode["SteaneCode"]["CheckMatrix"], Mesh -> All, ImageSize -> 260]
```

The block structure is visible: three generators with support only in the X half, three only in
the Z half. That is what makes it a CSS code, and it is a property the object will confirm
directly.

```wl
QECCode["SteaneCode"]["CSSQ"]
```

<!-- => True -->

## Codes and states

It is worth being precise about where this layer sits, because the engine already represents
stabilizer *states* and the two objects are easy to confuse.

A [PauliStabilizer]() holds $n$ independent commuting generators on $n$ qubits. That pins a single
state — there is exactly one joint $+1$ eigenvector — and the engine's job is to push it through
Clifford gates and measure it. A `QECCode` holds $m < n$ generators. That does not pin a state; it
pins a $2^{n-m}$-dimensional subspace, and the questions worth asking about it are different:
which operators act on the space without leaving it, how far apart its members are, what circuit
prepares one.

The two meet constantly, and deliberately in one direction: the code hands work to the engine
rather than reimplementing it. Completing a code's generators to a full set of $n$ pins one
fiducial codeword, and that codeword is an engine object:

```wl
Head[QECCode["SteaneCode"]["PauliStabilizer"]]
```

<!-- => PauliStabilizer -->

---

That is how an encoding circuit is validated, how the physical correction cycle runs, and how
syndromes are cross-checked: the answer computed here symplectically is compared against the
answer the engine measures on an actual state. The two agree because they share a row format, not
because they were written to agree.

The division of labour also runs the other way. Anything genuinely about a state at scale — long
Clifford circuits, measurement, entanglement entropy — belongs in the engine, where it has a
compiled bulk path; the encoder here delegates its gate application to exactly that path rather
than folding gates one at a time.

## Paulis and phases

Underneath the code object, a Pauli is a row of bits: $x_1 \ldots x_n$, then $z_1 \ldots z_n$,
then a phase. It is the same layout the stabilizer engine uses internally, so the two halves of
the framework speak one format.

A Pauli string converted to its row:

```wl
QECPauliVector["XZZXI"]
```

<!-- => {1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0} -->

---

One deliberate difference from the engine: the phase lives in $\mathbb{Z}_4$, meaning an overall
factor $i^e$, where the engine's tableau carries a single sign bit. The reason is that the Pauli
group is not closed under multiplication with a sign alone.

The product of X and Z is not Y:

```wl
QECPauliString[QECPauliProduct["X", "Z"]]
```

<!-- => "-iY" -->

---

A sign bit cannot express that factor, so a package that tracks signs only must either drop the
phase or refuse to multiply. Carrying $\mathbb{Z}_4$ costs one integer per Pauli and makes the
algebra closed:

```wl
QECPauliString[QECPauliProduct["XIX", "IZZ"]]
```

<!-- => "-iXZY" -->

Hermitian Paulis — every stabilizer generator, every error — have even phase, and half of it is
exactly the engine's sign bit, so nothing is lost in translation.

---

Commutation is the symplectic product of two rows, which is what makes it cheap:

```wl
Outer[QECPauliCommuteQ, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}, {"XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"}]
```

<!-- => 4x4 array of True -->

## Logical operators and distance

The logical information lives in the Paulis that commute with every check but are not themselves
in the check group — the normalizer modulo the stabilizer, $N(S) \setminus S$. Finding a
canonical pair is done on the standard form of the check matrix, following the closed formulas
of Gottesman's thesis, section 4.1.

The logical X of the five-qubit code, which is the operator that thesis reports:

```wl
QECCode["5QubitCode"]["LogicalX"]
```

<!-- => {"ZIIZX"} -->

---

And its logical Z:

```wl
QECCode["5QubitCode"]["LogicalZ"]
```

<!-- => {"ZZZZZ"} -->

---

The distance is the minimum weight of a logical operator, and it is the number that says how
many errors the code can correct: $\lfloor (d-1)/2 \rfloor$ of them.

```wl
QECCode["SteaneCode"]["Distance"]
```

<!-- => 3 -->

---

Computing it exactly is NP-hard in general, so the search is a search. Two things keep it honest
on the codes that fit. The logical operators from the standard form are themselves in
$N(S) \setminus S$, so their minimum weight is an upper bound and the scan never runs past it;
and the syndromes of a whole weight class are computed as one matrix product rather than one call
per candidate. The witness that realises the distance is available, so the answer can be checked
rather than trusted:

```wl
QECCode["5QubitCode"]["MinimumWeightLogical"]
```

<!-- => "ZIIZX" -->

On this code the witness is the logical X itself, which is a coincidence of the five-qubit code
being as tight as it is: there is no lighter operator anywhere in $N(S) \setminus S$ than the
canonical pair the standard form produces.

---

That witness is a logical operator, and its weight is the distance:

```wl
{QECPauliWeight[QECCode["5QubitCode"]["MinimumWeightLogical"]],
 QECCode["5QubitCode"]["LogicalPauliQ", QECCode["5QubitCode"]["MinimumWeightLogical"]]}
```

<!-- => {3, True} -->

## Syndromes and decoding

An error either commutes or anticommutes with each check. The pattern of anticommutations is the
syndrome, and it is the only thing a measurement is allowed to reveal.

That last clause is the whole trick, and it is worth spelling out. Measuring a check operator on a
codeword is a projective measurement, and on a state already in the code space it returns $+1$
with certainty and disturbs nothing. On a state that an error has pushed out, it returns $-1$ —
but the outcome depends only on whether the error anticommutes with *that* check, and every
logical operator commutes with every check by construction. So the measurement cannot distinguish
$|0_L\rangle$ from $|1_L\rangle$, or any superposition of them, no matter how many checks are
measured. The syndrome carries information about the error and none about the logical state.
Reading it is free.

What the syndrome does not do is identify the error. Two errors differing by a stabilizer produce
the same syndrome and need the same correction — they are genuinely the same problem, and a code
where many errors share a coset is called degenerate. Two errors differing by a *logical* operator
also produce the same syndrome and need different corrections, and no measurement can tell them
apart. That ambiguity is the entire difficulty of decoding, and the distance is the measure of how
much room the code leaves before it bites.

The syndrome of a single X error on the Steane code:

```wl
QECCode["SteaneCode"]["Syndrome", "XIIIIII"]
```

<!-- => {0, 0, 0, 0, 1, 0} -->

---

A sign on the error cannot change which checks it anticommutes with, so the syndrome ignores it:

```wl
QECCode["SteaneCode"]["Syndrome", "-XIIIIII"] === QECCode["SteaneCode"]["Syndrome", "XIIIIII"]
```

<!-- => True -->

---

Decoding is the inverse problem: given the syndrome, infer the error. The simplest honest answer
is a table from syndrome to a minimum-weight error producing it, built out to the weight the code
actually guarantees:

```wl
QECPauliString /@ QECCode["BitFlipCode"]["Decoder"]
```

<!-- => <|{0, 0} -> "III", {1, 0} -> "XII", {1, 1} -> "IXI", {0, 1} -> "IIX"|> -->

---

Running the full cycle — syndrome, inference, correction — reports what happened rather than a
bare success flag:

```wl
Dataset[QECCode["SteaneCode"]["CorrectionCycle", "YIIIIII"]]
```

---

The distinction the report draws matters. A residual that lands in the stabilizer group is a
correction that worked; one that lands in $N(S) \setminus S$ is a **silent logical error**, the
failure mode the whole subject exists to bound; and a syndrome the table cannot decode is a third
thing again. On the bit-flip code, a Z error is invisible to every check:

```wl
QECCode["BitFlipCode"]["CorrectionCycle", "ZII"]["Outcome"]
```

<!-- => "LogicalError" -->

## Encoding circuits

A code is not useful until something prepares a codeword. The encoder follows Procedure 6.6 of
Gottesman's 2026 book, with the gate-to-row-operation dictionary of its Table 6.1.

The construction is worth understanding, because it is the one place where a linear-algebra fact
turns directly into a circuit. Write the $m$ generators as the columns of a $2n \times m$ binary
matrix, the X parts stacked above the Z parts. Left-multiplying by a Clifford gate is then a *row
operation* on that matrix, and the dictionary is small and fixed: a Hadamard on qubit $q$ swaps
the $q$-th X row with the $q$-th Z row; a phase gate adds the X row into the Z row; a CNOT adds
one X row to another and one Z row back the other way; a controlled-Z crosses X into Z between two
qubits. Relabelling which product of generators counts as a generator is a *column* operation,
and costs no gate at all.

So the problem becomes Gaussian elimination with two kinds of move, and the target is the matrix
of the trivial code — every generator a single Z. Reduce to it, record the gates used, and the
encoder is that sequence run backwards, since every gate in the dictionary is its own inverse
except the phase gate, whose inverse is three of itself.

One step remains after inverting. Reducing and re-inflating gets the *group* right but not
necessarily the *signs*: the prepared state is stabilized by each generator up to a minus. A Pauli
correction fixes it, and finding it is again linear algebra — solve over GF(2) for the Pauli whose
symplectic products with the generators match the pattern of wrong signs, then apply it.

The book's own recipe for that fixup is to prepend an X on each qubit whose generator came out
negative; the package instead solves the general system, which handles the completed generator set
in one step. Both are correct and produce circuits of the same length, on every code tested.

The encoding circuit of the five-qubit code:

```wl
QECCode["5QubitCode"]["EncodingCircuit"]
```

---

The reduction can only ever emit the five gates of that dictionary, so the circuit's vocabulary is
fixed in advance — here, four of them:

```wl
DeleteDuplicates[First /@ QECCode["SteaneCode"]["EncodingGates"]]
```

<!-- => {"H", "CZ", "SWAP", "CNOT"} -->

---

The sign fixup is the one thing that can add to that vocabulary, and on codes that need it a
Pauli appears alongside the Cliffords:

```wl
DeleteDuplicates[First /@ QECCode["5QubitCode"]["EncodingGates"]]
```

<!-- => {"H", "CZ", "CNOT", "X"} -->

---

Reducing and inverting prepares a state stabilized by each generator up to sign; a Pauli fixup
solved over GF(2) corrects the signs that came out wrong. Whether the result is really a codeword
is checked against the engine, by measuring every check on the prepared state:

```wl
AllTrue[{"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"},
    QECCode[#]["EncodingCircuitValidQ"] &]
```

<!-- => True -->

## Building codes from other codes

Three constructions of Gottesman's thesis, section 3.5, plus the CSS construction of chapter 5 of
the book, all return the same code object, so everything above works on the result unchanged.

A CSS code is built from classical parity-check matrices: rows of one become X-type checks, rows
of the other Z-type, and they commute exactly when the matrices are orthogonal over GF(2). Fed the
classical Hamming matrix, the construction reproduces the Steane code:

```wl
QECCode["CSS", QECClassicalHammingMatrix[3]]["Parameters"]
```

<!-- => {7, 1, 3} -->

---

"Reproduces" deserves care: the generators come out in a different order, so the claim is about
the stabilizer *group*, checked in both directions and with signs:

```wl
Module[{css = QECCode["CSS", QECClassicalHammingMatrix[3]], steane = QECCode["SteaneCode"]},
    AllTrue[steane["GeneratorVectors"], css["StabilizerMemberQ", #] &] &&
    AllTrue[css["GeneratorVectors"], steane["StabilizerMemberQ", #] &]
]
```

<!-- => True -->

---

Concatenation re-encodes each qubit of an outer code with an inner one. The phase-flip code inside
the bit-flip code is the Shor code, which is how Shor found it:

```wl
QECConcatenate[QECCode["PhaseFlipCode"], QECCode["BitFlipCode"]]["Parameters"]
```

<!-- => {9, 1, 3} -->

---

Distances multiply only as a lower bound, not an equality. The bit-flip code has quantum distance
1, yet concatenating it under the five-qubit code buys distance 5, not 3, because it still corrects
bit flips at weight 3:

```wl
QECConcatenate[QECCode["5QubitCode"], QECCode["BitFlipCode"]]["Parameters"]
```

<!-- => {15, 1, 5} -->

---

Qubit removal is surgery: recombine the generators so that exactly one ends in X and one in Z at
the removed qubit, drop those two — they become the new logical pair — and truncate the rest. On
the five-qubit code it produces the $[[4,2,2]]$ code, with the generators Gottesman reports:

```wl
QECRemoveQubit[QECCode["5QubitCode"]]["Generators"]
```

<!-- => {"XZZX", "YXXY"} -->

## Noise

None of the above says anything about whether a code *works*, because a code only means something
relative to a noise model. The simplest honest one is code-capacity noise: each qubit independently
suffers X, Y or Z with fixed probabilities, between rounds, with perfect check measurements. That
last assumption is a real one, and the sections after the next give it up.

Depolarizing noise, with the rate left symbolic:

```wl
QECNoiseModel["Depolarizing", p]["Probabilities"]
```

<!-- => {1 - p, p/3, p/3, p/3} -->

---

Because the errors are independent and identically distributed, the probability of a specific
error depends only on how many X, Y and Z it carries — never on where they sit. That collapses a
sum over $4^n$ errors into a sum over the far smaller set of weight profiles, and it is what makes
everything in the next section exact rather than sampled:

```wl
QECNoiseModel["Depolarizing", p]["ErrorProbability", "XIIIIII"]
```

<!-- => ((1 - p)^6 p)/3 -->

---

The model is not a private random-number generator: it maps onto the engine's own channels, so
the probabilities can be checked against what the engine actually does to a stabilizer state.
Note that the conventions differ — the engine's depolarizing parameter gives each Pauli a quarter
of it, this model spreads its rate over three — so the constructor converts, and this is the test
that the conversion is right:

```wl
Module[{noise = QECNoiseModel["Depolarizing", p]},
    Simplify[Sort[noise["QuantumChannel", {1}][PauliStabilizer[1]][[All, 1]]] == Sort[noise["Probabilities"]]]
]
```

<!-- => True -->

---

For codes too large to enumerate, the model samples. A bit-flip channel can only produce X errors,
and does:

```wl
BlockRandom[SeedRandom[7]; QECNoiseModel["BitFlip", 1/3]["RandomErrors", 6, 5]]
```

## The logical error rate

This is the question the subject is about. Prepare a codeword, hit it with an error drawn from the
noise, measure the checks, decode, correct, and ask whether the logical qubit got flipped.

One observation makes it exact. Label every Pauli error by two things: which checks it
anticommutes with (its syndrome) and which logical operators it anticommutes with (its class).
Both are the same question, so both are one matrix product over GF(2). Two errors sharing a label
differ by a stabilizer, so the label's fibre *is* the coset; decoding is guessing the class half
from the syndrome half; and a correction succeeds exactly when the guess was right. No residual
algebra per trial — the label already decides it.

The exact logical error rate of the five-qubit code under depolarizing noise:

```wl
Expand[QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", p]]]
```

<!-- => 10 p^2 - (200 p^3)/9 + (160 p^4)/9 - (128 p^5)/27 -->

---

The leading term is not a fit. The five-qubit code is perfect, so every weight-one error is
corrected and every weight-two error fails; there are $9\binom{5}{2} = 90$ of the latter, each of
probability $(p/3)^2(1-p)^3$, giving exactly $10p^2$. The polynomial is checkable by hand.

The Steane code, which is not perfect, so some weight-two errors survive:

```wl
Normal[Series[QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", p]], {p, 0, 3}]]
```

<!-- => (49 p^2)/3 - 56 p^3 -->

---

The strongest check available on any of this is a closed form that owes the package nothing. The
repetition code under bit-flip noise fails exactly when a majority of its qubits flips, which is a
binomial tail:

```wl
Simplify[
    QECLogicalErrorRate[QECCode["Repetition", 5], QECNoiseModel["BitFlip", p]] ==
        Sum[Binomial[5, j] p^j (1 - p)^(5 - j), {j, 3, 5}]
]
```

<!-- => True -->

---

With an exact rate per code, the threshold — the rate below which adding qubits helps and above
which it hurts — stops being a fit to sampled curves and becomes a crossing of polynomials. For
repetition codes under bit flips it sits at exactly one half:

```wl
Plot[
    Evaluate[Table[QECLogicalErrorRate[QECCode["Repetition", n], QECNoiseModel["BitFlip", p]], {n, {3, 5, 7}}]],
    {p, 0, 1},
    PlotLegends -> {"n = 3", "n = 5", "n = 7"},
    AxesLabel -> {"physical error rate", "logical error rate"},
    ImageSize -> 380
]
```

Below one half the longer code is better and the gap widens with n; above it the order reverses.
That crossing is the whole argument for error correction, and here it is an exact statement about
three polynomials rather than a claim about three sampled curves.

---

For codes past the enumeration limit the same quantity is estimated by drawing errors, and the two
routes are independent implementations that check each other:

```wl
BlockRandom[SeedRandom[3];
    {N[QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", 1/20]]],
     QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Depolarizing", 1/20], 100000]}
]
```

---

The decoder is maximum likelihood: for each syndrome it picks the coset carrying the most
probability, not the lightest representative. The `"Decoder"` option chooses, and the difference
is not academic. Against the lookup table — the decoder anyone can afford on a large code — coset
weighing wins in two separate ways. On the degenerate nine-qubit Shor code, under ordinary
symmetric noise, because many light errors share a coset and the table sees only one of them:
$13p^2$ against $31p^2$. And on the Steane code under biased noise, as real hardware is, with Z
ten times likelier than X: $651p^2/25$ against $756p^2/25$. On a perfect code like the five-qubit
one there is nothing to gain and the two agree exactly.

The biased case:

```wl
Coefficient[
    Normal[Series[
        QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel[<|"X" -> p/10, "Y" -> p/10, "Z" -> p|>]],
        {p, 0, 2}]],
    p, 2]
```

<!-- => 651/25 -->

## Noisy syndrome extraction

Everything so far has assumed the checks are read out perfectly. That assumption is the one a
hardware person will not grant, and giving it up is what the rest of this note is about.

There are three levels, and they are the same object with different keys set. At **code
capacity** the data is hit between rounds and the checks are exact — that is everything above. At
the **phenomenological** level the checks additionally lie with some probability, which is the
first point at which repeating the measurement means anything. At **circuit level** there are no
"checks" at all: the extraction is written out as a subcircuit and every gate, reset and readout
in it can fail; faults attach to circuit *locations* with one rate per location type, which is a
Pauli specialisation of the basic model of Gottesman's §10.1.1.

Each level drops an assumption the one before it made, but they are not three estimates of one
number, and the middle one has a name in the literature that is worth carrying: Gottesman's §10.4
introduces the phenomenological model as the standard example of an **optimistic** assumption,
not of a more realistic one.

The level is derived from what you set, not declared:

```wl
{QECNoiseModel["Depolarizing", p]["Level"],
 QECNoiseModel["Depolarizing", p, "MeasurementError" -> q]["Level"],
 QECNoiseModel["Circuit", p]["Level"]}
```

<!-- => {"CodeCapacity", "Phenomenological", "Circuit"} -->

---

A circuit-level model puts a rate on each kind of location. One number spreads over all of them;
an Association sets them separately, and a missing key is zero, which is a clean way to ask what
one kind of failure does on its own:

```wl
QECNoiseModel["Circuit", <|"TwoQubit" -> p, "Measurement" -> q|>]["Rates"]
```

<!-- => <|"OneQubit" -> 0, "TwoQubit" -> p, "Measurement" -> q, "Reset" -> 0, "Idle" -> 0|> -->

Idling is modelled but stays at zero by default, and the next section says why that default is
optimistic rather than neutral.

## The extraction circuit

To measure a Pauli generator with an ancilla: rotate each data qubit in the generator's support
so that its letter there becomes $Z$, CNOT those data qubits into the ancilla, rotate back, and
measure the ancilla. The ancilla ends up holding the parity of the rotated qubits, which is the
generator's eigenvalue.

This is the circuit of Gottesman's §12.1.1, figure 12.1a, and that section's title is worth
quoting: *Non-Fault-Tolerant Measurement of Paulis*. The book's version puts the ancilla in
$|+\rangle$ and applies controlled-$P$ from ancilla to data; conjugating each data qubit by the
rotation that sends its letter to $Z$ turns that into a controlled-$Z$ ladder, and since
$H \cdot CZ \cdot H = CNOT$ the ancilla's Hadamards are absorbed into the ladder. Gate for gate it
is the same circuit with the basis change written out rather than folded into "controlled-$P$".

The rotation is $H$ for an $X$, the engine's `"V"` — the square root of $X$, which sends $Y$ to
$Z$ — for a $Y$, and nothing at all for a $Z$. A single generator carrying all three letters
shows every case at once:

```wl
QECCode[{"XYZ"}]["SyndromeCircuit"]["Instructions"]
```

<!-- => {{"R", 4}, {"H", 1}, {"V", 2}, {"CNOT", 1, 4}, {"CNOT", 2, 4}, {"CNOT", 3, 4}, {"H", 1}, {"Vdg", 2}, {"M", 4}} -->

---

The bit-flip code's generators are pure $Z$, so its circuit is nothing but CNOTs and
measurements — one ancilla per generator, reset before use:

```wl
QECCode["BitFlipCode"]["SyndromeCircuit"]["Instructions"]
```

<!-- => {{"R", 4}, {"CNOT", 1, 4}, {"CNOT", 2, 4}, {"M", 4}, {"R", 5}, {"CNOT", 2, 5}, {"CNOT", 3, 5}, {"M", 5}} -->

---

The ancillas are reset and reused between rounds, so a long experiment does not need more of
them. Five rounds of the Steane code run on the same thirteen qubits as one:

```wl
{QECCode["SteaneCode"]["SyndromeCircuit", 1]["Qubits"],
 QECCode["SteaneCode"]["SyndromeCircuit", 5]["Qubits"]}
```

<!-- => {13, 13} -->

---

What the whole thing costs, which is also the number of places a fault can happen:

```wl
QECCode["SteaneCode"]["SyndromeCircuit", 3]["GateCounts"]
```

<!-- => <|"R" -> 18, "H" -> 72, "CNOT" -> 72, "M" -> 18|> -->

---

Two limits are worth stating, and the first is not a rough edge. The ancillas are bare rather
than cat-state-verified, and that is the difference between fault tolerant and not, not a detail
alongside gate scheduling: §12.1.1 says of this exact circuit that "a single faulty gate can
cause multiple data errors, even if the fault doesn't directly affect the data block". The fix is
a verified ancilla feeding Shor (§12.2), Steane (§12.3) or Knill (§12.4) error correction. Using
the bare version anyway is a sourced choice rather than a shortcut — §12.5.1: "Frequently people
using surface codes don't even bother with Shor EC, instead using the non-FT measurement
technique of section 12.1.1 instead" — but that concession belongs to the surface-code setting
and is not a general licence.

The second: generators are extracted one after another rather than interleaved into parallel
layers, so the circuit is deeper than a real one, and the `"Idle"` rate is zero. Do not read that
as idling being undefined here. §15.5.2 resolves an ill-defined time step by *defining* one, the
longest gate, and charging the padding as storage error; and §15.5.1 shows that less parallelism
means *more* waiting, multiplying the storage rate, with the flat statement that "to have a
threshold, it is essential to do parallel gates, at least when the storage error rate $p_S$ is
non-zero". Zero idle noise is an optimistic simplification, and it is optimistic in exactly the
direction this circuit is already weakest.

Nothing in this circuit is taken on trust. Pushing a Pauli through the emitted instructions and
reading the ancilla measurements has to reproduce `code["Syndrome", …]`, which is a matrix
product that knows nothing about circuits; the same circuit is also run as an actual stabilizer
state on the engine and its checks read by expectation. Both agree, for every weight-one and
weight-two error on all five named codes.

## Detectors

Once a check can lie, a syndrome bit on its own is worthless: it can be wrong because the data
was hit or because the readout misreported, and one round cannot tell those apart. Repeating the
extraction is the old fix (§12.2.2). Decoding the *differences* rather than the syndromes is a
newer one, from the topological-code literature — Dennis, Kitaev, Landahl and Preskill,
*Topological quantum memory* ([arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143)),
where repeated noisy measurement turns the history into a lattice in one more dimension:

$$D_1 = S_1, \qquad D_r = S_r \oplus S_{r-1}, \qquad D_{R+1} = S_{\mathrm{final}} \oplus S_R$$

A detector is a parity that is deterministically zero when nothing goes wrong, so a fired
detector is unambiguous evidence of a fault. The last one closes the experiment against a
noiseless final readout, where $S_{\mathrm{final}}$ is the syndrome of whatever Pauli is left on
the data.

Two things follow from "deterministically zero" that are easy to skip past. It assumes the data
enters *inside* the code space — $D_1 = S_1$ is a detector only because the experiment prepares a
codeword first — so what is modelled here is a memory experiment, not an error-correction gadget
handed a block that already carries errors. And this is not the book's rule. Gottesman's answer to
an untrustworthy syndrome is repeat and *agree*: trust a run of $t+1$ consecutive identical
syndromes, since at least one of them was taken with no faults, discard the rest, and look that
single syndrome up (§12.2.2, §12.2.3). He raises the differencing strategy and sets it aside —
"it is difficult to analyze such strategies in the completely general case" — so what follows
computes exactly the thing the standard reference declines to analyse in general, and inherits
none of its fault-tolerance guarantees.

`QECDetectorModel` is the map from each individual fault to what it does: which detectors it
fires, which logical operators it flips, and with what probability. On the bit-flip code with one
round it is small enough to read in full — three data faults and two readout faults, four
detectors:

```wl
Module[{dem, bits},
    dem = QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p, "MeasurementError" -> q], 1];
    bits[v_] := StringJoin[ToString /@ v];
    Grid[
        Prepend[
            MapThread[
                {StringRiffle[ToString /@ #1, " "], bits[#2], bits[#3], #4} &,
                {dem["Locations"], dem["DetectorMatrix"], dem["ObservableMatrix"], dem["Probabilities"]}
            ],
            {"fault", "detectors", "observables", "probability"}
        ],
        Frame -> All,
        FrameStyle -> GrayLevel[0.85],
        Alignment -> Left,
        Background -> {None, {GrayLevel[0.94]}},
        ItemSize -> {{Automatic, 8, 8, Automatic}, Automatic}
    ]
]
```

---

That table is the point of the construction, and it is worth reading twice. A data error fires
its checks in round one and is *still there* at the final readout, so the difference detectors
stay quiet: one fired detector. A readout error fires its check in its own round and again in the
difference against the next one: two fired detectors, and nothing left on the data. The two
failures no longer look alike.

```wl
Module[{dem = QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p, "MeasurementError" -> q], 1]},
    {dem["DetectorMatrix"][[1]], dem["DetectorMatrix"][[4]]}
]
```

<!-- => {{1, 0, 0, 0}, {1, 0, 1, 0}} -->

---

The detector half of a data fault's row is exactly the syndrome computed the old way, which is
the cheapest check there is that the circuit measures what it claims to:

```wl
Module[{code = QECCode["SteaneCode"], dem},
    dem = QECDetectorModel[code, QECNoiseModel["BitFlip", p, "MeasurementError" -> q], 1];
    dem["DetectorMatrix"][[1, 1 ;; 6]] == code["Syndrome", "XIIIIII"]
]
```

<!-- => True -->

---

Some faults fire nothing at all and still move the logical qubit. Their existence is not a bug:
they are the circuit-level analogue of a logical operator. Under depolarizing noise the bit-flip
code has exactly three of them — the single $Z$ errors — which is the circuit-level way of saying
it has no phase protection:

```wl
QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["Depolarizing", p, "MeasurementError" -> q], 1]["UndetectableFaults"]
```

<!-- => {{"Data", 1, 1}, {"Data", 1, 2}, {"Data", 1, 3}} -->

---

Every detector's firing rate is available in closed form, which is the natural quantity to check
an external simulator against — it is basis-free and every detector is a separate test. With the
data noise switched off, each detector fires exactly when its own readout lies:

```wl
Simplify @ QECDetectorModel[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", 0, "MeasurementError" -> q], 1]["DetectorRates"]
```

<!-- => {q, q, q, q} -->

## The memory experiment

Prepare a codeword, run $r$ rounds of noisy extraction, read the data out perfectly at the end,
decode the whole detector history, and ask whether the logical qubit survived.

`"Rounds"` says how many and defaults to the code distance. That is the topological-memory
convention: fewer rounds than the distance protects the time direction less than the space ones,
which quietly caps the whole thing. It is deliberately not the book's criterion, and the two
should not be conflated — §12.2.2 asks instead for a run of $t+1$ *agreeing* syndrome
measurements, proves $(t+1)^2$ repetitions enough to find one, and explicitly rejects $2t+1$
repetitions with a majority vote, "since there are $2^{n-k}$ possible values of the syndrome …
There might not be a majority". That criterion buys a fault-tolerance guarantee; this default
buys a memory experiment whose time and space directions are equally protected.

The identity that ties all of this to the earlier sections: a phenomenological model with perfect
readout, run for one round, *is* code-capacity noise. Its rate had better be the same polynomial,
and it is:

```wl
Simplify[
    QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p, "MeasurementError" -> 0], "Rounds" -> 1] ==
        QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p]]
]
```

<!-- => True -->

That one line exercises the circuit emitter, the frame propagation, the detector definition, the
fault enumeration and the exact sum at once, against an answer already checked against a binomial
tail.

---

The rate is still exact, and now in two variables. Setting either to zero recovers a limit worth
knowing: $q = 0$ gives back $3p^2 - 2p^3$, and $p = 0$ gives $q^2$ — with no data errors at all, a
single readout lie is recognised by its two-detector signature and corrected, so failure needs
two of them:

```wl
Expand @ QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["BitFlip", p, "MeasurementError" -> q], "Rounds" -> 1]
```

<!-- => 3 p^2 - 2 p^3 + 6 p q - 12 p^2 q + 6 p^3 q + q^2 - 6 p q^2 + 9 p^2 q^2 - 4 p^3 q^2 -->

---

This is the part no sampling simulator can do, and it is worth saying how it is possible.
Enumerating fault configurations would mean $2^{\text{faults}}$ of them, which is hopeless. Instead
the *distribution over effects* is carried — a vector indexed by the possible detector-and-observable
patterns — and the locations are folded in one at a time,

$$P'[s] = P[s]\left(1 - \sum_i p_i\right) + \sum_i p_i \, P[s \oplus \sigma_i]$$

the inner sum running over the mutually exclusive outcomes of a single location, so a three-way
depolarizing channel is handled exactly rather than as three independent coins. The cost is
locations times state space rather than $2^{\text{faults}}$, and the state space —
$2^{(r+1)m + 2k}$ — is what bounds it, through `$QECExactDetectorLimit`.

Past that bound the same quantity is sampled, and the two routes check each other. The Steane
code at circuit level, exactly and then over two hundred thousand shots:

```wl
BlockRandom[SeedRandom[11];
    {N @ QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Circuit", 1/1000], "Rounds" -> 1],
     QECLogicalErrorRate[QECCode["SteaneCode"], QECNoiseModel["Circuit", 1/1000], 200000, "Rounds" -> 1]}
]
```

---

Because all three levels give polynomials, the price of honesty about noise is a picture rather
than an argument. The same code, the same rate $p$, three noise models:

```wl
Module[{bf = QECCode["BitFlipCode"], capacity, phenomenological, circuit},
    capacity = QECLogicalErrorRate[bf, QECNoiseModel["BitFlip", p]];
    phenomenological = QECLogicalErrorRate[bf, QECNoiseModel["BitFlip", p, "MeasurementError" -> p], "Rounds" -> 1];
    circuit = QECLogicalErrorRate[bf, QECNoiseModel["Circuit", p], "Rounds" -> 1];
    Plot[
        Evaluate[{capacity, phenomenological, circuit}], {p, 0, 1/5},
        PlotLegends -> {"code capacity", "phenomenological", "circuit"},
        AxesLabel -> {"physical error rate", "logical error rate"},
        ImageSize -> 380
    ]
]
```

Code capacity and the phenomenological level are both quadratic in $p$ and differ by a constant
factor — $3p^2$ against $10p^2$. The circuit-level curve is a different shape entirely, and the
next section is about why.

One caveat on that figure, and it is the book's: "A threshold in a phenomenological model should
not be confused with the threshold derived from a full circuit model; they are not at all
comparable" (§10.4). The three curves share an axis but not a quantity — the phenomenological $p$
is a rate per qubit per round, the circuit-level $p$ a rate per location, and the
phenomenological model ignores how many gates a syndrome bit costs. Read the picture as one code
losing ground as assumptions are dropped, not as three thresholds being compared.

## What circuit-level noise actually costs

Two facts fall out of these polynomials that are hard to state without them. Both are about the
extraction circuit rather than the code, and both are exact.

**Duality is a property of the codes, not of the circuits that measure them.** The bit-flip and
phase-flip codes are the same code with $X$ and $Z$ exchanged; both are $[[3,1,1]]$, and under
depolarizing noise at code capacity their rates are identical. At circuit level they are not:

```wl
Module[{lead = Normal[Series[#, {p, 0, 1}]] &},
    {lead @ QECLogicalErrorRate[QECCode["BitFlipCode"], QECNoiseModel["Circuit", p], "Rounds" -> 1],
     lead @ QECLogicalErrorRate[QECCode["PhaseFlipCode"], QECNoiseModel["Circuit", p], "Rounds" -> 1]}
]
```

<!-- => {32 p/15, 112 p/15} -->

---

The asymmetry is entirely the rotations. The phase-flip code's generators are $X$-type, so its
extraction needs eight $H$ gates that the bit-flip code's does not, and each is a fresh
depolarizing location sitting on a data qubit. Switch the one-qubit-gate noise off and the two
codes agree again, to the coefficient:

```wl
Module[{noise = QECNoiseModel["Circuit", <|"TwoQubit" -> p, "Measurement" -> p, "Reset" -> p|>],
        lead = Normal[Series[#, {p, 0, 1}]] &},
    {lead @ QECLogicalErrorRate[QECCode["BitFlipCode"], noise, "Rounds" -> 1],
     lead @ QECLogicalErrorRate[QECCode["PhaseFlipCode"], noise, "Rounds" -> 1]}
]
```

<!-- => {32 p/15, 32 p/15} -->

---

**This extraction circuit stops a distance-three code correcting single faults.** The five-qubit
code has distance 3, so at code capacity and at the phenomenological level its logical error rate
starts at $p^2$: any one error is corrected. With this circuit it starts at $p$:

```wl
Normal @ Series[
    QECLogicalErrorRate[QECCode["5QubitCode"], QECNoiseModel["Circuit", p], "Rounds" -> 1],
    {p, 0, 1}
]
```

<!-- => 29 p/5 -->

---

Note the phrasing: this is a defect of the *gadget*, not a property of circuit-level noise.
Gottesman's §10.2 exists precisely so that circuit noise need not cost the exponent — a gadget
satisfying its propagation properties keeps a distance-three code at $p^2$ under circuit noise.
Order $p$ is the signature of a construction that does not satisfy them.

The mechanism is the *hook error*, and it is worth being exact about, because the direction of the
CNOTs here is the mirror of figure 12.1b's and the obvious transcription is wrong. With the
ancilla as target, an $X$ on the ancilla reaches no data qubit at all — CNOT carries $X$ only from
control to target, so it sits there and flips one readout. The fault that spreads is a $Z$ (hence
also a $Y$): CNOT carries $Z$ from target back to control, the ancilla keeps its $Z$, and it lands
on every data qubit the ladder has **still to touch**. After the un-rotation the residual is
exactly the un-extracted tail of the generator — a mid-ladder $Z$ on the five-qubit code's
`XZZXI` leaves `IZZXI`, `IIZXI` or `IIIXI` — and, being a sub-product of the generator, it does
not fire that generator's own check in the round it happened. A $Z$ before the first CNOT leaves
the whole generator, which is a stabilizer and harmless; after the last, nothing. Only mid-ladder
hurts.

Count what that does to the decoder's job. At the phenomenological level every single fault has
its own detector signature, so no single fault is ever ambiguous. With this circuit the signatures
collide, and some of the colliding faults do different things to the logical qubit — at which
point no decoder reading detectors alone can tell them apart.

```wl
Module[{ambiguity},
    ambiguity[code_, noise_] := Module[{dem, d, o, grouped},
        dem = QECDetectorModel[code, noise, 1];
        d = dem["DetectorMatrix"]; o = dem["ObservableMatrix"];
        grouped = GroupBy[Range[Length[d]], d[[#]] &, DeleteDuplicates[o[[#]] & /@ #] &];
        <|"faults" -> Length[d], "distinct signatures" -> Length[grouped],
          "ambiguous signatures" -> Count[grouped, alternatives_ /; Length[alternatives] > 1]|>
    ];
    Dataset @ <|
        "phenomenological" -> ambiguity[QECCode["5QubitCode"], QECNoiseModel["Depolarizing", 1/1000, "MeasurementError" -> 1/1000]],
        "circuit" -> ambiguity[QECCode["5QubitCode"], QECNoiseModel["Circuit", 1/1000]]
    |>
]
```

Nineteen faults, nineteen distinct signatures, none ambiguous — against 288 faults collapsing
onto 71 signatures of which 29 carry conflicting logical effects. The code did not get worse; the
circuit measuring it is not fault tolerant.

This is the gap that verified ancillas exist to close, and the book spends a chapter on it: cat
states (§12.1.2–12.1.3) feeding Shor error correction (§12.2), Steane error correction (§12.3),
Knill error correction (§12.4), with §12.5 comparing their costs. One detail there is a useful
target: §12.3.3 argues Steane error correction needs no repeated syndrome measurement at all,
because its ancilla and measurement errors enter additively. None of these is implemented yet —
they are roadmap work — and the number above is what it costs not to have them.

## Handing it to Stim

Stim ([arXiv:2103.02202](https://arxiv.org/abs/2103.02202)) is the community's reference
stabilizer simulator, and the source of the detector error model as an explicit object. The point
of connecting to it rather than competing with it is that the two are good at different things: it
samples fast and large, this layer derives exactly and symbolically. Writing the memory experiment out in its syntax buys
an independent check on all of the above, and a path to the decoders built on top of it.

The emitted circuit encodes first — Stim starts in $|0\ldots0\rangle$, which is not a codeword for
a general code, and a memory experiment that does not start inside the code space has a
meaningless first round of detectors — then runs the noisy rounds, then a noiseless one:

```wl
Column @ Take[
    StringSplit[QECStimCircuit[QECCode["BitFlipCode"], QECNoiseModel["Circuit", 1/200], 2], "\n"],
    14
]
```

---

The structure scales as it should: one detector per check per round, plus the final noiseless
round, and a single logical observable. Three noisy rounds of the Steane code's six checks give
four layers of six:

```wl
Module[{source = QECStimCircuit[QECCode["SteaneCode"], QECNoiseModel["Circuit", 1/500], 3]},
    {StringCount[source, "DETECTOR"], StringCount[source, "OBSERVABLE_INCLUDE"]}
]
```

<!-- => {24, 1} -->

---

One asymmetry is deliberate. `QECLogicalErrorRate` asks whether the decoder got the residual
error's whole class right, covering logical $X$ and logical $Z$ damage at once, because that is
the right question for an unknown logical state. A Stim memory experiment is single-basis by
construction — prepare $|0_L\rangle$, keep it, measure $\bar{Z}$ — so the emitted circuit declares
one observable and `"Observable"` chooses which. The detector part is basis-independent, and that
is where the two models are compared.

That comparison has been run. Five cases — three at circuit level, two at phenomenological,
including the non-CSS five-qubit code whose $Y$ generators are the only thing exercising the
square-root-of-$X$ rotation — with 58 detector firing rates and 5 observable flip rates computed
exactly here and sampled over two million shots there. Everything agrees within 2.7 standard
errors, and PyMatching decodes every exported circuit unchanged, which is a second and structural
check: a malformed circuit does not produce a detector error model that decomposes into a matching
graph.

## Where this comes from

Every claim above is either computed by the package and checked against an independent oracle, or
taken from one of these. Section numbers are Gottesman's 2026 *Stabilizer Codes and Quantum Error
Correction* unless the thesis is named.

| Topic | Source |
|---|---|
| Stabilizer formalism, standard form, logical operators, distance | Thesis §3.2, §4.1 |
| Concatenation, qubit removal, pasting | Thesis §3.5 |
| The $[[n, n-2, 2]]$ family, the quantum Singleton bound | Thesis §8.1 |
| Repetition codes as neighbouring-pair checks | Thesis §2.2 |
| The encoder | Book §6.4.1, Procedure 6.6 with Table 6.1 |
| Circuit-level noise: locations, faults, one rate per location type | Book §10.1.1, Definitions 10.1–10.3 |
| What fault tolerance formally requires | Book §10.2 |
| The phenomenological model, and thresholds not being comparable across levels | Book §10.4 |
| The extraction circuit, and its non-fault-tolerance | Book §12.1.1, figures 12.1a and 12.1b |
| Measuring the whole syndrome and then repeating | Book §12.2.2, figures 12.6 and 12.7 |
| How many repetitions a fault-tolerance guarantee needs | Book §12.2.2 |
| Verified ancillas: cat states, Shor, Steane, Knill | Book §12.1.2–12.1.3, §12.2, §12.3, §12.4, §12.5 |
| Why surface-code practice uses the non-fault-tolerant circuit anyway | Book §12.5.1 |
| Reusing ancillas across rounds | Book §15.4, §15.4.3 |
| Serial extraction, waiting, and why idle noise matters more not less | Book §15.5.1, §15.5.2 |
| Detectors as differences of repeated noisy syndromes | Dennis, Kitaev, Landahl, Preskill, [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143) |
| The detector error model as an object, and the Stim syntax | Gidney, [arXiv:2103.02202](https://arxiv.org/abs/2103.02202) |
| Minimum-weight matching on a detector model | Higgott and Gidney, PyMatching |

Two places where this note is deliberately outside the sources, and says so rather than implying
cover it does not have. The exact circuit-level rate is a closed-form polynomial obtained by
folding fault locations into a distribution over effects; §10.4 lists rigorous proof, simulation,
analytical estimate and hybrid methods for getting at a threshold, and this is none of them. And
the decoder optimises over the whole spacetime fault history, which is the strategy §12.2.2 raises
and then sets aside as hard to analyse in general — so the numbers here are exact error rates for
a stated protocol, not fault-tolerance proofs.

## What is exact, and what is not

Five limits are worth stating plainly, because a result that looks exact and is not is worse than
an estimate.

**Enumeration is capped.** The code-capacity route walks all $4^n$ errors, so it is bounded by
`$QECExactEnumerationLimit`, which is $4^{10}$ by default. The nine-qubit Shor code takes about a
second; past ten qubits, pass a sample count instead and get an estimate with error bars you
control.

**The exact circuit-level route is capped by its state space, not its fault count.** Folding
locations into a distribution over effects costs $2^{(r+1)m + 2k}$ per step, bounded by
`$QECExactDetectorLimit`. That grows in the *rounds*, so the Steane code is exact at one round and
must be sampled at three. Sampling draws the true channels, not the independent-mechanism
approximation attached to the detector model's rows, so the sampled number is not the weaker
answer — only the slower-converging one.

**The decoder is chosen at a probe rate.** Maximum likelihood needs to compare probabilities, and
with a symbolic rate there is nothing to compare. The argmax is therefore resolved at a small
fixed rate, `$QECDecoderProbeRate`. Once chosen the decoder is a fixed function, and the
polynomial reported for it is exact at every p — but it is the polynomial for *that* decoder, not
for a decoder retuned at each p.

**The circuit-level decoder is minimum weight over fault subsets, up to a reach.** It is exact at
low weight and explodes combinatorially past it, so it stops and says so at
`$QECDecoderSubsetLimit`. A matching decoder plugs in at the detector model and returns the same
kind of answer; that is roadmap work, not a gap in the model.

**The extraction circuit is not fault tolerant, and is not claimed to be.** The ancillas are bare
rather than verified, which is what produces the hook errors measured two sections above. That is
an honest property of the circuit rather than an approximation: the rate reported is the true rate
*for this circuit*, and §12.1.1 is where to read what it costs.

**Zero idle noise is the default, and it is the one optimistic assumption left.** Idling is now
modelled rather than missing: the schedule is recovered from the instruction list by ASAP list
scheduling, a qubit no instruction of a time step touches is waiting, and each such slot carries
its three Paulis. But `"Idle"` still *defaults* to zero, and every circuit-level number in this
note is computed at that default. Generators are extracted sequentially, so the circuit waits more
than a parallel one would, not less; §15.5.1 quantifies the storage rate that serialisation
multiplies, and §15.5.2 says an ill-defined time step is to be defined and its padding charged,
never dropped. Set `"Idle"` and the numbers get worse, in the direction this circuit is already
weakest.

## How this is checked

Nothing above is trusted because the package computed it. Every fast routine is cross-checked
against a slow, obvious one, and the checks run as a suite of 473 tests.

| What | Checked against |
|---|---|
| Pauli algebra, products and phases | Dense $2^n \times 2^n$ matrices built from scratch |
| Row layout | The engine's own tableau, row for row |
| Syndromes | Expectation values measured on an actual stabilizer state |
| Encoding circuits | The state they prepare, checked check by check |
| Structural results | Gottesman's published values — $\bar{X} = ZIIZX$, the $[[4,2,2]]$ generators |
| CSS and concatenation | The textbook codes they must reproduce, as groups and with signs |
| Noise probabilities | The mixture the engine's channels actually return |
| Logical error rate | A closed-form binomial tail, and an independent sampled route |
| The extraction circuit | `code["Syndrome", …]`, and the circuit run as a state on the engine |
| Detectors and faults | Stim's own sampler, detector by detector, over two million shots |
| Circuit-level rate | The code-capacity polynomial it must reduce to, and a sampled route |
| Idle noise and the schedule | The Stim export, `DEPOLARIZE1` line for line against the model's slots |
| The exported circuit | PyMatching, which decodes it unchanged |

The habit is worth keeping when extending the layer: a result checked only against the code that
produced it is not checked.
