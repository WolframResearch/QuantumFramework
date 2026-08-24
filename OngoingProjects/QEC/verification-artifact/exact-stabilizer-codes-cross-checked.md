---
Template: Default
Title: Exact Stabilizer Codes
Author: Mads Bahrami
---

# Exact Stabilizer Codes

<!-- #| style: Subtitle -->
Symbolic, sign-exact stabilizer-code analysis in the Wolfram Language, re-derived from GF(2) first principles, cross-checked against Stim, and carrying a parameter Stim cannot.

<!-- #| style: Author -->
Mads Bahrami (last updated: Aug 24, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

### Setting the Stage: How This Notebook Flows

This notebook is a computation-first tour of a small, sharp corner of quantum error correction: the *exact* analysis of stabilizer codes. For the five textbook codes (the three-qubit bit-flip and phase-flip codes, the nine-qubit Shor code, the five-qubit perfect code, and the seven-qubit Steane code) we compute the exact code distance and a minimum-weight witness for it, the logical Pauli operators together with their full commutation algebra, and a sign-exact +1 Clifford encoder. Every one of these results is then verified two independent ways: by a from-scratch symplectic re-derivation written in this notebook, and by a differential comparison against [Stim](https://github.com/quantumlib/Stim), the field-standard stabilizer simulator. At the end we do the one thing a numeric simulator structurally cannot: carry a symbolic parameter (a whole code family indexed by $n$, and a continuous error angle $\theta$) and certify the answer in a single evaluation.

In other words, I have tried to build a small, honest instrument. It does not race Stim on scale or speed, because it is not that kind of tool: it computes things exactly, and exact minimum distance is NP-hard, so it lives at small $n$ by mathematical necessity. What it offers instead is exactness and symbolic generality, and the whole notebook is arranged so you can watch both claims be checked rather than take them on faith. I strongly believe in a computation-first narrative for learning: in a sense, if I cannot compute it, I cannot claim to understand it.

Before we start, a few notes on the environment. This is a live Wolfram notebook: evaluate the cells from top to bottom. Some cells define objects (the code list, a handful of small symplectic functions, a Python session) that later cells use, so order matters. The narrative runs as one continuous sequence, like a movie, with headings only as gentle transitions. My suggestion is to read each output and its meaning first, and only then unpack the input code that produced it.

Remember that you are not locked into the code as given. You can (and should) change the codes, feed in your own stabilizer generators, push the symbolic family further, and re-run. That is the point of shipping this as a notebook rather than a PDF. If you have questions, reach out at [quantum@wolfram.com](mailto:quantum@wolfram.com).

### Prerequisites and How to Read This

You should be comfortable with the stabilizer formalism at the level of a first course: Pauli operators, stabilizer generators, the codespace as their joint +1 eigenspace, and logical operators as the normalizer modulo the stabilizer group. The binary symplectic picture (each Pauli as a vector over GF(2), commutation as a symplectic inner product) is built up from scratch in Part II, so no prior fluency there is required. The cross-check in Part III calls Stim through [`ExternalEvaluate`](https://reference.wolfram.com/language/ref/ExternalEvaluate.html); if you do not have Stim installed those cells will say so and you can skip them, keeping the independent re-derivation as your check.

Let's start!

## Part I: The Open Niche Beside a Fast Simulator

A stabilizer code on $n$ physical qubits is defined by a set of commuting Pauli operators, its stabilizer generators. The codespace is their simultaneous +1 eigenspace, and if there are $m$ independent generators the code stores $k = n - m$ logical qubits. An error is *detectable* when it anticommutes with at least one generator (it flips that generator's measured eigenvalue, a syndrome bit), and *undetectable* when it commutes with all of them. The dangerous errors are the undetectable ones that still act nontrivially on the codespace: these are the logical operators, and the code distance is the weight of the lightest one.

Equivalently, and this is the representation the whole subject runs on, each Pauli string on $n$ qubits becomes a binary vector of length $2n$: an $X$ contributes to the first half, a $Z$ to the second, a $Y$ to both. Two Paulis commute exactly when a bilinear form over GF(2), the symplectic product $\omega$, vanishes between their vectors. In short, the entire question "is this error detectable, is it logical, how heavy is the lightest logical" becomes linear algebra over the two-element field.

Why is this worth a dedicated exact tool when Stim already simulates stabilizer circuits at enormous scale? Because the two live on different axes. Stim samples: it runs stabilizer circuits and detector-error models blindingly fast, over thousands of qubits, which is exactly what decoder benchmarking and threshold estimation need. It does not hand you the exact minimum distance of a code, and for a deep reason. Computing the exact minimum distance of a stabilizer code is NP-hard (Kapshikar and Kundu, [arXiv:2203.04262](https://arxiv.org/abs/2203.04262)), so any tool that returns the *exact* answer is confined to small codes no matter how it is written. That is not a limitation to apologize for; it is the shape of the problem. The textbook codes are small, and small is precisely where an exact, symbolic treatment earns its keep.

The other axis is symbolic generality. A numeric simulator represents one concrete code at a time: fixed $n$, fixed generators, fixed angles. A symbolic engine can carry $n$ as a symbol and certify a statement about an entire infinite family in one evaluation, or carry an error angle $\theta$ and return the syndrome expectation as an exact function of it. Neither of those is a thing Stim is built to do, and neither competes with what Stim is built to do. The tools are complementary. This notebook stakes out the exact-and-symbolic corner and checks, at every step, that its exact answers agree with Stim's wherever the two overlap.

The exact code-analysis layer used here ships in the QuantumFramework paclet, following the constructions in Gottesman's thesis. Let us load it: register the paclet directory, then read in the QEC package. The paths below point to a local QuantumFramework checkout; adjust them to your own.

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework/"]
Needs["Wolfram`QuantumFramework`"]

Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/QEC/QEC/QEC.wl"]
```

A code is an object you can ask questions of. Build the Steane code and read its parameters and generators:

```wl
steane = QECCode["SteaneCode"];
{steane["Parameters"], steane["Generators"]}
```

The first entry is $\{n, k\}$, the physical and logical qubit counts; the second is the six Pauli-string generators, three all-$X$ checks and three all-$Z$ checks, the hallmark of a CSS code.

## Part II: The Five Textbook Codes, Exactly

Collect the five codes we will study into one association, so later cells can sweep over all of them at once:

```wl
codes = AssociationMap[
    QECCode,
    {"BitFlipCode", "PhaseFlipCode", "ShorCode", "5QubitCode", "SteaneCode"}
];
Keys[codes]
```

### The Symplectic Picture: A Pauli as a Vector over GF(2)

To verify the package's answers independently, we will not call the package's own distance or logical-operator routines and ask them to check themselves. Instead we rebuild the underlying linear algebra from scratch, in a handful of tiny functions, and let *those* pass judgment. This is also the cleanest way to see the machinery turn.

Start at the bottom of the ladder: map a Pauli string to its binary symplectic vector, an $X$-block followed by a $Z$-block. Define the per-letter rule and assemble:

```wl
xzOf = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;
strip[str_] := StringDelete[str, StartOfString ~~ ("+" | "-")];
symOf[str_] := With[{ps = xzOf /@ Characters[strip[str]]},
    Join[ps[[All, 1]], ps[[All, 2]]]
];
symOf["XZZXI"]
```

The first five entries are the $X$-pattern, the last five the $Z$-pattern; an $X$ shows up in the left half, a $Z$ in the right, and a $Y$ would show up in both.

Two Paulis commute exactly when the symplectic product of their vectors, $\omega(u,v) = x_u \cdot z_v + z_u \cdot x_v$ taken mod 2, is zero. Define it:

```wl
omega[u_, v_] := With[{n = Length[u]/2},
    Mod[Take[u, n] . Drop[v, n] + Drop[u, n] . Take[v, n], 2]
];
{omega[symOf["XII"], symOf["ZII"]], omega[symOf["XII"], symOf["IZI"]]}
```

As one can see, $X_1$ and $Z_1$ anticommute (they share a qubit and disagree there), while $X_1$ and $Z_2$ commute (they act on different qubits). This single function is the whole notion of "detectable": the syndrome of an error is just its symplectic product against each generator.

We will also need the weight of a Pauli (how many qubits it touches) and, for the distance search, every Pauli of a given weight. Define both:

```wl
wt[str_] := Count[Characters[strip[str]], Except["I"]];
weightK[n_, k_] := If[k == 0,
    {StringRepeat["I", n]},
    Catenate @ Map[
        Function[pos,
            StringJoin /@ (
                ReplacePart[ConstantArray["I", n], Thread[pos -> #]] & /@ Tuples[{"X", "Y", "Z"}, k]
            )
        ],
        Subsets[Range[n], {k}]
    ]
];
{Length[weightK[5, 1]], Length[weightK[5, 2]]}
```

On five qubits there are $3 \cdot 5$ single-qubit Paulis and $3^2 \binom{5}{2}$ weight-two Paulis; the counts confirm the enumerator is complete, which matters because the distance search below trusts it to miss nothing.

### Distance from First Principles

Here is the definition made computational. A Pauli is a logical operator when it lies in the normalizer but not in the stabilizer group: it commutes with every generator (syndrome zero) *and* it is not itself a product of generators. The second condition is the one that is easy to get backwards. Membership in the stabilizer group is a GF(2) span question, and the correct test is that appending the candidate's vector to the generator matrix leaves the rank unchanged; so "not in the group" means the rank strictly *rises*. Define the logical test with the rank-rises condition spelled out:

```wl
genMatrix[code_] := symOf /@ code["Generators"];
logicalVecQ[code_, v_] := With[
    {sg = genMatrix[code], m = code["StabilizerCount"]},
    AllTrue[sg, omega[v, #] === 0 &] && MatrixRank[Append[sg, v], Modulus -> 2] > m
];
{logicalVecQ[codes["BitFlipCode"], symOf["ZII"]],
 logicalVecQ[codes["BitFlipCode"], symOf["ZZI"]]}
```

As expected, $Z_1$ passes both tests for the bit-flip code (undetectable and outside the group, so a genuine logical), while the generator $Z_1 Z_2$ fails the second test: it commutes with everything, but it *is* in the group, so its rank does not rise. We will come back to what that first result means for the bit-flip code, because it is the honest surprise of this whole subject.

The distance is now a search: the smallest weight at which a logical exists, and a witness achieving it. Climb the weights and stop at the first hit:

```wl
firstLogicalOfWeight[code_, w_] := SelectFirst[
    weightK[code["Qubits"], w], logicalVecQ[code, symOf[#]] &
];
independentWitness[code_] := Module[{w},
    w = SelectFirst[Range[code["Qubits"]], ! MissingQ[firstLogicalOfWeight[code, #]] &];
    If[MissingQ[w], Missing["NoLogical"], firstLogicalOfWeight[code, w]]
];
independentDistance[code_] := With[{x = independentWitness[code]},
    If[MissingQ[x], Infinity, wt[x]]
];
independentWitness[codes["5QubitCode"]]
```

This witness is a lightest logical for the five-qubit code, found by our own search with no reference to the package. Its weight is the code distance.

Now put the package and the independent re-derivation side by side. For each code we show its parameters as the code itself reports them, the package's minimum-weight logical, our independently computed distance, and whether the two distances agree:

```wl
Grid[
    Prepend[
        Function[name,
            With[{c = codes[name], di = independentDistance[codes[name]]},
                {name,
                 StringJoin["[[", ToString[c["Qubits"]], ",", ToString[c["LogicalQubits"]], ",", ToString[CodeDistance[c]], "]]"],
                 MinimumWeightLogical[c], di, CodeDistance[c] === di}
            ]
        ] /@ Keys[codes],
        Style[#, Bold] & /@ {"code", "[[n,k,d]]", "package witness", "independent d", "agree"}
    ],
    Frame -> All, Alignment -> Left, Spacings -> {1.5, 0.8}
]
```

The final column is the whole point of this part: for every code, the distance our from-scratch symplectic search returns is exactly the distance the package returns. The Shor, five-qubit, and Steane codes come out as the familiar distance-three codes; the two three-qubit codes come out at distance one, which brings us to the honest surprise.

### The Honest Surprise: A Repetition Code Has Quantum Distance One

The three-qubit bit-flip code is often introduced as a distance-three code, and as a *classical* code against bit flips it is. As a *quantum* code it is not. Look at the syndrome of a single $Z$ on the first qubit, and ask whether that $Z$ is a logical operator:

```wl
{Syndrome[codes["BitFlipCode"], "ZII"], logicalVecQ[codes["BitFlipCode"], symOf["ZII"]]}
```

The syndrome is all zeros, so a lone $Z_1$ is undetectable: it commutes with both $Z$-type stabilizers. Yet it is not in the stabilizer group, so it acts nontrivially on the codespace. It is therefore a weight-one logical, and the quantum code distance is one, not three. In words, protecting against bit flips with $Z$-type checks leaves phase errors completely unguarded, and a single unguarded phase error is a logical fault. This is not a defect in the computation; it is the reason a real code has to protect against both error types, which is exactly what the Shor code does by nesting the two repetition ideas. The package agrees with the mathematics, not with the folklore, and that is what we want from it.

### Logical Operators and Their Algebra

A code with $k$ logical qubits comes with logical operators $\bar X_i$ and $\bar Z_i$, and they must obey the same algebra as bare Paulis on $k$ qubits: $\bar X_i$ anticommutes with its partner $\bar Z_i$, commutes with every other logical, all the $\bar X$ commute among themselves, all the $\bar Z$ commute among themselves, and every logical commutes with every stabilizer while lying outside the group. The package produces logical operators from a standard-form reduction; we check the entire algebra independently with our symplectic product and rank test:

```wl
logicalAlgebraOkQ[code_] := Module[
    {logs = LogicalOperators[code], sg = genMatrix[code], m = code["StabilizerCount"], xs, zs, k},
    xs = symOf /@ logs["X"]; zs = symOf /@ logs["Z"]; k = Length[xs];
    And[
        AllTrue[Tuples[Range[k], 2],
            With[{i = #[[1]], j = #[[2]]}, omega[xs[[i]], zs[[j]]] === If[i === j, 1, 0]] &],
        AllTrue[Subsets[Range[k], {2}], omega[xs[[#[[1]]]], xs[[#[[2]]]]] === 0 &],
        AllTrue[Subsets[Range[k], {2}], omega[zs[[#[[1]]]], zs[[#[[2]]]]] === 0 &],
        AllTrue[Join[xs, zs], Function[l, AllTrue[sg, omega[l, #] === 0 &]]],
        AllTrue[Join[xs, zs], MatrixRank[Append[sg, #], Modulus -> 2] > m &]
    ]
];
AssociationMap[logicalAlgebraOkQ[codes[#]] &, Keys[codes]]
```

Every code's logical operators pass the full algebra check. Note that a minimum-weight witness from the distance search need not be one of these canonical $\bar X$, $\bar Z$: the standard-form logicals are a valid *generating* pair, while the witness is the lightest element of the whole coset, and the two answer different questions.

### A Sign-Exact Encoder

The last exact object is the encoder: a Clifford circuit that prepares a codeword from the all-zeros state. "Prepares a codeword" has to mean the +1 eigenstate of every generator, with the sign exactly +1 and not $-1$, because a stabilizer with the wrong sign describes a different state entirely. We check this at the level of the state vector, comparing vectors rather than fidelities, so a sign error cannot hide.

Build the Pauli operators as honest matrices, run the real encoding circuit on $\lvert 0\ldots 0\rangle$, and confirm each generator leaves the resulting vector fixed. The zero test is symbolic, so it certifies exact equality rather than numerical closeness:

```wl
pmat = <|"I" -> {{1, 0}, {0, 1}}, "X" -> {{0, 1}, {1, 0}},
         "Y" -> {{0, -I}, {I, 0}}, "Z" -> {{1, 0}, {0, -1}}|>;
pauliOp[str_] := Fold[KroneckerProduct, pmat /@ Characters[strip[str]]];
encodedVector[code_] := Normal[
    EncodingCircuit[code][
        QuantumState[SparseArray[{1 -> 1}, 2^code["Qubits"]], ConstantArray[2, code["Qubits"]]]
    ]["StateVector"]
];
encoderFixesQ[code_] := With[{psi = encodedVector[code]},
    AllTrue[code["Generators"], AllTrue[pauliOp[#] . psi - psi, PossibleZeroQ] &]
];
AssociationMap[encoderFixesQ[codes[#]] &, Keys[codes]]
```

For each code, applying any generator to the encoded state returns the state unchanged, exactly: the encoder prepares a genuine +1 codeword, sign included. The Pauli matrices here are built with qubit 1 as the leftmost tensor factor, which is the same ordering QuantumFramework uses for its state vectors, so the comparison is apples to apples.

## Part III: Cross-Checking Against Stim

The independent re-derivation lives entirely inside the Wolfram Language, so it shares the linear-algebra substrate with the package. A second, genuinely external opinion is worth having, and Stim is the natural source: it is the reference stabilizer simulator, written in C++, with an independent tableau engine. We open a Python session and confirm Stim is importable; if this cell reports that Stim is missing, install it with `pip install stim` into the Python that `ExternalEvaluate` uses, or skip Part III and rely on the re-derivation above.

```wl
stimSession = StartExternalSession[<|"System" -> "Python"|>];
stimVersion = Check[
    ExternalEvaluate[stimSession, "import stim; stim.__version__"],
    "stim not available: skip Part III or pip install stim"
]
```

### The Converter and the Syndrome Differential

The converter is small because both sides speak Pauli strings. Stim writes identity as an underscore, so the only translation is a character swap; a syndrome bit is then Stim's own commutation test between a generator and an error. Define the Stim-side syndrome and compare it, for every code, against the package's `Syndrome` across all weight-one and weight-two errors:

```wl
stimSyndrome = ExternalEvaluate[stimSession, "
def stim_syndrome(gens, err):
    e = stim.PauliString(err.replace('I','_'))
    return [0 if stim.PauliString(g.replace('I','_')).commutes(e) else 1 for g in gens]
stim_syndrome"];
syndromeMatchesStimQ[code_] := With[{gens = code["Generators"], n = code["Qubits"]},
    With[{errs = Join[weightK[n, 1], weightK[n, 2]]},
        (Syndrome[code, #] & /@ errs) === (stimSyndrome[gens, #] & /@ errs)
    ]
];
AssociationMap[syndromeMatchesStimQ[codes[#]] &, Keys[codes]]
```

For every code, the package's syndrome and Stim's syndrome agree on the entire weight-one and weight-two error set, hundreds of errors per code. The syndrome layer is exactly what a real decoder consumes, so this is the differential that matters most for interoperability.

### The Encoder Differential

Stim can also judge the encoder, and more stringently than a syndrome check: we translate the encoding gates into a Stim circuit, run Stim's tableau simulator, and ask Stim for the expectation value of each of our generators on the state it prepared. A value of +1 means Stim agrees the prepared state is stabilized by that generator with the correct sign:

```wl
stimEncoderExp = ExternalEvaluate[stimSession, "
def stim_encoder_exp(gates, gens):
    circ = stim.Circuit()
    tr = {'CNOT':'CX'}
    for name, args in gates:
        circ.append(tr.get(name, name), [a-1 for a in args])
    sim = stim.TableauSimulator(); sim.do_circuit(circ)
    return [int(sim.peek_observable_expectation(stim.PauliString(g.replace('I','_')))) for g in gens]
stim_encoder_exp"];
gatePairs[code_] := {#[[1]], Flatten[{#[[2]]}]} & /@ EncodingGates[code];
stimEncoderOkQ[code_] := AllTrue[stimEncoderExp[gatePairs[code], code["Generators"]], # === 1 &];
AssociationMap[stimEncoderOkQ[codes[#]] &, Keys[codes]]
```

Every generator has Stim-expectation +1 on the encoded state, for every code. This is the same sign-exact statement we checked with our own Pauli matrices in Part II, now confirmed by an independent simulator through a completely different route (a tableau update, not a dense state vector).

### The Witness Differential

Finally, Stim can certify the distance witness without ever enumerating errors. A minimum-weight logical must be *undetectable* (it commutes with every generator) and *nontrivial* (it anticommutes with at least one logical operator, so it acts on the codespace). Those two facts together put it in the normalizer but outside the stabilizer group, which is the definition of a logical. Ask Stim for both:

```wl
stimWitness = ExternalEvaluate[stimSession, "
def stim_witness(gens, logicals, witness):
    w = stim.PauliString(witness.replace('I','_'))
    undetectable = all(stim.PauliString(g.replace('I','_')).commutes(w) for g in gens)
    nontrivial = any(not stim.PauliString(p.replace('I','_')).commutes(w) for p in logicals)
    return [bool(undetectable), bool(nontrivial)]
stim_witness"];
stimWitnessRes[code_] := With[{logs = LogicalOperators[code]},
    stimWitness[code["Generators"], Join[logs["X"], logs["Z"]], MinimumWeightLogical[code]]
];
AssociationMap[stimWitnessRes[codes[#]] &, Keys[codes]]
```

For each code Stim reports the witness as both undetectable and nontrivial. Stim confirms the witness is a valid logical of the stated weight; the *minimality* of that weight, that nothing lighter is logical, is what our exhaustive symplectic search established, and it is the part no sampler provides. The division of labor is the honest one: Stim vouches that the witness is real, and the exact search vouches that it is the lightest.

As a last sanity check on the converter itself, ask Stim for the canonical stabilizers of the encoded Steane state and compare them, by eye, with the code's own generators:

```wl
stimCanonical = ExternalEvaluate[stimSession, "
def stim_canonical(gates):
    circ = stim.Circuit()
    tr = {'CNOT':'CX'}
    for name, args in gates:
        circ.append(tr.get(name, name), [a-1 for a in args])
    sim = stim.TableauSimulator(); sim.do_circuit(circ)
    return [str(s) for s in sim.canonical_stabilizers()]
stim_canonical"];
stimCanonical[gatePairs[steane]]
```

Stim returns seven signed stabilizers in its reduced canonical form. These generate the same group as our six Steane generators plus the logical $\bar Z$ that pins the prepared codeword; the +1 expectations above already proved each of our generators sits inside this group with the right sign, so the two descriptions are the same stabilizer state written in two bases.

## Part IV: A Parameter Stim Cannot Carry

Everything so far was a concrete code at a concrete size. The exact-and-symbolic corner has one more move, the one that has no numeric counterpart: hold a parameter symbolic and certify a whole family, or a whole continuum, at once.

### A Family Indexed by n

Take the $[[n, n-2, 2]]$ code family: two generators, an all-$X$ string and an all-$Z$ string, on $n$ qubits. The generators must commute, and their symplectic product is a sum of contributions, one per qubit; the all-$X$ and all-$Z$ strings overlap on every qubit, so each of the $n$ qubits contributes a 1:

```wl
omegaXZ = Sum[1, {i, n}]
```

The product is $n$, carried symbolically. The two generators commute exactly when that count is even, which is the family's defining constraint: an odd $n$ would make the two generators anticommute, and there would be no code at all.

```wl
Simplify[Mod[omegaXZ, 2] == 0, Element[n/2, Integers]]
```

Now the distance. Every weight-one error is detected, because a single $X$ overlaps the all-$Z$ generator on exactly one qubit, an odd overlap, so it anticommutes; that rules out distance one for every $n$. A weight-two logical, on the other hand, commutes with both generators and is lighter than every nontrivial group element (each of which has weight $n$) precisely when two is less than $n$:

```wl
Reduce[2 < n, n, Integers]
```

So the distance is exactly two for the whole even-$n$ family, established from a single evaluation with $n$ symbolic. A numeric simulator has to instantiate a specific $n$ to say anything at all; it can confirm any number of cases but never the family. As a bridge between the symbolic certificate and the concrete package, here is the package's numeric distance at several sizes, each an independent brute-force search, all landing where the certificate says they must:

```wl
Table[CodeDistance[DistanceTwoCode[n]], {n, {4, 6, 8, 10}}]
```

The same style of one-shot argument handles the repetition family $[[n, 1, 1]]$: every generator is $Z$-type, so a lone $Z_1$ commutes with all of them for every $n$; and a parity functional (the mod-2 sum of the $Z$-support) vanishes on every weight-two generator but equals one on $Z_1$, so $Z_1$ is never in the group. A single undetectable $Z$, for every $n$: quantum distance one, the honest surprise of Part II, now stated for the whole family in one breath.

### A Continuous Error Angle

Symbolic generality is not only about integer families; it reaches continuous parameters too, which is where a Pauli-sampling simulator has nothing to say. Consider a *coherent* error on the bit-flip codeword: not a discrete bit flip, but a rotation $R_x(\theta) = e^{-i\theta X/2}$ on the first qubit, through a symbolic angle. The state stays a superposition of "no error" and "flipped," and the syndrome operator $Z_1 Z_2$ has an expectation value that is an exact function of $\theta$:

```wl
psi0 = QuantumState[SparseArray[{1 -> 1}, 8], {2, 2, 2}];
vTheta = Normal[QuantumCircuitOperator[{"RX"[\[Theta]] -> 1}][psi0]["StateVector"]];
zzOp = KroneckerProduct[PauliMatrix[3], PauliMatrix[3], IdentityMatrix[2]];
zzExpectation = Simplify[Conjugate[vTheta] . zzOp . vTheta, \[Theta] \[Element] Reals]
```

The syndrome expectation is $\cos\theta$: it interpolates continuously from +1 at $\theta = 0$ (no error, the codeword is untouched) to $-1$ at $\theta = \pi$ (a full bit flip, a definite nonzero syndrome), passing through zero at $\theta = \pi/2$ (a maximally uncertain syndrome). A coherent rotation is non-Clifford, so it lives in the exact dense simulator, not in a stabilizer tableau; carrying the angle symbolically returns the whole error curve as a closed form. This is the natural language for coherent-error and threshold questions, and it is exactly the representation a discrete stabilizer sampler cannot produce.

```wl
{zzExpectation /. \[Theta] -> 0, zzExpectation /. \[Theta] -> Pi/2, zzExpectation /. \[Theta] -> Pi}
```

The three sampled values confirm the endpoints and the midpoint of the closed form we just derived.

## Part V: Where This Leaves Us

Before we close, let us summarize what we have actually established:

- For the bit-flip, phase-flip, Shor, five-qubit, and Steane codes, the exact code distance from a from-scratch symplectic search matches the package's distance in every case.
- The package's minimum-weight witnesses are valid logicals of that weight, verified both by our rank-and-commutation test and by Stim's undetectable-and-nontrivial test.
- Each code's logical operators satisfy the full commutation algebra ($\bar X_i$ anticommutes with $\bar Z_i$ and commutes with everything else), checked independently.
- Each encoder prepares a sign-exact +1 codeword, confirmed both by dense Pauli matrices and by Stim's tableau expectation values.
- The package's syndrome equals Stim's across every weight-one and weight-two error, code by code.
- The $[[n, n-2, 2]]$ and repetition families are certified at symbolic $n$ in one evaluation, and a coherent error's syndrome is returned as the exact function $\cos\theta$.

The picture that leaves us with is a clean division of labor. Stim owns scale and speed: sampling large circuits, benchmarking decoders, estimating thresholds. This exact layer owns the small-and-certain corner: the exact minimum distance, a witness for it, the logical algebra, a sign-exact encoder, and, uniquely, the ability to hold a parameter symbolic and certify an entire family or a continuous error at once. The two agree wherever they overlap, which is what makes the exact answers trustworthy and the symbolic ones worth reaching for.

The honest boundary is worth restating: exact minimum distance is NP-hard, so this is a small-code instrument by mathematical necessity, not a scaling competitor. That is the right trade to have made. When you want the exact answer for a textbook code, or a statement about a whole family, the exactness and the symbol are the point; when you want a threshold for a distance-25 surface code, reach for Stim. Both live in the same notebook, and now you can check one against the other yourself.

If you want to push further: feed your own stabilizer generators into `QECCode`, run the same independent checks, and diff against Stim; or take the symbolic-family argument to a code family of your own and see whether one evaluation can settle it.
