---
Template: Default
Title: Stabilizer Codes, Cross-Checked Against Stim
Author: Mads Bahrami
---

# Stabilizer Codes, Cross-Checked Against Stim

<!-- #| style: Subtitle -->
The shipped QuantumFramework stabilizer engine agrees with Stim where they overlap; the Wolfram Language around it keeps an unresolved measurement outcome as an algebraic symbol, derives a code family's condition at symbolic $n$, and returns a continuous error's syndrome in closed form.

<!-- #| style: Author -->
Mads Bahrami (last updated: Aug 25, 2026)

<!-- #| style: Affiliation -->
Wolfram Research Inc, USA

### Setting the Stage: How This Notebook Flows

This notebook is a computation-first tour of the stabilizer layer that ships in the QuantumFramework paclet, checked against [Stim](https://github.com/quantumlib/Stim), the field-standard stabilizer simulator. Part I frames the niche. Part II takes the shipped `PauliStabilizer` engine on the standard small codes (the three-qubit bit-flip and phase-flip codes, the nine-qubit Shor code, the five-qubit perfect code, and the seven-qubit Steane code) and checks, against Stim, that its stabilizer generators form a consistent group and its syndromes agree across every one- and two-qubit error, while a state-vector test confirms every codeword carries the correct sign; it then exercises a symbolic measurement that keeps an unresolved outcome as an algebraic symbol, which Stim has no counterpart for. Part III builds an exact code analysis from GF(2) first principles and computes the exact distance, a minimum-weight witness, and the logical algebra, cross-checking the witness against Stim. Part IV carries a symbolic parameter through the algebra: it derives a code family's existence condition at symbolic $n$, and returns a coherent error's syndrome as a closed-form function of a continuous angle. Part V draws the boundary.

The division of labor is clean. Stim owns scale and speed: sampling large circuits, benchmarking decoders, estimating thresholds over thousands of qubits. This layer owns the exact-and-symbolic corner: the exact minimum distance (which is NP-hard in the worst case, so an exhaustive search like the one here is a small-$n$ instrument), a sign-exact stabilizer state, and a parameter kept symbolic where the algebra allows. The two agree wherever they overlap, and that agreement is what makes the shared states and syndromes trustworthy; the distance's minimality rests on the exhaustive search, and the symbolic results on the algebra. Neither competes with the other.

Before we start, a few notes on the environment. This is a live Wolfram notebook: evaluate the cells from top to bottom. Some cells define objects (the code table, a handful of small symplectic functions, a Python session running Stim) that later cells use, so order matters. My suggestion is to read each output and its meaning first, and only then unpack the input code that produced it. You are not locked into the codes as given: feed in your own stabilizer generators, push the symbolic family further, and re-run. If you have questions, reach out at [quantum@wolfram.com](mailto:quantum@wolfram.com).

### Prerequisites and How to Read This

You should be comfortable with the stabilizer formalism at the level of a first course: Pauli operators, stabilizer generators, the codespace as their joint $+1$ eigenspace, and logical operators as the normalizer modulo the stabilizer group. The binary symplectic picture (each Pauli as a vector over GF(2), commutation as a symplectic inner product) is built up from scratch in Part III, so no prior fluency there is required. The Stim cross-check calls Python through [`ExternalEvaluate`](https://reference.wolfram.com/language/ref/ExternalEvaluate.html); if you do not have Stim installed, those cells say so and you can skip them, keeping the shipped-engine outputs and the from-scratch re-derivation as your checks.

Let's start!

## Part I: A Shipped Engine and an Open Niche Beside It

A stabilizer code on $n$ physical qubits is defined by a set of commuting Pauli operators, its stabilizer generators. The codespace is their simultaneous $+1$ eigenspace, and if there are $m$ independent generators the code stores $k = n - m$ logical qubits. An error is *detectable* when it anticommutes with at least one generator (it flips that generator's measured eigenvalue, a syndrome bit), and *undetectable* when it commutes with all of them. The dangerous errors are the undetectable ones that still act nontrivially on the codespace: these are the logical operators, and the code distance is the weight of the lightest one.

QuantumFramework ships a compiled stabilizer engine, `PauliStabilizer`, built on the Aaronson-Gottesman tableau. It carries a stabilizer state as $O(n^2)$ bits rather than $2^n$ amplitudes, applies Clifford gates and measurements in polynomial time, and already comes with the standard codes as named states. That engine is fast and scales; this notebook is not about its speed. It is about two things speed does not settle: does it get the *signs* and *syndromes* exactly right, checked against an independent engine, and can the surrounding language say something a stabilizer sampler cannot say at all.

Why is exactness worth isolating when Stim samples so well? Because computing the exact minimum distance of a stabilizer code is NP-hard in the worst case (Kapshikar and Kundu, [arXiv:2203.04262](https://arxiv.org/abs/2203.04262)): there is no known efficient general method, and the exhaustive weight-search this notebook uses is a small-code instrument by construction. That is a property of this algorithm together with the worst-case hardness, not a claim that every large code is out of reach; structured families and smarter exact solvers go further. The textbook codes are simply where an exhaustive exact treatment earns its keep. And why symbolic generality? Because a numeric simulator represents one concrete code at a time: fixed $n$, fixed generators, fixed angles. A symbolic engine can carry $n$ as a symbol and settle a family-level statement that a simulator tied to one concrete $n$ cannot reach, or carry an error angle $\theta$ and return a syndrome expectation as an exact function of it. Neither is a thing Stim is built to do.

Load the framework from a local checkout (adjust the path to your own), and build the Steane code as a shipped stabilizer state:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]
```

```wl
steane = PauliStabilizer["SteaneCode"];
steane["Stabilizers"]
```

The object is the codeword $\lvert 0_L\rangle$, carried as its stabilizer group. The first three generators are pure $X$-checks and the next three pure $Z$-checks, the defining CSS structure of the Steane code; the seventh, the all-$Z$ string, is the logical $\bar Z$ that pins this particular codeword inside the two-dimensional code space. Reading the parity checks apart from the logical is the one bookkeeping step the rest of Part II rests on.

## Part II: The Shipped Stabilizer Layer, Cross-Checked Against Stim

Collect the five codes into one table. The five-qubit, Steane, and Shor codes ship as named `PauliStabilizer` states; the two three-qubit repetition codes are not named states, so we build them directly from their generators. For each we record the shipped state, its parity checks (the generators other than the logical), and the logical $\bar X$, $\bar Z$ pair we will use later:

```wl
codes = <|
   "bit-flip" -> <|"state" -> PauliStabilizer[{"ZZI", "IZZ", "ZZZ"}],
     "checks" -> {"ZZI", "IZZ"}, "logicals" -> {"XXX", "ZII"}|>,
   "phase-flip" -> <|"state" -> PauliStabilizer[{"XXI", "IXX", "XXX"}],
     "checks" -> {"XXI", "IXX"}, "logicals" -> {"XII", "ZZZ"}|>,
   "5-qubit" -> <|"state" -> PauliStabilizer["5QubitCode"],
     "checks" -> Most[PauliStabilizer["5QubitCode"]["Stabilizers"]], "logicals" -> {"XXXXX", "ZZZZZ"}|>,
   "Steane" -> <|"state" -> PauliStabilizer["SteaneCode"],
     "checks" -> Most[PauliStabilizer["SteaneCode"]["Stabilizers"]], "logicals" -> {"XXXXXXX", "ZZZZZZZ"}|>,
   "Shor" -> <|"state" -> PauliStabilizer["9QubitCode"],
     "checks" -> Most[PauliStabilizer["9QubitCode"]["Stabilizers"]], "logicals" -> {"XXXXXXXXX", "ZZZZZZZZZ"}|>
   |>;
Keys[codes]
```

### A Sign-Exact +1 Codeword

The first thing to get right is the sign. "The state is a codeword" has to mean the $+1$ eigenstate of every generator, with the sign exactly $+1$ and not $-1$, because a stabilizer with the wrong sign describes a different state entirely. We check this at the level of the state vector, comparing vectors rather than fidelities. The choice matters: a fidelity $\lvert\langle\psi\rvert G\psi\rangle\rvert$ is $1$ whether $G$ fixes $\lvert\psi\rangle$ or sends it to $-\lvert\psi\rangle$, because it quotients out the very sign we are testing, so comparing the vectors is what keeps that sign visible. The shipped state has exact algebraic amplitudes (the Steane codeword's nonzero entries are $1/(2\sqrt 2)$, for instance), so the comparison is a symbolic zero test, not a numerical tolerance.

Build the Pauli operators as dense matrices (qubit 1 as the leftmost tensor factor, the ordering `PauliStabilizer` uses), read the state vector, and confirm each generator, sign included, leaves it fixed exactly:

```wl
strip[s_] := StringDelete[s, StartOfString ~~ ("+" | "-")];
sgn[s_] := If[StringStartsQ[s, "-"], -1, 1];
pauliOp[s_] := sgn[s] Fold[KroneckerProduct,
   PauliMatrix /@ (Characters[strip[s]] /. {"I" -> 0, "X" -> 1, "Y" -> 2, "Z" -> 3})];
denseFixesQ[c_] := With[{psi = Normal[c["state"]["State"]["StateVector"]]},
   AllTrue[c["state"]["Stabilizers"], AllTrue[pauliOp[#] . psi - psi, PossibleZeroQ] &]];
denseFixesQ /@ codes
```

For every code, applying any generator to the codeword returns the codeword unchanged, exactly: the shipped state is a $+1$ codeword, sign included. This is an internal-consistency check rather than an independent reconstruction, and worth reading as such: the state vector is materialized from the same tableau whose generators it tests, so what it certifies is that QuantumFramework's own state vector and stabilizer signs agree exactly. The externally independent checks come next, when Stim weighs in on the group and the syndromes, and in Part III, where a from-scratch search recomputes the distance. The same exact amplitudes make normalization a symbolic identity, not a numerical near-miss, so each codeword has unit length as an exact statement:

```wl
AllTrue[codes, With[{psi = Normal[#["state"]["State"]["StateVector"]]}, Simplify[Conjugate[psi] . psi] === 1] &]
```

Normalization holds as an exact identity, not a value merely close to one. Stim now enters as a second, independent engine, to confirm that these same generators form a consistent stabilizer group and read every syndrome the way the shipped engine does.

### Opening a Stim Session

The check so far lives entirely inside the Wolfram Language. A second, genuinely external opinion is worth having, and Stim is the natural source: the reference stabilizer simulator, written in C++, with an independent tableau engine. Open a Python session and confirm Stim is importable; if this reports Stim missing, install it with `pip install stim` into the Python that `ExternalEvaluate` uses, or skip the Stim cells and rely on the in-language checks.

```wl
stimSession = StartExternalSession[<|"System" -> "Python"|>];
stimVersion = Check[
   ExternalEvaluate[stimSession, "import stim; stim.__version__"],
   "stim not available: skip the Stim cells"]
```

Stim reads identity as an underscore, so the only translation between the two sides is a character swap. Define that once and reuse it for the differentials below:

```wl
toStim[s_] := StringReplace[s, "I" -> "_"];
```

### Stim Accepts the Generators, and the Syndromes Agree

Two differentials carry Part II. The first hands Stim our stabilizer generators and asks whether they even form a consistent stabilizer group: `stim.Tableau.from_stabilizers` builds a tableau only if the generators are independent and mutually commuting, and it hands back its own canonical form of the same group. This confirms the generators are a valid, full-rank stabilizer set; it takes the strings as given rather than reconstructing the state independently, so the load-bearing external check is the syndrome differential that follows.

```wl
stimGroupSize = ExternalEvaluate[stimSession, "
def group_size(gens):
    t = stim.Tableau.from_stabilizers([stim.PauliString(g) for g in gens], allow_redundant=False, allow_underconstrained=False)
    return len(t.to_stabilizers(canonicalize=True))
group_size"];
stimAcceptsQ[c_] := With[{gens = toStim /@ c["state"]["Stabilizers"]},
   stimGroupSize[gens] === Length[gens]];
stimAcceptsQ /@ codes
```

The second differential is the one that matters most for interoperability, because the syndrome is exactly what a real decoder consumes. For an error $E$ on a codeword, the syndrome bit of a check $g$ is $0$ when $E$ commutes with $g$ and $1$ when it anticommutes. On the Wolfram side we apply the error to the shipped state and read each check's expectation value; on the Stim side it is a commutation test. Compare them, for every code, across all weight-one and weight-two errors:

```wl
weightK[n_, k_] := If[k == 0, {StringRepeat["I", n]},
   Catenate@Map[Function[pos,
      StringJoin /@ (ReplacePart[ConstantArray["I", n], Thread[pos -> #]] & /@ Tuples[{"X", "Y", "Z"}, k])],
     Subsets[Range[n], {k}]]];
applyErr[ps_, e_] := Fold[
   With[{c = StringTake[e, {#2}]}, If[c === "I", #1, #1[c, #2]]] &,
   ps, Range[StringLength[e]]];
ourSyndrome[c_, e_] := (1 - (applyErr[c["state"], e]["Expectation", #] & /@ c["checks"]))/2;
stimSyndrome = ExternalEvaluate[stimSession, "
def stim_syndrome(gens, err):
    e = stim.PauliString(err)
    return [0 if stim.PauliString(g).commutes(e) else 1 for g in gens]
stim_syndrome"];
syndromeMatchesStimQ[c_] := With[{errs = Join[weightK[c["state"]["Qubits"], 1], weightK[c["state"]["Qubits"], 2]]},
   (ourSyndrome[c, #] & /@ errs) === (stimSyndrome[toStim /@ c["checks"], toStim[#]] & /@ errs)];
syndromeMatchesStimQ /@ codes
```

For every code, the shipped engine's syndrome and Stim's syndrome agree on the entire weight-one and weight-two error set, which is hundreds of errors for the larger codes. Combined with the sign-exact codeword above and Stim's acceptance of the generators, the shipped stabilizer layer and Stim describe the same states and read the same syndromes, wherever the two overlap.

### The One Thing Stim Cannot Do: A Symbolic Measurement

Everything so far Stim can also do, and does. Here is something with no Stim counterpart: a measurement outcome that is a coin flip need not be drawn as a number; it can be carried as an algebraic *symbol*. `PauliStabilizer`'s `"SymbolicMeasure"` writes each unresolved random outcome as a fresh formal bit $s_k$ inside a stabilizer sign, so the whole state stays a single object and the measurement is deferred, not sampled. Measure qubit 1 of a Bell pair symbolically and read the generators:

```wl
sm = PauliStabilizer[2][{"H", "CNOT"}]["SymbolicMeasure", 1];
sm["Stabilizers"]
```

The sign of the first generator now carries a formal bit $s$ (displayed as $s_1$ here; the subscript is a running session counter, so a re-run allocates $s_2$, $s_3$, and so on), $+1$ at $s = 0$ and $-1$ at $s = 1$: the formal bit *is* the unmade measurement. Extract whichever bit was allocated and substitute both values to reproduce the two post-measurement branches exactly:

```wl
s = First@Cases[sm["Phase"], \[FormalS][_], Infinity];
{sm["SubstituteOutcomes", {s -> 0}]["Stabilizers"],
 sm["SubstituteOutcomes", {s -> 1}]["Stabilizers"]}
```

The deeper point is that a formal bit is spent only where the randomness is genuine. Prepare two qubits in $\lvert +\rangle$, so each carries a fair-coin $Z$-outcome, and copy their parity onto a third with two `CNOT`s. The prepared state carries a conserved $Z$-parity, the stabilizer $Z_1 Z_2 Z_3$:

```wl
parityCode = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
parityCode["Stabilizers"]
```

The last generator, $Z_1 Z_2 Z_3$, is the conservation law: because it stabilizes the state, the three $Z$-outcomes cannot be independent. Measuring all three symbolically therefore allocates only two formal bits, one per genuinely random qubit, not three:

```wl
smPar = parityCode["SymbolicMeasure", {1, 2, 3}];
freeBits = DeleteDuplicates@Cases[smPar["Phase"], \[FormalS][_], Infinity]
```

Two bits, not three. The third outcome is not discarded: the engine stores it as an explicit function of those two free bits, in the `"Outcomes"` map of the measured object.

```wl
smPar["Outcomes"]
```

That is the engine's closed-form record of the third qubit's outcome as a Boolean function of the two free bits, written in an unreduced form; over $\{0, 1\}$ it reduces to $m_1 \oplus m_2$. This is the object no sampler carries: the dependency held as an algebraic expression rather than a statistic. Substituting the two free bits over all assignments makes the same relation legible as a table:

```wl
readBits[ps_] := First@Keys@Select[
    AssociationThread[Tuples[{0, 1}, ps["Qubits"]], Normal[ps["State"]["StateVector"]]], # =!= 0 &];
readBits[smPar["SubstituteOutcomes", Thread[freeBits -> #]]] & /@ Tuples[{0, 1}, Length[freeBits]]
```

The third bit is the sum of the first two in every row, consistent with the conserved $Z_1 Z_2 Z_3$: the engine spent a formal bit only on the free coins and stored the determined one as the function above, rather than allocating a third random bit. That parity is the syndrome of the simplest possible code, a single $Z$-check, and the same relation, written for faults rather than measurements, is what a decoder consumes.

The difference from a sampler here is specific, and worth stating precisely rather than broadly. Stim propagates discrete outcome and fault dependencies very efficiently, through Pauli-frame sampling, sweep bits, and detector-error models that map error mechanisms to detector symptoms and logical frame changes; it does not re-traverse the circuit for every shot. What it does not do is keep a dependency as an *algebraic* object. Here the relation is held as the `"Outcomes"` symbol shown above, something you can substitute, simplify, or reason about in closed form, not only sample. That algebraic form, not the propagation of the dependency, is the part with no Stim counterpart.

## Part III: An Independent Exact Analysis, Cross-Checked Against Stim

The shipped engine reports syndromes and states. It does not, on its own, hand you the exact code distance or a lightest logical operator, and neither does Stim. Part III builds that exact analysis from the ground up, in a few small functions over GF(2), and then checks its answers both against the definition and against Stim. Rebuilding the linear algebra from scratch (rather than calling a packaged distance routine to check itself) is also the cleanest way to watch the machinery turn.

### The Symplectic Picture: A Pauli as a Vector over GF(2)

Each Pauli string on $n$ qubits becomes a binary vector of length $2n$: an $X$ contributes to the first half, a $Z$ to the second, a $Y$ to both. Two Paulis commute exactly when a bilinear form over GF(2), the symplectic product $\omega$, vanishes between their vectors. Define the per-letter rule, the vector, and the product:

```wl
xzOf = <|"I" -> {0, 0}, "X" -> {1, 0}, "Y" -> {1, 1}, "Z" -> {0, 1}|>;
symOf[s_] := With[{p = xzOf /@ Characters[strip[s]]}, Join[p[[All, 1]], p[[All, 2]]]];
omega[u_, v_] := With[{n = Length[u]/2},
   Mod[Take[u, n] . Drop[v, n] + Drop[u, n] . Take[v, n], 2]];
{omega[symOf["XII"], symOf["ZII"]], omega[symOf["XII"], symOf["IZI"]]}
```

As one can see, $X_1$ and $Z_1$ anticommute (they share a qubit and disagree there), while $X_1$ and $Z_2$ commute (they act on different qubits). This single function is the whole notion of "detectable": the syndrome of an error is just its symplectic product against each generator.

### Distance from First Principles

A Pauli is a logical operator when it lies in the normalizer but not in the stabilizer group: it commutes with every generator (zero syndrome) *and* it is not itself a product of generators. The second condition is the one that is easy to get backwards. Membership in the group is a GF(2) span question, and the correct test is that appending the candidate's vector to the check matrix leaves the rank unchanged; so "not in the group" means the rank strictly *rises*. Spell that out as `logicalVecQ`, then make the distance a weight-increasing search that stops at the first logical it finds. The search builds the check matrix and its rank once and reuses them across every candidate, rather than rebuilding them at each step:

```wl
rank2[rows_] := MatrixRank[rows, Modulus -> 2];
logicalVecQ[checks_, v_] := With[{sg = symOf /@ checks},
   With[{m = rank2[sg]}, AllTrue[sg, omega[v, #] === 0 &] && rank2[Append[sg, v]] > m]];
distanceWitness[checks_] := With[{sg = symOf /@ checks, n = StringLength[strip[First[checks]]]},
   With[{m = rank2[sg]},
    With[{logicalQ = Function[v, AllTrue[sg, omega[v, #] === 0 &] && rank2[Append[sg, v]] > m]},
     Catch[
      Scan[Function[w,
        With[{cand = SelectFirst[weightK[n, w], logicalQ[symOf[#]] &]},
         If[! MissingQ[cand], Throw[{w, cand}]]]],
       Range[n]];
      Missing[]]]]];
distanceWitness[codes["5-qubit"]["checks"]]
```

The output is a weight and a witness at that weight: a lightest logical for the five-qubit code, found by our own search with no reference to any distance routine. Its weight is the code distance. Now put the whole table together, showing each code's parameters, the independent distance, and the witness the search returned:

```wl
Grid[
   Prepend[
    Function[name,
     With[{c = codes[name], n = codes[name]["state"]["Qubits"], dw = distanceWitness[codes[name]["checks"]]},
      {name, StringJoin["[[", ToString[n], ",", ToString[n - rank2[symOf /@ c["checks"]]], ",", ToString[dw[[1]]], "]]"], dw[[2]]}]] /@ Keys[codes],
    Style[#, Bold] & /@ {"code", "[[n,k,d]]", "witness"}],
   Frame -> All, Alignment -> Left, Spacings -> {2, 0.8}]
```

The Shor, five-qubit, and Steane codes come out as the familiar distance-three codes. The two three-qubit codes come out at distance one, the result worth pausing on.

### The Surprise: A Repetition Code Has Quantum Distance One

The three-qubit bit-flip code is often introduced as a distance-three code, and as a *classical* code against bit flips it is. As a *quantum* code it is not. Look at the syndrome of a single $Z$ on the first qubit, and ask whether that $Z$ is a logical operator:

```wl
{AllTrue[codes["bit-flip"]["checks"], omega[symOf["ZII"], symOf[#]] === 0 &],
 logicalVecQ[codes["bit-flip"]["checks"], symOf["ZII"]]}
```

The lone $Z_1$ commutes with both $Z$-type checks, so it is undetectable, yet it is not in the stabilizer group, so it acts nontrivially on the codespace. It is therefore a weight-one logical, and the quantum code distance is one. In words, protecting against bit flips with $Z$-type checks leaves phase errors completely unguarded, and a single unguarded phase error is a logical fault. This is not a defect in the computation; it is the reason a real code protects against both error types, which is exactly what the Shor code does by nesting the two repetition ideas. The analysis reports the mathematics, not the folklore.

### Logical Operators and Their Algebra

A code with one logical qubit comes with a logical pair $\bar X$, $\bar Z$ that must obey the algebra of a bare qubit: $\bar X$ anticommutes with $\bar Z$, both commute with every stabilizer, and neither lies in the group. Check the whole algebra independently with the symplectic product and the rank test, for every code:

```wl
logicalAlgebraQ[c_] := With[{sg = symOf /@ c["checks"], xbar = c["logicals"][[1]], zbar = c["logicals"][[2]]},
   With[{m = rank2[sg]},
    And[omega[symOf[xbar], symOf[zbar]] === 1,
     AllTrue[c["checks"], omega[symOf[xbar], symOf[#]] === 0 &],
     AllTrue[c["checks"], omega[symOf[zbar], symOf[#]] === 0 &],
     rank2[Append[sg, symOf[xbar]]] > m, rank2[Append[sg, symOf[zbar]]] > m]]];
logicalAlgebraQ /@ codes
```

Every code's logical pair passes the full algebra check. Note that the minimum-weight witness from the distance search need not be one of these canonical operators: the logical pair is a valid *generating* set, while the witness is the lightest element of the whole coset, and the two answer different questions.

### Stim Confirms the Witness

Stim can certify the distance witness without enumerating a single error. A minimum-weight logical must be *undetectable* (it commutes with every check) and *nontrivial* (it anticommutes with at least one logical operator, so it acts on the codespace). Those two facts together place it in the normalizer but outside the stabilizer group, which is the definition of a logical. Ask Stim for both:

```wl
stimWitness = ExternalEvaluate[stimSession, "
def stim_witness(gens, logs, w):
    W = stim.PauliString(w)
    undetectable = all(stim.PauliString(g).commutes(W) for g in gens)
    nontrivial = any(not stim.PauliString(p).commutes(W) for p in logs)
    return [bool(undetectable), bool(nontrivial)]
stim_witness"];
stimWitnessQ[c_] := With[{w = distanceWitness[c["checks"]][[2]]},
   stimWitness[toStim /@ c["checks"], toStim /@ c["logicals"], toStim[w]] === {True, True}];
stimWitnessQ /@ codes
```

For each code Stim reports the witness as both undetectable and nontrivial. Stim vouches that the witness is a valid logical; the *minimality* of its weight, that nothing lighter is logical, is what our exhaustive symplectic search established, and it is the part no sampler provides. The division of labor is clear: Stim confirms the witness is a real logical, and the exact search confirms it is the lightest.

## Part IV: A Parameter Stim Cannot Carry

Everything so far was a concrete code at a concrete size. The exact-and-symbolic corner has one more move, the one with no numeric counterpart: hold a parameter symbolic and reason about a whole family, or a whole continuum, rather than one instance at a time.

### A Family Indexed by n (even n ≥ 4)

Take the $[[n, n-2, 2]]$ code family: two generators, an all-$X$ string and an all-$Z$ string, on $n$ qubits. It is defined for even $n \geq 4$, and both halves of that domain fall out of the symbolic analysis below. Start with the symplectic product of the two generators: it is a sum of per-qubit overlaps, and the all-$X$ and all-$Z$ strings meet on every qubit, so each of the $n$ qubits contributes a $1$ and the product is $n$:

```wl
omegaXZ = Sum[1, {i, n}]
```

Carried symbolically, the product is $n$. The generators commute when it is even and anticommute when it is odd. Substitute an even index $n = 2k$ and an odd index $n = 2k+1$ and reduce the commutator mod 2 to derive both cases at once, in the integer $k$:

```wl
Simplify[{Mod[omegaXZ /. n -> 2 k, 2], Mod[omegaXZ /. n -> 2 k + 1, 2]}, k \[Element] Integers]
```

The two parities are $0$ and $1$: the two generators commute exactly on the even integers, so an odd $n$ is not a code at all. That is the first half of the domain, derived for the whole family in $k$ rather than checked case by case.

The second half is the distance, a size-independent argument the from-scratch search then confirms at sample sizes. Build the family's checks at any size, and put a candidate weight-two operator, $X_1 X_2$, to the same logical test the concrete codes used:

```wl
familyChecks[n_] := {StringRepeat["X", n], StringRepeat["Z", n]};
logicalVecQ[familyChecks[6], symOf["XXIIII"]]
```

$X_1 X_2$ is a logical: it meets the all-$Z$ generator on two qubits (an even overlap, so it commutes), meets the all-$X$ generator trivially, and (for $n > 2$) is not a product of the generators. That argument is size-independent, even though the check above ran at a representative $n$: $X_1 X_2$ has weight two, while the only nontrivial group elements are $X^{\otimes n}$, $Z^{\otimes n}$, and $Y^{\otimes n}$, each of weight $n$. So for even $n \geq 4$ the lightest logical has weight two: no weight-one error is undetectable, and a weight-two logical exists. The lower bound $n > 2$ is not decoration; at $n = 2$ the construction degenerates, with $k = 0$ and $X_1 X_2 = X^{\otimes 2}$ equal to the generator itself, so there is no logical to weigh. Run the search at several sizes to confirm the distance lands where the argument says:

```wl
Table[distanceWitness[familyChecks[n]][[1]], {n, {4, 6, 8, 10}}]
```

A numeric simulator has to instantiate a specific $n$ to say anything at all; it can confirm any number of cases but never the family. The repetition family $[[n, 1, 1]]$ yields to the same treatment. Its checks are the neighboring $Z$-pairs, which generate the even-weight $Z$-strings, so a lone $Z_1$ commutes with all of them yet is not their product: an undetectable weight-one logical, for every $n$.

```wl
repChecks[n_] := StringJoin /@ (ReplacePart[ConstantArray["I", n], {# -> "Z", # + 1 -> "Z"}] & /@ Range[n - 1]);
distanceWitness[repChecks[7]]
```

A single undetectable $Z$, of weight one: the repetition code's quantum distance is one for the whole family, the surprise of Part III stated once for every $n$.

### A Continuous Error Angle

Symbolic generality reaches continuous parameters too. Consider a *coherent* error on the bit-flip codeword: not a discrete bit flip, but a rotation $R_x(\theta) = e^{-i\theta X/2}$ on the first qubit, through a symbolic angle. The state stays a superposition of "no error" and "flipped," and the syndrome operator $Z_1 Z_2$ has an expectation value that is an exact function of $\theta$:

```wl
psi0 = QuantumState[SparseArray[{1 -> 1}, 8], {2, 2, 2}];
vTheta = Normal[QuantumCircuitOperator[{"RX"[\[Theta]] -> 1}][psi0]["StateVector"]];
zzOp = KroneckerProduct[PauliMatrix[3], PauliMatrix[3], IdentityMatrix[2]];
zzExpectation = Simplify[Conjugate[vTheta] . zzOp . vTheta, \[Theta] \[Element] Reals]
```

The rotation is unitary, so the state it produces stays normalized at every angle, and the expectation above is read in a properly normalized state:

```wl
Simplify[Conjugate[vTheta] . vTheta, \[Theta] \[Element] Reals]
```

The syndrome expectation itself is $\cos\theta$: it interpolates continuously from $+1$ at $\theta = 0$ (no error, the codeword is untouched) to $-1$ at $\theta = \pi$ (a full bit flip, a definite nonzero syndrome), passing through zero at $\theta = \pi/2$ (a maximally uncertain syndrome). A coherent rotation is non-Clifford for generic $\theta$ (it happens to be Clifford at the special angles $0$, $\pi/2$, $\pi$ sampled just below), so at a symbolic angle it lives in the exact dense simulator, not in a stabilizer tableau; carrying the angle symbolically returns the whole error curve as a closed form, the syndrome as a function of the continuous angle rather than a value estimated at one fixed angle.

```wl
{zzExpectation /. \[Theta] -> 0, zzExpectation /. \[Theta] -> Pi/2, zzExpectation /. \[Theta] -> Pi}
```

The three sampled values confirm the endpoints and the midpoint of the closed form we just derived.

The endpoints are the whole story only at large angle. For a *weak* coherent error the quantity that matters is the probability that the syndrome fires, $(1 - \langle Z_1 Z_2\rangle)/2$, and its behavior as the rotation turns on is the leading term of its series in $\theta$:

```wl
Normal@Series[(1 - zzExpectation)/2, {\[Theta], 0, 2}]
```

To leading order the syndrome fires with probability $\theta^2/4$: the error rate is quadratic in the rotation angle, vanishing faster than linearly as $\theta \to 0$. That is the kind of input a threshold calculation consumes, delivered here as a closed form in the angle rather than a number at a fixed rate.

One qualification keeps the example honest. The observable here, the $Z_1 Z_2$ syndrome, sees only the flip *probability*: a stochastic bit-flip channel with $p = \sin^2(\theta/2)$ gives $\langle Z_1 Z_2\rangle = 1 - 2p = \cos\theta$, the identical curve. So this demonstrates a symbolic continuous-parameter calculation, not an observable that tells a coherent rotation apart from incoherent noise; distinguishing those would take an interference-sensitive measurement, such as reading the qubit in a rotated basis. What is genuinely beyond a Pauli sampler is the closed form itself, the syndrome returned as a function of $\theta$ rather than a statistic estimated at each fixed angle.

## Part V: Where This Leaves Us

Before we close, a plain summary of what has been established:

- For every one of the five codes, applying any shipped generator to the codeword leaves it fixed exactly, sign included; Stim accepts the same generators as a consistent stabilizer group; and the shipped syndrome equals Stim's across every weight-one and weight-two error.
- `PauliStabilizer`'s symbolic measurement carries a random outcome as a formal bit, and a forced outcome comes back as an explicit algebraic function of the free bits, which a sampler propagates as concrete frames but does not retain in closed form.
- An independent GF(2) search returns the exact distance and a minimum-weight witness for each code; the logical operators satisfy the full commutation algebra; and Stim confirms each witness is a valid logical, while the exhaustive search establishes it is the lightest.
- The $[[n, n-2, 2]]$ code condition is derived at symbolic $n$, its distance argued size-independently and confirmed at sample sizes; the repetition family's distance is one; and a coherent error's syndrome is returned as the exact function $\cos\theta$, with a small-angle syndrome rate of $\theta^2/4$.

The picture is a division of labor. Stim owns scale and speed: sampling large circuits, benchmarking decoders, estimating thresholds. This layer, checked here against Stim, gets the states, signs, and syndromes exactly right and drops into the rest of a quantum-mechanics toolkit; and the Wolfram Language around it holds a parameter symbolic where the algebra allows, deriving a family's code condition at symbolic $n$ and returning a coherent error's syndrome as a closed form in the angle. The two agree wherever they overlap, which is what makes the shared states and syndromes trustworthy; the distance's minimality rests on the exhaustive search, not on Stim, and the symbolic results on the algebra.

The boundary is worth naming plainly: the exhaustive distance search is a small-code instrument by construction, and exact minimum distance is NP-hard in the worst case, so this is not a scaling competitor. When you want the exact answer for a textbook code or a statement about a family, that exactness is what this buys; when you want a threshold for a distance-25 surface code, reach for Stim. Both live in the same notebook, and you can now check one against the other yourself. To push further, adapt the from-scratch functions to your own stabilizer generators (they assume a well-formed, independent set of checks) and diff against Stim, or carry the symbolic-family argument to a code family of your own.
