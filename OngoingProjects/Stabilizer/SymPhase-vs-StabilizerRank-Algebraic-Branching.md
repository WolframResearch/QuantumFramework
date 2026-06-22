---
Title: SymPhase versus stabilizer rank: when does "branching becomes algebra" transfer from measurement to the T gate?
Author: investigation in the Wolfram QuantumFramework
Description: A computation-first study of whether the SymPhase symbolic-sign trick (Fang and Ying, arXiv:2311.03906) that replaces measurement branching with an algebraic bit-polynomial transfers to non-Clifford (T-gate) branching in a StabilizerFrame, and where it fundamentally cannot. Every claim is verified with wolframscript against the working-tree paclet.
---

# SymPhase versus stabilizer rank: algebraic branching for measurement and for the T gate

## The question

The stabilizer tutorial `Stabilizer-Formalism-By-Computation.md` exhibits two kinds of branching and treats them very differently.

A **measurement** of an undetermined qubit forks the history into two equiprobable outcomes. The SymPhase technique of Fang and Ying replaces that fork with algebra: each unresolved random outcome is carried as a fresh formal bit $s_k \in \mathbb{F}_2$ written into a stabilizer sign, the circuit is traversed once, and all $2^r$ branches are served by substitution. The state stays a single `PauliStabilizer`.

A **non-Clifford $T$ gate** also branches: the framework turns $T\lvert\psi\rangle$ into a `StabilizerFrame`, a superposition of two stabilizer states, and $k$ successive $T$ gates give up to $2^k$ components, the *stabilizer rank*.

The natural question, asked by a reader who has just seen SymPhase defeat measurement branching: *can the same "branching becomes algebra" move be applied to the $T$-gate branching, carrying the magic as a symbolic phase instead of forking into $2^k$ components?*

This report answers it, with computation rather than recall. The short answer:

> **The trick transfers exactly once, and not where you would guess.** It does **not** fail for the reason one expects (the two $T$-components having different tableaux); in the framework's representation they share a tableau, just like a measurement's two branches. It fails for a deeper reason: a symbolic sign encodes a classical **disjunction** over a fixed tableau (one valid stabilizer state per assignment of the bits), whereas a $T$ gate produces a coherent **superposition** with complex weights that must be summed, a state that is *no* single stabilizer state. The one genuine transfer is to **diagonal** Clifford+$T$ circuits, where the magic really does become an algebraic object: a degree-$\le 3$ phase polynomial over $\mathbb{F}_2$ valued in $\mathbb{Z}_8$, the exact $\mathbb{Z}_8$ analogue of SymPhase's degree-$1$ $\mathbb{Z}_2$ sign polynomial, and the rank stays $1$. The moment a Hadamard is interleaved between two $T$ gates the phase polynomial can no longer stay diagonal, the frame branches as $2^k$, and the cost reappears as a $\#\mathrm{P}$-hard Gauss sum at contraction time. The exponential is **relocated, not removed**, which is exactly what Clifford+$T$ universality demands.

All experiments below were run against the working-tree paclet at

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

(paclet `Wolfram/QuantumFramework 2.0.0`). The driver scripts are in `symphase-vs-rank-scripts/` next to this file.

---

## 1. What makes SymPhase work: signs are a $\mathbb{Z}_2$, degree-$1$ algebra over a frozen tableau

SymPhase rests on two facts about the Aaronson-Gottesman tableau (Fang and Ying, arXiv:2311.03906, Facts 1 and 2).

**Fact 1.** The Pauli gates $X, Y, Z$ affect **only** the sign columns $\mathbf{R}, \bar{\mathbf{R}}$ of the tableau, never the symplectic $\mathbf{X}, \mathbf{Z}$ part. In the framework this is literally how the gate updates are written: `GateUpdates.m` implements $Z_j$ as "XOR the phase with the $X$-row, tableau unchanged."

**Fact 2.** The control flow of the A-G algorithm (which row anticommutes with a measurement, which row operations fire) depends only on $\mathbf{X}, \mathbf{Z}$, never on the signs. The signs only decide measurement *outcomes*.

Together these say: all $2^r$ branches of a sequence of measurements share **one** $\mathbf{X}, \mathbf{Z}$ tableau and differ only in the sign column, and that sign column is an $\mathbb{F}_2$-linear (XOR, degree-$1$) form in the fresh bits. We can watch both facts directly.

```wl
base = PauliStabilizer[3][{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 3}];
{base["Tableau"] === base["X", 1]["Tableau"],   (* X: tableau frozen *)
 base["Tableau"] === base["Z", 2]["Tableau"],   (* Z: tableau frozen *)
 base["Tableau"] === base["Y", 1]["Tableau"],   (* Y: tableau frozen *)
 base["Phase"]   =!= base["X", 1]["Phase"]}     (* but the phase moved *)
```

returns `{True, True, True, True}`: the three Pauli gates leave the symplectic tableau bit-identical and move only the phase.

The framework's `"SymbolicMeasure"` is the direct kernel realization (`Stabilizer/SymbolicMeasure.m`). A random (anticommuting) measurement allocates one fresh symbol `\[FormalS][k]` and writes it into a **single** sign position of the **shared** tableau:

```wl
sm = PauliStabilizer[2][{"H", "CNOT"}]["SymbolicMeasure", 1];
b0 = sm["SubstituteOutcomes", {\[FormalS][1] -> 0}];
b1 = sm["SubstituteOutcomes", {\[FormalS][1] -> 1}];
{sm["Tableau"] === b0["Tableau"] === b1["Tableau"],  (* one tableau *)
 b0["Signs"] =!= b1["Signs"]}                        (* two sign choices *)
```

returns `{True, True}`. The two outcome branches are the physical post-measurement states $\lvert 00\rangle$ and $\lvert 11\rangle$, one valid stabilizer state per value of the formal bit. Fact 2 in action: a state and its sign-flipped sibling measure to the **same** tableau update,

```wl
psA = PauliStabilizer[2][{"H", "CNOT"}];
psB = psA["X", 1];                                   (* same tableau, flipped signs *)
mA = psA["SymbolicMeasure", 1]; mB = psB["SymbolicMeasure", 1];
mA["Tableau"] === mB["Tableau"]                      (* True *)
```

Finally the outcome correlations come out as **degree-$1$** relations. For the parity gadget (two $\lvert+\rangle$ qubits whose parity is copied onto a third), measuring all three symbolically spends two free bits and leaves the conserved $ZZZ$ generator untouched:

```wl
parityCheck = PauliStabilizer[3][{"H" -> 1, "H" -> 2, "CNOT" -> {1, 3}, "CNOT" -> {2, 3}}];
smPar = parityCheck["SymbolicMeasure", {1, 2, 3}];
smPar["Phase"]
```

gives `{0, 0, 0, \[FormalS][i], \[FormalS][j], 0}`, two fresh bits in two signs with the third generator bare (the subscripts $i, j$ are a running session counter, so their exact values depend on how many symbols were already allocated). Tracing the four assignments of those two bits reproduces $m_3 = m_1 \oplus m_2$ throughout, i.e. $m_1 \oplus m_2 \oplus m_3 = 0$, a linear (degree-$1$) relation over $\mathbb{F}_2$. Crucially, substituting the bits yields **one basis state**, never a superposition:

```wl
free = DeleteDuplicates @ Cases[smPar["Phase"], _\[FormalS], Infinity];
oneShot = smPar["SubstituteOutcomes", Thread[free -> {1, 0}]];
Count[Chop[Normal[oneShot["State"]["StateVector"]]], x_ /; Abs[x] > 10^-9]   (* 1 *)
```

This is the structure that the question proposes to reuse: **a single frozen tableau whose sign column carries a low-degree $\mathbb{F}_2$ polynomial in formal bits, and a classical disjunction over those bits.**

---

## 2. The $T$ gate: same tableau, but a coherent complex sum, not a disjunction

The framework's `"T"` and `"P"[\[Theta]]` (`Stabilizer/GateUpdates.m`) turn a `PauliStabilizer` into a `StabilizerFrame` of two components,
$$
T\lvert s\rangle = \tfrac{1 + e^{i\pi/4}}{2}\,\lvert s\rangle \;+\; \tfrac{1 - e^{i\pi/4}}{2}\,Z_q\lvert s\rangle ,
$$
with the second component built by applying $Z_q$ to the first. Build $T\lvert+\rangle$ and inspect it:

```wl
frame = PauliStabilizer[1][{"H", "T"}];   (* T|+> *)
{frame["Length"], frame["Coefficients"]}
```

returns `{2, {(1 + E^(I/4 Pi))/2, (1 - E^(I/4 Pi))/2}}`: two components, with **complex** weights $(1 \pm e^{i\pi/4})/2$, not signs.

Here is the first surprise, and the place where the naive expectation is wrong. One imagines the two stabilizer components have *different* $\mathbf{X}, \mathbf{Z}$ tableaux, so that the single-frozen-tableau picture obviously cannot hold both. It is not so. Because the framework decomposes $T\lvert s\rangle$ using $Z_q\lvert s\rangle$, and $Z_q$ is a Pauli (Fact 1 again: it freezes the tableau), the two components **share a tableau** and differ only in a sign, exactly like a measurement's two branches:

```wl
comps = frame["Stabilizers"];
{comps[[1]]["Tableau"] === comps[[2]]["Tableau"],     (* True: one shared tableau *)
 comps[[1]]["Signs"]   =!= comps[[2]]["Signs"],       (* True: only signs differ *)
 {comps[[1]]["Stabilizers"], comps[[2]]["Stabilizers"]}}   (* {+X} vs {-X} *)
```

returns `{True, True, {{X}, {-X}}}`. Structurally, $T\lvert+\rangle$ looks *identical* to a symbolic measurement: one tableau, two sign patterns. This is no accident of $\lvert+\rangle$; the framework maintains exactly this invariant for every gate-built frame through its "relating Paulis" bookkeeping (`StabilizerFrame.m`), so **all** $2^k$ components of a $k$-$T$ frame share the reference tableau and differ only by relating Paulis, which are sign-like.

So the obstruction is **not** "different tableaux." The obstruction is the operation that combines the branches. For a measurement the two branches are a classical **disjunction**: pick one, weighted by probability, and the result is a physical state. For the $T$ gate the two components are a coherent **superposition**: they must be **added** with their complex weights, and the result is a state that is *no* single stabilizer state.

```wl
Count[Chop[Normal[frame["StateVector"]]], x_ /; Abs[x] > 10^-9]   (* 2: a superposition *)
```

The coherent sum $\tfrac{1+e^{i\pi/4}}{2}\lvert+\rangle + \tfrac{1-e^{i\pi/4}}{2}\lvert-\rangle$ equals $T\lvert+\rangle = (\lvert0\rangle + e^{i\pi/4}\lvert1\rangle)/\sqrt2$ up to a global phase,

```wl
dense = QuantumCircuitOperator[{"H" -> 1, "T" -> 1}][QuantumState["0"], Method -> "Schrodinger"];
QuantumSimilarity[frame["State"], dense]   (* 1 *)
```

but it is **neither** component. No assignment of a sign reaches it, because $T\lvert+\rangle$ is not a stabilizer state at all: no signed or phased Pauli stabilizes it.

```wl
vT = Chop[Normal[frame["StateVector"]]];
Or @@ Flatten @ Table[Chop[(s P) . vT - vT] === {0, 0},
   {P, {PauliMatrix[1], PauliMatrix[2], PauliMatrix[3]}}, {s, {1, -1, I, -I}}]   (* False *)
```

returns `False`: $T\lvert+\rangle$ has stabilizer rank $\ge 2$. This is the precise content of the disanalogy. A `PauliStabilizer`, with concrete *or* symbolic signs, is **always exactly one** stabilizer state for each assignment of its formal bits (the stabilizer group has $n$ generators, which pins the state uniquely up to global phase). A symbolic sign therefore ranges over a set of stabilizer states; it can express "this branch **or** that branch," never "this branch **plus** that branch."

| | measurement / Pauli fault | $T$ gate |
|---|---|---|
| tableau | one, frozen | one, frozen (relating Paulis) |
| sign carrier | $s_k \in \mathbb{F}_2$ in the sign column | complex weight $(1\pm e^{i\pi/4})/2$ |
| combine branches by | disjunction (OR, probabilities) | superposition (PLUS, complex amplitudes) |
| substitute a bit gives | one physical stabilizer state | one *summand*, not the state |
| is the result a stabilizer state? | yes, each branch | no (rank $\ge 2$) |

### 2a. Why a $\mathbb{Z}_8$ symbolic phase on one tableau cannot rescue it

One might hope to upgrade the sign column from $\mathbb{Z}_2$ ($\pm1$) to $\mathbb{Z}_8$ (eighth roots of unity) and absorb the $e^{i\pi/4}$ into a single tableau. This fails at the most elementary level: a sign $\pm1$ keeps the generator an **involution**, $(\pm P)^2 = I$, so $\pm P$ has a $+1$ eigenspace and defines a stabilizer state. A primitive eighth-root phase does not:

```wl
zeta = Exp[I Pi/4];
{Chop[(zeta PauliMatrix[1]) . (zeta PauliMatrix[1])],   (* {{I,0},{0,I}} = i I, NOT identity *)
 Chop[Eigenvalues[zeta PauliMatrix[1]]]}                (* {-e^{i pi/4}, e^{i pi/4}}, no +1 *)
```

$\big(e^{i\pi/4}X\big)^2 = iI \ne I$, and its eigenvalues are $\pm e^{i\pi/4}$, neither equal to $1$. There is **no** state stabilized by $e^{i\pi/4}X$. The magic eighth root simply cannot live in a sign: it breaks the involution condition that makes a tableau a state. It can only live as a relative amplitude in a superposition (the `StabilizerFrame`), or, in one special regime, as a phase polynomial on a computational-basis support, which is the next section.

---

## 3. The one genuine transfer: diagonal Clifford+$T$ is a $\mathbb{Z}_8$ phase polynomial, rank $1$

There is a regime where the reader's intuition is exactly right and the magic does become algebra. Restrict to **diagonal** gates: $T$, $S$, $Z$, $CZ$, $CCZ$, that is, Clifford+$T$ with no Hadamard. Acting on $\lvert+\rangle^{\otimes n}$, a diagonal gate multiplies each computational-basis amplitude by a phase, and the phases compose into a single **phase polynomial** $\varphi: \mathbb{F}_2^n \to \mathbb{Z}_8$:
$$
U\lvert+\rangle^{\otimes n} \;=\; \frac{1}{\sqrt{2^n}}\sum_{x \in \mathbb{F}_2^n} \omega^{\varphi(x)}\,\lvert x\rangle, \qquad \omega = e^{i\pi/4},
$$
where $\varphi$ has **degree at most $3$** over $\mathbb{F}_2$ (the CNOT-dihedral normal form: $T_j$ contributes a linear term $x_j$, $CZ_{ij}$ a quadratic term $4 x_i x_j$, $CCZ_{ijk}$ a cubic term $4 x_i x_j x_k$). Take $U = T_1\,T_2\,CZ_{12}$, for which $\varphi(x) = x_1 + x_2 + 4 x_1 x_2 \pmod 8$:

```wl
w = Exp[I Pi/4];
denseDiag = QuantumCircuitOperator[{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}][
   QuantumState["00"], Method -> "Schrodinger"];
phi[x1_, x2_] := Mod[x1 + x2 + 4 x1 x2, 8];
vPoly = Flatten @ Table[(1/2) w^phi[x1, x2], {x1, 0, 1}, {x2, 0, 1}];
Max[Abs[N[denseDiag["StateVector"] - vPoly]]] < 10^-12    (* True *)
```

The single phase polynomial reproduces the exact state. The cubic term is reached by a $CCZ$ on $\lvert+\rangle^{\otimes 3}$, $\varphi(x) = 4 x_1 x_2 x_3$:

```wl
ccz = DiagonalMatrix[{1, 1, 1, 1, 1, 1, 1, -1}];
vCCZ = Flatten @ Table[(1/Sqrt[8]) w^Mod[4 x1 x2 x3, 8], {x1, 0, 1}, {x2, 0, 1}, {x3, 0, 1}];
Max[Abs[N[ccz . (ConstantArray[1, 8]/Sqrt[8]) - vCCZ]]] < 10^-12   (* True *)
```

This is the precise $\mathbb{Z}_8$ analogue of SymPhase. SymPhase carries a degree-$1$ polynomial over $\mathbb{F}_2$ in the **sign** column ($\mathbb{Z}_2$-valued); a diagonal Clifford+$T$ circuit carries a degree-$\le 3$ polynomial over $\mathbb{F}_2$ in the **phase** ($\mathbb{Z}_8$-valued). Same move, one ring up and two degrees up:

| | measurement / Pauli (SymPhase) | diagonal Clifford+$T$ |
|---|---|---|
| carrier | sign column $\mathbf{R}$ | phase of each amplitude |
| ring | $\mathbb{Z}_2$ (bits $\pm1$) | $\mathbb{Z}_8$ (eighth roots) |
| polynomial degree over $\mathbb{F}_2$ | $1$ (affine XOR) | $\le 3$ (CNOT-dihedral) |
| support | the stabilizer's affine subspace | unchanged affine subspace |
| branching | none: one tableau | none: one phase-polynomial state |
| stays rank $1$ | yes | yes |

In this regime the magic is genuinely "absorbed into phases, no branching," exactly as the question hoped. The state is still non-stabilizer (the eighth roots are real magic), but it is stored as **one** algebraic object, a support plus a low-degree polynomial, rather than as a tree.

**Honest caveat about the framework.** The framework's `StabilizerFrame` does **not** exploit this. Its `"T"` update unconditionally doubles the component count, even when the circuit is diagonal and the true object is rank $1$:

```wl
qfFrame = PauliStabilizer[2][{"H" -> 1, "H" -> 2, "T" -> 1, "T" -> 2, "CZ" -> {1, 2}}];
{qfFrame["Length"], QuantumSimilarity[qfFrame["State"], denseDiag]}   (* {4, 1} *)
```

The frame is **exact** (similarity $1$) but carries $4$ components for a state a phase polynomial describes in one. So the "rank stays $1$" statement is a property of the phase-polynomial representation, an algorithm the framework could add, not of the current `StabilizerFrame`. A phase-polynomial / CNOT-dihedral simulator is the right tool for the diagonal fragment, and it is a clean potential extension.

---

## 4. The obstruction: a Hadamard between two $T$ gates, and the exponential reappears as a $\#\mathrm{P}$-hard sum

The diagonal fragment is not universal and not even closed under Clifford: it has no Hadamard. As soon as an $H$ is interleaved between $T$ gates the phase polynomial cannot stay diagonal, because $H$ Fourier-transforms the computational basis and the next $T$ then acts in the dual basis. The framework's frame branches once per $T$, to $2^k$:

```wl
stabRank[x_] := If[MatchQ[x, _StabilizerFrame], x["Length"], 1];
interleaved[k_] := PauliStabilizer[1][Prepend[Flatten[Table[{"T", "H"}, k]], "H"]];
Table[{k, stabRank[interleaved[k]]}, {k, 0, 8}]
```

returns $\{1, 2, 4, 8, 16, 32, 64, 128, 256\}$: exactly $2^k$.

What happened to the algebra? It did not vanish; it moved into the **contraction**. An interior Hadamard introduces a summed $\mathbb{F}_2$ variable: $H\lvert x\rangle = \tfrac{1}{\sqrt2}\sum_a (-1)^{ax}\lvert a\rangle$. So the amplitude of $H\,T\,H\lvert0\rangle$ is an internal Gauss sum,
$$
\langle x \lvert H T H \rvert 0\rangle = \frac{1}{2}\sum_{a \in \mathbb{F}_2} \omega^{a}\,(-1)^{ax},
$$
which we can confirm term by term:

```wl
ampHTH[x_] := (1/2) Sum[w^a (-1)^(a x), {a, 0, 1}];
denseHTH = QuantumCircuitOperator[{"H" -> 1, "T" -> 1, "H" -> 1}][QuantumState["0"], Method -> "Schrodinger"];
Max[Abs[N[denseHTH["StateVector"] - {ampHTH[0], ampHTH[1]}]]] < 10^-12   (* True *)
```

Each interior Hadamard doubles the number of summed terms, so a circuit with $d$ interior $H$-layers carries a $2^d$-term sum per amplitude. The framework's frame is the state-vector face of the same exponential: materializing `frame["StateVector"]` sums all $2^k$ components (`StabilizerFrame.m` `Total`s the relating-Pauli contributions). The cost is exponential in the $T$-count, and we can watch it, clearing the system cache each row:

```wl
benchRow[k_] := Module[{f = interleaved[k]},
  {k, stabRank[f], ByteCount[f], First @ AbsoluteTiming[f["StateVector"];]}];
Table[ClearSystemCache[]; benchRow[k], {k, 2, 12, 2}]
```

| $k$ ($T$-count) | rank $=2^k$ | `ByteCount` | materialize (ms) |
|---|---|---|---|
| 2 | 4 | 9 128 | 2 |
| 4 | 16 | 37 352 | 2 |
| 6 | 64 | 151 368 | 5 |
| 8 | 256 | 605 560 | 17 |
| 10 | 1024 | 2 426 200 | 62 |
| 12 | 4096 | 9 709 624 | 257 |

(`ByteCount` is exact; the times are representative single-run figures on an M2 Pro and vary run to run.) `ByteCount` and time both scale as $2^k$. This is not a defect to be optimized away: it is forced. An amplitude of a general Clifford+$T$ circuit is a normalized exponential (Gauss) sum over $\mathbb{F}_2$ of a degree-$\le 3$ phase polynomial mod $8$, and evaluating such sums exactly is $\#\mathrm{P}$-hard in general. Equivalently, exact strong simulation of a universal gate set (and Clifford+$T$ is universal) would place $\mathrm{BQP}$ inside $\mathrm{P}$ and collapse the polynomial hierarchy (via $\mathrm{PostBQP} = \mathrm{PP}$). No Gottesman-Knill-style polynomial extension can exist past Clifford; the magic that the phase polynomial absorbed for diagonal circuits becomes, once Hadamards interleave $T$ gates, a sum no algebra can collapse. The exponential is **relocated** from the branching tree into the Gauss sum, never deleted.

This also reframes the disjunction-versus-superposition point of Section 2 at the level of complexity. SymPhase samples a probability distribution: substitute random bits, read off one branch, and the $2^r$ tree is a red herring because the branches are mutually exclusive classical events. A $T$-rich amplitude is a coherent sum of $2^k$ interfering terms, and interference is precisely what a single forward-and-substitute pass cannot capture: there is no "one branch" to sample; every term contributes to every amplitude at once.

---

## 5. The partial wins, and what they cost

The verdict is not "nothing transfers." Three genuine footholds remain, all paid for in the rank.

**Diagonal fragment (Section 3).** For Clifford+$T$ with no $H$ between $T$ gates, the phase-polynomial representation is rank $1$ and exact: storage and single-amplitude readout are polynomial. This is the true algebraic analogue of SymPhase and is a clean candidate feature for the framework.

**Bounded $T$-count / bounded rank.** When the magic is small the frame is exactly the right object: faithful and cheap. The cost is the rank, and the framework keeps every term, $2^t$, with no compression toward the true stabilizer rank (which is far smaller: $\chi(\lvert T\rangle^{\otimes t}) = t$ for $2 \le t \le 4$ and $\le t+1$ for $t \le 5$, Bravyi et al., arXiv:1808.00128; asymptotically $\chi_t \le 2^{0.47 t}$ and approximate rank $2^{\gamma t}$ with $\gamma = -2\log_2\cos(\pi/8) \approx 0.23$, Bravyi and Gosset, arXiv:1601.07601):

```wl
exactRankT = <|0 -> 1, 1 -> 1, 2 -> 2, 3 -> 3, 4 -> 4, 5 -> 5|>;
Table[{t, 2^t, stabRank[PauliStabilizer[1][Join[{"H"}, Table["T", t]]]], exactRankT[t]}, {t, 0, 5}]
```

| $t$ | $2^t$ naive | QF frame `Length` | true exact rank $\chi$ |
|---|---|---|---|
| 0 | 1 | 1 | 1 |
| 1 | 2 | 2 | 1 |
| 2 | 4 | 4 | 2 |
| 3 | 8 | 8 | 3 |
| 4 | 16 | 16 | 4 |
| 5 | 32 | 32 | 5 |

The framework stores $2^t$ where $t$ would do (or $\approx 2^{0.4 t}$ for generic magic): exact and faithful, but not a scalable magic-state simulator. Compressing to the true rank is the Bravyi-Gosset program (sparsification, CH-form, stabilizer extent, arXiv:1808.00128), which the frame does not yet implement.

**Local-observable expectations.** A single expectation $\langle\psi\lvert P\rvert\psi\rangle$ on a rank-$\chi$ frame is a sum of $\chi^2$ stabilizer inner products, each $O(n)$, so it is polynomial in $n$ for fixed rank, no full $2^n$ contraction:

```wl
fr = PauliStabilizer[1][{"H", "T"}];   (* rank 2 *)
{fr["Length"], fr["Length"]^2,
 N @ Re[Conjugate[#] . (PauliMatrix[2] . #)] &[Normal[fr["StateVector"]]]}
```

returns `{2, 4, 0.7071...}`: rank $\chi = 2$, $\chi^2 = 4$ inner-product terms, and the expected $\langle Y\rangle = +1/\sqrt2$. For a *local* observable one need only keep the $T$ gates inside its backward light cone, so the relevant rank is $2^{(\text{cone }T\text{-count})}$, often far below $2^t$, the same light-cone economy the tensor-network engine exploits for Clifford-dominated circuits. And for **sampling** (weak simulation) rather than exact amplitudes, the approximate-rank methods of Bravyi and Gosset run in $2^{0.23 t}$, a genuine improvement the frame could host but does not.

---

## 6. Verdict

**(a) Why symbolic signs transfer for measurement and Pauli but not for the generic $T$.** SymPhase works because a measurement (and a Pauli fault) is a classical **disjunction** over a single frozen $\mathbf{X}, \mathbf{Z}$ tableau: Facts 1 and 2 guarantee that all branches share the tableau and differ only in an $\mathbb{F}_2$-linear sign, and a `PauliStabilizer` with symbolic signs is exactly one valid stabilizer state per bit-assignment. The $T$ gate produces a coherent **superposition** with complex weights $(1\pm e^{i\pi/4})/2$; the result is no stabilizer state (rank $\ge 2$). The framework's two $T$-components, perhaps unexpectedly, **also** share a tableau, so the obstruction is not "different tableaux" but "PLUS versus OR": a symbolic sign can encode "this branch or that," never "this branch plus that," and a fractional ($\mathbb{Z}_8$) phase cannot even sit in a sign because it breaks the involution $(\pm P)^2 = I$ that makes a tableau a state.

**(b) The one genuine transfer.** For **diagonal** Clifford+$T$ (CNOT-dihedral: $T$, $S$, $CZ$, $CCZ$, no $H$) the magic is exactly an algebraic object: a degree-$\le 3$ phase polynomial over $\mathbb{F}_2$ valued in $\mathbb{Z}_8$, the $\mathbb{Z}_8$ analogue of SymPhase's degree-$1$ $\mathbb{Z}_2$ sign polynomial, with the rank staying $1$. This is verified exactly for $T_1 T_2 CZ_{12}$ and for the cubic $CCZ$. The framework does not currently use it (its frame branches anyway), so this is a real, clean extension opportunity, not a present capability.

**(c) The fundamental obstruction.** Interleaving a Hadamard between two $T$ gates forbids the diagonal form: the frame branches as $2^k$, and extracting an amplitude or probability becomes an exponential Gauss sum over $\mathbb{F}_2$ that is $\#\mathrm{P}$-hard to evaluate exactly. The exponential is relocated from the branching tree into the contraction, not removed. This is mandatory: Clifford+$T$ is universal, so any polynomial Gottesman-Knill-style extension that captured generic $T$ branching algebraically would collapse the polynomial hierarchy.

**(d) The partial wins and their price.** The diagonal fragment (rank $1$), bounded-$T$ / bounded-rank circuits (cost $=$ rank, but the framework keeps $2^t$ rather than the true $\approx 2^{0.4 t}$), local-observable expectations ($\chi^2$ inner products, light-cone-limited), and approximate sampling ($2^{0.23 t}$) are all real footholds, each paying the exponential in the one currency that cannot be avoided: the number of stabilizer terms.

In one sentence: SymPhase moves measurement randomness from a tree into a degree-$1$ $\mathbb{F}_2$ sign polynomial because the branches are a classical disjunction over a frozen tableau; the same move reaches non-Clifford magic only for diagonal circuits, as a degree-$\le 3$ $\mathbb{Z}_8$ phase polynomial, and stops dead at the first Hadamard-between-$T$s, where coherent interference turns the would-be algebra into a $\#\mathrm{P}$-hard sum that Clifford+$T$ universality guarantees no algebra can collapse.

---

## References

- Y. Fang and M. Ying, *SymPhase: Phase Symbolization for Fast Simulation of Stabilizer Circuits*, arXiv:2311.03906 (2023). Facts 1 and 2; symbolic $\mathbb{F}_2$ sign polynomials; sampling as bit-matrix multiplication.
- S. Bravyi and D. Gosset, *Improved classical simulation of quantum circuits dominated by Clifford gates*, Phys. Rev. Lett. 116, 250501 (2016), arXiv:1601.07601. Stabilizer rank $\chi_t$, runtime $2^{0.23 t}$, $\gamma = -2\log_2\cos(\pi/8)$.
- S. Bravyi, D. Browne, P. Calpin, E. Campbell, D. Gosset, M. Howard, *Simulation of quantum circuits by low-rank stabilizer decompositions*, Quantum 3, 181 (2019), arXiv:1808.00128. Stabilizer extent, CH-form, exact and approximate rank bounds.
- M. Howard and E. Campbell, *Application of a resource theory for magic states to fault-tolerant quantum computing*, Phys. Rev. Lett. 118, 090501 (2017), arXiv:1609.07488. Robustness of magic; $\lvert H\rangle = (\lvert0\rangle + e^{i\pi/4}\lvert1\rangle)/\sqrt2$; quasiprobability / sum-over-Cliffords.
- The CNOT-dihedral / phase-polynomial normal form (diagonal Clifford+$T$ as a degree-$\le 3$ polynomial mod $8$): M. Amy, D. Maslov, M. Mosca, and related work; A. Montanaro, *Quantum circuits and low-degree polynomials over $\mathbb{F}_2$*, arXiv:1607.08473 (2016).
- H. J. García and I. L. Markov, *Simulation of Quantum Circuits via Stabilizer Frames*, arXiv:1712.03554 (2015). The `StabilizerFrame` construction the framework follows.
- S. Aaronson and D. Gottesman, *Improved simulation of stabilizer circuits*, Phys. Rev. A 70, 052328 (2004), arXiv:quant-ph/0406196. The tableau and its $O(n^2)$ measurement update.

Kernel anchors: `QuantumFramework/Kernel/Stabilizer/SymbolicMeasure.m` (`freshOutcome`, the `"Outcomes"` map, `"SymbolicMeasure"` / `"SubstituteOutcomes"` / `"SampleOutcomes"`), `Stabilizer/StabilizerFrame.m` and `Stabilizer/GateUpdates.m` (the `"P"[\[Theta]]` / `"T"` frame doubling and relating Paulis). Verified against the working tree, paclet `Wolfram/QuantumFramework 2.0.0`. Driver scripts: `symphase-vs-rank-scripts/exp1-facts-and-disanalogy.wls`, `exp2-phasepoly-and-blowup.wls`, `exp3-analogy-table.wls`.
