# The non-Clifford boundary: how QF's `StabilizerFrame` blows up on T / arbitrary rotations

**What this is.** A verified, run-on-the-live-kernel study of the claim that the stabilizer
formalism is "Clifford only," and that adding non-Clifford gates to an otherwise-Clifford circuit
makes QuantumFramework's (QF) `StabilizerFrame` representation blow up. The headline result is
confirmed exactly: the moment a $T$ gate (or any single-qubit phase rotation $P[\theta]$) lands on a
stabilizer state, QF replaces the $O(n^2)$ tableau with a `StabilizerFrame` whose **component count
doubles per non-Clifford gate**, giving a literal $2^k$ for $k$ such gates. It is fine for a handful
of $T$ gates and exponential for $T$-rich circuits. Along the way two things that the brief did not
anticipate were found and are reported honestly: (1) arbitrary rotations `RY[theta]`/`RZ[theta]` are
not represented at all (they return `$Failed`, not a frame), and (2) the frame's dense
materialization *was* numerically wrong for genuine superpositions, with the one test guarding it too
weak to notice. Defect (2) is now **fixed in the kernel** (Section 4): a gate-built frame carries a
relating Pauli per component, so the dense readout materializes one reference component coherently and
obtains the rest by the actual Pauli operators that relate them, recovering the relative phases. The
fix is verified amplitude-exact against the dense `Method -> "Schrodinger"` state across 216 random
Clifford+$T$ circuits, and the guarding test is strengthened to a `Max`-based amplitude assertion that
fails on the old code and passes on the new.

**Companion.** This extends the stabilizer caveat in
[`../Measurement Dilation/Computing-Probabilities-And-Expectations.md`](../Measurement%20Dilation/Computing-Probabilities-And-Expectations.md)
(Section 4), which says the stabilizer route is the win for Clifford + Pauli measurements but "cannot
represent the `RY[theta]` benchmark circuit." This document is the empirical proof of that sentence,
plus the $T$-count scaling law, the mechanism, and the materialization caveat. It is the non-Clifford
companion to the benchmark-backed
[`../Platform Comparison/QF-Stabilizer-vs-Packages.md`](../Platform%20Comparison/QF-Stabilizer-vs-Packages.md)
(which covers the Clifford regime) and the tutorial
[`Stabilizer-Formalism-By-Computation.md`](Stabilizer-Formalism-By-Computation.md).

**Provenance.** Every QF claim was run on the live working-tree paclet (loaded with
`PacletDirectoryLoad`, not the lagging installed copy) and every number below was captured from that
run. Kernel anchor: live repo HEAD `a601a187`. The non-Clifford blow-up findings (Sections 1-3) were
first run at `f370f568`; the Section 4 materialization fix landed in `61bc0e39`
("Fix StabilizerFrame dense materialization to be phase-coherent"), which modified
`Stabilizer/StabilizerFrame.m`, `Stabilizer/GateUpdates.m`, and `Tests/Stabilizer/Formulas.wlt`. All
`file:line` cites below were re-verified against the live kernel at this HEAD: the `"P"`/distribute/
materialization handlers in those two files shifted with the fix (so Section 5 cites them by rule name
rather than line), while unchanged files (`PauliStabilizer.m`, `Conversions.m`, `InnerProduct.m`) keep
their line numbers. Machine: Apple M2 Pro (Darwin 23.2.0), Wolfram Language 15.0.0. `ClearSystemCache[]`
before every timing; no `Quiet`. Repro scripts: [`scripts/`](scripts/) (one per experiment, echoed
inline in the appendix).

---

## 0. The one-paragraph answer

A Clifford circuit applied with `Method -> "Stabilizer"` stays a single `PauliStabilizer` (an
Aaronson-Gottesman tableau, $O(n^2)$ storage). Insert one $T$ or one $P[\theta]$ and QF returns a
`StabilizerFrame`: a closure object $\sum_i c_i |s_i\rangle$ that stores a list of
$\{\text{coefficient}, \text{PauliStabilizer}\}$ components. A non-Clifford diagonal gate is a
two-term combination $P[\theta] = \tfrac{1+c}{2}\,\mathbb{1} + \tfrac{1-c}{2}\,Z$, so each gate maps
every component to **two**, and the implementation is naive: it Flattens the doubling with no merging,
deduplication, or rank reduction. The measured law is therefore exactly $2^k$ components for $k$
non-Clifford gates, confirmed bit-for-bit out to $k=18$ ($262144$ components, $596$ MB, $8.2$ s).
Arbitrary non-diagonal rotations `RY[theta]` (and even the diagonal `RZ[theta]`) are not handled at
all: they reach the tableau as a generic $R$ rotation, miss every gate rule, and return `$Failed`
with one `PauliStabilizer::nonclifford` message. This is exactly why the measurement-dilation report's
`RY[theta]` benchmark cannot use the stabilizer route. One caution remains: QF's `"P"[\theta]` frame
handler uses a half-angle convention that disagrees with the gate-level `QuantumOperator["P"[\theta]]`.
A second defect, the frame's dense `["StateVector"]` being physically wrong for any frame with $\ge 2$
distinct components, has been **fixed in the kernel** (Section 4): each gate-built component now carries
the Pauli that relates it to a single reference component, so the dense sum is phase-coherent.

---

## 1. The dichotomy: Clifford stays a tableau, one non-Clifford gate changes the head

Circuit backbone (pure Clifford): `{"H"->1, "CNOT"->{1,2}, "S"->2, "H"->3, "CNOT"->{2,3}}`. Applied
via `Method -> "Stabilizer"` it returns a single `PauliStabilizer` with `StabilizerStateQ -> True`.
Now append one non-Clifford gate (generic numeric angle $\theta = 0.7$ where applicable) and record
what comes back:

| Appended gate | Result head | Components | Message |
|---|---|---:|---|
| (none, Clifford only) | `PauliStabilizer` | 1 | none; `StabilizerStateQ` $\to$ `True` |
| `"T" -> 1` | `StabilizerFrame` | 2 | none |
| `"P"[0.7] -> 1` | `StabilizerFrame` | 2 | none |
| `"RZ"[0.7] -> 1` | `$Failed` | n/a | `PauliStabilizer::nonclifford` on `R[0.7, {Z -> {1}}]` |
| `"RY"[0.7] -> 1` | `$Failed` | n/a | `PauliStabilizer::nonclifford` on `R[0.7, {Y -> {1}}]` |

So the dichotomy is real but finer than "Clifford vs not":

- **Diagonal phase gates `"T"`, `"T"`$^\dagger$, `"P"[\theta]`** are lifted to a `StabilizerFrame`.
  The $T$ frame's two components are honest `PauliStabilizer`s with coefficients exactly
  $\left\{\tfrac{1+e^{i\pi/4}}{2},\ \tfrac{1-e^{i\pi/4}}{2}\right\}$ (verified to symbolic zero).
- **Everything else non-Clifford, including the diagonal `"RZ"[\theta]`, returns `$Failed`.** QF's
  `QuantumShortcut` lowers `"RZ"[\theta]` and `"RY"[\theta]` to a generic rotation
  `R[\theta, {Z|Y -> {q}}]`, which matches no tableau gate rule, so the dispatcher emits one
  `PauliStabilizer::nonclifford` and aborts. The stabilizer route never sees that `RZ` equals
  $P[\theta]$ up to a global phase, so it does not represent it even though the frame structurally
  could.

**Gate-convention caveat (verified).** The frame's `"P"[\theta]` handler is *not* the same gate as
QF's `QuantumOperator["P"[\theta]]`:

$$\texttt{QuantumOperator["P"[}\theta\texttt{]]} = \mathrm{diag}(1, e^{i\theta}), \qquad
\text{frame } \texttt{"P"[}\theta\texttt{]} \to \mathrm{diag}(1, e^{i\theta/2}) = P[\theta/2].$$

It uses $c = e^{i\theta/2}$ internally (Section 4a). The two agree for `"T"` (both
$\mathrm{diag}(1, e^{i\pi/4})$, because `"T"` hard-routes to `"P"[\pi/2]`), but a generic `"P"[\theta]`
gate applied through the stabilizer path realizes the *half-angle* gate. Use `"T"`/`"T"`$^\dagger$ for
unambiguous results, or double the angle if you intend $P[\theta]$.

---

## 2. Headline scaling: components $= 2^k$ exactly

Fix a 4-qubit Clifford backbone (Hadamard wall + CNOT ladder, 7 gates), insert $k$ $T$ gates cycling
over the four qubits, apply with `Method -> "Stabilizer"`, and measure the `StabilizerFrame`
component count (`["Length"]`), wall time, and `ByteCount`:

| $k$ | head | components | $2^k$ | $=2^k$ | $t$ (ms) | MB |
|---:|---|---:|---:|:--:|---:|---:|
| 0 | `PauliStabilizer` | 1 | 1 | yes | 406* | 0.00 |
| 1 | `StabilizerFrame` | 2 | 2 | yes | 12.6 | 0.00 |
| 4 | `StabilizerFrame` | 16 | 16 | yes | 13.0 | 0.04 |
| 8 | `StabilizerFrame` | 256 | 256 | yes | 21.6 | 0.58 |
| 10 | `StabilizerFrame` | 1024 | 1024 | yes | 47.7 | 2.33 |
| 12 | `StabilizerFrame` | 4096 | 4096 | yes | 145.6 | 9.31 |
| 14 | `StabilizerFrame` | 16384 | 16384 | yes | 530.1 | 37.25 |
| 16 | `StabilizerFrame` | 65536 | 65536 | yes | 2073.9 | 149.0 |
| 17 | `StabilizerFrame` | 131072 | 131072 | yes | 4158.3 | 298.0 |
| 18 | `StabilizerFrame` | 262144 | 262144 | yes | 8204.4 | 596.0 |

(*$k=0$ is the all-Clifford compiled-fold fast path plus first-call warmup; for $k \ge 1$ a single
non-Clifford gate disables that path and forces the per-gate object fold.)

(Storage/timing note. The numbers above predate the Section 4 materialization fix. Carrying one
relating Pauli per component costs about $+15\%$ `ByteCount` and $+8$-$12\%$ build time, measured at
the same $4$-qubit backbone against the fixed kernel: $k=12$ goes $9.31 \to 10.72$ MB and
$145.6 \to 156.9$ ms, $k=14$ goes $37.25 \to 42.87$ MB and $530.1 \to 592.2$ ms. The $2^k$ component
law and the exponential scaling are unchanged; the relating Pauli is a fixed $O(n)$ per-component
overhead.)

Readings:

- **The law is $2^k$ exactly, with no exceptions.** Every row satisfies `components === 2^k`; QF does
  no merging or deduplication (Section 3 makes this dramatic).
- **Both memory and time double per gate past $k \approx 10$.** Storage is flat per component
  ($\approx 2.4$ KB each: $596\ \text{MB}/262144$), so the curve is the pure exponential in component
  count. The exponential is independent of $n$ (the backbone size sets the per-component tableau cost,
  not the growth rate).
- The crossover from "negligible" to "painful" sits near $k \approx 14$-$16$ on this machine: $k=16$
  already costs $2$ s and $149$ MB, $k=18$ costs $8$ s and $0.6$ GB. A few $T$ gates: free. A
  $T$-rich circuit: hopeless.

---

## 3. The implementation is naive: $2^k$ even when the true state is rank 1

The doubling is unconditional. It does not look at whether the resulting components are distinct, nor
whether the state is still a single stabilizer state. Three demonstrations:

**(a) Arbitrary rotation `P[theta]` blows up identically to `T`.** Replacing the $k$ $T$ gates with
$k$ generic `"P"[0.7]` gates gives the same $2^k$ ladder ($1, 2, 4, 8, \ldots, 256$ for
$k = 0\ldots 8$). This is the direct tie-back to the measurement-dilation report: the "arbitrary
single-qubit rotation" that the frame *does* represent (a diagonal phase) suffers the identical
blow-up; the non-diagonal `RY[theta]` it cannot represent at all (Section 1).

**(b1) $k$ $T$ gates on $|0\rangle$: $2^k$ identical components for a rank-1 state.** Since
$T|0\rangle = |0\rangle$, the true state after any number of $T$ gates on $|0\rangle$ is just
$|0\rangle$ (a single stabilizer state). QF nonetheless produces $2^k$ components, and they are all
the **same** tableau:

| $k$ | components | distinct tableaux | materialized state $= |0\rangle$ |
|---:|---:|---:|:--:|
| 1 | 2 | 1 | yes |
| 2 | 4 | 1 | yes |
| 4 | 16 | 1 | yes |
| 8 | 256 | 1 | yes |

256 stored copies of one tableau for a state that is one basis vector. No rank reduction whatsoever.

**(b2) Period-8 return: $T^8|+\rangle = |+\rangle$ with 256 components.** $T^8 = \mathbb{1}$, so the
true state cycles back to $|+\rangle$ every 8 gates, yet the component count keeps doubling. The
overlap $|\langle + | \psi \rangle|^2$ tracks the genuine $T^k|+\rangle$ physics exactly while the
storage explodes:

| $k$ | components | $\lvert\langle+|\psi\rangle\rvert^2$ | exact $T^k|+\rangle$ |
|---:|---:|---:|---|
| 1 | 2 | 0.853553 | $\cos^2(\pi/8)$ |
| 2 | 4 | 0.500000 | $S|+\rangle$ |
| 4 | 16 | 0.000000 | $T^4 = Z,\ Z|+\rangle = |-\rangle$ |
| 7 | 128 | 0.853553 | $T^7|+\rangle$ |
| 8 | 256 | 1.000000 | $T^8 = \mathbb{1},\ |+\rangle$ again |

The physics is right; the representation is profligate.

---

## 4. Correctness anchor: the decomposition is exact, the dense readout is now coherent (fixed)

### 4a. The per-gate decomposition rule is an exact matrix identity

The mechanism (Section 5) is the operator decomposition the frame applies for a phase-gate spec
`"P"[\theta]`, with $c = e^{i\theta/2}$:

$$\frac{1+c}{2}\,\mathbb{1} + \frac{1-c}{2}\,Z \;=\; \mathrm{diag}(1, c) \;=\; P[\theta/2],$$

verified on the kernel:

- $c = e^{i\pi/4}$ (the value `"T"` uses via `"P"[\pi/2]`): $\tfrac{1+c}{2}\mathbb{1} +
  \tfrac{1-c}{2}Z = \mathrm{diag}(1, e^{i\pi/4}) = T$ exactly (residual $= 0$).
- For a generic `"P"[\theta]` the same identity gives $\mathrm{diag}(1, e^{i\theta/2})$, the
  half-angle gate, quantitatively reproducing the convention mismatch of Section 1
  ($[2,2]$ entry $0.9394 + 0.3429 i = e^{i\cdot 0.35}$ from the frame versus
  $0.7648 + 0.6442 i = e^{i\cdot 0.7}$ from `QuantumOperator["P"[0.7]]`).

This identity is materialization-independent: it certifies that the frame's *symbolic content*
(coefficients and tableaux, both verified correct per gate in Section 1) is the right stabilizer
decomposition. That is the part of the correctness anchor that holds.

### 4b. The dense materialization was wrong for genuine superpositions (now fixed)

The part that did **not** hold was `fr["StateVector"]` / `fr["State"]`. Before the fix, materializing
the frame and comparing to the dense `Method -> "Schrodinger"` statevector (fidelity
$|\langle \text{frame} | \text{dense} \rangle|^2$, global-phase blind) gave:

| $n$ | $k$ | components | fidelity(frame, dense), **pre-fix** |
|---:|---:|---:|---:|
| any | 0 | 1 | **1.000** |
| 1 | 1 | 2 | 0.500 |
| 2 | 2 | 4 | 0.250 |
| 3 | 1 | 2 | 0.500 |
| 3 | 2 | 4 | 0.250 |
| 3 | 3 | 8 | 0.125 |

Fidelity was $1$ **only** for single-component (rank-1) frames; for any genuine superposition the dense
readout was a different state, with fidelity degrading as components multiply.

**Root cause.** `StabilizerFrame["StateVector"]` summed each component's
`PauliStabilizer["State"]["StateVector"]`. That per-component materialization
(`Stabilizer/Conversions.m:69-80`) builds the stabilizer-group $+1$ eigenvector as
`Normalize @ Total @ NullSpace[...]`, whose overall phase is an arbitrary `NullSpace` representative,
and applies a recovered `"GlobalPhase"` only if that key is present. The key is set **only** by the
tomography constructors `PauliStabilizer[qs]` / `[qo]` (`Conversions.m:62-64`), never by the gate
updates that build frame components. So the components carried mutually incoherent global phases, and the
coherent sum $\sum_i c_i |s_i\rangle$ mixed them. Concretely, for $T|+\rangle$ the second component is
$Z|+\rangle = |-\rangle$, but `ps["Z",1]["State"]` returns $\{-1,1\}/\sqrt2 = -|-\rangle$; the spurious
$-1$ flipped the relative phase and the frame yielded $\{e^{i\pi/4}, 1\}/\sqrt2$ instead of the correct
$\{1, e^{i\pi/4}\}/\sqrt2$ (fidelity $0.5$). This was the
`PauliStabilizer[...][gates...]["State"]` "correct up to a global phase" contract
(`Stabilizer/GateUpdates.m`) leaking into a context where the relative phases are physical: the same
tableau states $|s_i\rangle$ with the right $c_i$ can describe different physical states (they differ by
per-component signs), so the abstract frame data $\{c_i, \text{ps}_i\}$ alone does not fix the state. A
materialization that re-derives each component's phase from the final tableau is therefore
fundamentally insufficient (an independent check confirmed it: a destabilizer-based phase recovery
matched the dense state on only $86/216$ random Clifford+$T$ circuits).

**The fix (in kernel).** The missing phase information is supplied at construction, not recovered at
read-out. A key structural fact makes this cheap: a $Z$ update leaves the AG tableau unchanged and a
Clifford update transforms it identically for every component, so **all components of a gate-built
frame share one tableau**, differing only in signs (verified across all $216$ random circuits). Each
component therefore equals a definite Pauli operator applied to a single reference component
($|s_i\rangle = P_i\,|s_1\rangle$). A gate-built `StabilizerFrame` now carries that relating Pauli
$P_i$ per component (encoded as $x/z$ bits plus a phase in $\{1,-1,i,-i\}$, with $P_1$ the identity).
Two $O(n)$ primitives maintain it: conjugation by a Clifford gate, $P_i \to G P_i G^\dagger$ (acting
only on the one or two qubits the gate touches, precomputed as a lookup), and left-multiplication by
$Z_q$ at a frame doubling. Materialization then computes the reference vector once,
$v_1 = $ `canonVec(ps_1)`, and forms $\sum_i c_i\, P_i\, v_1$, which recovers all relative phases up to
one overall global phase (physically irrelevant). Frames *without* relating Paulis (hand-built via
`StabilizerFrame[{{c, ps}, ...}]`, or produced by `Plus`, where components may not even share a
tableau) keep the previous independent per-component readout: it is correct per tableau, with the
relative phase unspecified, which is the inherent ambiguity of a hand-built frame.

After the fix, `fr["StateVector"]` is **amplitude-exact** against the dense state up to global phase
across $216$ random Clifford+$T$ circuits ($n=1\ldots3$, max amplitude error $0$), and $T|+\rangle$
returns exactly $\{1, e^{i\pi/4}\}/\sqrt2$. Generic `"P"[\theta]` is also materialized coherently; it
still realizes the documented half-angle gate $\mathrm{diag}(1, e^{i\theta/2})$, i.e. frame
`"P"[\theta]` matches the dense `"P"[\theta/2]` exactly (Section 1/4a caveat unchanged). Cites:
`Stabilizer/StabilizerFrame.m` (relating-Pauli helpers, branched materialization, gate distribute,
doubling, `Times`), `Stabilizer/GateUpdates.m` (initial $P[\theta]$ split sets the identity and $Z_j$
relating Paulis).

### 4c. The guarding test was too weak; now strengthened

`Tests/Stabilizer/Formulas.wlt` (`S11-T-on-Plus-equals-T-State`) intended to verify
$T|+\rangle = \{1, e^{i\pi/4}\}/\sqrt2$, but its check was `Total[diff] == 0` (the sum of the difference
vector). For this state the pre-fix error was a swap of the two amplitudes, whose differences are
$\{a, -a\}$, so `Total` was $0$ and the test passed even though `Max[Abs[diff]] \approx 0.46$ and the
fidelity was $0.5$. The companion tests `S11-T-on-Zero` and `S11-T-Tdag` passed legitimately, but only
because those target states ($|0\rangle$ and $|+\rangle$) are symmetric under the swap. The three S11
frame-state tests are now strict amplitude-wise assertions, `Chop[Max[Abs[N[... - expected]]]] == 0`,
and a new test `S11-Frame-Materialization-Matches-Dense` compares the frame readout for an entangled,
post-$T$-Clifford, two-$T$ circuit against the dense state up to one global phase. Verified
fail-before / pass-after: with the kernel reverted, `S11-T-on-Plus` and the new test both **fail**
($131/133$); with the fix, the full `Formulas.wlt` is green ($133/133$), as is the rest of the
stabilizer suite ($907$ tests across `Formulas`, `PauliStabilizer`, `Roundtrips`, `HybridInterop`,
`CliffordChannel`, `Correctness_TextbookResults`, `AuditMatrix`, `Connections`).

### 4d. Consequence for numeric readout

Single-state observables route through the same `StateVector`, so they are now correct: $\langle Y
\rangle$ on $T|+\rangle$ comes out $+0.707$ (it was $-0.707$ before the fix), and any
$\langle \psi | O | \psi \rangle$ from one coherent frame is exact because the overall phase cancels.
The frame inner product (`Stabilizer/InnerProduct.m:111-116`) still uses `StateVector`: its **magnitude
is correct**, but its **phase between two distinct states remains gauge-limited**, because each state's
overall phase is an independent arbitrary `canonVec` representative and the abstract stabilizer data
does not fix it (this is the same `PauliStabilizer` global-phase contract; e.g.
$\langle T_+|S_+\rangle$ has the right modulus $0.924$ but its phase carries the relative gauge of the
two `canonVec`s). **Net: the `StabilizerFrame` is now a correct dense-readout object for a single
state** (state vector, expectations) and remains a correct symbolic/counting object; the only residual
numeric caveat is the phase of an inner product *between* two separately-gauged stabilizer states. It
is still not a *scalable* stabilizer-rank simulator (no rank compression; Section 3), so for large
$T$-counts use `Method -> "Schrodinger"` or the tensor-network path.

---

## 5. Mechanism and kernel citations

The whole path, end to end:

1. **Dispatcher.** `PauliStabilizerApply[qco, state]` (`Stabilizer/PauliStabilizer.m:123-171`) is the
   `Method -> "Stabilizer"` entry point. It tries a compiled all-Clifford fold and, finding any
   non-encodable gate, falls through to a per-gate `Fold`. For each gate it inspects the result
   (`:150-160`): a `PauliStabilizer` (Clifford) is carried; a `StabilizerFrame` (the $P/T$ boundary)
   is carried; a `Plus` of stabilizer terms is carried; anything else triggers one
   `Message[PauliStabilizer::nonclifford, gate]` and returns `$Failed` (`:160`). This is the
   `RZ`/`RY` exit.

2. **Lift.** A bare `PauliStabilizer` meeting a diagonal phase gate
   (`Stabilizer/GateUpdates.m`, the `ps["P"[phase], j]` rule) returns a two-component frame, now with
   the per-component relating Paulis attached (identity for the surviving branch, $Z_j$ for the
   $Z$-inserted branch):
   ```
   ps["P"[phase], j] := With[{c = Exp[I phase/2], n = ps["Qubits"]},
       StabilizerFrame[<|
           "Components" -> {{(1 + c)/2, ps}, {(1 - c)/2, ps["Z", j]}},
           "Paulis" -> {identityPauli[n], Z_j pauli} |>]]
   ps["T", j]            := ps["P"[Pi/2], j]
   ps[SuperDagger["T"], j] := ps["P"[-Pi/2], j]
   ```
   This is the $c = e^{i\theta/2}$ (half-angle) convention of Section 1/4a.

3. **The doubling.** Once the state is a `StabilizerFrame`, a further phase gate
   (`Stabilizer/StabilizerFrame.m`, the `f["P"[phase], q]` rule) maps **every** component to two, and
   the new $Z$-inserted sub-component gets the relating Pauli $Z_q\,P_i$:
   ```
   f["P"[phase], q] := With[{c = Exp[I phase/2]},
       StabilizerFrame[<|
           "Components" -> Flatten[MapThread[
               {{(1 + c)/2 #1[[1]], #1[[2]]}, {(1 - c)/2 #1[[1]], #1[[2]]["Z", q]}} &,
               {Components, Paulis}], 1],
           "Paulis" -> Flatten[MapThread[{#2, timesZ[#2, q]} &, {Components, Paulis}], 1] |>]]
   ```
   The two entries per component are the literal source of $|\text{components}| \to
   2\,|\text{components}|$. There is no `DeleteDuplicates`, no coefficient merging, no rank reduction.
   Clifford gates, by contrast, distribute over components without growth, conjugating each relating
   Pauli by the gate (`Stabilizer/StabilizerFrame.m`, the `f[gate, args]` rule).

4. **Materialization.** `f["StateVector"]` (`StabilizerFrame.m`) now branches: a frame with a
   `"Paulis"` key materializes the reference component once,
   $v_1 = $ `ps_1["State"]["StateVector"]` (`Conversions.m:69-80`), and forms $\sum_i c_i\,P_i\,v_1$
   (each $P_i$ applied as an $O(2^n)$ phase-and-permutation, no $2^n \times 2^n$ matrix), which is
   phase-coherent (Section 4b). A frame without the key keeps the independent per-component sum.

5. **The cheap Clifford-regime baseline** is `stabilizerExpectation`
   (`Stabilizer/InnerProduct.m:156-193`): an $O(n^2)$ anticommutation check plus an $O(n^3)$
   $\mathbb{F}_2$ `LinearSolve`, for a *pure* `PauliStabilizer` only. It does not apply to a frame;
   once you are in the frame, that poly-time path is gone.

### Physics: the Bravyi-Gosset picture, and where QF sits in it

A $T$ gate applied to a stabilizer state is a coherent superposition of two stabilizer states,
$T_q|s\rangle = \tfrac{1+e^{i\pi/4}}{2}|s\rangle + \tfrac{1-e^{i\pi/4}}{2}Z_q|s\rangle$. Each
non-Clifford gate therefore multiplies the *stabilizer rank* (the number of stabilizer states needed
to write the state) by at most $2$, so a naive decomposition of a $t$-doped circuit has $2^t$ terms.
This is the stabilizer-decomposition / stabilizer-rank simulation idea (Bravyi-Gosset). The point of
that line of work is that the rank can be made far smaller than $2^t$: the stabilizer extent and rank
of $t$ magic states scale as roughly $\chi \sim 2^{\alpha t}$ with $\alpha \approx 0.396$ (and
better with approximate decompositions and term recombination), so $t = 50$ $T$ gates need
$\sim 2^{20}$ terms, not $2^{50}$.

QF's `StabilizerFrame` implements the *naive* end of this spectrum: it realizes the per-gate doubling
faithfully (Sections 1-2) but performs **no** rank compression, as Section 3 proves directly (it keeps
$2^k$ identical components for a rank-1 state). The file header (`StabilizerFrame.m:18-19`) cites
Garcia-Markov "Simulation of Quantum Circuits via Stabilizer Frames" (Quipu, arXiv:1712.03554), whose
contribution is precisely the frame *merging* that bounds the term count; that compression is not yet
implemented here. So relative to the literature QF is exponentially worse than an optimized
stabilizer-rank simulator in term count. Its dense state-vector readout is now correct (Section 4); the
only residual numeric caveat is the phase of an inner product between two separately-gauged stabilizer
states. It is best read today as a structural object that exhibits the boundary and gives a correct
dense readout for one state, not as a scalable magic-state simulator.

---

## 6. Caveats and provenance

- **Anchor drift.** Live HEAD `a601a187` (Sections 1-3 first run at `f370f568`; the Section 4 fix is
  `61bc0e39`, which changed `Stabilizer/StabilizerFrame.m` and `Stabilizer/GateUpdates.m`). Cites were
  re-verified against the live kernel at this HEAD; re-verify again if the repo moves (per the QF audit
  ritual).
- **Timings** are single-run wall times on one M2 Pro after `ClearSystemCache[]`; treat them as
  order-of-magnitude. The qualitative law (components $= 2^k$ exactly; memory and time $\sim 2^k$) is
  exact and reproduced across $k$.
- **Kernel changes made (Section 4 fix).** The materialization defect and its masking test are now
  fixed: `Stabilizer/StabilizerFrame.m` gained the relating-Pauli bookkeeping, a branched coherent
  `["StateVector"]`, and Pauli-tracking in the gate-distribute / doubling / `Times` rules;
  `Stabilizer/GateUpdates.m`'s initial $P[\theta]$ split now attaches the identity and $Z_j$ relating
  Paulis; `Tests/Stabilizer/Formulas.wlt` strengthened the three S11 frame-state tests to `Max`-based
  amplitude assertions and added `S11-Frame-Materialization-Matches-Dense`. Verified fail-before /
  pass-after and the full stabilizer suite green ($907$ tests). Remaining (not addressed here)
  kernel-side follow-ups: route `"RZ"[\theta]` to the `"P"` frame handler (with the angle convention
  reconciled) so the diagonal rotation is representable; track a global phase on `PauliStabilizer`
  gate updates so the inner product *between* distinct stabilizer states is phase-correct (not just
  modulus-correct); and, for real utility, add frame merging (the Quipu compression) to bound the term
  count below $2^k$.
- **`PackageScope` predicate trap.** `StabilizerFrameQ` and `PauliStabilizerQ` are `PackageScope`, so
  bare calls after `Needs` stay unevaluated (silent false negatives). Classify via `Head[...]` and the
  public `StabilizerStateQ`, as the scripts do.

---

## Appendix: repro scripts

Saved under [`scripts/`](scripts/); each loads the working tree with `PacletDirectoryLoad`, uses no
`Quiet`/`Print`, and ends with a result Association (echo with
`wolframscript -code 'Print[Get["scripts/<name>.wls"]]'`).

- [`sf_probe1_dichotomy.wls`](scripts/sf_probe1_dichotomy.wls) - Section 1: Clifford vs one
  non-Clifford gate; the `P[\theta]` vs gate-level convention; `RZ`/`RY` `$Failed`.
- [`sf_probe2_scaling.wls`](scripts/sf_probe2_scaling.wls) - Section 2: components / time / memory vs
  $k$ (default $k = 0\ldots16$; the report also ran $k = 17, 18$).
- [`sf_probe3_naive.wls`](scripts/sf_probe3_naive.wls) - Section 3: `P[\theta]` blow-up; $2^k$
  identical components on $|0\rangle$; period-8 $T^8|+\rangle = |+\rangle$.
- [`sf_probe4_correctness.wls`](scripts/sf_probe4_correctness.wls) - Section 4: matrix-identity
  decomposition; the `Total`-vs-`Max` test-masking demonstration. (The frame-vs-dense fidelity table it
  prints documented the *pre-fix* behavior; against the fixed kernel the fidelities are now $1$.)
- [`sf_probe5_corollaries.wls`](scripts/sf_probe5_corollaries.wls) - Section 4d (pre-fix snapshot):
  the then-wrong $\langle Y \rangle$ and the inner-product gauge issue.
- [`sf_probe6_fix_verification.wls`](scripts/sf_probe6_fix_verification.wls) - Section 4b/4c fix
  verification: $T|+\rangle$ exact; $216$ random Clifford+$T$ circuits amplitude-exact vs dense;
  $\langle Y\rangle = +0.7071$; coherent half-angle `"P"[\theta]`; manual frames unchanged.
