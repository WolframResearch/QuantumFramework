# Claim 6 for QF-Showcase.md: coherent measurement records and feedforward. Full implementation handoff

**Purpose of this file.** Self-contained instructions for adding the measurement-feedforward material to `OngoingProjects/Why QF/QF-Showcase.md`, plus the designs for the examples that were deliberately deferred to other documents. A cold session should be able to execute this file top to bottom with no memory of the conversation that produced it.

**Provenance.** Investigation of `QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1}, {2}, {3}, "CNOT" -> {-1, 2}}]`, 2026-06-12, at working-tree HEAD `6f953ad5`. Companion analysis doc: `OngoingProjects/Why QF/QF-Feedforward-Pipeline.md` (capability matrix, platform comparison, full example ladder). Verification probes from that session: `/tmp/qf-dilation-probe{1..6}.wls` (tmp files; re-create from the code blocks below if gone). Everything marked ✅ below was executed and its stated output captured; everything marked ⚠️ is designed but not yet run.

**How to verify anything here.** Run via `wolframscript -file ...` (never the WolframLanguageEvaluator MCP), loading the working tree, not the installed paclet:

```wolfram
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

No `Quiet`, no `Print` in showcase cells; end each Input with the expression to show.

---

## 1. The decision (and why)

**Add to QF-Showcase.md: a new Claim 6**, not a subsection of Claim 2.

> **Claim 6: Measurement is unitary physics inside the model: records are quantum wires, and feedforward is a controlled gate.**

Content, in order:
1. Two-sentence mechanics of pointer wires, hung on the seed example `"CNOT" -> {-1, 2}` (§3 below).
2. **Teleportation as the deferred-measurement identity** (§5, ✅ verified). The recognizable two-cell on-ramp.
3. **The QEC cycle with a coherent error** (§6, ✅ verified). The payoff: digitization of errors as an executed symbolic identity, post-recovery fidelity $\equiv 1$ for every error angle.
4. One honest platform paragraph (§8): classical feedforward exists everywhere; coherent records exist only here.
5. Optional third beat if the section can afford it: **magic-state injection** (§7, ⚠️ designed).

**Deferred to other documents** (do NOT put in QF-Showcase.md): constant-depth dynamic circuits (§9.1), the MBQC ladder (§9.2), and Wigner's friend (§9.3, flagship of a separate foundations showcase).

Why this split: QF-Showcase.md's house style is "every example ends in an exact symbolic statement other platforms cannot produce" (the Rabi propagator, $T_2 \le 2T_1$ as a theorem, the QAOA maximum provably equal to the true cut). The QEC example is the one feedforward candidate whose conclusion is a theorem rather than a demo. Teleportation alone reads trivial; as the definition-bearing on-ramp it is the recognizable rung. Wigner's friend changes the document's register (foundations no-go staging vs. working-physicist closed forms) and deserves its own document, the way the open-systems examples got `QF-OpenSystems-Showcase.md`.

Placement: between the current Claim 2 and Claim 3 (it extends Claim 2's "one object model" scope and inherits Claim 1's symbolics), renumbering the later claims, **or** appended as Claim 6 after Claim 5 to avoid renumbering. Appending is the low-churn choice; the claim list at the top of the document must be updated either way, and the "Where This Leaves Us" closing paragraph should gain one sentence.

---

## 2. Background: what QF's measurement actually is (for the prose)

QF implements a projective measurement as the von Neumann dilation: the isometry

$$
V \;:\; |\psi\rangle_S \;\longmapsto\; \sum_k \bigl(P_k|\psi\rangle\bigr)_S \otimes |k\rangle_A ,
$$

where the appended pointer qudit $A$ is the measuring device, its basis labelled by eigenvalues. No branch is discarded; "collapse" is what the system looks like with the pointer traced out. The cost is one extra qudit per measurement ($2^n \to 2^{n+k}$ after $k$ measurements); the payoff is that the record is a first-class quantum wire.

Kernel anchors (verified live at `6f953ad5`):
- Dilation construction: `QuantumMeasurementOperator/Properties.m:206-276` (the `"SuperOperator"` property).
- Pointer-wire assignment in circuits: `getNext` in `QuantumCircuitOperator/Properties.m:117-124`. Records get non-positive wire numbers top-down: first measurement $\to 0$, second $\to -1$, third $\to -2$.
- Longer cost discussion: `OngoingProjects/Measurement Dilation/measurement-hilbert-space-blowup.md`.

---

## 3. The seed example, decoded (use as the two-sentence opener) ✅

```wolfram
qc = QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1}, {2}, {3}, "CNOT" -> {-1, 2}}];
qc["Diagram"]
```

Decoding: `{1}, {2}, {3}` are three separate computational measurements; their records land on wires $0, -1, -2$. So `"CNOT" -> {-1, 2}` is a CNOT **controlled by the measurement record of qubit 2**, targeting qubit 2: an active reset, returning qubit 2 to $|0\rangle$ in every branch.

Verified facts to quote if needed:
- `qc[]` returns a `QuantumMeasurement` whose outcome distribution is exactly $p(x) = \tfrac14\cos^2(\pi x/8)$, $x = 0,\dots,7$, i.e. $|\langle x|\,\mathrm{QFT}_8\,|\mathrm{GHZ}\rangle|^2$; the conditional gate leaves the record statistics untouched, as it must.
- `qc["QuantumOperator"]` compiles to a single `QuantumMeasurementOperator` with output order $\{-2,-1,0,1,2,3\}$: one 6-wire isometry for the whole experiment.

---

## 4. Traps the implementer must respect (all found and verified this session)

1. **Degenerate observables: use the projector-list constructor, never the observable form.**
   `QuantumMeasurementOperator[QuantumOperator["ZZ"], {1,2}]` is rank-1 in the eigenbasis: it produces a **dim-4** pointer (one slot per eigen*vector*, eigenvalues sorted ascending, $-1$ slots before $+1$ slots) and collapses superpositions *inside* a degenerate eigenspace. It is not a parity measurement. The correct syndrome measurement is the projector list:
   ```wolfram
   QuantumMeasurementOperator[{Pplus, Pminus}, {1, 2}]
   ```
   which gives a dim-2 pointer (Type "POVM"; outcome 1 $\leftrightarrow$ first listed operator) and preserves in-eigenspace coherence (verified: $a|00\rangle + b|11\rangle$ survives with amplitudes intact).
2. **Pointer wires cannot be re-measured.** Measurement targets must be positive (`targetQ`); `QuantumMeasurementOperator[{0}]` returns a `Failure`, and `{0}` inside a circuit does not become a measurement. Gates on pointer wires are fine. Record marginals in any basis are read off the dilated state via `QuantumPartialTrace`.
3. **`Eigenorder` bookkeeping disagreement (cosmetic here).** `qc["Eigenorder"]` drops a pointer wire from its list as soon as a downstream unitary touches it (`Properties.m:326-327`), but the compiled operator and the `QuantumMeasurement` result still count all pointers as eigen-registers. Do not introspect `qc["Eigenorder"]` programmatically in showcase cells.
4. **Mixed-polarity conditioning call form.** Control-on-1 list and control-on-0 list:
   ```wolfram
   QuantumOperator["C"[QuantumOperator["X", {1}], {0}, {-1}]]   (* X on wire 1 if wire 0 reads 1 AND wire -1 reads 0 *)
   ```
   Named shortcuts also work on pointer wires: `"CNOT" -> {-1, 3}`, `"CZ" -> {0, 3}`, `"Toffoli" -> {0, -1, 2}` (controls first, target last).
5. **Result head varies.** If every pointer wire is left untouched after its measurement, application returns a `QuantumMeasurement`; once a unitary touches a pointer (teleportation, QEC), it returns a `QuantumState` on the dilated register. Write cells so either head is handled, or pin with `If[Head[...] === QuantumMeasurement, #["State"], #]&`.
6. **Wire layout of the dilated state.** The output order is the sorted wire list (e.g. $\{-1, 0, 1, 2, 3\}$), so system wire $q$ sits at position `Position[Sort[order], q]` in the state's qudit list. Compute positions, never hardcode.
7. **Cost honesty.** Each measurement adds a qudit. Keep showcase examples at $\le 8$ total wires and say so in one clause ("the price of keeping every branch alive").

---

## 5. Section beat 1: teleportation as the deferred-measurement identity ✅

**Punchline for the prose:** output fidelity is exactly $1$, and there is *no classical communication step anywhere in the code*: the corrections $X^{m_2} Z^{m_1}$ are ordinary controlled gates whose controls are the measurement records. Nielsen-Chuang §4.4's deferred-measurement principle is the primitive here, not a compiler pass.

Verified code (fidelity $1.0$ on a Haar-random input state):

```wolfram
in = QuantumState["RandomPure"[1]];
tele = QuantumCircuitOperator[{
    "H" -> 2, "CNOT" -> {2, 3},        (* Bell pair on wires 2,3 *)
    "CNOT" -> {1, 2}, "H" -> 1,        (* rotate 1,2 into the Bell basis *)
    {1}, {2},                          (* records: wire 0 and wire -1 *)
    "CNOT" -> {-1, 3},                 (* X^{m2} on wire 3 *)
    "CZ"   -> {0, 3}                   (* Z^{m1} on wire 3 *)
}];
out = tele[QuantumTensorProduct[in, QuantumState["00"]]];   (* a QuantumState on 5 wires *)
```

Readout cell (wire 3 is the last qudit of the sorted dilated register $\{-1,0,1,2,3\}$):

```wolfram
QuantumSimilarity[QuantumPartialTrace[out, {1, 2, 3, 4}], in, "Fidelity"]
```

Verified output: `1.` (to machine precision; measured $0.9999999999999993$).

**Upgrade worth one verification cell (⚠️ not yet run):** replace the random input by the symbolic `QuantumState[{Cos[\[Alpha]], Sin[\[Alpha]] Exp[I \[Beta]]}]` and `FullSimplify` the output qudit's Bloch vector against the input's, under real assumptions. Symbolic flow through measurement + feedforward is verified in general (the QEC example below ran symbolically), so risk is low; if simplification stalls, ship the random-input version, whose prose punchline is already strong. Suggested check:

```wolfram
psi = QuantumState[{Cos[\[Alpha]], Sin[\[Alpha]] Exp[I \[Beta]]}];
out = tele[QuantumTensorProduct[psi, QuantumState["00"]]];
FullSimplify[ComplexExpand[
   QuantumPartialTrace[out, {1, 2, 3, 4}]["BlochVector"] - psi["BlochVector"]],
  {\[Alpha] \[Element] Reals, \[Beta] \[Element] Reals}]
```
Expected: `{0, 0, 0}`.

Diagram cell (`tele["Diagram"]`) is worth including: the records visibly feed the controlled gates.

---

## 6. Section beat 2 (the payoff): QEC with a coherent error, digitization as a theorem ✅

**Punchline for the prose:** apply a *coherent* error $e^{-i\theta X_2}$ (not a Pauli, not stochastic) to the 3-qubit repetition code, extract the two stabilizer syndromes $Z_1Z_2$ and $Z_2Z_3$ with genuine degenerate (coherence-preserving) measurements, recover with three record-controlled gates, and the codeword fidelity is exactly $1$ **for every $\theta$, symbolically**. The syndrome measurement digitizes the continuous error into $\{\mathbb{1}, X_2\}$ with weights $\cos^2\theta$ and $\sin^2\theta$, and the feedforward fixes both branches. The digitization-of-errors argument at the heart of fault tolerance (Shor, PRA 52, R2493 (1995); the standard argument in Nielsen-Chuang §10.2), normally stated in words, executed as an identity.

### 6.1 Verified numeric version (θ = 0.7, fidelity 1.0000000000000004)

```wolfram
(* degenerate stabilizer projectors: P± = (1 ± Z⊗Z)/2, as a projector-list QMO (dim-2 pointer) *)
id2 = IdentityMatrix[2];
zz  = KroneckerProduct[PauliMatrix[3], PauliMatrix[3]];
syn[q1_, q2_] := QuantumMeasurementOperator[
   {QuantumOperator[(KroneckerProduct[id2, id2] + zz)/2, {q1, q2}],
    QuantumOperator[(KroneckerProduct[id2, id2] - zz)/2, {q1, q2}]}, {q1, q2}];

a = Cos[0.3]; b = Sin[0.3];
in   = QuantumState[{a, b}];
psiL = QuantumState[a UnitVector[8, 1] + b UnitVector[8, 8]];   (* a|000> + b|111> *)

qec = QuantumCircuitOperator[{
    "CNOT" -> {1, 2}, "CNOT" -> {1, 3},                       (* encode *)
    QuantumOperator["RX"[2 \[Theta]], {2}] /. \[Theta] -> 0.7, (* coherent error exp(-i θ X₂) *)
    syn[1, 2],                                                  (* record on wire 0 *)
    syn[2, 3],                                                  (* record on wire -1 *)
    QuantumOperator["C"[QuantumOperator["X", {1}], {0}, {-1}]], (* X₁ if syndrome (1,0) *)
    QuantumOperator["C"[QuantumOperator["X", {2}], {0, -1}, {}]],(* X₂ if syndrome (1,1) *)
    QuantumOperator["C"[QuantumOperator["X", {3}], {-1}, {0}]]  (* X₃ if syndrome (0,1) *)
}];
st = qec[QuantumTensorProduct[in, QuantumState["00"]]];
```

Readout (system wires 1,2,3 sit at positions 3,4,5 of the sorted dilated register $\{-1,0,1,2,3\}$; the original verification used ancilla-based extraction with the same recovery and layout $\{-1,0,1,\dots,5\}$, positions computed not hardcoded):

```wolfram
rho = QuantumPartialTrace[st, Complement[Range[st["Qudits"]], {3, 4, 5}]];
QuantumSimilarity[rho, psiL, "Fidelity"]
```

Verified output: `1.0000000000000004`.

> ⚠️ **One re-verification required before publishing:** the session's fidelity-$1$ run used **ancilla-based** syndrome extraction (CNOTs onto explicit ancilla qubits 4 and 5, then `{4}, {5}`; probe 4A). The projector-list `syn` measurement was verified separately to be degenerate and coherence-preserving (probe 5A). The combination above (projector-list syndromes feeding the same recovery) is the cleaner showcase form but was not run end-to-end. Run it once; if anything is off (e.g. outcome-slot polarity), fall back to the verified ancilla version:
> ```wolfram
> "CNOT" -> {1, 4}, "CNOT" -> {2, 4},   (* ancilla 4 = Z1Z2 parity *)
> "CNOT" -> {2, 5}, "CNOT" -> {3, 5},   (* ancilla 5 = Z2Z3 parity *)
> {4}, {5},                              (* records: wires 0 and -1 *)
> ```
> with the identical three recovery gates. That version's fidelity-1 output is captured.
> Polarity note for the projector-list version: outcome slot 1 corresponds to the **first** listed projector ($P_+$), slot 2 to $P_-$; the recovery's control-on-1 means "second outcome fired", i.e. eigenvalue $-1$. If recovery misfires, swap $P_+ \leftrightarrow P_-$ in `syn` rather than rewriting the controls.

### 6.2 Verified symbolic version (the actual theorem)

Symbolic parameters flow through the whole cycle. Verified: with input $\{\alpha, \beta\}$ and error angle $\theta$ symbolic, the syndrome distribution comes out $\propto \cos^2\theta$ on "no error" and $\propto \sin^2\theta$ on "$X_2$", independent of $\alpha, \beta$ (probe 4C, ancilla-extraction form; simplification needs `Assuming[{\[Theta] \[Element] Reals, ...}, ...]` plus normalization of $|\alpha|^2 + |\beta|^2 = 1$). For the showcase, the strongest closing cell is the symbolic fidelity:

```wolfram
(* same circuit, θ symbolic, input {α, β} symbolic *)
Assuming[{\[Theta] \[Element] Reals, \[Alpha] \[Element] Reals, \[Beta] \[Element] Reals,
    \[Alpha]^2 + \[Beta]^2 == 1},
  FullSimplify[ComplexExpand[1 - QuantumDistance[rhoSym, psiLSym, "Fidelity"]]]]
```

Expected: `1` identically. ⚠️ The symbolic *fidelity* cell specifically was not run (only the symbolic syndrome distribution was); if `FullSimplify` stalls on the mixed-state fidelity, present the symbolic syndrome weights $\{\cos^2\theta, \sin^2\theta\}$ plus the numeric fidelity sweep `Table[..., {θ, 0., 1.5, 0.1}]` returning all 1.s, which makes the same point.

### 6.3 What the prose must name

- The error is **not** in the code's error model a priori: it is a continuous rotation. Digitization is the content.
- The syndrome measurement is **degenerate by construction** (projector list): it learns the eigenvalue, never the eigenvector, so the encoded superposition survives. Contrast in one clause with the observable form (trap §4.1), which would destroy the code: the difference between the two constructors *is* the physics of syndrome extraction.
- Three recovery gates cover the three single-qubit error locations: $(1,0) \to X_1$, $(1,1) \to X_2$, $(0,1) \to X_3$.
- Cost clause: 7 wires total (3 data + 2 ancillas + 2 records in the ancilla version; 3 data + 2 records in the projector version).

---

## 7. Optional third beat: magic-state injection (⚠️ designed, not run)

**Punchline:** the $T$-gate gadget, the reason feedforward is the engine of fault tolerance, compiles to exactly $T$. Prove-by-compiling, same register as the QAOA "provably equals the cut".

```wolfram
(* resource state |T> = T|+> on wire 2; data on wire 1 *)
tState = QuantumOperator["T"][QuantumState["+"]];
gadget = QuantumCircuitOperator[{
    tState -> {2},
    "CNOT" -> {2, 1},          (* note: ancilla controls data in the standard injection *)
    {1},                       (* measure data; record on wire 0 *)
    QuantumOperator["C"[QuantumOperator["S", {2}], {0}]]   (* S correction if record = 1 *)
}];
```

Verification target: `gadget[QuantumTensorProduct[psi, ...]]` reduced on the output wire equals `QuantumOperator["T"][psi]` for symbolic `psi`, or `gadget["QuantumOperator"]` restricted to the data wire equals $T$ up to the dilation registers.

Risks (why this is ⚠️): the injection circuit has two common conventions (CNOT direction, which wire is measured, whether the output lands on the ancilla wire); expect one iteration to get the right one. Phase convention: compare up to global phase (`Abs[Tr[U†V]]`), per the KAK lesson in `audit/mistakes.md`. Literature: Gottesman-Chuang, *Nature* 402, 390 (1999); Bravyi-Kitaev, PRA 71, 022316 (2005).

Include only if the Claim 6 section can afford a third example without bloating; teleportation + QEC already carry the claim.

---

## 8. The honesty paragraph (platform comparison, one paragraph in the showcase)

Required content, in the document's voice:

- Classical feedforward is table stakes: Qiskit dynamic circuits (`if_test` on clbits, runs on IBM hardware), Cirq's classically controlled operations, Stim's `CX rec[-1] 0`, Q#'s control flow.
- PennyLane goes one step further internally: `qml.defer_measurements` performs **exactly this dilation** (one ancilla per measurement, exponential overhead acknowledged in its docs) to support `qml.cond` on simulators. But it is a hidden compiler transform: the ancilla is unaddressable, the joint state is not an object, and the only operation on a record is classical conditioning.
- What exists only here: the record is a quantum wire you can *control on* (this section), *rotate and entangle*, *un-measure* (the adjoint dilation restores the input exactly, relative phase included; verified, see §9.3), and read jointly with the system; and the whole experiment, measurement and feedback included, is one symbolic object (the QEC theorem above).
- The cost sentence: each measurement adds a qudit, $2^n \to 2^{n+k}$; that is the price of keeping every branch alive and precisely what the capabilities above are made of.

References for the writer: [PennyLane dynamic circuits docs](https://docs.pennylane.ai/en/stable/introduction/dynamic_quantum_circuits.html), [qml.defer_measurements](https://docs.pennylane.ai/en/stable/code/api/pennylane.defer_measurements.html), [qml.cond](https://docs.pennylane.ai/en/stable/code/api/pennylane.cond.html).

---

## 9. Deferred examples (for other documents; designs preserved)

### 9.1 Constant-depth long-range entanglement with feedforward (⚠️ designed)
Miniature of Bäumer et al., *PRX Quantum* 5, 030339 (2024) ([arXiv:2308.13065](https://arxiv.org/abs/2308.13065)): teleport a CNOT across a chain by consuming intermediate Bell pairs, measuring the middle qubits, feeding the records forward; dynamic depth is constant while unitary depth grows with distance. In QF the identity "dynamic circuit $=$ long-range CNOT" is provable by compiling. Keep to 5-7 system qubits (records double the register). Destination: a standalone "dynamic circuits, coherently" note or a follow-up section of the platform-comparison doc.

### 9.2 MBQC ladder (⚠️ designed)
A 1D cluster-state wire executing a programmed single-qubit rotation: measure each site in a rotated basis (QF: arbitrary-basis `QuantumMeasurementOperator["X", {q}]`-style or explicit basis QMOs), correct byproduct Paulis with record-controlled gates; output deterministic with corrections, maximally mixed without, symbolically in the measurement angle. **Constraint discovered:** measurement *bases* cannot be adaptive on records in QF (only unitary corrections can); use the one-bit-teleportation ladder formulation, which needs only byproduct corrections. Literature: Raussendorf-Briegel, PRL 86, 5188 (2001); Raussendorf-Browne-Briegel, PRA 68, 022312 (2003). Destination: tech-note / tutorial material.

### 9.3 Wigner's friend, executable (⚠️ designed; primitives ✅ verified)
Flagship of a **separate foundations showcase** (sibling of `QF-OpenSystems-Showcase.md`). The friend's measurement is the dilation isometry kept coherent; Wigner can (i) read the record, (ii) **reverse the friend's measurement**, (iii) measure system+friend jointly in an entangled basis. Verified primitive for (ii), exact including relative phase:

```wolfram
m  = QuantumMeasurementOperator[{1}];
V  = m["SuperOperator"];                       (* the dilation isometry; V†V = I verified *)
in = QuantumState[{Cos[0.4], Sin[0.4] Exp[I 0.9]}];
QuantumCircuitOperator[{m, V["Dagger"]}][in]   (* returns exactly `in`, phase included *)
```

Also verified: system-apparatus entanglement entropy of measure-after-Hadamard is exactly 1 bit (`QuantumEntanglementMonotone[st, {{1},{2}}, "EntanglementEntropy"]` on the dilated state). Design constraint: pointer wires cannot be re-measured (trap §4.2), so Wigner's own final measurement must target positive wires; model the friend's lab on explicit positive wires where Wigner needs to measure it, or read Wigner's statistics from the final dilated state. Literature: Wigner (1961); Deutsch, IJTP 24, 1 (1985); Frauchiger-Renner, *Nat. Commun.* 9, 3711 (2018); Bong et al., *Nat. Phys.* 16, 1199 (2020).

### Dropped (do not resurrect without new reasons)
- **Iterative/Kitaev phase estimation with qubit reuse**: the hardware win is reusing one ancilla; QF's registers accumulate instead, so the example showcases the cost, not the feature.
- **Delayed-choice eraser via pointer re-measurement**: blocked by trap §4.2; with the marker on a positive wire it collapses to a standard two-qubit demo every platform has.

---

## 10. Editing checklist for QF-Showcase.md

1. Add Claim 6 to the numbered list in the intro (one sentence: *measurement is unitary physics inside the model; records are quantum wires and feedforward is a controlled gate*).
2. Insert the Claim 6 section (seed mechanics, teleportation, QEC, honesty paragraph), following the house style: every cell followed by outcome prose with the result inline (per the `feedback_outcome_prose_inline_math` convention: "The fidelity returns as exactly $1$...", never a bare displayed pseudo-output).
3. Re-run **every** new cell with `wolframscript` against the working tree before committing the text; update stated outputs to the captured ones. Outstanding pre-publication runs: §5 symbolic teleportation, §6.1 projector-list end-to-end, §6.2 symbolic fidelity, §7 if included.
4. Extend "Where This Leaves Us" with one sentence (e.g. "...and a measurement became a wire: teleportation ran without a classical channel, and a coherent error was digitized and corrected with the codeword fidelity returning as the identity, symbolically").
5. Track status in `OngoingProjects/Why QF/QF-Showcase-Pipeline.md` (add a Claim 6 row) and update the `project_qf_feedforward_showcase` memory entry when the section lands.
6. Style rules in force: no em dashes anywhere; all math as proper TeX (`$...$` / `$$...$$`); no `Quiet`; no `Print`; comments label physics, not history.
