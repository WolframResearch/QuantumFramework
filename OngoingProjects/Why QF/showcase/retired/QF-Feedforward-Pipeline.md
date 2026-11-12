# Coherent measurement records and feedforward in QF: audit + showcase pipeline

**Status:** investigation complete, examples designed; 3 of 6 verified on live kernel. Started 2026-06-12.
**Anchor:** working-tree HEAD `6f953ad5` (4 commits past audit anchor `3225a899`; those commits touch only Stabilizer/QuantumEvolve, and every mechanic below was re-verified live via `wolframscript`, probes in `/tmp/qf-dilation-probe{1..6}.wls`).

---

## 1. The seed example, audited

```wolfram
QuantumCircuitOperator[{"GHZ", "Fourier"[3], {1}, {2}, {3}, "CNOT" -> {-1, 2}}]
```

### What it does, mechanically (all verified live)

1. `"GHZ"` prepares $(|000\rangle + |111\rangle)/\sqrt{2}$, `"Fourier"[3]` applies the 8-dimensional QFT.
2. Each `{q}` is a *separate* computational measurement of qubit $q$. QF realizes it as the von Neumann dilation isometry
   $$
   V \;:\; |\psi\rangle_S \mapsto \sum_k \bigl(P_k|\psi\rangle\bigr)_S \otimes |k\rangle_A ,
   $$
   appending one **pointer qudit** per measurement. Pointer wires are assigned non-positive indices top-down by `getNext` (`QuantumCircuitOperator/Properties.m:117-124`): the first measurement's record sits on wire $0$, the second on $-1$, the third on $-2$.
3. `"CNOT" -> {-1, 2}` is therefore a CNOT **controlled by the measurement record of qubit 2**, targeting system qubit 2. Since in each branch the post-measurement qubit 2 equals its record, this is an *active reset*: it returns qubit 2 to $|0\rangle$ deterministically, conditioned on the outcome. Feedforward, expressed as an ordinary controlled unitary.
4. Verified output: `qc[]` returns a `QuantumMeasurement` whose outcome distribution is exactly $p(x) = \tfrac14\cos^2(\pi x/8)$, $x = 0,\dots,7$, i.e. $|\langle x|\,\mathrm{QFT}_8\,|\mathrm{GHZ}\rangle|^2$; the conditional CNOT does not disturb the record statistics, as it must not.
5. The compiled circuit is a single `QuantumMeasurementOperator` with output order $\{-2,-1,0,1,2,3\}$: a 6-qubit isometry for a 3-qubit experiment. The cost model is the known dilation blow-up, $2^{n} \to 2^{n+k}$ after $k$ measurements (see `OngoingProjects/Measurement Dilation/measurement-hilbert-space-blowup.md`).

### Bookkeeping subtlety found

`qc["Eigenorder"]` returns `{0, -2}`, not `{0,-1,-2}`: the fold at `Properties.m:326-327` *removes* a pointer wire from the eigen-register list as soon as a downstream unitary touches it. Physically right (a record that gets acted on coherently is no longer a passive classical outcome), but the compiled `QuantumMeasurementOperator` and the `QuantumMeasurement` result still report 3 eigenqudits. The two bookkeepers disagree; harmless for these examples but worth knowing.

---

## 2. Verified capability matrix (live kernel, 2026-06-12)

| capability | status | evidence |
|---|---|---|
| Gates controlled on pointer wires (`"CNOT" -> {-1, q}`, `"CZ" -> {0, q}`, `"Toffoli" -> {0,-1,q}`) | ✅ | probes 2, 4 |
| Generic control-1/control-0 mix on pointers: `QuantumOperator["C"[QuantumOperator["X",{1}], {0}, {-1}]]` | ✅ | probe 4 |
| Teleportation with quantum-controlled corrections, no classical conditioning | ✅ fidelity $1.0$ on a Haar-random input | probe 2 |
| Full QEC cycle (encode, coherent error, syndrome, 3-branch feedforward recovery) | ✅ fidelity $\equiv 1$ | probe 4 |
| True **degenerate** (eigenvalue-only) projective measurement via projector list: `QuantumMeasurementOperator[{P_+, P_-}, {1,2}]` gives a dim-2 pointer and preserves in-eigenspace coherence | ✅ | probe 5A |
| Symbolic parameters through measurement + feedforward (error angle $\theta$, amplitudes $\alpha,\beta$) | ✅ syndrome probs $\propto \cos^2\theta,\ \sin^2\theta$ | probe 4C |
| Coherent operations *on the record* (H on wire 0) | ✅ | probes 2C, 5B |
| System-apparatus entanglement entropy from the dilated state | ✅ exactly 1 bit for an H-then-measure | probe 5C |
| **Un-measuring**: `{qmo, qmo["SuperOperator"]["Dagger"]}` restores the input exactly, relative phase included | ✅ | probe 6 |
| Re-measuring a pointer wire (`QuantumMeasurementOperator[{0}]` or `{0}` in a circuit) | ❌ targets must be positive (`targetQ`); read record marginals from the dilated state via `QuantumPartialTrace` instead | probe 2C |
| `QuantumMeasurementOperator[QuantumOperator["ZZ"], {1,2}]` as a syndrome measurement | ⚠️ NOT degenerate: dim-4 pointer, one slot per eigen**vector**; it collapses in-eigenspace superpositions. Use the projector-list form for stabilizer/syndrome measurements | probe 3A, 5A |

---

## 3. How unique is this, honestly (platform comparison)

Every serious platform now has *classical* feedforward, so the showcase must not claim "conditional operations on measurement results" as such:

- **Qiskit**: dynamic circuits (`if_test` on clbits), runs on IBM hardware.
- **Cirq**: `ClassicallyControlledOperation` with sympy conditions.
- **Stim**: measurement-record feedback `CX rec[-1] 0` (Clifford only).
- **Q#**: full classical control flow on results.
- **PennyLane**: `qml.measure` + `qml.cond`; on statevector devices this is implemented by `qml.defer_measurements`, which performs **exactly QF's dilation** (one ancilla per measurement, exponential overhead acknowledged in their docs). But it is a hidden compiler transform: the ancilla is inaccessible, the joint state is not an object, and the only thing you can do with a record is classically condition on it.
- **QuTiP**: measurement functions, no circuit-level record object.

**The defensible unique claims for QF** (each backed by a probe above):

1. **The record is a first-class quantum wire.** Conditioning is an ordinary controlled unitary; the deferred-measurement principle is the *primitive*, not a compile pass. Nielsen-Chuang §4.4 becomes a one-line circuit element.
2. **Coherent access to the record after the fact**: rotate it, entangle it, control on it, un-measure it (verified exact), or compute system-apparatus entanglement and mutual information. Classical-bit platforms cannot express any of these.
3. **Measurement generality**: arbitrary bases, observables, POVMs, projector lists, qudit pointers of any dimension; mid-circuit, in one object model.
4. **Symbolic end-to-end**: outcome distributions and post-feedforward fidelities as closed forms in circuit parameters. No mainstream platform does this at all.
5. **One object**: the whole experiment (prep, measure, feedforward) compiles to a single isometry whose matrix you can inspect, conjugate, and verify identities with.

The cost trade-off must be stated alongside: each measurement adds a qudit, $2^{n}\!\to\!2^{n+k}$, so examples live at $n \lesssim 8$. That is the price of keeping every branch alive, and it is exactly what enables points 1-5.

---

## 4. Showcase ladder (recognizable → consequential), PhD-level

Per-example: what it shows about the *object model*, target audience, literature anchor, verification status.

### 4.1 Teleportation as the deferred-measurement identity ✅ verified
Bell measurement, then `"CNOT" -> {-1, 3}`, `"CZ" -> {0, 3}`. Output qubit fidelity $1.0$ with a random input, with **no classical communication step anywhere in the code**: the circuit *is* the deferred-measurement form of teleportation, and compiling it proves the identity $X^{m_2} Z^{m_1}$-correction $=$ controlled gates off the records.
- Object-model property: pointer wires + controlled unitaries; `qc["QuantumOperator"]` as an identity-prover.
- Audience: anyone; the recognizable rung.
- Literature: Nielsen-Chuang §4.4 (deferred measurement); Gottesman-Chuang, *Nature* 402, 390 (1999) for where this leads.

### 4.2 Full QEC cycle with a coherent error: digitization made visible ✅ verified
3-qubit bit-flip code, error $e^{-i\theta X_2}$ (not a Pauli!), syndrome extraction, recovery as three feedforward gates with mixed 1/0-controls on the two syndrome records:
```wolfram
QuantumCircuitOperator[{
  "CNOT" -> {1,2}, "CNOT" -> {1,3},
  QuantumOperator["RX"[2 θ], {2}],
  QuantumMeasurementOperator[{Pplus12, Pminus12}, {1,2}],   (* dim-2 syndrome pointer *)
  QuantumMeasurementOperator[{Pplus23, Pminus23}, {2,3}],
  QuantumOperator["C"[QuantumOperator["X",{1}], {0}, {-1}]],
  QuantumOperator["C"[QuantumOperator["X",{2}], {0,-1}, {}]],
  QuantumOperator["C"[QuantumOperator["X",{3}], {-1}, {0}]]
}]
```
Post-recovery codeword fidelity is **exactly 1 for every $\theta$**, and symbolically so: the syndrome measurement digitizes the continuous error into $\{\mathbb{1}, X_2\}$ with weights $\cos^2\theta, \sin^2\theta$, and the recovery fixes both branches. This is the standard fault-tolerance digitization argument, normally stated in words, here executed as one symbolic object. The ancilla-free degenerate stabilizer measurement (projector list, coherence-preserving, probe 5A) is itself a differentiator.
- Object-model property: projector-list QMOs, multi-pointer mixed controls, symbolic parameters.
- Audience: QEC/FT people; the consequential rung.
- Literature: Shor PRA 52, R2493 (1995); the digitization-of-errors argument (e.g. Nielsen-Chuang §10.2); deferred-measurement QEC.

### 4.3 Gate teleportation / magic-state injection (designed, unverified)
Inject $T$ via $|T\rangle = T|+\rangle$: CNOT data→ancilla, measure ancilla, apply $S$ to data conditioned on the record (`QuantumOperator["C"[QuantumOperator["S",{1}], {0}]]`). Then the *whole gadget compiled* equals $T$ exactly: `gadget["QuantumOperator"]` vs `QuantumOperator["T"]`, up to the dilation registers. The cleanest entry point to "why feedforward is the engine of fault tolerance."
- Literature: Gottesman-Chuang, *Nature* 402, 390 (1999); Bravyi-Kitaev, PRA 71, 022316 (2005).
- Risk: low; same verified primitives as 4.1.

### 4.4 Constant-depth long-range entanglement with feedforward (designed, unverified)
Reproduce the dynamic-circuits headline result in miniature: teleport a CNOT across a chain by consuming intermediate Bell pairs, measuring the middles, and feeding the records forward; the dynamic circuit has constant depth where the unitary circuit's depth grows with distance. In QF the identity "dynamic circuit $=$ long-range CNOT" is *provable* by compiling. Timely: this is what IBM markets as dynamic circuits, executed there with real-time classical electronics, expressible here as one coherent object.
- Literature: Bäumer et al., *PRX Quantum* 5, 030339 (2024) (CNOT teleportation over 101 qubits via 99 fed-forward measurements).
- Risk: moderate; register count grows with chain length ($n+k$ qubits), keep to 5-7 system qubits.

### 4.5 Measurement-based quantum computing, determinism from feedforward (designed, unverified)
A 1D cluster-state wire executing a programmed single-qubit rotation: each step measures in a rotated basis (QF: arbitrary-basis QMO) and the byproduct Pauli is corrected by a record-controlled gate. Show the output is **deterministic** with the corrections and a 50/50 mess without, symbolically in the measurement angle.
- Object-model property: arbitrary-basis mid-circuit measurement + adaptive corrections, the two halves of MBQC.
- Literature: Raussendorf-Briegel, PRL 86, 5188 (2001); Raussendorf-Browne-Briegel, PRA 68, 022312 (2003).
- Risk: moderate; adaptive *bases* (not just byproduct corrections) cannot depend on records in QF, only unitary corrections can. The one-bit-teleportation ladder formulation avoids basis adaptivity; use it.

### 4.6 Wigner's friend, executable (designed; primitives verified)
The friend's measurement is the dilation unitary $V$, kept coherent. Wigner can then (i) read the record (the classical story), (ii) **reverse the friend's measurement** ($V^\dagger$, verified exact in probe 6), or (iii) measure system+friend jointly in an entangled basis. Extended versions walk straight into the no-go-theorem scenarios. No classical-record platform can even *write down* option (ii) or (iii); QF runs them. The honest framing of the user's original instinct: QF measurement "explained properly" is precisely the observer-as-quantum-system model.
- Audience: foundations; the flagship "only here" rung.
- Literature: Wigner (1961); Deutsch, IJTP 24, 1 (1985) (in-principle reversal); Frauchiger-Renner, *Nat. Commun.* 9, 3711 (2018); Bong et al., *Nat. Phys.* 16, 1199 (2020) (local friendliness).
- Caveat to design around: pointer wires can't be re-measured, so Wigner's own final measurement must target positive wires; model the friend's lab on explicit positive wires and use the dilation for the friend's record only where it is acted on, or read Wigner's statistics from the final state.

### Dropped candidates
- **Iterative (Kitaev) phase estimation with qubit reuse**: the hardware win is reusing one ancilla; QF's registers accumulate instead, so the example would showcase the *cost*, not the feature.
- **Delayed-choice quantum eraser via pointer re-measurement**: blocked by the no-measurement-on-pointer-wires limitation; the eraser works fine with the marker on a positive wire, but then it is a standard two-qubit demo every platform has.

---

## 5. Assumptions surfaced

1. All examples stay at $\le 8$ total wires (system + pointers); the dilation cost is acknowledged in-text, not hidden.
2. "Unique" is claimed only for the coherent-record capabilities (§3 items 1-5), never for feedforward per se.
3. Degenerate stabilizer measurements always go through the projector-list constructor; the observable form over-measures (probe 3A) and must not appear in QEC examples.
4. Eigenvalue-to-pointer-slot ordering is ascending by eigenvalue (probe 3C: $-1$ before $+1$); control polarities in recovery circuits depend on it and must be re-checked per example.
5. The `Eigenorder` vs compiled-eigenqudits disagreement (§1) is cosmetic for these examples; revisit if an example reads `qc["Eigenorder"]` programmatically.

## 6. Next steps

- Verify 4.3 and 4.4 numerically (same primitives as the verified ones; low risk).
- Prototype 4.5 as a 3-rung one-bit-teleportation ladder with symbolic angle.
- Design the 4.6 notebook narrative (the no-go scenario needs careful staging to stay honest).
- Then promote the survivors into `QF-Showcase.md` / a notebook, with the §3 comparison table as the framing.
