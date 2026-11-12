# Talk plan: collapse, revival and Wigner negativity, live in QF

Built from a full read of every notebook source delivered in the last two months (40 files, listed at the end) and two kernel probes run on 2026-10-28 against QF HEAD. Assumes a 45-minute slot; the segment times scale.

## One decision before anything else

The abstract asks whether the cavity field's Wigner function has negative regions at the moment the atomic oscillation returns. On 2026-10-05 the N17 essay dropped that exact framing ("revival vs negativity at matched times") as loosely posed and replaced it with: given only the atomic record z(t), how much negativity is guaranteed? The plan below keeps the abstract's question as the hook, answers it by computation, and then uses N17's result as the payoff: the atom's record cannot see the negativity it is being asked about. If you would rather end on the yes/no alone, drop segment 12.

## Every claim in the abstract already has a QF cell behind it

| Abstract claim | Where the cell lives | Status |
|---|---|---|
| Lossy atom-field state, zero-loss check against exact JC | N17 essay (QuantumEvolve with jump `a`, inversion vs Σ pₙ cos(2√(n+1)τ)) | delivered, evaluated 2026-10-05 |
| Exact JC solution as a QF object | QO book ch4: interaction-picture `MatrixExp`, g and t symbolic, cutoff 50 | delivered, untimed |
| QAOA cost landscape in closed form | QF-Showcase Claim 1 (K4, `Maximize`); UMD demo (5-ring) | probe: triangle landscape in 2.3 s, correct at γ=β=0 |
| Two coupled spins, entropy vs time | QF-Showcase Claim 1 (Heisenberg, `QuantumEvolve` + `QuantumPartialTrace`) | probe: binary entropy of cos²2t in 0.17 s |
| Same objects reach Fock space | Showcase coherent-state evolution; QO ch2; boson-sampling Part I; interferometer showcase | delivered |
| Channels with a surviving parameter | Showcase amplitude damping and symbolic Lindblad; QI `"Loss"[η]` | delivered |

What does NOT back the abstract: the open-systems catalog (20 entries), the continuous-QEC essay, the syndrome derivation, record analysis and feedback, and the two-examples essay are all plain WL matrices with no QF object. They make good physics but cannot be shown as evidence that QF composes. The only QF notebook in that folder family is the Heron dissipative-entanglement one.

## Running order

### A. The question (3 min)
Pose the 1980 result and the question. Show the lossless inversion collapsing and returning (QO ch4 plot or N17 `same-record.png`), precomputed.

### B. Symbolic parameters survive composition (12 min)
1. **A qubit family keeps its symbols.** Showcase: `QuantumState[{Cos[α], Sin[α] Exp[I β]}]["BlochVector"]`. Instant.
2. **Two coupled spins.** Showcase cells `hheis`, `psit = QuantumEvolve[hheis, QuantumState["01"], t]`, `s = FullSimplify[QuantumPartialTrace[psit, {2}]["VonNeumannEntropy"], ...]`. Then differentiate `s` and solve for the time of maximal entanglement: the abstract's "differentiate" in one line.
3. **QAOA landscape as a formula.** Showcase K4 build (`stateprep`, `qcost`, `mixer`, `meanvalue = oneLayer /* hcost /* oneLayer["Dagger"]`, `landscape`, `Maximize`). Close with one numeric spot-check of `landscape` against `meanvalue` at random angles. Fallback if K4 FullSimplify is slow: the triangle version from the probe.
4. **A channel is a formula too.** Showcase `QuantumChannel["AmplitudeDamping"[γ]][psi]["Purity"]`, then the symbolic qubit Lindblad for every initial ρ and both rates, ending in `Reduce` giving T2 ≤ 2T1. This is the bridge to loss.

### C. The same objects in Fock space (8 min)
5. **Coherent state under the oscillator**, α, ω, t symbolic: `uosc @ CoherentState[] == E^(-I ω t/2) CoherentState[][\[FormalAlpha] E^(-I ω t)]` (Showcase, QO ch2). Say "closed form at fixed cutoff": the norm is the truncated exponential series.
6. **A circuit of Fock-space operators with a symbolic angle.** Either the HOM coincidence with polarization overlap χ (boson-sampling Part I, χ→0 and χ→π/2 limits) or the MZI with `PhaseShiftOperator[φ]` (QO ch6). Use `Method -> "Recurrence"`; the default `"MatrixExp"` beam splitter with too few levels freezes |1,1⟩.
7. **Loss with symbolic η.** Interferometer showcase §12: `QuantumInterferometer[{"Loss"[η]}, "Parameters" -> {η}]`, the binomial loss law for symbolic n and k, lossy HOM survival η².
8. **Negativity in one picture.** Cat-State-Wigner: even and odd cats, origin value locked to parity (+1/π, -1/π), coherent state never negative.

### D. One atom, one mode, photon loss (15 min)
9. **One excitation, exactly solvable with loss.** Damped-JC: the 2×2 no-jump generator's eigenvalues give strong coupling, crossover and Purcell; QF's `"Hamiltonian"[H, {σminus, a}, {γ, κ}]` evolution lands on the exact curves. This is "test numerics against a formula" with loss already in.
10. **Coherent field: build and check.** N17 cells: `h = sigmaMinus[SuperDagger[a]] + SuperDagger[sigmaMinus][a]`, excitation number commutes with `h`, the jump lowers it by one. Lossless record against the exact sum `design . photonProbabilities`. This is the abstract's zero-loss check.
11. **Turn on loss and look at the field.** N17 `trajectories` at a few κ/g, field Wigner snapshots near half the revival and at the revival, negative volume vs time. Answer the abstract's question from the plot, not from memory. Pick κ/g small: QO ch8's exercise (γ1 = γ2 = 0.5 at g = 10) wipes the oscillations out entirely.
12. **The twist.** The phase-randomized field gives the identical atomic record at every loss rate, yet its negativity is only the ring part. So the atom certifies only the phase-averaged negative volume; the rest, most of it near half the revival, is invisible to the record (N17 final section).

### E. Where it stops (3 min)
Closed forms hold while the reduced block stays small (Showcase "Where This Stops"); the lossy JC field is numeric. Dense Lindblad is the slow path. A measurement record is a random draw, so it is never symbolic. Name these before the audience does.

## Not in the talk
The named-circuit tour (it was the Fujitsu talk). The TransmonLab notebooks (not QF; internal contradictions in Case 18). DRAG: every result is numeric, and the only closed forms are hand-derived in QolabHardwareCases.md, so it supports none of the symbolic claims. The QEC tech notes (QECCore is unshipped). The open-systems essays, except as an answer to "why not QF there?".

## Preparation
1. Assemble one talk `.md` (md2nb) from the cells above; run every `wl` cell through `wolframscript` on the talk machine with QF HEAD and record its time.
2. Precompute and load from disk: N17 trajectories, the ch4 symbolic `U` at cutoff 50, the K4 `landscape` if it exceeds a few seconds. `WignerRepresentation` costs about a second per state; use N17's angular-harmonic Wigner (milliseconds per state) for the negative-volume curve and QF's representation only for the cross-check and snapshots.
3. Check before stage: the Showcase `QuantumOperator["ZZ"/2, #]` (a string divided by 2); QO ch5's `QuantumState["1", atomBasis]` is |g⟩, not the |e⟩ the text says; QF-Showcase is stamped QF 2.0.0 while the folder runs 2.1.x.
4. Keep off screen: any cell with `Quiet`; the default `QuantumDistance` "Fidelity" (it is a distance); `"Probability"` on truncated states (it renormalizes and hides lost norm); one-element Kraus lists with an explicit order (the channel lands on the wrong wire).
5. Do not say on stage: "QF checks physics" (it checks shapes; normalization and trace preservation are not enforced); "the formal semantics covers Fock space" (it does not).

## Sources read in full
N17, Damped-Jaynes-Cummings, Cat-State-Wigner, Lossy-NOON, Driven-Dissipative-Kerr, QF-Showcase, Why QF README, QF-Numerics-20-Questions; TransmonGates V1-V3, QolabHardwareCases, TransmonGates revision notes, StatePreparation_LBC, Case18; open-system-simulation-catalog, three-qubit continuous QEC, syndrome populations, pure-dephasing clock, two examples from minimal assumptions, record analysis and feedback, dissipative-entanglement Heron, EQTC2026 abstract; interferometer showcase, boson-sampling, why boson sampling is hard, atom-cavity mean field, named-circuits catalogue, Fujitsu talk; StabilizerCodes, FaultTolerantGadgets, QECCoreInternals, QF formal semantics, structured arrays, continuous-space Solutions 1-3, UMD QC-for-engineers, railway-track distance; Quantum Optics book chapters 1-8.
