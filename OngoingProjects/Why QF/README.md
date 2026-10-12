# Why QF: folder index

A map of what lives in this folder, grouped by effort and tagged by status. This file is the index; keep it in sync when files are added, moved, or change status.

## Orientation

Two efforts share this folder, and they take different positions on what makes QuantumFramework worth choosing.

- **Numerics-priorities research (Oct 2026): the current direction.** The case that the persuasive pitch is *verified numerics against incumbents' fast numerics*: match QuTiP, Qiskit, Stim and the rest on their own published benchmarks, then win on convergence and precision certification where their own trackers record failures. It measures QF's numerics honestly, including real defects.
- **The showcase (Jun to Aug 2026): earlier, finished.** The case that QF is distinctive because its objects are *exact and symbolic*. `QF-Showcase.md` is the flagship; it is QA'd and has open items (below).

If the two seem to disagree about whether "symbolic" is the edge, the October report is the newer position: it explicitly reframes the showcase's positioning and warns against an audience that "already files QF under symbolic work."

## Status tags

- **Active**: current, still being developed.
- **Research**: evidence base for a deliverable.
- **Deliverable**: a finished or near-finished artifact meant for others to read.
- **Process log**: decision or QA diary; provenance, not a plan.
- **QA**: an audit of a deliverable.
- **Held**: parked and verified, waiting for a future use.
- **Draft (pending)**: a proposed change not yet integrated; a decision is open.
- **Spent**: served its purpose; the work it fed has landed elsewhere.
- **Asset**: non-prose supporting file (figures, scripts, built notebooks).

## Numerics-priorities research (Oct 2026, current)

- **`reports/Quantum software numerics priorities.md`** (_Research_). The strategic case: compete on accuracy and convergence where incumbents chase throughput. Ranks twelve problem classes and names the top four to lead with (driven-dissipative Kerr, multilevel transmon with DRAG, small codes beyond Pauli noise, cavity-QED spectra).
- **`research_notes/Quantum software numerics priorities/`** (_Research_). Six primary-source notes behind the report, one per domain: `open_systems_libraries.md`, `gate_model_sdks.md`, `bosonic_cv_photonics.md`, `qec_and_tensor_networks.md`, `hardware_domain_tools.md`, `pain_points_and_benchmarks.md`.
- **`QF-Numeric-Inventory.md`** (_Active_). What QF can do numerically today, measured rather than claimed, scored against the report's ranking. Lists the blocking defects (dense Lindblad 30 to 250 times slower than incumbents; subsystem jump operators fail silently; non-prime Fock cutoffs segfault). Probes in `numeric-inventory-probes/`.
- **`QF-Capability-Gap-Research-2026-10-03.md`** (_Research_, untracked). Five scientifically useful calculations QF supports beyond stock Qiskit: symbolic logical-error channels, antibunching interference, Fubini-Study metric rank, Franck-Condon matrix elements, Bayesian entanglement posteriors.
- **`N17-Atomic-Revival-and-Wigner-Negativity.md`** (_Active_, draft essay, untracked). Computation-first essay asking how much cavity Wigner negativity the atomic population record can guarantee: the record sees only the photon-number distribution, so the guarantee is the phase-averaged field's negative volume. First concrete instance of the capability-gap thesis. Figures in `img/n17/`.

## Showcase deliverables (Jun to Aug 2026)

- **`QF-Showcase.md`** (_Deliverable_, done with open items). The seven-claim verified showcase of QF's symbolic and object-model capabilities. Open items: a citation to a non-existent companion `QF-Stabilizer-vs-Packages.md` (line 727) that backs the Stim numbers; the "Where This Stops" Abel-Ruffini boundary contradicts the Claim 2 five-site ring (a 5x5 operator with a closed-form spectrum); version provenance drift (header 2.0.0 vs a 2.1.x tree). See `QF-Showcase-Assessment.md`.
- **`QF-OpenSystems-Showcase.md`** (_Deliverable_, done). Six open-system worked examples positioned against QuTiP; each wins on the form of the answer, not speed.
- **`QF-Webinar-From-Qubits-to-Fock-Space.md`** (_Deliverable_, done). A compact webinar essay distilled from the showcase: seven transitions across circuits, open systems, basis changes, tensor networks and optics.

## Showcase process and QA

- **`QF-Showcase-Pipeline.md`** (_Process log_). The 105-row decision log for how the showcase examples were chosen and revised. Provenance, not a current plan.
- **`QF-Showcase-Assessment.md`** (_QA_, partly superseded). The August adversarial audit, findings F1 to F20. F3 and F4 (the string-shorthand "hangs") were later found non-reproducible on 2.1.0; F2, F5, F6, F15 and the overclaim findings are still open.
- **`QF-OpenSystems-Pipeline.md`** (_Process log_). Decision log for the open-systems showcase; also records kernel bugs worth upstreaming (symbolic Kossakowski rates, matrix-rate Sqrt, dropped DSolve options).

## Held, drafts, and spent

- **`QF-Showcase-Qudit-Extras.md`** (_Held_). Three verified qutrit and qudit examples lifted out of Claim 2 to keep it on the object model (discrete-Wigner magic, KCBS contextuality, single-qudit Grover). Parked for a future higher-dimensional section.
- **`QF-Showcase-TN-LightCone-Draft.md`** (_Draft (pending)_). Proposes replacing the Claim 4 Grover tensor-network example with a light-cone magnetization win. Not integrated; its "automatic locality" framing is deflated by the inventory (the cone is built by hand).
- **`QF-Showcase-Claim6-Feedforward-Implementation.md`** (_Spent_). Cold-start handoff to add the feedforward material; it landed as Claim 3 of the showcase. Residual value only in its deferred designs (constant-depth entanglement, MBQC, Wigner's friend).
- **`QF-Feedforward-Pipeline.md`** (_Spent_). The predecessor investigation to the handoff above; fully absorbed into it and the main pipeline.

## Assets

- **`QF-Showcase.nb`**, **`QF-Webinar-From-Qubits-to-Fock-Space.nb`** (_Asset_). Built notebooks generated from the matching `.md` sources; do not hand-edit.
- **`verify-showcase.wls`** (_Asset_). The showcase verifier script (an earlier environment; see the Assessment for provenance caveats).
- **`numeric-inventory-probes/`** (_Asset_). Runnable WL probes behind `QF-Numeric-Inventory.md`, plus `results.log` and `runprobes.sh`.
- **`img/showcase/`**, **`img/n17/`** (_Asset_). Figures for the showcase and the N17 essay.
