# Why QF: folder index

A map of what lives in this folder, grouped by effort and tagged by status. This file is the index; keep it in sync when files are added, moved, or change status.

```
Why QF/
├── img/          figures used by showcase/QF-Showcase.md
├── showcase/     the symbolic showcase (Jun to Aug 2026)
│   └── retired/  handoffs whose content has landed in the showcase
├── numerics/     the numerics-priorities research (Oct 2026)
│   ├── reports/
│   ├── research_notes/
│   ├── probes/
│   └── explicit-examples/
└── talks/        talk decks, scripts and plans (Oct to Nov 2026)
    ├── modeling-batches/
    └── QF-sweep/  (not tracked in git)
```

## Orientation

Two arguments for QuantumFramework live here, plus the talks that draw on both.

- **The showcase (Jun to Aug 2026): earlier, finished.** QF is distinctive because its objects are *exact and symbolic*. `showcase/QF-Showcase.md` is the flagship; it is QA'd and has open items (below).
- **Numerics-priorities research (Oct 2026): the newer position.** The persuasive pitch is *verified numerics against incumbents' fast numerics*: match QuTiP, Qiskit, Stim and the rest on their own published benchmarks, then win on convergence and precision certification where their own trackers record failures. It explicitly reframes the showcase's positioning and warns against an audience that "already files QF under symbolic work."
- **Talks (Oct to Nov 2026).** Slide decks built with md2nb from the `.md` files beside them, for the EU tour 2026, the WTC keynote, and general physicist audiences.

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

## numerics/

- **`reports/Quantum software numerics priorities.md`** (_Research_). The strategic case: compete on accuracy and convergence where incumbents chase throughput. Ranks twelve problem classes and names the top four to lead with (driven-dissipative Kerr, multilevel transmon with DRAG, small codes beyond Pauli noise, cavity-QED spectra).
- **`reports/Lanyon AI profile and quantum fit.md`** (_Research_). Profile of Lanyon AI; its finding is that Lanyon proves stencil algebra, not convergence.
- **`research_notes/`** (_Research_). Primary-source notes behind the two reports, one folder per report.
- **`QF-Numeric-Inventory.md`** (_Active_). What QF can do numerically today, measured rather than claimed, scored against the report's ranking, with the blocking defects it found. Probes in `probes/`.
- **`QF-Numerics-20-Questions.md`** (_Research_). Twenty quantum questions chosen to develop with QF, scored for payoff, numerical substance and scope.
- **`QF-Capability-Gap-Research-2026-10-03.md`** (_Research_). Five scientifically useful calculations QF supports beyond stock Qiskit: symbolic logical-error channels, antibunching interference, Fubini-Study metric rank, Franck-Condon matrix elements, Bayesian entanglement posteriors.
- **`explicit-examples/`** (_Active_). Worked numeric essays (cat-state Wigner function, damped Jaynes-Cummings, driven-dissipative Kerr, lossy NOON metrology, N17 atomic revival and Wigner negativity), each with a `backbone-*.md` statement layer and its graph, plus `TransmonExample/` (transmon gates in three versions and the Qolab hardware cases). The N17 essay links figures in `img/n17/`, which does not exist yet.
- **`probes/`** (_Asset_). Runnable WL probes behind `QF-Numeric-Inventory.md` (`runprobes.sh`, `results.log`), and the short-pulse comparison against QuTiP (`shortpulse-head-to-head.md` with its scripts and tables; the raw `shortpulse_data*/` folders are not tracked).

## showcase/

- **`QF-Showcase.md`** (_Deliverable_, done with open items). The seven-claim verified showcase of QF's symbolic and object-model capabilities. Open items: a citation to a non-existent companion `QF-Stabilizer-vs-Packages.md` that backs the Stim numbers; the "Where This Stops" Abel-Ruffini boundary contradicts the Claim 2 five-site ring (a 5x5 operator with a closed-form spectrum); version provenance drift. See `QF-Showcase-Assessment.md`.
- **`QF-OpenSystems-Showcase.md`** (_Deliverable_, done). Six open-system worked examples positioned against QuTiP; each wins on the form of the answer, not speed.
- **`QF-Webinar-From-Qubits-to-Fock-Space.md`** (_Deliverable_, done). A compact webinar essay distilled from the showcase: seven transitions across circuits, open systems, basis changes, tensor networks and optics.
- **`QF-Showcase-Pipeline.md`** (_Process log_). The decision log for how the showcase examples were chosen and revised.
- **`QF-Showcase-Assessment.md`** (_QA_, partly superseded). The adversarial audit, findings F1 to F20; its status line records which findings were rechecked and on which commit.
- **`QF-OpenSystems-Pipeline.md`** (_Process log_). Decision log for the open-systems showcase; also records kernel bugs worth upstreaming.
- **`QF-Showcase-Qudit-Extras.md`** (_Held_). Three verified qutrit and qudit examples lifted out of Claim 2 (discrete-Wigner magic, KCBS contextuality, single-qudit Grover).
- **`QF-Showcase-TN-LightCone-Draft.md`** (_Draft (pending)_). Proposes replacing the Claim 4 Grover tensor-network example with a light-cone magnetization win. Not integrated; the numeric inventory found the cone is built by hand, which weakens its "automatic locality" framing.
- **`verify-showcase.wls`** (_Asset_). The showcase verifier script.
- **`QF-Showcase.nb`**, **`QF-Webinar-From-Qubits-to-Fock-Space.nb`** (_Asset_). Built from the matching `.md` sources; do not hand-edit.
- **`retired/`** (_Spent_). `QF-Showcase-Claim6-Feedforward-Implementation.md`, the handoff that landed as Claim 3 of the showcase (its deferred designs, constant-depth entanglement, MBQC and Wigner's friend, are the only residual value), and `QF-Feedforward-Pipeline.md`, its predecessor investigation.

## talks/

Each deck `X.md` builds to `X.nb` with MarkdownToNotebook; the build line is the first comment in the `.md`.

- **EU tour 2026** (_Active_). `EU-Tour-2026-Showcase.md` (the deck), `EU-Tour-2026-Core-Selection.md` (the 33 examples that cover every stop), `EU-Tour-2026-Running-Order.md` (slides per stop).
- **Modeling quantum physics** (_Active_). `Modeling-Slides.md` (the deck); `Modeling-Audit-Opinion.md` and `Momentum-Composition-Audit.md` (audits of its slides, with the two Momentum figures); `TN-Slide-Candidates.md` (tensor-network slides not yet merged); `modeling-batches/` (Batch-A to Batch-V reviews and `Findings-outside-the-deck.md`).
- **Continuous measurement** (_Active_). `Continuous-Measurement-Slides.md` and two talk scripts, full length and 20 minutes.
- **QF sweep over the decks**. `QF-issues-from-slides.md` (_Active_, temporary): QF issues found while revising the Modeling and Continuous Measurement decks, to be filed as QF issues. `QF-sweep/` holds the batch prompts, proposals and backups of that sweep; it is not tracked in git.
- **General physicist talk**. Three decks, `Talk-Slides.md`, `Talk-Slides-v3.md` (loads `TalkSetup.wl`) and `Talk-Slides-50min.md`, with `Talk-Outline-50min.md`, `Talk-Topic-Summaries.md`, `Talk-Source-Index.md`, `Talk-Catalogue-Unique-WL-QF-Examples.md`, `Catalogue-Examples.md` and `Talk-Plan-JC-Revival-and-Wigner-Negativity.md`.
- **WTC keynote**. `WTC-Keynote-Examples.md`: bosonic algebra, stabilizers and interferometers.
