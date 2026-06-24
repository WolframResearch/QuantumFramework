# Stabilizer: project map and status

This folder is the home of two things that grew together and are now tangled in one flat
directory:

1. The **QuantumFramework stabilizer subsystem** itself: the shipped paclet code
   (`PauliStabilizer`, `StabilizerFrame`, `CliffordChannel`, `GraphState`,
   `LocalComplement` and their methods, living in `QuantumFramework/Kernel/Stabilizer/`),
   together with its design audits, API reference, test battery, and the pedagogical
   tutorial built on top of it.
2. An **open research program** asking how far the symbolic, exact-arithmetic strengths of
   the Wolfram Language can push *past* the Clifford boundary: stabilizer frames for magic,
   the SymPhase symbolic-sign idea, and the symbolic path-sum simulator for Clifford+$T$.

Concept first: the stabilizer formalism describes a quantum state not by its $2^n$
amplitudes but by the $n$ commuting Pauli operators that fix it, an $O(n^2)$ description that
makes thousand-qubit Clifford circuits, error-correcting codes, and graph states tractable.
The shipped subsystem implements that (the Gottesman-Knill world) and is feature-complete and
benchmarked against Stim and QuantumClifford.jl. The research program is about the part the
tableau *cannot* hold: a single $T$ gate turns one Pauli into a sum, magic enters, and the
cost reappears as a $\#\mathrm{P}$-hard Gauss sum. The open question is whether carrying that
sum *symbolically* (a phase polynomial over $\mathbb{F}_2$ valued in $\mathbb{Z}_8$, reduced
by rewriting rather than enumerated) buys anything no numeric simulator has. The honest
verdict so far: it does not break the wall, but exact **symbolic-parameter** path sums are a
genuine WL-only capability, and the Gröbner/variety point count is an unbuilt realization of
an open problem.

---

## The sub-projects and their status

| # | Project | What it is | Status |
|---|---------|-----------|--------|
| 1 | **Kernel subsystem** | The 16-file `Kernel/Stabilizer/` paclet code: tableau states, frames, Clifford channels, graph states, symbolic + Pauli-string measurement, random Clifford, inner products, entropy. | **Shipped.** Phases 1-8 done; 5 public symbols + method-grade ops; ~965 stabilizer tests green. 3 roadmap items still open (optimal synthesis, VOP tracking, `CliffordTableau`/auto rank decomposition). |
| 2 | **Performance rework** | Bit-packed gates + a Compile-to-C generator-major gate fold; packed measurement; exact-rational Mallows sampler. | **Shipped to kernel.** 266-1095x speedups; beats QuantumClifford.jl at $n=1000$, within 2.5-5.3x of Stim. Full report lives in `../Platform Comparison/`. |
| 3 | **Tutorial** | `Stabilizer-Formalism-By-Computation.md/.nb`: a computation-first tour, single qubit to 1000-qubit GHZ, codes, graph states, magic, channels, with an honest Stim / QC.jl comparison. | **Committed, in polish.** `.md` + built `.nb`; both modified in the working tree. |
| 4 | **Formula verification** | `Formula_Test/`: 155 strict `VerificationTest`s mapping 1:1 to formulas distilled from 44 papers; surfaced and fixed 3 real kernel bugs (F1/F2/F3). | **Green (155/155).** Mirrored into `Tests/Stabilizer/Formulas.wlt`. |
| 5 | **Non-Clifford blow-up** | `StabilizerFrame-NonClifford-Blowup.md`: verifies the exact $2^k$ frame blow-up per $k$ non-Clifford gates, and a phase-coherent dense-materialization **fix**. | **Verified; fix shipped to kernel** (`61bc0e39`). `RY`/`RZ` still return `$Failed`; half-angle convention caveat noted. |
| 6 | **Frame observables** | `StabilizerFrame-Observables-Design.md` + `stabilizerframe-observables-scripts/`: design + out-of-tree prototype for native `<P>`/Born/inner-product on a frame at cost $\chi^2\mathrm{poly}(n)$ instead of $2^n$. | **Design + verified prototype; kernel change deferred.** Not yet committed. |
| 7 | **SymPhase vs rank** | `SymPhase-vs-StabilizerRank-Algebraic-Branching.md`: where the symbolic-sign trick transfers from measurement to $T$-gate branching (answer: only the diagonal Clifford+$T$ fragment). | **Analysis complete; no kernel change.** |
| 8 | **Beyond SymPhase / magic** | `Beyond-SymPhase-Symbolic-Algebra-For-Magic.md` + `pathsum-build/`: the literature landscape + a working WL **symbolic path-sum simulator** whose Hadamard-elimination folds all Clifford branching to closed form and leaves the magic as a residue. | **Prototype built + verified** against dense amplitudes. Staged build plan (rank cap, phase-poly compression, full confluent rewriter) not yet in kernel. |
| 9 | **Novelty re-examination** | `WL-Symbolic-Novelty-Reexamination.md` + `wl-symbolic-toolkit-docs/` (91 WL doc pages) + `novelty-demos.wls`: a skeptical audit of whether a WL-unique method was missed. | **Verdict delivered.** No $\#\mathrm{P}$ breakthrough; WL-unique capability = free-symbol path sums; Gröbner/variety direction = unbuilt realization of an Amy-Stinchcombe open problem. Not committed. |
| 10 | **Doc-page drafts** | `kernel-subsystem/doc-page-drafts/`: 10 draft reference-page `.nb`s for the public symbols, plus the `build_*.wl` generators. | **Draft, ready to promote.** The 3 kernel bugs (`A.11/A.12/A.13`) that blocked them are now fixed; pages can move to `QuantumFramework/Documentation/`. |
| 11 | **Literature library** | `literature/`: `paper-bibliography.md` + `paper-fetch-report.md` + `tex/` (44 paper sources) + `pathsum-literature/` (14 path-sum papers) + `External Packages/` (Stim, QC.jl clones). | **Reference inputs.** `tex/`, `External Packages/`, and `pathsum-literature/` are all gitignored. |

**Design-driver documents** (inputs that shaped project 1, now mostly historical, all under
`kernel-subsystem/`): `package-design-synthesis.md` (the 28-paper distillation that set the API),
`external-packages-audit.md` (Stim + QC.jl as upgrade material), `paulistabilizer-source-audit.md`
(pre-rebuild audit, marked historical), and the `api.md` / `roadmap.md` / `postmortem.md` triple.

---

## Where to start reading

- **"What can the subsystem do?"** -> `kernel-subsystem/api.md`, then `tutorial/`.
- **"What is left to build?"** -> `kernel-subsystem/roadmap.md` (open items A.5, A.6, A.8; deferred B.1-B.13).
- **"What is the research frontier?"** -> `magic-and-pathsums/WL-Symbolic-Novelty-Reexamination.md` (the verdict),
  then `magic-and-pathsums/Beyond-SymPhase-Symbolic-Algebra-For-Magic.md` and `magic-and-pathsums/pathsum-build/PathSum-WL-Note.md`.
- **"Is it correct?"** -> `Formula_Test/` and `Tests/Stabilizer/` in the paclet.

---

## Folder layout

Five concern-based groups (the folder was previously flat, ~20 top-level files plus subdirectories
mixing all five). `Formula_Test/` keeps its name and place because `Tests/Stabilizer/` references it
by path; the gitignored corpora (`tex/`, `External Packages/`, `pathsum-literature/`) live under
`literature/` and stay ignored.

```
Stabilizer/
├── README.md                       # this map
├── .gitignore                      # ignores tex/, External Packages/, eprint/, pathsum-literature/
│
├── kernel-subsystem/               # docs + audits for the shipped paclet code
│   ├── api.md  roadmap.md  postmortem.md
│   ├── package-design-synthesis.md
│   ├── external-packages-audit.md
│   ├── paulistabilizer-source-audit.md         (historical, pre-rebuild)
│   └── doc-page-drafts/                         # 10 reference-page .nb + their build_*.wl generators
│
├── tutorial/
│   └── Stabilizer-Formalism-By-Computation.{md,nb}
│
├── Formula_Test/                   # 155-test formula battery (referenced by Tests/Stabilizer/)
│
├── magic-and-pathsums/             # the OPEN research cluster (beyond Clifford)
│   ├── StabilizerFrame-NonClifford-Blowup.md
│   ├── StabilizerFrame-Observables-Design.md
│   ├── SymPhase-vs-StabilizerRank-Algebraic-Branching.md
│   ├── Beyond-SymPhase-Symbolic-Algebra-For-Magic.md
│   ├── WL-Symbolic-Novelty-Reexamination.md   novelty-demos.wls
│   ├── pathsum-build/              # the WL path-sum prototype + benchmarks + tests
│   ├── wl-symbolic-toolkit-docs/   # 91 WL doc pages + the toolkit-mapping note
│   └── scripts/  symphase-vs-rank-scripts/  beyond-symphase-scripts/  stabilizerframe-observables-scripts/
│
└── literature/                     # reference inputs (mostly gitignored)
    ├── paper-bibliography.md  paper-fetch-report.md
    ├── tex/                  (gitignored)
    ├── pathsum-literature/   (gitignored)
    └── External Packages/    (gitignored)
```

Each research report keeps its own `*-scripts/` folder as a sibling inside `magic-and-pathsums/`,
so every report-to-script link still resolves. Still untracked and worth committing once settled:
the new research outputs (`StabilizerFrame-Observables-Design.md`, `WL-Symbolic-Novelty-Reexamination.md`,
`novelty-demos.wls`, `pathsum-build/`).
