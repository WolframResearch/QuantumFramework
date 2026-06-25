# Pipeline & Brief: "Quantum in Discrete Space" QF comprehensive-exam manual

This file is the reusable brief and the build tracker. The deliverable is
`Quantum-in-Discrete-Space-QF-Comprehensive-Exam.md`.

## DRAFT BRIEF

- **Goal:** A self-learning manual structured as a comprehensive exam. The reader, working
  through it with a kernel open, should be able to *do* each task in QuantumFramework (QF):
  prepare states, build operators/channels, measure, quantify entanglement, simulate open
  systems, run the standard algorithms. After each question the reader can reproduce the
  worked QF solution and see the verified output.
- **Reader:** Advanced BSc / early MSc student in quantum physics. Strong linear algebra and
  quantum mechanics; little or no QF experience. Wants "how do I compute X", not a proof course.
- **Claim (spine):** Every concept in finite-dimensional quantum theory has a short, exact,
  executable QF realization; if you can compute it you understand it.
- **Genre/length:** One large Markdown file, ~80 question/answer items, each = clearly stated
  question, concept gloss, paste-safe QF code, embedded verified output, one-line interpretation.
- **Scope.** IN: finite-dimensional ("discrete space") quantum theory: qubits and qudits,
  observables, measurement (projective + POVM), composite systems and entanglement, mixed
  states/channels, multi-qubit gates and circuits, stabilizer formalism + QEC, discrete phase
  space, the standard algorithms, and information-theoretic limits. Plus a short CV/Fock coda
  using the SecondQuantization sub-package. OUT: full quantum field theory, hardware-vendor
  plumbing (IBM/AWS job submission), deep performance tuning (those live in other docs).
- **Voice:** learning-by-computing (concept -> computation -> interpretation, short paragraphs,
  first-person plural, every assertion earned by a cell) layered with quantum-review pedagogy
  (physically motivated, operationally precise, careful about phase/gauge/basis subtleties).
- **Verification:** every code cell is run in a live Wolfram kernel against the QF working tree
  (`PacletDirectoryLoad`); the embedded output is the real returned value. Final pass = `/wl-verify`.
- **Hard constraints (user standing rules):** no em dash anywhere; valid delimited TeX for all
  math in prose; code is paste-safe text WL (NO front-end-only `^\[Dagger]` superscripts or
  `|psi>`/`<psi|` glyphs inside code; use `qs["Dagger"]`, `["Scalar"]`, etc.).

## Pressure-test (skeptic's pass)

- *Weakest assumption:* that QF has a clean one-liner for every item. Mitigation: a few items
  (HHL, Simon, QAOA, mana) are heavier; each is verified live and trimmed to what actually runs,
  with the honest caveat stated rather than faked.
- *Strongest objection:* "this is just a function catalog." Answer: each item leads with the
  physics question and ends with interpretation; the QF call is the middle, not the point.
- *Scope vs length:* 80 items is large but each is short; the exam framing justifies breadth.

## Anchor / provenance

- QF working tree HEAD `d0e71690`; audit anchor `b4fcdb31` (only delta is `SymbolicMeasure.m`,
  no `PacletInfo` symbol change). Paclet `Wolfram/QuantumFramework 2.0.0`.
- Paste-safe idioms verified live before authoring: bra `qs["Dagger"]`; 0-qudit scalar
  `["Scalar"]`/`["Number"]`; expectation `QuantumMeasurementOperator[A][qs]["Mean"]` or
  `Tr[A["Matrix"].qs["DensityMatrix"]]`; `qs["VonNeumannEntropy"]` returns a `Quantity[_,"Bits"]`.

## Taxonomy (target ~80 items, A-M)

- **A. Discrete state space & the qubit** (6): amplitudes/Born, normalization & global phase,
  Bloch vector, Bloch-angle prep, overlap & |overlap|^2, qudit prep + uniform superposition.
- **B. Observables & single-qubit operators** (7): Pauli algebra, expectation (4 ways),
  eigen-decomposition, rotations + U3, unitary/Hermitian checks, ZYZ decomposition, commutator+uncertainty.
- **C. Measurement** (6): projective outcomes+probs, non-computational basis, post-measurement
  states, POVM/SIC, finite-shot counts, mid-circuit measurement + feedforward.
- **D. Composite systems & entanglement** (7): tensor product, Bell/GHZ/W, partial trace,
  entanglement test (PPT), entropy/concurrence/negativity, Schmidt, reduced-state no-signalling.
- **E. Mixed states & density matrices** (5): ensemble + purity, von Neumann entropy + caveat,
  Werner threshold, fidelity + trace distance, purification.
- **F. Multi-qubit gates, circuits, universality** (6): CNOT/CZ/SWAP/Toffoli/controlled-U,
  compose + draw, depth/layers, overall unitary, KAK decomposition, universality demo.
- **G. Qudits & higher-dim discrete** (5): clock/shift Weyl-Heisenberg, qutrit gate + SUM,
  spin-j operators, MUB/SIC, KCBS contextuality.
- **H. Channels & open systems** (7): Kraus channel, Choi + complete positivity, named noise,
  Stinespring, Lindblad (T1/T2), Liouvillian + steady state, entanglement sudden death.
- **I. Stabilizer formalism & QEC** (7): tableau, Clifford updates, large Clifford circuit,
  graph/cluster state, 3-qubit code + syndrome, 5-qubit/Steane stabilizers, magic / T gate.
- **J. Discrete phase space** (3): discrete Wigner function, Wigner negativity, mana.
- **K. Quantum algorithms** (12): Deutsch-Jozsa, Bernstein-Vazirani, Simon, Grover, QFT,
  phase estimation, HHL/linear-solve, VQE, QAOA, teleportation, superdense coding, ent. swapping.
- **L. Foundations & info limits** (4): no-cloning, CHSH Tsirelson bound, Holevo, GHZ paradox.
- **M. CV/Fock coda** (5): Fock + ladder operators, coherent state, displacement/squeezing,
  beam splitter / Hong-Ou-Mandel, Jaynes-Cummings.

## Build status

| Part | Items | Harness | Verified | Prose written |
|------|-------|---------|----------|---------------|
| A | 6 | yes | yes | yes |
| B | 7 | yes | yes | yes |
| C | 6 | yes | yes | yes |
| D | 7 | yes | yes | yes |
| E | 5 | yes | yes | yes |
| F | 6 | yes | yes | yes |
| G | 5 | yes | yes | yes |
| H | 7 | yes | yes | yes |
| I | 7 | yes | yes | yes |
| J | 3 | yes | yes | yes |
| K | 12 | yes | yes | yes |
| L | 4 | yes | yes | yes |
| M | 5 | yes | yes | yes |

All 80 items authored with live-verified embedded outputs. **/wl-verify: OPEN ISSUES 0**
(2 rounds; round 1 reproduced all ~80 cells and refuted only F5's abbreviated KAK label;
round 2 confirmed the F5 correction `{{GlobalPhase,1},{U,9},{CX,2}}` and H5/H6/I7 regressions).
Deliverable: `Quantum-in-Discrete-Space-QF-Comprehensive-Exam.md`. Embedded output blocks were
then dropped on request (81 removed). **Then all answers were rewritten** in the minimal WL+QF
idiom of `introduction-to-QIS-revised.md`: most answers are now direct Wolfram Language on matrices
and vectors (`PauliMatrix`, `TensorProduct`, `KroneckerProduct`, `MatrixExp`, `DSolveValue`,
`FourierMatrix`, `Eigenvalues`, `==`/`FullSimplify`), with QuantumFramework kept only where it is
the most minimal tool: stabilizer formalism (Part I), discrete Wigner/mana (Part J), circuit
depth/KAK/decompose (F3/F5/F6), and many-wire circuit protocols (K3, K10-K12). A reusable toolkit
(Pauli vector, blochVector, densityMatrix; reduceA/reduceB/partialT/vne) is defined once and reused.
Full manual re-extracted and run top-to-bottom against the working tree: 84/84 cells, exit 0, empty
stderr. Per-part rewrite harnesses at `verify/rw-A.wls`..`rw-LM.wls`. Now ~1613 lines.
