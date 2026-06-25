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

**v4 (curriculum question-list revision, adversarial-draft Revise track):** the user raised concerns
about the *list of questions* itself and asked, as a professor wanting students to learn finite-dim
quantum in WL at BSc+MSc level, to revise it for complete topic coverage. Extracted all 80 questions
to a new sibling artifact `Question-List.md`, then ran the adversarial-draft loop on the curriculum.
Decisions (AskUserQuestion): drop the bosonic Part M (it is infinite-dimensional, violating the one
hard constraint; spin coherent states + discrete Wigner are the finite-dim analogs); target a
comprehensive ~130-160 set; tag every item [BSc]/[MSc]; quantum-review rigor lens. Rebuilt as **171
questions across 21 parts (71 BSc / 100 MSc)** with a coverage table. Two fresh-context review rounds:
round 1 found 38 issues (8 physics mis-statements: HHL inverts not applies, Gleason needs d>=3,
five-qubit code is non-CSS / saturates the quantum Hamming bound, Mermin-Peres = state-INDEPENDENT KS
contextuality, Berry phase as discrete Bargmann invariant, bound entanglement needs 3x3/2x4, CKW =
tangle, Stern-Gerlach = spin-projection; plus redundancies, 8 canonical gaps now filled: perturbation
theory both kinds, variational principle, interaction picture + optical Bloch, symmetry/conservation
block-diagonalization, identical-particle exchange symmetry, magic-state distillation, decoherence/
pointer basis, even-d Wigner, KS uncolorable ray set; replaced a bad Holevo-vs-measurement-problem
item with PBR; split TP/Stinespring and parameter-shift/barren-plateau; PPT re-tagged MSc; Part 8
reordered BSc-first). Round 2 found only 2 (Fermi golden rule needs the dense-band limit, not a finite
level set; minor adiabatic shorthand), both fixed. Converged.

**v4b (quantum-review clarity + granularity pass):** audited all questions for clarity/
unambiguousness (quantum-review skill). 12 IMPRECISE prompts disambiguated (1.3 "which system's
density matrix", 9.2 "read the Bloch ball", 10.7 "the quantities above" + which capacity, 15.6
"T gate leaves the Clifford class", 16.5/16.8 vague verbs, 21.9 PBR made computational, etc.; no
physics errors). Then, per user policy "split distinct-concept pairs list-wide" (keep single-thread
bundles like the inner-product chain, instances like the four noise channels, and build-and-verify
pairs), split 8 genuinely fused questions: 1.1 (state rep | Born rule), 1.2 (normalization | global
phase), 2.11 (PT series | avoided crossing), 4.3 (driven Rabi | free Larmor), 14.2 (SUM gate |
qudit Fourier), 16.6 (code distance | logical operators), 19.6 (amplitude amplification |
estimation), 21.1 (no-cloning | no-deleting). Renumbered all affected parts. **Now 179 questions /
21 parts / 76 BSc / 103 MSc** + coverage table; sequential numbering verified, no dashes. (My rough
"~200-210" split estimate was high; most multi-output questions are single concepts the policy
keeps, so the principled count is 179.)

**v4c (completeness + dependency-ordering audit, adversarial-draft + fresh reviewer):** user required
(1) concept completeness and (2) a valid prerequisite chain (no question needs a concept introduced
by a LATER question). Fresh-context audit found 10 forward-reference violations (the keystone:
density operator defined in old Part 9 but used from Part 1; same for tensor product, multi-qubit
gates, projective measurement, commutator, vN/Shannon entropy, channels-before-info) and 12
completeness gaps. Restructured into a dependency-correct **25-part** order: hoisted density operators
(new Part 3), tensor product + partial trace (Part 4), elementary circuits (Part 10) to the front;
moved the commutator into Part 2; basic projective measurement (Part 5) before spin; channels (Part
13) before quantum-information/entropy (Part 14). Added the 12 missing concepts (pure-vs-mixed test,
maximally mixed state, convexity/extreme points, ensemble ambiguity, Gibbs state, Uhlmann fidelity,
no-broadcasting, Zeno, time-energy/Mandelstam-Tamm, three-picture equivalence, quantum trajectories,
Hilbert-Schmidt inner product). **Now 194 questions / 25 parts / 83 BSc / 111 MSc** + dependency-
ordered coverage table. Self-run fresh dependency re-read (Agent tool was down) caught one
reorder-induced forward reference: the Zeno effect needs unitary evolution (Part 7), so moved it from
Part 5 to 7.15; confirmed variational/perturbation items (2.11-2.13) only need H as a Hermitian
operator (Part 2), not dynamics. Numbering verified sequential 1..N in all 25 parts, no dashes.

NOTE: the worked-answer manual
`Quantum-in-Discrete-Space-QF-Comprehensive-Exam.md` still reflects the OLD 84-cell A-M set; if it is
rebuilt to the new curriculum, Part M must be dropped and the ~90 new finite-dim items authored.

**v5 (two-answer WL+QF manual, NEW deliverable):** user asked for a NEW md giving each Question-List
question TWO worked answers, (1) native-WL only and (2) QuantumFramework, no hardcoding, QF answers
using property downvalues. Read the two reference notebooks fully (introduction-to-QIS-revised twin =
native-WL idiom; Introduction-to-quantum-computing twin via NotebookToMarkdown = QF idiom) + qf skill.
Deliverable `Quantum-in-Finite-Dimensions-WL-and-QF-Answers.md`. **Part 1 pilot (9 Q x 2 answers)
DONE + verified**: extract-run 20/20 cells exit 0 zero messages, WL===QF on every item, and
/wl-verify wl-verifier OPEN ISSUES 0 (one round, converged). Reusable QF idioms confirmed live:
qs["StateVector"]//Normal, ["ProbabilitiesList"], ["Norm"], ["Normalize"], (bra["Dagger"]@ket)
["Scalar"] for inner product, ["Formula"] (Dirac sum, Head RawBoxes), ["BlochVector"]/
["BlochSphericalCoordinates"], ["Dimension"] to drive dof. TRAP reconfirmed: QuantumDistance
[.,.,"Fidelity"] is the DISTANCE 1-Sqrt[F] (0 for global-phase-equal states), NOT fidelity; pilot
avoids it.

**Settled answer-style rules (apply to all parts), after user feedback:** (1) SYMBOLIC wherever
possible: general amplitudes {\[Alpha],\[Beta]}, angles {\[Theta],\[Phi]} (the symbolic input also
makes every operation "bite" in full generality, satisfying the earlier pick-a-non-degenerate-example
rule). (2) NO helper functions: write the computation inline (e.g. the explicit Bloch formula
{2 Re[Conjugate[a]b], 2 Im[Conjugate[a]b], Abs[a]^2-Abs[b]^2}) because the code IS the explanation.
(3) Prefer WL BUILT-INS that shorten the answer: ToSphericalCoordinates (= CoordinateTransform but
shorter) for Cartesian->spherical. (4) Avoid hardcoding even minor cases: Normalize@ConstantArray[1,d]
not ConstantArray[1,d]/Sqrt[d]. (5) Use QF NAMED/built-in features over manual reconstruction:
QuantumState["UniformSuperposition",d] not ConstantArray normalized. (6) COMPUTE FROM THE
DEFINITION, do not write a derived closed-form and call it the definition: e.g. the Bloch vector is
DEFINED as the Pauli expectation values <psi|sigma_j|psi> (compute Table[Conjugate[psi].PauliMatrix[j]
.psi,{j,3}]); the component triple {2Re[Conj[a]b],2Im[Conj[a]b],Abs[a]^2-Abs[b]^2} is the RESULT of
that evaluation, shown in prose, never presented as the definition. (User caught this on 1.8.) Verified-live QF symbolic facts:
QuantumState[{a,b}] accepts symbolic; ["ProbabilitiesList"] auto-normalizes (match with
Abs[Normalize[{a,b}]]^2 on the WL side); ["Norm"]/["Normalize"] symbolic; == returns True symbolically
for global-phase-equal states (=== False); (bra@ket)["Scalar"] gives the RAW (unnormalized) inner
product a2 Conj[a1]+b2 Conj[b1]; ["BlochVector"] & ["BlochSphericalCoordinates"] FullSimplify (with
0<theta<Pi && -Pi<phi<Pi) to {Sin t Cos p,...} and {1,theta,phi}. Part 1 rewritten symbolic, /wl-verify
OPEN ISSUES 0 (3 verifier rounds total across the pilot's revisions). FORMAT SETTLED. Next: Parts 2-25
(~185 questions x 2) on these rules.

**v3 (blind-exam pass, user request "drop outputs"):** the file already had no fenced output
blocks; the user meant the result-revealing prose. Stripped the explicit computed values (numbers,
matrices, True/False) from every interpretation paragraph across all 13 parts, and trimmed the
"-> value" result annotations in code comments to operation-describing ones. KEPT as physics:
definitions/setups (normalization, the chosen hidden string s=101, the Tsirelson target bound, the
GHZ-paradox explanation, the g^(2) bunching/Poissonian/antibunching classifications) and all
qualitative interpretation. Now reads as a true blind exam: question + concept + runnable code +
physics interpretation, with the numeric answer found only by running the cell. Re-run clean 84/84,
exit 0, zero messages.
