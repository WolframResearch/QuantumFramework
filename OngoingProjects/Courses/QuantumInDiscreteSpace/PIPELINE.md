# Pipeline & Brief: "Quantum in Discrete Space" QF comprehensive-exam manual

This file is the reusable brief and the build tracker. Original deliverable:
`Quantum-in-Discrete-Space-QF-Comprehensive-Exam.md`. Spun-off deliverables in this folder:
`Question-List.md` (the curriculum) and `Quantum-in-Finite-Dimensions-WL-and-QF-Answers.md`(+`.nb`)
(the two-answer key). The section just below distills the DNA the latter is built on.

## Answer-key design DNA (foundational principles, dialectically distilled; govern Parts 3-25 and future WL teaching-docs)

Distilled by a two-agent dialectic (one agent boiled each session directive to an example-independent
principle, the other generalized it, named its failure mode, and transferred it to later Parts). The
through-line is **faithfulness + generality + transparency + verification**. Each rule is stated
*independently of the example that produced it* so it carries to any Part or topic; that carrying is
itself the last principle. **Format premise:** every question gets a native-WL answer and a QF answer,
agreeing by *exact* equality (WL===QF), proven in-kernel, so the bare linear algebra and the object
model stand side by side. Each principle below names the failure mode it prevents.

**Content and pedagogy**
1. **Reader computes, never reads the conclusion.** Withhold the worded result; the learner runs or
   derives it. Test: *could the reader have written this sentence before evaluating the cell?* If yes,
   it is a spoiler, cut it. (Failure: spoiler prose, the key answers itself.)
2. **Derive from definitions; never restate the result as the premise.** The definition is the
   starting point; the formula/value is the computed conclusion. (Failure: question-begging, premise = conclusion.)
3. **Maximize generality: symbolic over numeric** wherever the kernel allows ($\alpha,\beta,\theta,\phi$,
   a generic Hermitian $\{\{a_1,a_2+ia_3\},\dots\}$); a number is an incidental instance of the general
   case. (Failure: premature numerics, a fact stated too narrowly to generalize.)
4. **Discriminating examples only.** An example is valid iff *toggling the operation under study
   changes the rendered output*; otherwise it is degenerate (complex amplitudes so conjugation bites,
   unequal weights so the Born rule discriminates, off-axis states, non-commuting operators).
   (Failure: degenerate witness, the example is a silent no-op.)
5. **Atomicity: one concept per question, one computation per cell.** A single mathematical object (a
   vector, a distribution, one boolean, a tight derivation chain like ip/|ip|/|ip|^2) stays whole; a
   list of *independent* computations splits into separate captioned cells. (Failure: concept bundling, a miss cannot be localized.)
6. **Physicist register.** Open with the physics (what the operator *does* to the state); gloss or
   banish software-engineering jargon. (Failure: register mismatch.)

**Code and tool craft (trust the object)**
7. **Minimal idiomatic tool** = the lightest construct whose *own representation still exhibits the
   property under study*; escalate to a heavier abstraction only when it is itself the minimal faithful
   form (a `QuantumChannel` for a channel, but a plain `QuantumOperator`, not a channel, for a unitary).
   (Failure: abstraction inflation, ceremony obscures the physics.)
8. **Code is the explanation: explicit and inline, no hidden helper definitions.** If a step matters it
   is visible at the call site. (Failure: off-screen magic, the load-bearing step is invisible.)
9. **Let the abstraction do its work.** Compute from the object's own properties/behavior, never a
   hand-rolled proxy or hardcoded literal (`Normalize` not `/Sqrt[d]`; named states; `state == state`
   not a `ProbabilitiesList` comparison; `EigenvalueDecomposition`/`PauliDecompose`). (Failure: proxy drift, the stand-in and the object diverge unnoticed.)
10. **Show the object, not a matrix of it.** Render the object (assign without `;` so the summary box
    displays); a matrix/amplitude-vector is one view read from it. Stay operator-level (`==`, `@`,
    `Commutator`). (Failure: flattening, the rich structure that is the point is thrown away.)
11. **Well-posed and precisely named.** Unambiguous referents ("the partial trace over $B$"); name a
    property in prose as `"X"` or `obj["X"]`, NEVER a bare bracket `["X"]` (strict, qf skill 0.5).
    (Failure: ambiguous referent.)
12. **Valid delimited TeX, no em/en dash (build correctness, not style).** Malformed `$...$` fails to
    render; Unicode math in prose is wrong by policy. (Failure: render rot, the published artifact silently breaks.)

**The Output Boundary** (reconciles 1 with 8/10): withhold the *worded answer*; render the *object the
code produces*. The summary box, tableau, plot, or entropy value is the computation the reader is meant
to obtain, not a spoiler of it. "Drop outputs" means drop worded answers, never rendered objects.

**Curriculum (the question list)**
13. **Complete and dependency-ordered.** Cover the whole declared scope; every prerequisite precedes
    its first use (no forward reference); one concept per question; tag [BSc]/[MSc]. (Failure: scope hole or forward-reference.)

**Verification and process**
14. **Run the user's literal first; verify, don't recall.** Test the exact supplied token/value before
    "correcting" it; never overwrite supplied ground truth from memory; aim the first check at the
    load-bearing premise, before asserting. (Failure: confabulated correction, e.g. AntiCommutator over the real Anticommutator.)
15. **Fix the cause, not the symptom.** When something was missed or wrong, diagnose the gap-generating
    cause so it cannot recur. (Failure: whack-a-mole patching.)
16. **Fresh-context adversarial review before "done."** Hand finished WL work to a clean-context
    reviewer that reproduces every claim in a kernel (`wl-verifier` / `/wl-verify` / `quantum-review`).
    (Failure: self-certification, the author who erred grades the work.)
17. **Single source, standard pipeline.** The literate `.md` (with a `Template:` key, here `Default`)
    is canonical; regenerate the `.nb` via `MarkdownToNotebook[src,out,"Evaluate"->False]` (md2nb),
    never hand-edit it; read it back to confirm structure and clean re-evaluation. (Failure: source desync.)
18. **Spin off out-of-scope work; don't gold-plate inline.** A kernel bug or adjacent topic surfaced
    mid-Part becomes a separate task, not a detour. (Failure: scope creep.)

**Carrier**
19. **Principles transfer, instances don't.** Apply every rule above list-wide and to later Parts
    without being re-prompted; a one-off correction is a standing policy in disguise. This rule
    transports 1-18 to Parts 3-25 untouched. (Failure: instance myopia, re-litigating the same rule per item.)

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
OPEN ISSUES 0 (3 verifier rounds total across the pilot's revisions). FORMAT SETTLED.

**Part 2 (Observables, spectra, approximation; 13 Q x 2 answers) DONE.** Applied all Part-1 lessons.
Verified-live QF operator API: operator algebra `a @ b` = product, `+/-` combine, `["Matrix"]//Normal`;
symbolic QuantumOperator accepted; `["HermitianQ"]`, `["Eigenvalues"]`, `["Eigenvectors"]`,
`["Eigensystem"]`, `["PauliDecompose"]` (assoc); NO `"Variance"` measurement prop (compute <A^2>-<A>^2);
multi-qubit operator = `QuantumTensorProduct[qoA, qoB]` (the `\[CircleTimes]`/`QuantumOperator[{"X","X"}]`
forms do NOT work). Answers: 2.1 commutator/anticommutator (explicit a.b-b.a, QF @ algebra); 2.2 Levi-
Civita Outer (user's idiom); 2.3 HermitianQ + real-spectrum proof; 2.4 spectral decomp reconstruction;
2.5 functional calculus FROM the spectral definition + ==MatrixFunction check; 2.6 <sigma_z>=cos th,
Var=sin^2 th; 2.7 resolution of identity; 2.8 XX/ZZ commute, Bell simultaneous eigenbasis; 2.9 Z(x)I,
I(x)Z CSCO joint spectrum; 2.10 Hilbert-Schmidt c_j=Tr[sigma_j A]/2 vs QF PauliDecompose; 2.11 Rayleigh-
Ritz (Minimize reaches E0=-Sqrt2); 2.12 RS-PT 1st {w,-w} + 2nd order; 2.13 avoided crossing gap
2 Sqrt[Delta^2+lambda^2]. Full file (Parts 1+2) extract-run = 46 cells exit 0 ZERO messages, WL===QF on
all. psmallmatrix->pmatrix for GitHub KaTeX. 2.1 also shows the built-ins {Commutator[a,b,Dot],
Anticommutator[a,b,Dot]} per user; NOTE the System symbol is **Anticommutator** (lowercase c) and
Commutator; **AntiCommutator (capital C) does NOT exist** (I wrongly "corrected" the user's correct
lowercase spelling from memory -> verify-don't-recall failure; the kernel and Names["*ommutator*"] settled
it). wl-verifier (fresh-context adversarial) OPEN ISSUES 0 (confirmed all 13 WL===QF, the Anticommutator
built-in, the physics, no helpers/hardcoding/Fidelity, Part 1 unregressed). PARTS 1+2 DONE & VERIFIED.

**Style rule (6b), QF results are OBJECTS not matrices (user feedback on 2.1):** a QF answer should
return/keep the rich object (QuantumOperator, QuantumState), NOT collapse it to `["Matrix"]//Normal`
as the deliverable. Show a matrix/vector only as ONE representation and say the object is richer (it
acts on a target state, carries its qubit Order, answers ["HermitianQ"]/["Eigenvalues"]/["Dagger"],
changes basis, composes). 2.1 QF rewritten to return {Commutator[a,b], a@b+b@a} (two QuantumOperators)
+ a richness cell (c=[X,Y]=2iZ: c[|0>]=2i|0>, HermitianQ False, Eigenvalues {2i,-2i} = anti-Hermitian,
Matrix=2i sigma_z shown LAST). Intro updated with the principle. 2.8/2.9 commute test changed from
matrix extraction to operator equality `a @ b == b @ a`.
**Verified-live QF operator facts:** `Commutator[qoA,qoB]` REDUCES to a QuantumOperator (QF teaches the
System `Commutator` its operator product); `Anticommutator[qoA,qoB]` does NOT reduce (Head Plus, stays
`qo**qo+qo**qo` with NonCommutativeMultiply which QF never defines) -> write the anticommutator as
`a @ b + b @ a`. QF operator `==` is a meaningful exact equality: `a@b==b@a` is True for commuting,
False otherwise (use it for commute tests, not matrix compares). QF Dot (`.`) is NOT defined for
operators (stays symbolic); use `@` for the operator product. 47-cell full-file extract-run exit 0 zero
messages; wl-verifier OPEN ISSUES 0 (confirmed the object-returns, the Commutator/Anticommutator
asymmetry, operator==, no regression).

**Style rule (7), DERIVE the answer, never HARDCODE the thing the question asks you to find (user
flagged 2.8).** 2.8 asks to FIND the simultaneous eigenbasis; the old answer hardcoded the Bell basis
(bell={{1,0,0,1},...}) and just checked it. Fixed: since X(x)X and Z(x)Z are each degenerate (eigs +-1
doubled), DIAGONALIZE a generic non-degenerate combination `Eigenvectors[a + Pi b]` -> the eigenvectors
ARE the simultaneous eigenbasis (the Bell states fall out, derived not assumed). QF object-native check:
each found QuantumState v satisfies `a[v] == v` (state == is up-to-phase, valid here because eigs +-1
are unit-modulus; returns False for a non-eigenstate, so it discriminates). 2.2 QF also switched to the
user's fully operator-level form: With[{sigma=...},{com=Outer[Commutator[#1,#2]&,sigma,sigma,1]},
Table[com[[i,j]]==2 I Sum[LeviCivita[[i,j,k]] sigma[[k]],{k,3}],...]] -> 3x3 all True, no ["Matrix"].
(3-arg With sets a global `com` via inner Set; harmless, only used in that cell.) wl-verifier OPEN
ISSUES 0 (no hardcoded bell in code; eigenstate test has teeth; no regression).
**2.4 now uses EigenvalueDecomposition (user suggestion, doc read + kernel-verified).** EVD facts
(verified live, not recalled): `EigenvalueDecomposition[m]` yields `{S, D}` for a diagonalizable square
matrix, eigenvectors as COLUMNS of S, eigenvalues on Diagonal[D], and `m == S . D . Inverse[S]` (uses
Inverse[S], NOT S-dagger: it does not normalize eigenvectors, so S is a general similarity, not unitary;
for Hermitian m, normalizing S to unitary U gives the spectral/projector form U D U-dagger = sum a|a><a|).
EVD does NOT accept a QuantumOperator (errors `::matsq` "not a square matrix") -> QF side stays object-
native: qo["Eigenvalues"], qo["Eigenvectors"], qo["Diagonalize"] (the latter returns a QuantumOperator
= the operator in its own eigenbasis, diagonal with eigenvalues on it = the spectral form as an OBJECT).
WL EVD eigenvalues === QF Eigenvalues (as a set). wl-verifier OPEN ISSUES 0; 47-cell file clean.
**NB build (md2nb):** built `Quantum-in-Finite-Dimensions-WL-and-QF-Answers.nb` from the .md twin via
`MarkdownToNotebook[src, out, "Evaluate" -> False]` (local function at GitHub/MarkdownToNotebook).
**Template = Default** (user switched from TechNote). Default builder (`defaultNotebook`) maps
`$headingStyleMap = <|1->Title, 2->Section, 3->Subsection, 4->Subsubsection|>`, takes the Title from an
H1 `#` in the BODY (NOT frontmatter Title, unlike TechNote), and adds NO template slots (no
MoreAbout/Keywords/Metadata/Categorization ceremony) -> `StyleDefinitions -> Default.nb`. So for Default
the source needs frontmatter `Template: Default` + a restored H1 title (TechNote took Title from
frontmatter and had no level-1). Rebuilt: 145 KB; 1 Title (from H1), 3 Sections, 22 Subsections, 47
Input, 93 Text, 187 InlineFormula, FormalA 3/3. Evaluate->False = Input cells only (no Output), which matches the source design (reader
runs cells) and dodges released-vs-working-tree eval. Built nb 150 KB; structure verified by reading it
back: 1 Title, 3 Sections (Setup, Part 1, Part 2 = the 3 `##`), 22 Subsections (the 22 `###` questions),
47 Input, 93 Text, 187 InlineFormula (the `$...$` math), 0 empty `\[Null]` mathboxes (no dropped LaTeX).
FormalA preserved 3/3 (raw-nb grep == source). The 47 reparsed Input cells from the BUILT nb run clean
top-to-bottom against the working tree (exit 0, zero msgs) = md2nb code-parse is faithful. md2nb GOTCHA
noted: NotebookToMarkdown renders `\[FormalA]` as its display glyph `a` and `$A$` as the math-italic
char in the twin (cosmetic round-trip rendering; the .nb itself preserves `\[FormalA]`). Rebuild the nb
whenever the .md changes (it is a build output, not hand-edited).
**One-computation-per-cell + narrative + show-objects pass (user).** Split every bundled-independent-
computation list into separate ```wl cells, each preceded by a concise narrative; kept single objects
whole (Bloch vector, prob dist, single boolean, the inner-product triple ip/|ip|/|ip|^2, an Assoc of
coeffs, a Table). QF objects SHOWN (assign w/o `;` so the summary box displays, then extract). 1.2 got
`// FullSimplify` (worded precisely: it reduces to a single ratio with alpha Conjugate[alpha] in the
numerator, not literal Abs^2). Answer file now **93 cells** (was 50). **Strict QF prose rule added to
the qf SKILL.md (section 0.5):** in prose never write a property as a bare attached bracket
`["StateVector"]`; use the quoted name `"StateVector"` or a full object call `qs["StateVector"]`. Fixed
the one doc violation (1.1 prose). Full audit via quantum-review + a fresh-context reviewer (4 axes:
correctness / one-per-cell / narrative-before-each / no-bare-property-prose): A1=2 (both in 2.8, now
split), A2=0, A3=0, A4=0; physics re-derived clean. Source extract-run 93 cells exit 0 zero msgs; nb
rebuilt (Default template, Evaluate->False) Input=93/Section=3/Subsection=22/Title=1, FormalA preserved,
all 93 reparsed cells re-eval clean. (TRAP noted: don't name a check-script var `t` while evaluating the
answer cells, 2.6/2.11 use `t` as the trial angle -> Minimize::ivar collision; false alarm, not a cell bug.)
**Part 3 (States as density operators, 3.1-3.8) BUILT to the DNA + 3-adversary loop.** Also folded the
external dof report (`/private/tmp/qubit-dof-two-variables.md`) into 1.5: replaced the hardcoded
`2d-2` assertion with the report's **Jacobian-rank** computation (rank 3 with global phase, 2 without
= the qubit's 2 physical dof); same idiom reused in 3.5 (qubit Hermitian-unit-trace rank 3 = d^2-1).
Verified-live QF density-operator facts (probe-first, per DNA 14): ensemble mixture = **sum of
`["DensityMatrix"]`** (a weighted sum of QuantumStateS is a coherent SUPERPOSITION, not a mixture, no
named/assoc mixture ctor exists); `["Purity"]`, `["PureStateQ"]`, `["BlochVector"]` (works on mixed),
`["VonNeumannEntropy"]` (Quantity in Bits); expectation on a mixed state = `QuantumMeasurementOperator[
QuantumOperator["Z"]][state]["Mean"]` (the OPERATOR, not the string "Z" which gives a different value;
no `["Variance"]` -> use <Z^2>-<Z>^2; symbolic Bloch needs `Simplify[...,-1<=rz<=1]` since Mean goes via
Born probs). Bloch map is AFFINE (convexity). 3 adversaries run in parallel: pipeline-enforcer (4: 3.7
WL/QF were DIFFERENT mixtures so WL===QF was coincidental -> made WL use the same |0>,|R> pair; 3.2 and
3.7 bundled independent computations -> split; 3.4 QF object suppressed with `;` -> shown), minimal-code
(5 verified: `Through[{Re,Im}]&/@` -> built-in `ReIm`; `Table[PauliMatrix[j],{j,3}]` -> Listable
`PauliMatrix[{1,2,3}]`), quantum-correctness wl-verifier (1: 3.1 prose falsely said summing the
{|0>,|+>} kets gives |+> -> reworded). All applied. Source extract-run **126 cells exit 0 zero msgs**,
WL===QF on all Part 3 (3.7 now genuine 7/9=7/9). nb regenerated (Default, Evaluate->False) Input=126
Section=4 Subsection=30; nb-cell re-eval confirmed clean (background; Part 3 symbolic QMO/FullSimplify
is slow, ~exceeds 2min foreground).
Next: Parts 4-25 (same probe-first + 3-adversary loop, one part per turn).

**v3 (blind-exam pass, user request "drop outputs"):** the file already had no fenced output
blocks; the user meant the result-revealing prose. Stripped the explicit computed values (numbers,
matrices, True/False) from every interpretation paragraph across all 13 parts, and trimmed the
"-> value" result annotations in code comments to operation-describing ones. KEPT as physics:
definitions/setups (normalization, the chosen hidden string s=101, the Tsirelson target bound, the
GHZ-paradox explanation, the g^(2) bunching/Poissonian/antibunching classifications) and all
qualitative interpretation. Now reads as a true blind exam: question + concept + runnable code +
physics interpretation, with the numeric answer found only by running the cell. Re-run clean 84/84,
exit 0, zero messages.
