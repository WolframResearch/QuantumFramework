# QEC Roadmap — September 2026

*Companion to `QEC-Development-Plan.md`, which stays as the original statement of
intent. That plan described the code layer as a
pre-QA prototype with 179 tests; it has since been rebuilt and items (a) through (d)
below are done. This revision records where that leaves us and what the remaining work
actually is, mapped against Gottesman's* Stabilizer Codes and Quantum Error Correction
*(2026), which is the reference the layer is built and cited against.*

Quantum error correction protects one *logical* qubit by spreading it across many
*physical* qubits, in a subspace chosen so that the errors we care about push the state
out of that subspace in a way we can detect by measuring a fixed set of parity checks,
and then undo. Everything below is organized around that one idea: build the code, read
out the error's fingerprint, guess the error, undo it, and ask whether the logical qubit
actually survived the noise.

---

## 1. Where we are now

### The stabilizer engine (shipped in the paclet: solid)

Unchanged from the first version of this plan, and still the sound base: `Kernel/Stabilizer/`
holds an n-qubit stabilizer state compactly, pushes it through Clifford gates in the
Heisenberg picture, measures any qubit or check operator, and reads expectations, inner
products and entanglement entropy. Cross-checked against Stim and QuantumClifford.jl.

Two defects found while building on it, both worth knowing:

- **`ps["M", q]` returns an outcome that depends on which generating set the state was
  built from.** Two `PauliStabilizer` objects that are the same physical state give
  different measurement outcomes. Documented with a minimal reproduction in
  `EngineMeasurementBug.md`; the QEC layer reads deterministic outcomes with
  `(1 - ps["Expectation", p])/2` instead and never calls `["M", q]` on a
  generator-constructed state. Mads' commits `f2a887ed` and `4d190ce8` harden exactly
  that code path and may have fixed it — re-test before removing the workaround.
- The summary box of a symbolic `QECNoiseModel` failed to render its chart, because
  `expr /. _Symbol -> 0.3` rewrites heads and `List` is itself a Symbol. Fixed.

### The code layer (rebuilt, in `OngoingProjects/QEC/`)

The layer is no longer a prototype. It is a functional-idiom package of 16 files with
**469 tests**, built on GF(2) symplectic rows with Z4 phases rather than strings, and
every exported name is prefixed `QEC`. What it does, in the order the physics builds up:

**The code object.** `QECCode` from generators, from a check matrix, or by name. Check
matrix, standard form, logical operators, exact distance with a minimum-weight witness,
`[[n,k,d]]`, encoding circuit, syndromes, CSS / concatenation / qubit-removal
constructions, and the 3/5/7/9-qubit, repetition and Hamming-based families.

**Two decoders.** Minimum-weight over a lookup table, and maximum likelihood over
cosets. The coset label — syndrome bits and logical-class bits in one GF(2) matrix
product — is what makes both cheap.

**Three noise levels in one object.** `QECNoiseModel` at code-capacity, phenomenological
(readout can lie), and circuit level (four rates: one-qubit gate, two-qubit gate, reset,
measurement, plus idle). Rates may be **symbolic**, and that is the point.

**The syndrome-extraction circuit, and honest noise on it.** One bare ancilla per check
— explicitly the *non*-fault-tolerant construction of the book's sec. 12.1.1 — plus a
Pauli frame propagator, a detector error model mapping every circuit fault to the
detectors it fires and observables it flips, a layered schedule so that waiting is a
fault location of its own, exact and sampled logical error rates, and a Stim export.

**Documentation.** Two tech notes built from literate markdown through
MarkdownToNotebook: `StabilizerCodes` (the layer from the outside) and
`QECCoreInternals` (the circuit-level modules function by function, 80 evaluated cells).

### Results worth showing

- **Logical error rates come back as exact polynomials**, not sampled estimates:
  `3p² - 2p³` for the bit-flip code at code capacity, and a closed-form rational
  function in `p` at circuit level.
- **A consistency identity ties the levels together.** A phenomenological model with
  perfect readout over one round *is* code-capacity noise, so two entirely different
  code paths — enumeration over cosets, and a fold over detector effects — must give
  the same polynomial. They do.
- **The bit-flip / phase-flip duality survives one level and breaks at the next**:
  `32p/15` against `112p/15` at circuit level, entirely because the X-type checks need
  eight extra Hadamards. Switch off one-qubit-gate noise and they agree again. That is a
  statement about a *circuit*, not a code, and it is awkward to make any other way.
- **The five-qubit code drops from `p²` to `29p/5`** at circuit level, and the layer says
  why: 288 circuit faults collapse onto 71 detector signatures, of which 29 carry
  conflicting logical effects.
- **Cross-checked against Stim**: 58 detector firing rates and 5 observable flip rates
  across five cases, computed exactly here and sampled there, all agreeing within 2.7σ
  of two million shots, with PyMatching decoding every exported circuit unchanged.

### Coverage against the book

The framing matters, because the two are not after the same thing. **The book proves
fault-tolerance properties of gadgets; this layer computes exact error rates of stated
protocols.** So the layer has none of the book's guarantees — no ECRP, ECCP, GPP or GCP
— and the book does not compute what the layer computes: sec. 12.2.2 raises decoding the
whole spacetime history and sets it aside as *"difficult to analyze ... in the completely
general case"*.

| have | book |
|---|---|
| stabilizer codes, logical operators, exact distance | ch. 3–5 |
| syndromes, cosets, two decoders | ch. 4 |
| encoding circuits (non-FT, which is what ch. 6 gives) | ch. 6 |
| CSS, concatenation, qubit removal | ch. 8, sec. 9.1 |
| locations, faults, one rate per location type | sec. 10.1.1, Def 10.1–10.3 |
| error propagation through Cliffords | sec. 10.1.2 |
| the Pauli frame | sec. 12.5.2 |
| non-FT Pauli measurement | sec. 12.1.1, fig. 12.1a |
| cat states and their verification | sec. 12.1.2–12.1.3 |
| ancilla reuse by reset, residual as a preparation error | sec. 15.4, 15.4.3 |
| waiting as a fault location; time steps defined, not dropped | sec. 10.1.1 Def 10.1; sec. 15.5.2 |
| detectors, and post-selection accounted for honestly | *not the book* — DKLP, Gidney |

| do not have | book |
|---|---|
| FT Pauli measurement (repetition and majority) | sec. 12.1.4–12.1.5, Thm 12.1 |
| Shor / Steane / Knill error correction | sec. 12.2 / 12.3 / 12.4 |
| GPP, GCP, ECRP, ECCP as checkable properties | sec. 10.2 |
| transversal gates | ch. 11 |
| FT preparation of encoded states | sec. 13.1 |
| gate teleportation, Clifford hierarchy, magic states, distillation | ch. 13 |
| adversarial noise, exRecs, benign/malignant sets, the \*-decoder | sec. 14.1–14.4 |
| level reduction, the threshold theorem, pseudothresholds | sec. 14.6–14.7 |
| Solovay–Kitaev, short-range gates, slow measurement, leakage, biased and non-Markovian noise | ch. 15 |
| surface and toric codes | ch. 17 |

The one-line summary: **we can build a code, put honest circuit-level noise on it, and get
an exact answer for whether the logical qubit survives — but we cannot yet build a gadget
that is fault tolerant, and therefore cannot yet reproduce the book's threshold results.**

---

## 2. What to develop

Items (a) through (d) of the first version — the code core, the decoder, the noise
object, and the logical error rate — are done. What follows replaces them.

### The critical path: a fault-tolerant simulation

Everything the book proves rests on one object: `FT(C)`, the fault-tolerant simulation of
an ideal circuit `C` (Def 10.6). Building it is the gate to chapters 12, 13 and 14, and
its dependencies are forced, because `FT(C)` needs a gadget for each kind of location —
preparation, gate, wait, measurement — and every one of them needs error correction,
which needs fault-tolerant measurement, which needs verified cat states.

| phase | deliver | book | status |
|---|---|---|---|
| **A** | cat states, verified, with a **derived** check set | sec. 12.1.2–12.1.3 | **done** |
| **B** | `QECPauliMeasurement` — cat measurement, `2t+1` repetitions, majority | sec. 12.1.4–12.1.5, Thm 12.1 | next |
| **C** | `QECErrorCorrection` — Steane EC | sec. 12.3 | |
| **D** | multi-block addressing: a logical qubit is a *block* | — | |
| **E** | transversal gates for the 7-qubit code | ch. 11, sec. 11.3 | |
| **F** | **`QECFaultTolerant[circuit, code]`** — the assembler | Def 10.6, ch. 14 | |
| **G** | exRecs, benign/malignant sets, the rigorous threshold bound | sec. 14.2–14.7 | |

Three decisions already taken, and the reasons, so they do not get relitigated:

- **The 7-qubit Steane code and Clifford gates**, not full generality. It is the book's
  own worked example, and sec. 14.5.3 computes `A = 735`, `Q = 63p` and
  `p_T > 3.6 × 10⁻⁶` for exactly that protocol — so the numbers can be *reproduced and
  validated against*, which no other choice offers.
- **Steane EC first.** CSS only, but the `e + f + g` additivity means no repetition logic
  is needed, and the detector model is already linear over GF(2); its ancillas double as
  preparation gadgets (sec. 13.1.2); and `codeCSSQ` already gates the precondition.
- **Corrections stay in the Pauli frame.** No classically-conditioned instructions are
  needed: for a Clifford circuit with Pauli errors, tracking the correction is
  equivalent to applying it (sec. 12.5.2), which is what real experiments and Stim both
  do. This removed a whole subsystem from the critical path.

**The milestone worth aiming at.** Phases B + C + D give `FT` of a *memory* experiment
without needing E — prepare, wait, measure — and that is enough to compare fault-tolerant
extraction against the bare-ancilla circuit we have today, and **see the `p²` restored
against the current `O(p)`**, measured in our own layer. That single plot is the evidence
that the hook error is a defect of the gadget and not a property of circuit-level noise.

### Standalone, and cheap: things with no missing dependencies

These do not wait on the critical path and each is days rather than months.

**(e) `QECPseudothreshold` — the threshold, honestly (sec. 14.7.6).** The one item of the
first plan still outstanding, and it is now nearly free. Because the detector model
returns an *exact polynomial* in `p`, the level-0/level-1 crossing is a polynomial root:
`Root` or `NSolve`, arbitrary precision, no Monte Carlo, no fit, reproducible bit for
bit. Sweeping the direction vector over the four circuit-level rates traces a **threshold
surface** by root-finding along rays. Two honesty requirements the book is emphatic
about: a pseudothreshold is *not* a threshold — sec. 14.7.6's own example lies
**outside** the true surface — and a threshold is a property of (code, gadget set, noise
model, decoder, level), never of a code alone. The return value should carry all of them,
with no bare-number path.

**`QECWeightEnumerator`, signed (sec. 7.4.1).** `A(x) = Σ_{M∈S} phase(M) x^{wt M}`.
Pure GF(2), enumerable for the codes we care about, and it reproduces the `1/6`
acceptance of the 5-qubit magic-state protocol analytically. Highest value per unit
effort in the whole list. The sign matters: our codes come from check matrices, which
discard it.

**`QECCliffordHierarchyLevel` (Def 13.1).** Conjugate the `2n` Pauli generators and
classify the image. Easy for `k ≤ 3`, which is all that is needed, and it unlocks the
vocabulary of ch. 13.

**A truly layered emitter (sec. 15.5.1).** The schedule is *discovered* from a
sequential instruction list today; emitting in layers would make idle noise realistic
rather than merely chargeable, and shortens the circuit. Note the current default
`Idle -> 0` is optimistic rather than neutral: serial extraction is the arrangement with
the *most* waiting.

**(f) The families that scale (ch. 17).** `QECSurfaceCode[d]`, `QECToricCode[...]`,
`QECBivariateBicycleCode[...]`, returning the same code object so distance, decoding and
rates work unchanged. Needed before any threshold estimate over a code *family*, and the
natural consumer of the PyMatching path.

### (g) Longer term: computing on the protected qubit

Chapter 13 in full — gate teleportation, magic states, distillation — and then resource
estimation. The five subsystems that gate it, in the order they bite:

1. **verified ancillas** (done), which was the gate to all of ch. 12;
2. **multi-block addressing** — Thm 13.2 needs `2m` blocks for a gate touching `m`;
3. **transversal-gate machinery** — which Cliffords are transversal on a code, and with
   what logical action;
4. **a Clifford frame, not just a Pauli frame** — a `𝒞₃` correction is a Clifford and
   cannot be absorbed into a Pauli frame;
5. **the `[[15,1,3]]` code** for distillation, which needs Reed–Muller — absent from the
   repo, and its generator matrix is in sec. 11.5.1, not ch. 13.

And a boundary to respect: **magic states are not stabilizer states.** `R_{π/8}|+⟩` and
`|R⟩` cannot be held in a tableau, so the DEM and the Stim export cannot see them. The
honest architecture is what real pipelines do: keep the magic state as an *input
assumption* with an injection error rate, export only the Clifford remainder, and verify
the injection gadget separately in the dense simulator at physical size.

### The one thing that is genuinely ours

Unchanged, and now demonstrated rather than asserted: everything computed exactly on a
small or symbolic code is checked against Stim, and the move Stim structurally cannot
make is to carry a symbolic parameter and answer for a whole family at once. `32p/15`
against `112p/15`, and `288 → 71 → 29`, are results of that kind. Stim simulates fast and
large; we derive and certify exactly and symbolically.

---

## 3. What we should not build

Not a fast large-scale simulator (Stim owns that), not a production or real-time decoder
(PyMatching and Riverlane own that), not lattice-surgery compilation at industrial scale
(Qiskit and tket own that), not GPU decoding. Connect to those; do not rebuild them.

Two additions learned since:

- **Not a unitary synthesizer.** Ch. 14 gives universality by a density theorem (Thm 6.9)
  and no synthesis algorithm; `𝒞ₖ` is finite. Promising a `QECSynthesize` would be
  promising something the theory here does not provide.
- **Not the proof machinery.** The \*-decoder's Choi–Jamiołkowski / Uhlmann argument, the
  `ℰ` map of Def 14.8, the r-filters of Def 10.7 and the graphical decoder-pushing
  calculus are proof devices with no computational content. They should never appear in
  the paclet.

---

## 4. Order of work

1. **Finish the critical path to a memory experiment**: phases B, C, D. This is the block
   that makes the layer able to say something the book would recognize, and it ends in
   the `p²`-versus-`O(p)` comparison.
2. **In parallel, land the cheap standalone items**: the signed weight enumerator, the
   Clifford hierarchy test, and `QECPseudothreshold`. None of them waits on anything, and
   `QECPseudothreshold` closes the last open item of the original plan.
3. **Then phases E and F**: transversal gates and `FT(C)`, which is the object ch. 14 is
   written about.
4. **Then the scalable families (f)** and the empirical thresholds over a family that
   they enable.
5. **Then phase G**: exRecs and malignant-set counting, and the rigorous threshold bound
   — with the caveat that it will land 1–3 orders of magnitude below any simulated
   number, which is what the book's own figures show (`3.6 × 10⁻⁶` crude,
   `2.7 × 10⁻⁵` with careful counting, against `0.7%` simulated for the surface code).
6. **Ch. 13 and resource counting (g)**, longer term.

Two standing practices: check every exact result against Stim as you go, and put the
weight on the exactness-and-symbolic angle wherever it is genuinely ours.

### Housekeeping

- Migrate `QECCore/*.wl` into `QuantumFramework/Kernel/QEC/`, then the guide page and the
  symbol reference pages; harvest `Courses/QuantumErrorCorrection` and the scattered
  5/7/9 notebooks. Note that folder is now untracked here and lives in a private repo.
- Re-test `EngineMeasurementBug.md` against Mads' hardening commits and retire the
  workaround if it is fixed.
- A recurring trap worth a lint rule: a `Table` with two or more iterators whose body
  returns a list needs `Flatten[..., 2]`, not `Catenate`. It has silently produced a
  malformed result three times in this package.

---

## References

- Gottesman, *Stabilizer Codes and Quantum Error Correction* (2026 draft), and the 1997
  thesis, [arXiv:quant-ph/9705052](https://arxiv.org/abs/quant-ph/9705052): the
  formalism, the gadget constructions, and the threshold theorem this plan is mapped
  against. Sections cited throughout the layer's file headers and the tech notes.
- Aaronson & Gottesman, *Improved simulation of stabilizer circuits*, Phys. Rev. A **70**,
  052328 (2004), [arXiv:quant-ph/0406196](https://arxiv.org/abs/quant-ph/0406196): the
  tableau engine underneath.
- Dennis, Kitaev, Landahl & Preskill, *Topological quantum memory*,
  [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143): decoding the
  spacetime history, which is where detectors come from rather than from the book.
- Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum **5**, 497 (2021),
  [arXiv:2103.02202](https://arxiv.org/abs/2103.02202); Higgott & Gidney, sparse-blossom
  matching (PyMatching).
- Google Quantum AI, below-threshold surface code (Willow),
  [Nature 2024](https://www.nature.com/articles/s41586-024-08449-y).
- Bravyi et al., bivariate-bicycle qLDPC codes / the Gross code [[144,12,12]],
  [arXiv:2308.07915](https://arxiv.org/abs/2308.07915).
- Kapshikar & Kundu, minimum distance of stabilizer codes is NP-hard,
  [arXiv:2203.04262](https://arxiv.org/abs/2203.04262): why the exact distance search is
  a small-code instrument by construction.
