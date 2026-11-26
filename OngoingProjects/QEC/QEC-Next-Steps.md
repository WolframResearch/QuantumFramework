# QEC: next steps

*November 2026. What is left to develop in the QEC layer (`OngoingProjects/QEC/`), in the order
recommended, with what each step delivers, what it rests on and how it will be checked.*

---

## Where we are

The layer covers Gottesman's book from chapter 10 to chapter 14 and the memory side of
chapter 17:

- codes as objects, with distance, logical operators, encoders, weight enumerators and the
  Clifford hierarchy;
- noise at three levels (code capacity, phenomenological, circuit), the detector error model,
  and exact or sampled logical error rates, cross-checked against Stim;
- decoders as objects: minimum weight, maximum likelihood, the detector model's own, and a
  native minimum-weight perfect matching that reproduces PyMatching;
- the surface and toric families, `QECPseudothreshold`, and `QECThreshold` both as a family
  crossing and as the concatenation lower bound;
- the fault-tolerant gadgets (verified cats, Steane EC with verified ancillas, transversal
  gates, `FT(C)`), extended rectangles with malignant-pair counting, and the π/8 gate by
  magic-state injection.

831 tests in 31 files, all green against the current engine.

---

## 0. Speed: done in November 2026

The circuit-level computation of `surface-code.nb` (exact rate of the d = 3 surface code,
one noisy round, CNOT depolarizing and readout errors, two decoders) now runs on the layer
alone:

| step | surface-code.nb | the layer |
|---|---|---|
| fault effects | 192 `PauliStabilizer` runs, 3.3 s | `QECDetectorModel`, 0.08 s, the same 352 faults |
| exact rate, single-fault decoder | 1 s | 0.15 s, the same polynomial |
| exact rate, best decoder | pairs decoder | maximum likelihood, 0.3 s, lower rate |

What changed:

- **Exact rates through the Walsh-Hadamard transform** (`Transform.wl`), generalising the
  notebook's trick to any detector model, heralds included. Same polynomial as the old fold,
  50-200x faster; `$QECExactDetectorLimit` goes from 2^16 to 2^24.
- **`"Decoder" -> "MaximumLikelihood"`** for detector models: the optimal decoder, from the
  exact distribution.
- **`"Basis" -> "X" | "Z"`**: a memory in one basis, which halves the bits; it puts the exact
  rate of d = 5 (one round) and d = 3 (five rounds) within reach, about 20 s each (2^25
  patterns, one past the default bound, so `$QECExactDetectorLimit` is raised for them).
- **Bit-parallel frame propagation**: every fault carried through the circuit at once. Building
  the model of d = 7 over seven rounds went from 128 s to 10 s.
- `"Decoder"` with an unknown value is now refused with a message, at the circuit and
  phenomenological levels.
- **The framework's spelling throughout**: parameters inside the name (`QECCode["Surface"[3]]`,
  `QECNoiseModel["Depolarizing"[p]]`, `QECThreshold["Surface"[{3, 5, 7}], ...]`), bare names as
  defaults, `$QECCodeNames` and `$QECNoiseNames`, `Failure["InvalidName" | "InvalidArguments", ...]`
  for bad constructors, and circuits as gate -> qubits rules, as in `QuantumCircuitOperator`. The
  old forms are refused with a message showing the new one.
- **Closer still to the framework**: circuits use the framework's own gate names ("Reset",
  "Measurement", SuperDagger["S"], heralds as labelled QuantumMeasurementOperators), so
  QuantumCircuitOperator accepts them as they are; dec[syndrome] and code[QECPauli[...]];
  Equal by content for codes, noise and Paulis; QuantumTensorProduct for codes; the
  PauliStabilizer names "7QubitCode" and "9QubitCode"; dem["Stim"].
- **Objects like the framework's** (QECCore/Objects.wl): a predicate per head (QECCodeQ, ...)
  backed by System`Private`HoldSetValid, obj["ValidQ"] and obj["AllProperties"], ::undefprop
  with the call left unevaluated for an unknown property, a Failure tagged with the message's
  name for every refusal that has a message, and memoisation that obeys
  $QuantumFrameworkPropCache.
- **Retired names removed**: Deprecated.wl is gone (the ten old functions, code["StimCircuit"]
  and dem["StimString"]); an old property name gets ::undefprop like any unknown one. The last
  silent $Failed on public paths (matching on a code it cannot pair, a bad noise for the logical
  channel, a threshold with no buildable location) are Failures with their message; the bare
  $Failed left are internal sentinels, all checked with FailureQ.
- Still to align: usage messages in Usage.m and reference pages, at the paclet migration
  (they live in QuantumFramework/Kernel/).

Where it stops now: exact up to about 2^26 patterns (memory, ~4 GB); sampling with matching
handles d = 7 over seven rounds (20 000 shots in 7 s); building the model of d = 9 over nine
rounds takes 45 s, which is the next thing to speed up (the per-fault bookkeeping after the
propagation, which can be done on the whole matrix at once).

---

## 1. Small fixes and housekeeping

Short, and worth doing before anything new.

| item | what to do | size |
|---|---|---|
| Silent fallback in `"Decoder"` | done for detector models (circuit and phenomenological levels); the code-capacity path still falls back to its default | a few lines and one test |
| `EngineMeasurementBug.md` | the engine now gives the right deterministic outcome on generator-built states (re-tested November 2026); mark the note resolved and update the roadmap | documents only |
| Clifford conjugation proposal | `QF-Clifford-Conjugation-Proposal.md` is still open on the engine side (`CliffordChannel[QuantumOperator[...]]` and `cc["Conjugate", P]` do not exist); remind Mads, or close it | a conversation |
| Decisions A–D | the four provisional decisions of the API redesign (`QEC-API-Redesign-Plan.md`) are waiting for review | a review |
| `QEC-String-Inputs.pdf` | decide whether the reference of string inputs belongs in the repo | none |

---

## 2. A standard surface-code circuit

**Why.** The circuit-level threshold of the surface code comes out at about 0.27 % here,
against the 0.5–1 % usually published. One known reason is the extraction: an X check is
measured by rotating its data qubits with H before and after the CNOTs, which puts extra
fault locations on the data every round. The standard circuit prepares the X-check ancilla in
|+⟩, measures it in X, and leaves the data alone; it also runs all checks of a round in four
layers of CNOTs.

**What.** An `"Extraction"` choice (for example `"Standard"`) that emits that circuit, used by
the detector model, the Stim writer and the rates, with `code["CheckOrder"]` deciding the
order of the four layers.

**How it is checked.** Detector rates against Stim on the new circuit (within a few σ, as for
the existing cases); the threshold of the surface family at circuit level, d = 3, 5, 7,
compared with the published range.

**Size.** Small to medium. It touches `Circuit.wl` and `Measurement.wl`; the rest of the
layer consumes instruction lists and needs no change.

---

## 3. qLDPC codes and a decoder for them

**Why.** Bivariate bicycle codes (IBM's [[144, 12, 12]] "gross" code among them) are the
current alternative to the surface code, and they have been on the roadmap since the start.
Matching cannot decode them, since one error trips more than two checks.

**What.**

- `QECCode["BivariateBicycle", l, m, A, B]` (or by name for the published members), built as
  a CSS code from two commuting matrices;
- a belief-propagation decoder with ordered-statistics post-processing (BP-OSD), as a
  `QECDecoder` method usable at code capacity and on the detector model.

**How it is checked.** The published parameters [[n, k, d]] of the listed members (k by rank;
d by search where it is feasible, by an upper-bound witness where it is not); BP-OSD against
maximum likelihood on codes small enough for both; logical error rates against the published
curves.

**Size.** Medium for the codes, large for BP-OSD. The family is a few dozen lines in
`Families.wl`; the decoder is a new file.

---

## 4. Moving into the paclet

**Why.** Today the layer is loaded with a separate `Get`, and `?QECCode` or F1 shows nothing.
This is what stands between the work and a release.

**What.**

- move `QECCore/*.wl` into `QuantumFramework/Kernel/QEC/` and register the contexts in
  `PacletInfo.wl`; the tests move to `Tests/QEC/` and the local runner goes away;
- one reference page per public function (23) and a guide page; the tech notes move to the
  paclet's documentation tree and go through `DocumentationBuild`;
- harvest the older scattered QEC notebooks (Courses, the 5/7/9-qubit examples).

**How it is checked.** The paclet's own test battery, the documentation build, and the tech
notes evaluating in a fresh kernel from the installed paclet.

**Size.** Medium, mostly mechanical. **Waits on a decision**: nothing goes into
`QuantumFramework/Kernel/` until it is approved.

---

## 5. The rest of the book

In rough order of value.

| topic | book | size | what it adds |
|---|---|---|---|
| Shor and Knill error correction | sec. 12.2, 12.4 | medium | fault-tolerant EC for codes that are not CSS (the 5-qubit code has none today) |
| GPP, GCP, ECRP, ECCP as checkable properties | sec. 10.2 | medium | each fault-tolerance property tested on its own, rather than through the residual weight |
| `FT(C)` with T gates | sec. 13.3 | large | universal circuits end to end; the S correction depends on an outcome, so it needs a Clifford frame, not a Pauli frame |
| Magic-state distillation | sec. 13.5 | large | the T gate as it is done in practice |
| Malignant triples, codes with t > 1 | sec. 14.5 | medium | a slightly better bound for the 7-qubit code; a bound at all for codes correcting two errors |
| The \*-decoder and level reduction, built | sec. 14.4, 14.6 | large | level 2 of concatenation as a circuit, not only as a consequence |
| Noise beyond Pauli | ch. 15 | varies | leakage, biased noise, slow measurement, nearest-neighbour gates |
| Lattice surgery and logical gates on the surface code | sec. 17.3 | large | the surface code as a computer, not only a memory |

---

## 6. Faster matching

The native matcher is exact but builds a complete graph on the lit detectors and an all-pairs
distance matrix; it is comfortable up to d ≈ 7. Sparse blossom (or union-find, which is
approximate but near-linear) would take it to the sizes where PyMatching is used today. Until
then, `QECDecoder[code, f]` wraps PyMatching through `ExternalEvaluate`.

---

## Recommended order

1. The small fixes of section 1.
2. The standard surface-code circuit (section 2): little work, and a number comparable with
   the literature.
3. qLDPC codes and BP-OSD (section 3): the most current part of the field, and on the roadmap
   since the start.
4. The move into the paclet (section 4), as soon as it is approved.
5. The rest of the book (section 5) and faster matching (section 6), by interest.

Standing rules, unchanged: nothing goes into `QuantumFramework/Kernel/` without approval; the
`.wlt` suite runs on bare Wolfram; documents are in English; every new number is checked
against an independent source (Stim, PyMatching, the book, or a direct simulation).
