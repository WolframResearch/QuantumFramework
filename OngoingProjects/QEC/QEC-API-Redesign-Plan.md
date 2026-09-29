# QEC API Redesign: Plan of Work

*Response to `QEC-API-Audit-and-Redesign.md`. The audit's observations are accepted; this
document says in what order they are carried out, what each step touches, how each one is
verified, and the four design decisions that should be settled before the third step
starts. Nothing here changes the mathematics: the 593 tests are the invariant every step
is measured against. (All eight steps are done; the suite now stands at 685.)*

## 0. What this plan accepts, and the four constraints it works under

**Accepted in full.** Every observation of the audit is answered by a step below; §1 is the
traceability table. Two of them were verified against the source before planning: the
`"Decoder"` duplicate is real (the key appears in both `$codeDerivedProperties` and
`$codeParametrizedProperties`, so the concatenated `"Properties"` list shows it twice), and
the `Quiet` at `Measurement.wl:212` is real and does contradict the house rule stated in
`GF2.wl`. The inventory is also correct: 28 `PackageExport`, 21 implementation files plus
the loader, 19 `.wlt` files plus the runner.

Four constraints shape *how* the work is done. They are not objections to the audit; three
of them are the audit's own reasoning applied to places where it stops short.

**C1 — The mathematics does not change.** Every step ends with 593 green from a terminal
`wolframscript` run. A step that requires a test to change must say why in its commit
message; "the test now expects the new name" is a reason, "the number moved" is not.

**C2 — The exact symbolic rate is the one capability nobody else has, and it survives.**
The layer returns `3p² − 2p³`, `32p/15` against `112p/15`, a closed-form rational function
at circuit level — computed by enumeration over cosets and by a fold over detector effects,
with `p` left symbolic. The audit's §4.2 makes the rate a functional of the logical channel,
and that is compatible with this, *provided the logical channel is built the right way*.
There are two ways to build it, and only one keeps the rate symbolic:

- **By composing physical channels**, `D ∘ (Σ R_s ∘ M_s) ∘ N ∘ E` taken literally. The noise
  on the 7-qubit code is then a `4⁷ × 4⁷` superoperator, about 2.7×10⁸ entries, each a
  polynomial in `p`. That does not fit in memory, and every product grows the expressions.
  It works at `n = 3` and dies long before `n = 7`.
- **From the engine the layer already has.** For Pauli noise on a stabilizer code the logical
  channel is *always* a Pauli channel on the `k` logical qubits: each physical error `E` has a
  syndrome `s`, the recovery `R_s` is a Pauli, so the residue `R_s·E` lies in the normalizer —
  a logical Pauli times a stabilizer — and its whole effect on the encoded state is that
  logical Pauli. So the channel is `4ᵏ` numbers, the probability of each logical residue, and
  those are exactly the coset probabilities the maximum-likelihood decoder already computes,
  exactly and symbolically.

Checked before writing this down: the Steane code under depolarizing noise, `p` symbolic,
comes out as a logical Pauli channel of four exact degree-7 polynomials in 0.07 s, and
`1 − q_Ī` agrees with `QECLogicalErrorRate` identically (the difference simplifies to `0`).
So **`code["LogicalChannel", noise]` is built from the coset / detector-model engine as a
logical Pauli channel, never by composing `4ⁿ` physical channels, and the rate is a
functional of it with nothing lost.** The one real limit is non-Pauli noise (amplitude
damping, say): there the residue is no longer a logical Pauli, the channel is no longer a
Pauli channel, the coset engine does not apply, and the only exact route is the dense one,
for small codes only. The audit's §3.4 ("the operational objects are the interface and the
semantics, not the computational representation") is the same argument; this constraint
makes it concrete for the channel the rate is read from.

**C3 — The public surface is a face, not an engine.** The layer runs on raw symplectic rows
`{x₁…xₙ, z₁…zₙ, e}` in the frame propagator, the detector model and the label matrix, and
those stay raw. `QECPauli` wraps a row; it does not replace it. Same rule as C2, applied to
the Pauli layer.

**C4 — Renames land in one commit, not in a trickle.** Every rename ripples into 593 tests,
three tech notes (~250 evaluated cells) and two rebuilt notebooks. Ten small rename commits
cost ten notebook rebuilds and ten chances to leave a stale hint. Step 5 is therefore one
step, done once.

---

## 1. Traceability: every observation to a step

| # | Audit observation | § | Action | Step | Risk | Status |
|---|---|---|---|---|---|---|
| 1 | `QECCode["Properties"]` lists `"Decoder"` twice | 1 | de-duplicate the listing | 1 | none | done |
| 2 | `Measurement.wl:212` uses a forbidden `Quiet` | 1 | replace with a `QECPauliQ` guard | 1 | none | done |
| 3 | `QECClassicalHammingMatrix` is a classical helper in a quantum namespace | 2.2, 5 | internalize to `PackageScope` | **5** (moved, see below) | low | done: `PackageScope`, alias warns |
| 4 | gadgets return flat instruction lists, not QF objects | 4.2, 5 | shared circuit trait; `"QuantumCircuitOperator"`, `"Diagram"` | 2 | low | done: `Gadget.wl` |
| 5 | `QECSyndromeCircuit` should expose its measurement | 4.2 | `"SyndromeMeasurement"` -> `QuantumMeasurementOperator` | 3 | medium | done: `sc["SyndromeMeasurement"]` |
| 6 | a code should present the operational objects | 3.3, 5 | `"Encoder"`, `"Codespace"`, `"SyndromeMeasurement"`, `"Recovery"`, `"LogicalChannel"` | 3 | medium | done: plus `"Codewords"`, `"KnillLaflammeMatrix"`, `"CorrectableQ"` |
| 7 | noise is a channel and should be consumable as one | 4.2, 5 | general Pauli `QuantumChannel`; accept a channel wherever a noise model is accepted | 3 | medium | done: `QECNoiseModel[qc]`, `noise["QuantumChannel"]` for every Pauli model |
| 8 | `QECLogicalErrorRate` is polymorphic in its return | 4.2, 5 | one return type; report moves to the DEM; `All` for the full report | 4 | medium | done: `[…, All]`, `dem["LogicalErrorRate" | "Acceptance" | "Failure"]` |
| 9 | there is no first-class decoder object | 2.3, 5 | add `QECDecoder` over `QECDetectorModel` | 4 | low | done: four constructors incl. the function seam |
| 10 | 7 `QECPauli*` verbs should be one object | 2.2, 5 | `QECPauli` (+ `QECPauliQ` kept as the guard) | 5 | low, wide | done: `QECPauli`, six aliases |
| 11 | `QECFaultTolerant` is an adjective | 2.2 | rename to a noun (decision D) | 5 | low, wide | done: `QECFaultTolerantCircuit` |
| 12 | `QECStimCircuit` returns a String | 2.2 | `QECStim` + `dem["StimString"]` | 5 | low | done |
| 13 | `QECCodeCatalog` is a free function | 5 | `QECCode["Catalog"]` | 5 | low | done |
| 14 | the code -> code constructions are morphisms | 3.1, 5 | see decision B | 5 | low | done: functions kept, properties added |
| 15 | `QECTransversalGate` hard-codes gate matrices | 5 | read them from `QuantumOperator`; propose the QF primitive | 6 | low | done: read from the engine; proposal in `QF-Clifford-Conjugation-Proposal.md` |
| 16 | operator-algebra QEC as the frame for the bosonic branch | 3.5 | cross-check against the now-existing `Bosonic-QEC-Plan.md` | 3 (design only) | none | done: same names adopted (step 3 note) |

---

## 2. The steps

### Step 0 — Baseline (half a day)

Branch off `main`. Record, from a terminal run rather than the MCP evaluator: `593/593
GREEN`, the three notebooks built and balanced (40 / 91 / 116 cells), and the reference
values the audit itself used (Steane `{7,1,3}`; bit-flip depolarizing rate
`(2p(9−9p+4p²))/9`; `QECPauliProduct["X","Z"]` = `−iY`; transversal S on Steane = `Sdg`).
Those five numbers are the regression suite for the whole redesign: they are re-checked at
the end of every step and must not move.

### Step 1 — The three one-liners (1 hour)

Behaviour-preserving, no decision required, do first.

- **De-duplicate `"Decoder"`.** The two entries are two different call shapes
  (`code["Decoder"]`, the table; `code["Decoder", noise]`, the coset map), so the fix is in
  the listing, not the dispatch: list it once, and make `"Properties"` report the
  parametrized ones distinguishably (this is also the moment to decide whether
  `"Properties"` should return a flat list or an association of the three groups it already
  maintains internally).
- **Replace the `Quiet`.** `v = Quiet[Check[QECPauliVector[p], $Failed]]` becomes a
  `QECPauliQ` guard with the message the house rule wants — which is what `QECPauliQ` was
  written for.
- ~~**Internalize `QECClassicalHammingMatrix`**~~ — **moved to step 5**, on evidence. It has
  one caller inside the package (`Families.wl`), but six uses in the tests and two evaluated
  cells in `StabilizerCodes.md`, where it is the argument of the CSS constructor. So it is
  not a one-liner: it is a surface change with documentation reach, which is what step 5 is
  for, and folding it in there costs one notebook rebuild instead of two. The doc cells
  should be rewritten with the literal 3x7 parity-check matrix rather than another helper —
  showing the matrix is better documentation of what the CSS constructor takes.

*Verify:* 593 green; new tests asserting the `QECPauliMeasurement::pauli` message now
surfaces where the `Quiet` used to swallow it. **Surface: unchanged at 28** (the one symbol
this step was going to remove moved to step 5).

**Done.** The guard also tightened two latent holes the `Quiet` was hiding: the old size
condition (`Length[v] < 2n`) accepted a Pauli on *more* qubits than the code, and the
`MatchQ[v, {(0|1)..}]` test rejected a phased Pauli by accident rather than by decision.
Both are now explicit conditions with the same message, and each has a test:
`QEC-Measure-refuses-what-is-not-a-Pauli`, `-a-Pauli-of-the-wrong-size`,
`-a-Pauli-carrying-a-phase`. `QECCode["Properties"]` went from 34 entries with `"Decoder"`
twice to 33 with it once. **596/596 GREEN.**

### Step 2 — The gadget trait, and the circuit bridge (1–2 days)

The single highest-value addition, and it is isolated. The template already exists in the
code: `QECCode["EncodingCircuit"]` is one line,
`QuantumCircuitOperator[codeEncodingGates[a]]`. What is missing is the same line for the
seven circuit-bearing heads, off a shared substrate rather than seven copies.

- Factor the circuit properties (`"Instructions"`, `"Depth"`, `"InstructionCount"`,
  `"GateCounts"`, `"Qubits"`) into one trait used by `QECSyndromeCircuit`, `QECCatState`,
  `QECPauliMeasurement`, `QECErrorCorrection`, `QECTransversalGate`, `QECFaultTolerant`
  and the encoder.
- Add `["QuantumCircuitOperator"]` and `["Diagram"]` to all of them.
- Accept a `QuantumCircuitOperator` on the way in, wherever an instruction list is accepted.

Two real problems to solve rather than discover later. **`instructionEngineGates` drops `R`
and `M`** — it was written for unitary validation only, and a drawn circuit that silently
loses its resets and measurements is worse than no drawing; the bridge needs them as
`QuantumMeasurementOperator` and reset. And **`Sdg`/`Vdg` expand to three gates** for the
engine; a faithful drawing shows three boxes where the instruction list has one. Decide
once: draw the instruction, annotate the expansion.

*Verify:* for each gadget, the circuit operator applied to the encoded state must reproduce
the engine result the existing tests already assert (`ErrorCorrection.wlt`,
`Transversal.wlt`, `FaultTolerant.wlt` all have such a test). Any divergence is a bridge
bug, caught immediately.

**Done**, in `QECCore/Gadget.wl` (the trait and the bridge) and `Tests/Gadget.wlt` (10
tests). What the framework turned out to accept settled the two open questions: `"Reset" -> q`
and `"Measurement" -> q` exist, so `R` and `M` survive the translation; `"M" -> q` is read by
`QuantumCircuitOperator` as the *identity*, which is exactly the silent lie this bridge had
to avoid; `QuantumOperator["S", {q}]["Dagger"]` is one box with the right matrix, so `Sdg`
and `Vdg` draw as themselves rather than as three gates; and a herald measurement carries
its own `"MH"` label, so a reader can see which readouts reject the shot. `ftLocations` now
delegates to `gadgetLocations`, so "a location is an instruction or a wait" is written once.
**606/606 GREEN.**

### Step 3 — The operational objects (1 week; the core of the redesign)

This is where the design earns its keep, and where decisions A and C must already be
settled.

- **`code["Encoder"]` -> `QuantumOperator`**, the isometry `V`, computed from
  `codeEncodingGates` (a Clifford tableau / gate list) and materialized densely only on
  request.
- **`code["Codespace"]` / `["Projector"]` -> `QuantumOperator`**, `P = V·V†`.
- **`code["SyndromeMeasurement"]` -> `QuantumMeasurementOperator`**, the instrument `{M_s}`;
  also exposed from `QECSyndromeCircuit` (observation 5).
- **`code["Recovery", s]` -> `QuantumChannel`**, from the decoder table already computed.
- **`code["LogicalChannel", noise]` -> `QuantumChannel`**, the composed `L`, built from the
  coset / detector-model engine as a logical Pauli channel of `4ᵏ` probabilities (constraint
  C2) — symbolic, and cheap at every size the exact rate already handles.
- **Noise as a channel.** `noise["QuantumChannel"]` today returns `Missing` for anything but
  a *named* channel — `noiseChannelName` covers depolarizing / bit-flip / phase-flip and
  nothing else. Build the general Pauli channel from the four probabilities, then accept a
  bare `QuantumChannel` wherever a `QECNoiseModel` is accepted.
- **A size guard.** These properties materialize `4ⁿ`-sized objects. Add
  `$QECDenseQubitLimit` alongside the three existing limits, with a message that names the
  cheap alternative (the detector model) rather than just refusing.

*Verify:* three new tests that are physics, not plumbing. The **Knill–Laflamme condition**
`P Eᵢ†Eⱼ P = λᵢⱼ P` computed from `["Codespace"]` and a noise channel, on the 3-qubit code,
must agree with the layer's own correctability verdict. The **logical channel of the
bit-flip code at code capacity** must reproduce `3p² − 2p³` — the same polynomial the
enumeration route gives, which is the cross-check that the new face and the old engine
agree. And `["Encoder"]` applied to `|0…0⟩` must give the state `code["State"]` already
returns.

*Design note, no code:* cross-check §3.5 of the audit (operator-algebra QEC as the frame
for subsystem, gauge and bosonic codes) against `Bosonic-QEC-Plan.md`, which did not exist
when the audit was written and does now. If the bosonic plan's object is the same
interface, say so in one paragraph; if it is not, that is a finding worth reporting back.

**Done**, in `QECCore/Operational.wl` and `Tests/Operational.wlt` (30 tests), with the
noise side in `Noise.wl`. What exists now:

| Property | Returns | Dense? |
|---|---|---|
| `code["Encoder"]` | `QuantumOperator`, the isometry `V` (k → n qubits) | yes |
| `code["Codewords"]` | the `2ᵏ` codewords as `QuantumState`s, in logical basis order | yes |
| `code["Codespace"]` | `QuantumOperator`, `P = V·V†` | yes |
| `code["SyndromeMeasurement"]` | `QuantumMeasurementOperator`, one projector per syndrome; outcome `i` is syndrome `IntegerDigits[i−1, 2, m]`. Also `QECSyndromeCircuit[…]["SyndromeMeasurement"]` | yes |
| `code["Recovery", s]` | `QuantumChannel` of the decoder's correction for `s` | yes |
| `code["LogicalPauliProbabilities", noise]` | `<\|"I" -> q_I, "X" -> …\|>`, the `4ᵏ` weights of the logical channel, exact and symbolic | no |
| `code["LogicalChannel", noise]` | `QuantumChannel` on the `k` logical qubits built from those weights | no |
| `code["KnillLaflammeMatrix", errors]` | the matrix `h_ab`, or `Missing["NotCorrectable", …]` naming the pair whose product is a logical | no |
| `code["CorrectableQ", errors]` | `True` / `False`; `errors` is a list or `"Weight"[t]` | no |

`noise` may be a `QECNoiseModel` or a one-qubit Pauli `QuantumChannel` everywhere, including
`QECLogicalErrorRate`. `QECNoiseModel[qc]` reads a channel back (Kraus → Pauli weights, then
a superoperator check) and refuses a non-Pauli one (`::notpauli`, e.g. amplitude damping);
`noise["QuantumChannel"]` now works for every Pauli model, not only the named ones. The
dense objects stop at `$QECDenseQubitLimit` (default 10) with a message naming the cheap
route. The logical channel is code-capacity only: at circuit level it returns
`Missing["NotAvailable", "Circuit"]` and points to the detector model (step 4 is where that
side gets its single return type).

The verifications the step asked for, as they came out: the logical channel's
`1 − q_I` equals `QECLogicalErrorRate` identically on the Steane code (symbolic, 0.07 s);
bit flip gives `q_X = (3 − 2p)p²`, i.e. the rate `3p² − 2p³`; the Knill–Laflamme matrix is
computed symplectically and checked against its dense definition `⟨W_i|E_a†E_b|W_j⟩` on
the 3-qubit code, and reports weight-1 correctable / weight-2 not on Steane, weight-1
correctable on the 5-qubit code. (`"Projector"` was not added as a second name for
`"Codespace"`: one name per object. The check "`Encoder` on `|0…0⟩` equals `code["State"]`"
was replaced by the projector cross-check, since `code["State"]` is only defined up to the
logical freedom when `k > 0`.) **636/636 GREEN.**

Two findings, recorded here rather than fixed at their source:

- *The encoding circuit's logical labelling is its own.* On the 5-qubit code
  `codeEncodingGates` sends `|0…0⟩` to the −1 eigenstate of the `Z̄` that
  `code["LogicalZ"]` reports. It is still a valid encoder, but an isometry built directly
  from it would relabel `|0_L⟩ ↔ |1_L⟩` against the code's own logical operators. The
  operational layer canonicalizes (`|0_L⟩` is moved into the +1 eigenspace of every `Z̄_j`,
  `|c_L⟩ = X̄^c|0_L⟩`), and a test pins `V†Z̄V = Z`, `V†X̄V = X` on five codes.
  `Encoder.wl` is untouched.
- *Engine: a one-element Kraus list with an order acts on the wrong wires.*
  `QuantumChannel[{X}, {1}]["Order"]` is `{{0}, {1}}` and `QuantumChannel[{P}, {1,2,3}]` is
  `{{0,1,2}, {1,2,3}}`: the environment wire takes the place of the first system wire, so
  a unitary channel given this way silently acts on the wrong qubits. The layer routes a
  single Kraus operator through `QuantumChannel[QuantumOperator[…]]` instead
  (`krausChannel`), and two tests pin it. Worth reporting upstream; nothing in the kernel
  was changed.

*Design note — the bosonic plan.* `Bosonic-QEC-Plan.md` and its Part A spec describe the
same frame the audit's §3.5 asks for: a code is a subspace (codewords / isometry), an
error set, a correctability question, and a logical channel. The qubit layer now answers
that frame with the same spellings the bosonic spec chose — `"Codewords"`,
`"KnillLaflammeMatrix"`, `"CorrectableQ"` — so the two branches share an interface rather
than two dialects. What differs is only how each is computed: symplectically here, by
truncated Fock-space linear algebra there. Nothing in the bosonic spec had to change.

### Step 4 — Return types, and the decoder (3–4 days)

- **De-polymorphize `QECLogicalErrorRate`.** One return type: a number or a symbolic
  polynomial, `Indeterminate` when acceptance is zero. The post-selected association moves
  onto the detector model as `dem["LogicalErrorRate"]`, `dem["Acceptance"]`,
  `dem["Failure"]` — where its siblings `"DetectorRates"` and `"ObservableRates"` already
  live — and the full report stays reachable as `QECLogicalErrorRate[code, noise, All]`, on
  the `QuantumLinearSolve[m, b, All]` precedent.
- **Add `QECDecoder`**, constructed from a detector model or from `[code, noise]`, with
  `dec["Decode", syndrome]`. The three decoders that exist inside today — minimum weight
  over the lookup table, maximum likelihood over cosets, and the DEM's lightest-set rule —
  become its built-ins, and the seam where PyMatching / BP-OSD plug in is documented rather
  than implied.

*Verify:* the `::noacceptance` path of `Herald.wlt` becomes an assertion about
`dem["Acceptance"]`; every rate test states which of the two shapes it expects. The
decoder's built-ins must reproduce, decision for decision, what the current code path
produces — including the herald-row filter, which is what moved the measured exponent from
1.73 to 1.965 and is easy to lose in a refactor.

**Surface: 27 -> 28** (one added object).

**Done**, in `QECCore/Decoder.wl` (new), `ErrorRate.wl`, `DetectorModel.wl`, and
`Tests/Decoder.wlt` (20 tests).

- `QECLogicalErrorRate[code, noise]` and `[code, noise, n]` now return **the rate only** —
  a number, or an exact expression, `Indeterminate` (with `::noacceptance`) when a herald
  rejects everything. `QECLogicalErrorRate[code, noise, All]` / `[code, noise, n, All]`
  give the report `<|"Rate", "Acceptance", "Failure"|>` (sampled:
  `<|"Rate", "Acceptance", "Accepted", "Shots"|>`), with the same keys whether or not
  anything was post-selected; without heralds the acceptance is 1 and the failure is the
  rate. The one test that read `["Rate"]` off the old association
  (`QEC-Extraction-restores-the-p-squared`) now reads the number directly.
- The detector model answers `dem["LogicalErrorRate"]`, `dem["LogicalErrorRate", All]`,
  `dem["LogicalErrorRate", n]`, `dem["Acceptance"]`, `dem["Failure"]`; the exact report is
  computed once and shared (memoised, cleared by `QECClearCache`).
- `QECDecoder` has four constructors: `QECDecoder[code]` (minimum weight),
  `QECDecoder[code, noise]` (maximum likelihood at code capacity; the detector-model
  decoder at the other two levels, with `"Rounds"` and `"Extraction"`),
  `QECDecoder[dem]` (lightest fault set), and `QECDecoder[code, f]` — **the seam**: `f`
  maps a syndrome to a correction Pauli or `Missing[]`, so PyMatching or BP-OSD plug in
  through `ExternalEvaluate` without the layer knowing. `dec["Decode", s]`,
  `dec["Table"]`, `dec["Method"]`, `dec["Reach"]`. The rate accepts it as
  `"Decoder" -> dec`, so a custom decoder is scored by the same machinery as the built-ins.

Every built-in delegates to the engine function that already decoded; the tests compare
tables, not rates. The herald filter is pinned twice: the object's table equals the
engine's filtered table on the 5-qubit transversal extraction, **and** that table differs
from the one built with the heralds ignored, so the filter is shown to be doing something.
**656/656 GREEN** (see the final count in step 7).

### Step 5 — The renames, all at once (2–3 days)

One commit, one notebook rebuild, aliases for one release.

- **`QECPauli`** absorbs the seven verbs; `QECPauliQ` stays public as the guard. Internally
  the row functions keep their names in `PackageScope` (constraint C3). Measured cost: 73
  uses inside `QECCore/`, 135 in the tests, 32 in the docs — mechanical, but wide.
- **`QECFaultTolerant` -> a noun** (decision D).
- **`QECStimCircuit` -> `QECStim`**, plus `dem["StimString"]`.
- **`QECCodeCatalog` -> `QECCode["Catalog"]`.**
- **The code -> code constructions** (decision B).
- Deprecation aliases with a message, to be removed one release later.

*Verify:* 593 green with the new names; a test per alias asserting the old name still works
and warns. Rebuild the three notebooks at the end of the step, not during it.

**Surface: 28 -> about 16.**

**Done**, in one pass over `QECCore/`, `Tests/`, the three tutorial sources and
`StimCrossCheck/emit.wls`, with `QECCore/Deprecated.wl` (new) and `Tests/Surface.wlt` (29
tests, new).

| Before | After |
|---|---|
| `QECPauliVector`, `QECPauliString`, `QECPauliWeight`, `QECPauliPhase`, `QECPauliCommuteQ`, `QECPauliProduct` | `QECPauli[p]` with `"Vector"`, `"String"`, `"Weight"`, `"Phase"`, `"CommuteQ"`, `"Product"` (and `**`), plus `"Qubits"`, `"Support"`, `"HermitianQ"`, `"Matrix"`, `"QuantumOperator"`; `QECPauliQ` stays as the guard |
| `QECFaultTolerant` | `QECFaultTolerantCircuit` (decision D) |
| `QECStimCircuit` | `QECStim`, plus `dem["StimString"]` |
| `QECCodeCatalog[]` | `QECCode["Catalog"]` |
| `QECClassicalHammingMatrix` | `PackageScope`; `QECCode["Hamming", r]` builds the quantum code, and the doc cells write the matrix out |
| `QECConcatenate`, `QECPasteCodes`, `QECRemoveQubit` | unchanged (decision B), plus `code["Concatenate", inner]`, `code["Paste", other, s1, s2]`, `code["RemoveQubit"]`, `code["RemoveQubit", q]` |

How the rename was kept honest: the six verbs became the `PackageScope` row functions
`pauliVector` … `pauliPhase` under their old definitions (constraint C3), so every internal
call site changed name and nothing else; the messages moved to `QECPauli::invalid`,
`::badrow`, `::size`. The test files call the row functions through `Symbol[...]` the way
they already reach other internals, so no expected value changed. The retired ten names
live in `Deprecated.wl`: each returns exactly what it returned and warns once per session
(`General::qecdeprecated`); `Surface.wlt` checks value and warning for all ten, and that the
second use is silent. Decision C held: `code["Generators"]` still returns strings.

Public functions after the step: 19 (plus 4 `$` limits), from 26 — not the audit's "about
16", because steps 3–4 added `QECDecoder` and the step keeps the three constructions as
functions. **685/685 GREEN.**

### Step 6 — Transversal matrices from `QuantumOperator`, and the QF gap (2 days)

`Transversal.wl` carries its own `$transversalGateMatrix` table. The tests already pin it
against `QuantumOperator[name]["MatrixRepresentation"]`, so reading it from the engine
directly is safe and removes a duplicated convention.

The audit's step 6 also asks for something QF does not have: a **phase-carrying Clifford
conjugation primitive**. That is the one genuine capability gap the layer works around —
`transversalAction` builds dense matrices and reads coefficients off traces because nothing
in the framework conjugates a Pauli through a Clifford *with its `ℤ₄` phase*. Propose it to
the QF side as a small, separable addition; if it lands, `transversalAction` retires.

**Done.** `$transversalGateMatrix` is now read from `QuantumOperator[name]` (and
`["Dagger"]` for `Sdg`, `Vdg`) at load time instead of transcribed; the 24 transversal tests
pass unchanged. The capability gap is written up as a proposal to the QF side in
`QF-Clifford-Conjugation-Proposal.md`: `CliffordChannel` already holds the phase-carrying
conjugation table for a Clifford unitary, and what is missing is a constructor from a Clifford
`QuantumOperator` and a `cc["Conjugate", P]` query. Nothing in the kernel was changed.

### Step 7 — Documentation and close (2 days)

Rewrite the affected sections of the three tech notes, `validate.wls`, `build.wls`, update
the roadmap, and re-walk the audit's own tables point by point so the response is
checkable. The `FaultTolerantGadgets` note is the most affected, since it is written
entirely around the six gadget objects.

**Done.** `StabilizerCodes` gained two sections, *The code as quantum operations* (encoder,
codewords in the code's own basis, syndrome instrument and recovery, the logical channel from
the coset engine and its agreement with the rate, noise read from a `QuantumChannel`,
Knill–Laflamme) and *Decoders as objects* (the built-ins, the function seam, the rate's `All`
report), plus the post-selected report where the transversal extraction is measured and six
new rows in its checks table. `QECCoreInternals` gained *Operational.wl*, *Decoder.wl* and
*Deprecated.wl*. All three notes use the new names throughout (`QECPauli`,
`QECFaultTolerantCircuit`, `QECStim`, the Hamming matrix written out). `validate.wls` runs
all 262 cells with the only non-results being the two deliberate refusals the notes show;
the three notebooks are rebuilt. The roadmap and the README carry the new counts. The audit's
observations are walked point by point in the status column of §1.

**Final: 685/685 GREEN**, 23 test files, on bare Wolfram. Nothing committed; nothing in
`QuantumFramework/Kernel/`.

---

## 3. Four decisions to settle before step 3 starts

*Status: work proceeds on the four recommendations below, adopted provisionally so steps 3–7
are not blocked. Any of them that review changes is revisited afterwards; each is isolated
enough (A and C in step 3, B and D in step 5) that reversing one does not touch the others.*

**A. Is the rate computed from the logical channel?** *Recommendation: yes, as the audit
proposes — with one condition on how that channel is built.* `code["LogicalChannel", noise]`
is constructed from the coset / detector-model engine as a logical Pauli channel (`4ᵏ`
exact probabilities), never by composing the `4ⁿ`-dimensional physical channels; the rate is
then `1 − q_Ī`, a functional of that channel, and stays an exact symbolic polynomial. The
Steane check in constraint C2 shows the two agree identically. Non-Pauli noise is the one case
that falls back to the dense route, and it should be refused beyond small codes rather than
attempted. *(An earlier revision of this plan recommended against the audit here, on the
grounds that a dense channel cannot carry the symbolic rate. That was the right worry aimed
at the wrong target: the dense construction is what cannot, and it is not the only one.)*

**B. Do `QECConcatenate`, `QECPasteCodes` and `QECRemoveQubit` become properties of
`QECCode`?** *Recommendation: keep them as free functions, add the properties as sugar.*
The audit's own argument is that they are **morphisms** (Cowtan–Burton), and a morphism
between two codes is naturally spelled as a function of two codes; `QECCode["Concatenate",
other]` privileges one argument over the other for no reason, and concatenation is not a
property of the inner code. Folding them onto the head does reduce the symbol count, which
is the audit's goal — the sugar gets that without the asymmetry.

**C. Does `code["Generators"]` start returning `QECPauli` objects instead of strings?**
*Recommendation: no, and say so explicitly.* It is the natural next step of observation 10
and it is a much wider break than the seven functions: everything that compares generators
against strings stops working, in the package, the tests and both older notes. Readable
string output is a feature. If it is wanted, it belongs in a later release with its own
migration, not inside step 5.

**D. `QECFaultTolerantCircuit` or `QECFaultTolerantProtocol`?** *Recommendation: `Circuit`,
for now.* The object is one assembled circuit plus its gadget table and overheads, and
`QECFaultTolerantCircuit` is parallel to `QuantumCircuitOperator`. `Protocol` is the right
noun for the Definition-10.5 object — a code plus a full set of gadgets — which the layer
does not yet have as a thing. Leave the name free for it.

---

## 4. Out of scope, stated so it is not assumed

- **The physics.** No result moves: `4 -> 1` residual weight, exponent `1.00 -> 1.97`,
  transversal S = logical `S†`, the three overheads. If one moves, the step is wrong.
- **Phase G** (exRecs, benign/malignant sets, the rigorous threshold bound). Orthogonal to
  every step here; it can proceed in parallel or after, and its only interaction is that it
  should be written against the post-step-4 return types if it comes later.
- **The migration into `QuantumFramework/Kernel/QEC/`.** It comes *after* this work, which
  is the audit's own timing argument: names and return types should be settled before
  reference pages freeze them.
- **The bosonic branch.** Its own plan, its own kickoff. This work touches it only through
  the design note in step 3.

---

## 5. Order of merit

| Step | What it buys | Risk | Effort | Blocked by |
|---|---|---|---|---|
| 0 | a baseline the rest is measured against | none | ½ day | — |
| 1 | the three defects the audit names | none | 1 hour | — |
| 2 | every gadget becomes drawable and composable | low | 1–2 days | — |
| 3 | the object-model alignment; the operational face | medium | ~1 week | decisions A, C |
| 4 | one return type; a first-class decoder | medium | 3–4 days | step 3 |
| 5 | 28 symbols -> ~16; the field's names | low but wide | 2–3 days | decisions B, D |
| 6 | one convention less; a QF gap named | low | 2 days | step 2 |
| 7 | the documentation matches the surface | none | 2 days | steps 1–6 |

Steps 0–2 are safe to start immediately: no decision depends on them and nothing they touch
is contested. Steps 3–6 are the ones the audit correctly says must be settled before the
migration, and they are also the ones where the four decisions above should be agreed in
writing first — half a day of agreement now against a week of rework later.

*Total: about three weeks of work, of which one week is the object model and three days are
mechanical renaming. The audit's judgement that this is the right moment — after the
fault-tolerance layer doubled the surface, before the paclet migration freezes it — is the
part of it that is least arguable.*
