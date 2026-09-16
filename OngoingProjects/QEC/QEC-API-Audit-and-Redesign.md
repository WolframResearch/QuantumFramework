# QEC Layer: API Audit and Redesign

*Audit of `OngoingProjects/QEC/QECCore/` (Maurice Engelhardt's functional-idiom
`Package[]` on top of the shipped stabilizer engine), against the three concerns
the review is about: too many functions, names that are not the best, and
uncertain signatures/downvalues. This is an update-and-redesign of a good first
implementation, not a rewrite. The starting point is sound: 593 tests pass GREEN,
the physics is cross-checked against dense matrices, the engine, and Stim, and the
functional core is real. The work below is about the public surface, not the
mathematics under it.*

## 0. What this rests on

- **Read in full**, one file at a time: the loader `QECCore.wl`, all 21
  implementation `.wl` files (`GF2`, `Pauli`, `Code`, `Structure`, `Syndrome`,
  `Encoder`, `Circuit`, `Constructions`, `Families`, `Noise`, `ErrorRate`,
  `DetectorModel`, `Cat`, `Memory`, `Stim`, `Cache`, `Register`, `Measurement`,
  `ErrorCorrection`, `Transversal`, `FaultTolerant`), all 20 `.wlt` test files plus
  the runner, and the four planning documents (`QEC-Development-Plan.md`,
  `QEC-Roadmap-2026-09.md`, `README.md`; `EngineMeasurementBug.md` read as
  historical, now fixed in the engine).
- **The 28 public exports and every property list were dumped from a live kernel**
  (WL 15.0, paclet 2.1.0 via `PacletDirectoryLoad`, then `Get` of the loader and
  `Needs["Wolfram`QuantumFramework`QEC`"]`), not inferred from source. Sanity checks
  agree with the code: Steane `{7,1,3}`; bit-flip depolarizing rate
  `(2p(9−9p+4p²))/9`; `QECPauliProduct["X","Z"]` gives the row `{1,1,3}` = −iY;
  the transversal S on Steane reports logical `"Sdg"`;
  `QECPauliMeasurement[...]["MeasurementCorrectQ"]` without EC is `False`.
- **The `Bosonic-QEC-Plan.md` named in the brief does not exist** in the tree or in
  git history. The forward-compatibility target below is taken from the brief's own
  description of it: a future bosonic code object sitting beside the qubit `QECCode`.

Commit history frames the concern precisely. The layer grew in five real commits:
`ebd05e73` (prototype, retired), `ea21d3b6` (the code-layer rebuild, plan items a-d),
`7ec522a2` (cat states, idle noise, post-selection), `187e4e7c` (depth = time steps),
and `4acb4be7` (2026-09-15, the fault-tolerance layer). That last commit added five
files (`ErrorCorrection`, `FaultTolerant`, `Measurement`, `Register`, `Transversal`)
and most of the recent symbols. The surface roughly doubled in one step, which is why
the naming and the count are worth settling now, before the migration into
`QuantumFramework/Kernel/QEC/` freezes it behind reference pages.

---

## 1. Inventory: every public symbol, by module

28 public exports: **10 objects** (a head wrapping one `Association`, with property
dispatch and a summary box), **8 free functions**, **7 Pauli utilities**, **3 global
limits**. Signatures and property lists below are the live ones.

### 1.1 The Pauli layer (`Pauli.wl`): 7 functions

Rows are `{x1..xn, z1..zn, e}` with `e` in Z4 (`i^e`), the same symplectic layout the
engine's `PauliRow` uses (`Kernel/Stabilizer/Conversions.m`); the tests pin row-for-row
agreement with `PauliStabilizer["Stabilizers"|"Destabilizers"|"Matrix"]`.

| Symbol | Signature | Returns |
|---|---|---|
| `QECPauliVector` | `[s_String]`, `[v:{__Integer}]`, `[{__String}]` | row; strings and rows both accepted; `$Failed` + `::invalid`/`::badrow` on bad input |
| `QECPauliString` | `[v:{__Integer}]`, `[s_String]`, `[{__List}]` | string, with leading `-`,`i`,`-i` for the phase |
| `QECPauliQ` | `[s_String]`, `[v:{__Integer}]`, `[_]` | Bool (structural, message-free) |
| `QECPauliWeight` | `[p]` | Integer (phase does not count) |
| `QECPauliCommuteQ` | `[p, q]` | Bool, or `$Failed` + `::size` |
| `QECPauliProduct` | `[ps__]` | row, phase carried in Z4; `$Failed` + `::size` |
| `QECPauliPhase` | `[p]` | `0|1|2|3` |

`PackageScope` here (correctly internal): `pauliQubits`, `symplecticPart`, `phasePart`,
`symplecticProduct`, `pauliIdentity`, `weightKVectors`, `weightOneVectors`,
`$pauliLetterXZ`.

### 1.2 The code object (`Code.wl` + `Structure` + `Syndrome` + `Encoder`): `QECCode`

Constructors:

```
QECCode[{"ZZI","IZZ"}]               (* from Pauli-string generators *)
QECCode[{row, row, ...}]             (* from Pauli rows *)
QECCode["SteaneCode"]                (* named, Families.wl *)
QECCode["Repetition", 5]             (* family member, Families.wl *)
QECCode["CSS", hx, hz]  QECCode["CSS", h]     (* construction, Constructions.wl *)
QECCode["Hamming", r]  QECCode["DistanceTwo", n]  QECCode["PhaseRepetition", n]
QECCode[ps_PauliStabilizer]          (* convert a stabilizer state to a code *)
```

Live `["Properties"]` (34 entries, three groups). Note **`"Decoder"` is listed twice**
(it appears in both the derived and the parametrized group); confirmed in the live dump.

- **Direct (8):** `CheckMatrix`, `Phases`, `Qubits`, `Generators`, `GeneratorVectors`,
  `Signs`, `StabilizerCount`, `LogicalQubits`.
- **Derived (18):** `Parameters` (`{n,k,d}`), `Distance`, `MinimumWeightLogical`,
  `LogicalOperators` (`<|"X"->…,"Z"->…|>`), `LogicalX`, `LogicalZ`, `StandardForm`,
  `CompletedGenerators`, `SyndromeTable`, `Decoder`, `PerfectQ`, `CSSQ`,
  `SyndromeCircuit`, `EncodingGates`, `EncodingCircuit` (a `QuantumCircuitOperator`),
  `EncodingCircuitValidQ`, `PauliStabilizer`, `State` (a `QuantumState`).
- **Parametrized (8):** `["Syndrome", err]`, `["Decode", syn]`, `["Decoder", w_Integer]`
  and `["Decoder", noise_QECNoiseModel]`, `["LogicalErrorRate", noise, (count)]`,
  `["CorrectionCycle", err]`, `["PhysicalCorrectionCycle", err]`,
  `["LogicalPauliQ", p]`, `["StabilizerMemberQ", p]`. Plus `["SyndromeCircuit", rounds]`
  and `["StimCircuit", …]` (added by `Stim.wl`).

Dispatch shape: one hard-coded downvalue per property, `QECCode[a_Association]["X"] := …`,
memoized on internal `PackageScope` symbols via self-blocking `f[a] := f[a] = …`, closed
by a catch-all `QECCode[a_Association][prop_String] := (Message[QECCode::noprop, prop];
Missing["NotFound", prop])`. There is **no `QECCodeQ` guard** on the pattern.

### 1.3 The gadgets and analysis objects: 8 more heads

Each is a head over one `Association` with property dispatch, a `::noprop` catch-all,
and an `ArrangeSummaryBox`. The six circuit-carrying ones share a substrate (see §2.1e).

| Head | Constructors | Options | `["Properties"]` count |
|---|---|---|---|
| `QECSyndromeCircuit` | `[code]`, `[code, rounds]` | none | 12 |
| `QECNoiseModel` | `["Depolarizing"\|"BitFlip"\|"PhaseFlip"\|"BitPhaseFlip", p]`, `[<\|"X"->…\|>]`, `[{pI,pX,pY,pZ}]`, `["Circuit", p\|rates]` | `"MeasurementError"->0` | 13 |
| `QECDetectorModel` | `[code, noise]`, `[code, noise, rounds]` | `"Extraction"->"BareAncilla"` | 21 |
| `QECCatState` | `[m]`, `[qubits, check]` | `"Pairs"->Automatic`, `"Repetitions"->1` | 13 |
| `QECRegister` | `[code]`, `[code, blocks]` | none | 16 |
| `QECPauliMeasurement` | `[code, P]`, `[code, P, reps]` | `"Pairs"`, `"CatRepetitions"`, `"ErrorCorrection"->None` | 25 |
| `QECErrorCorrection` | `[code]`, `[code, offset]` | `"Order"->"BitFlipFirst"` | 21 |
| `QECTransversalGate` | `[code, gate]`, `[code, All]` | none | 12 |
| `QECFaultTolerant` | `[circuit_List, code]` | `"Order"->"BitFlipFirst"` | 22 |

Notable parametrized properties: `QECRegister` carries `["BlockRange", b]`, `["Index", b, q]`,
`["Lift", instr, b]`, `["TransversalCNOT", c, t]`, `["LogicalAction", c, t]`;
`QECTransversalGate` carries `["LogicalGate"]`, `["LogicalAction"]`, `["StabilizerImages"]`;
`QECFaultTolerant` carries `["GadgetInstructions", kind]`, `["Overheads"]`, `["OpenAssumptions"]`.

### 1.4 Free functions (8)

| Symbol | Signature | Returns |
|---|---|---|
| `QECLogicalErrorRate` | `[code, noise, opts]`, `[code, noise, count_Integer, opts]` | **polymorphic**: symbolic polynomial (exact code-capacity), machine number (sampled/numeric), or `Association` `<\|"Rate","Acceptance",…\|>` (circuit level with heralds). Options `"Decoder"`, `"DecoderReach"`, `"Rounds"`, `"Extraction"`. |
| `QECConcatenate` | `[outer_QECCode, inner_QECCode]` | `QECCode` |
| `QECRemoveQubit` | `[code]`, `[code, q]` | `QECCode` |
| `QECPasteCodes` | `[c1, c2, solo1, solo2]` | `QECCode` |
| `QECStimCircuit` | `[code]`, `[code, noise]`, `[code, noise, rounds]` | **String** (Stim source). Option `"Observable"->"Z"`. |
| `QECCodeCatalog` | `[]` | `Dataset` of named codes/families with `{n,k,d}` |
| `QECClassicalHammingMatrix` | `[r_Integer]` | integer matrix (classical parity check) |
| `QECClearCache` | `[]` | Integer (count of dropped memoized rules) |

### 1.5 Global limits (3)

`$QECExactEnumerationLimit` (`4^10`, code-capacity enumeration cap),
`$QECExactDetectorLimit` (`2^16`, circuit-level state-space cap),
`$QECDecoderSubsetLimit` (`2·10^6`, fault-subset cap for the DEM decoder table).

### 1.6 Internal surface worth naming

The `PackageScope` layer is large and mostly healthy: GF(2) primitives (`gf2*`), the
code internals (`codeStandardForm`, `codeLogicalVectors`, `codeMinimumLogical`,
`codeStabilizerElement`, `codeDecoderToWeight`, `codeEncodingGates`, …), the frame
propagator (`framePropagate`, `circuitSchedule`, `circuitIdleSlots`), the detector-model
engine (`codeDetectorModel`, `demExactFailure`, `demSampledFailure`, `demDecoderTable`),
the extraction abstraction (`codeExtraction`, `$codeExtractions`, `recordSyndromes`),
and the FT internals (`transversalAction`, `ftAssemble`, `registerConjugate`). Two
internal facts matter for the redesign:

- **`transversalMatrix` hardcodes a gate-matrix table** `$transversalGateMatrix`
  (I, X, Y, Z, H, S, Sdg, V, Vdg, CNOT, CZ, SWAP) that a test keeps in sync with
  `QuantumOperator[name]["MatrixRepresentation"]`.
- **`$QECMemoisedFunctions` (in `Cache.wl`) is a hand-maintained list** of every
  memoized internal, including the FT-layer ones. `QECClearCache[]` walks it. A new
  memoized function that is not added here silently escapes the cache reset.

---

## 2. Audit against the three concerns

### 2.1 "We have too many functions"

The 10 object heads are **not** the problem: property dispatch over many separate
functions is exactly the QF house pattern (`QuantumOperator` alone exposes 209
properties), and each head is a genuinely distinct object with its own physics. The
189 properties spread across the ten objects are a feature, not bloat. The bloat is in
the flat function/utility surface around the objects.

**(a) The 7 `QECPauli*` utilities are a qubit-specific mini-package leaked into the
public API.** They operate on the internal `{x|z|e}` row representation, which is
deliberately the engine's own layout. A QF user does Pauli work with strings and
`PauliStabilizer`; these seven are the QEC layer's private arithmetic that happens to
be exported. They are also the one part of the public surface that is intrinsically
about qubits, which is the wrong thing to freeze into the top level right before a
bosonic sibling arrives (§4.5). **Recommendation: collapse the seven into one object,
`QECPauli`, with property dispatch, and demote the rest.** `QECPauli["XZZXI"]` becomes
the object; `p["Weight"]`, `p["Phase"]`, `p["String"]`, `p["Vector"]`,
`QECPauli["X"] @ QECPauli["Z"]` (or `p["Times", q]`), `p["CommutesWith", q]`. That is
7 top-level symbols down to 1, in the QF idiom, and it reads better than seven verbs.

**(b) Code-building constructions are split across two mechanisms, inconsistently.**
CSS lives on the head (`QECCode["CSS", hx, hz]`); concatenation, qubit removal, and
pasting are three standalone symbols (`QECConcatenate`, `QECRemoveQubit`,
`QECPasteCodes`). The roadmap gave the last three their `QEC` prefix precisely because
the prototype's bare `RemoveQubit`/`PasteCodes` collided; folding them onto the head
solves the collision the same way CSS already is solved. **Recommendation: make every
code→code construction a named constructor on `QECCode`** (`QECCode["Concatenate",
outer, inner]`, `QECCode["RemoveQubit", code, (q)]`, `QECCode["Paste", c1, c2, s1, s2]`),
matching the family constructors that are already there. That retires three top-level
symbols and makes "how do I build a code" one question with one answer.

**(c) `QECClassicalHammingMatrix` is a classical coding-theory tool in the quantum-QEC
namespace.** It builds a classical parity-check matrix and is used only inside
`QECCode["Hamming", r]`. It is neither quantum nor an object property. **Recommendation:
demote to `PackageScope`.** If a public classical-coding helper is ever wanted, it does
not belong under the `QEC` object prefix.

**(d) Two thin wrappers can move onto the object they describe.** `QECCodeCatalog[]`
is naturally `QECCode["Catalog"]` (or a `$QECCodeNames` list next to the family
constructors). `QECStimCircuit` is discussed under naming below; it can stay a function
but should also be reachable as a property.

**(e) Six heads re-implement one "instruction-carrying gadget" by hand.**
`QECSyndromeCircuit`, `QECCatState`, `QECPauliMeasurement`, `QECErrorCorrection`,
`QECTransversalGate`, and `QECFaultTolerant` each store `<|"Instructions"->…,
"Qubits"->…|>` and each independently defines `"Instructions"`, `"Qubits"`, `"Depth"`,
`"InstructionCount"`, `"GateCounts"`, `"Measurements"`/`"Heralds"`, and a nearly
identical `BarChart` summary box. This is not public bloat (the heads should stay
separate), but it is internal duplication that guarantees drift: today only some
gadgets have `"Heralds"`, only some have `"MeasurementCount"`, and none has a
`QuantumCircuitOperator` view. **Recommendation: one shared internal gadget trait**
(`qecGadgetProp[a, prop]` for the circuit-like properties, plus a shared box builder)
that every gadget head delegates to. Then a new circuit-like property is added once and
appears on all of them.

**Three global limits and `QECClearCache`** are legitimate (`$…Limit` globals are the
`$RecursionLimit` idiom; a cache-reset is a real maintenance need). Keep them. The one
improvement is on the cache mechanism itself (§2.3).

Net effect of (a)-(d): about 28 public symbols down to roughly 16, with no object head
lost and no capability removed.

### 2.2 "Function names are not the best"

Measured against the QF house style, where heads are nouns
(`QuantumState`, `QuantumOperator`, `QuantumChannel`, `QuantumCircuitOperator`,
`PauliStabilizer`, `CliffordChannel`, `StabilizerFrame`, `GraphState`):

- **`QECFaultTolerant` is the one head that is an adjective, not a noun.** It denotes
  `FT(C)`, the fault-tolerant simulation of a circuit. **Rename to a noun:
  `QECFaultTolerantCircuit`** (parallel to `QuantumCircuitOperator`), or `QECProtocol`
  if the intent is the broader Def-10.5 object. `QECFaultTolerant[circuit, code]`
  reads like a predicate; `QECFaultTolerantCircuit[circuit, code]` reads like the
  object it is.
- **`QECStimCircuit` returns a String, but "Circuit" names an object.** Every other
  `…Circuit` in the layer (`QECSyndromeCircuit`, `EncodingCircuit`) is a WL object.
  QF's own foreign-format exporter is `QuantumQASM`, a head, not a `…Circuit`. **Rename
  to `QECStim`** (parallel to `QuantumQASM`) so the name stops promising an object it
  does not return, and additionally expose it as `dem["StimCircuit"]` /
  `QECSyndromeCircuit[…]["Stim"]` for discoverability.
- **`QECCatState` names a gadget, not a state.** The object holds the preparation and
  verification instructions for a cat state, not the state. This is a mild mismatch;
  "cat state" is the accepted name for the GHZ ancilla, so the name is defensible, but
  the object is a `QECCatGadget` in kind. Low priority; keep unless the gadget family
  is renamed as a set.
- **`QECClassicalHammingMatrix`** is verbose and, once internal (§2.1c), moot.
- **`QECPauli*` (seven verbs)** collapse to one noun object `QECPauli` (§2.1a), which
  is both fewer symbols and a better name pattern.

Everything else is well named and idiomatic. `QECCode`, `QECNoiseModel`,
`QECDetectorModel`, `QECSyndromeCircuit`, `QECRegister`, `QECPauliMeasurement`,
`QECErrorCorrection`, `QECTransversalGate`, `QECLogicalErrorRate` are all clear nouns or
clear operations, and the property names inside the objects are consistent and readable
(`"Distance"`, `"LogicalOperators"`, `"DetectorMatrix"`, `"UndetectableFaults"`,
`"OpenAssumptions"`, `"Overheads"`). The property vocabulary is a strength.

### 2.3 "Signatures and downvalues are uncertain"

- **`QECLogicalErrorRate` has a polymorphic return that forces the caller to branch on
  `Head`.** It returns a symbolic polynomial (exact, code capacity), a machine number
  (sampled, or numeric exact), or an `Association` with keys `"Rate"`, `"Acceptance"`,
  `"Failure"` (exact, circuit level with heralds) or `"Rate"`, `"Acceptance"`,
  `"Accepted"`, `"Shots"` (sampled, post-selected). The `Herald.wlt` and `Extraction.wlt`
  suites lean on this and the design reason is honest (a conditional rate must carry its
  acceptance). But a function whose return type depends on the noise level and the
  presence of heralds is exactly the "uncertain signature" the review names. **Fix, the
  QF-idiomatic way: `QECLogicalErrorRate[code, noise, …]` always returns the bare rate**
  (Indeterminate when acceptance is 0), and the acceptance/failure decomposition becomes
  a property surface on `QECDetectorModel`, which already owns `"DetectorRates"`,
  `"ObservableRates"`, `"HeraldRates"`: add `dem["LogicalErrorRate"]`,
  `dem["Acceptance"]`, `dem["Failure"]`. The post-selected `Association` is a
  property-dispatch object waiting to be born; the detector model is its home. If a
  one-call full report is wanted, use the `QuantumLinearSolve[m, b, All]` precedent:
  `QECLogicalErrorRate[code, noise, All]` returns the full `Association`, the bare form
  returns the number. Either way the default return is one type.

- **`QECCode` (and the other heads) dispatch without a validity guard, with a
  hand-maintained property list.** The pattern is `QECCode[a_Association]["Prop"]` with
  no `QECCodeQ` predicate, so a malformed `QECCode[<|garbage|>]["Distance"]` still tries
  to dispatch, and the `"Properties"` list is a literal constant that has already drifted
  (`"Decoder"` appears twice). QF's own heads use `head[prop_?propQ, args] /;
  HeadQ[Unevaluated[head]]` delegating to an internal `HeadProp`, with uniform
  `undefprop`/`failprop` messages and a cache wrapper. **Fix: adopt that shape.** Add a
  `QECCodeQ` (and per-head `…Q`) predicate guard, route through a single internal
  `qecCodeProp[a, prop, args]`, keep the `::noprop` message but add the `::failprop`
  half, and derive `"Properties"` from the actual rule set (or at minimum de-duplicate
  it now). This removes the double `"Decoder"` and makes the surface self-describing.

- **One `Quiet` violates the package's own stated rule.** `Measurement.wl:212` reads
  `v = Quiet[Check[QECPauliVector[p], $Failed]]`, while `GF2.wl:67` states the house
  rule ("a message either matters or should not be raised") and the package is otherwise
  `Quiet`-free and `Print`-free (verified: one `Quiet`, zero `Print`/`Echo`). The intent
  is to catch a bad-Pauli input and re-report it as `QECPauliMeasurement::pauli`, but
  the message-free predicate for exactly this already exists: **guard on `QECPauliQ[p]`
  first, then convert.** Concrete one-line fix.

- **Argument order across the free functions is consistent and correct.**
  `[code, noise, count/rounds, opts]` holds for `QECLogicalErrorRate`,
  `QECDetectorModel`, and `QECStimCircuit`. Option handling uses `OptionsPattern` and
  `OptionValue` throughout. No reordering needed.

- **Idiom, calibrated.** The rebuild is functional-leaning: `With` dominates (101
  occurrences) and the code is built on `Map`/`Fold`/`Table`/`Association` with
  pattern-matched dispatch. It is not maximally idiomatic: `Module` is used 63 times to
  `Block`'s 1, and 28 `Do` / 2 `While` / 2 `AppendTo` remain. Most of that is
  concentrated in algorithms that are inherently imperative (GF(2) echelon reduction
  with pivoting and paired column swaps in `Structure.wl`/`Encoder.wl`, ASAP scheduling
  in `Circuit.wl`), where local mutable state is a defensible choice and `Internal`Bag`
  is already used in place of `AppendTo` in the hot spots. This is polish, not a defect.
  The one idiom item that is actually load-bearing is the `Quiet` above.

---

## 3. Interaction with QuantumFramework and the stabilizer formalism

The layer already composes with the engine in several clean places, and diverges from
it in a few that are worth naming.

### 3.1 Clean, keep and lean on

- **Row layout is the engine's.** `QECPauliVector` produces `Join[xbits, zbits, {e}]`,
  the `PauliRow` format, so a code's rows and a `PauliStabilizer`'s tableau are the same
  objects (`Pauli.wlt` checks it row for row). This is the single best integration
  decision in the layer.
- **Two-way traffic with `PauliStabilizer` and `QuantumState`.** `QECCode[ps]` ingests
  a stabilizer state; `code["PauliStabilizer"]` and `code["State"]` emit one and a
  `QuantumState`. The encoder runs through the engine's compiled bulk path
  `ps["ApplyCircuit", gates]` rather than gate-by-gate (about 25x faster, per the source
  note).
- **`code["EncodingCircuit"]` returns a real `QuantumCircuitOperator`.** A user can
  draw and run it.
- **`noise["QuantumChannel", qubits]` returns a `QuantumChannel`.** The depolarizing
  convention is converted correctly (the model's `p` spread over three Paulis maps to
  the engine's `q = 4p/3`); `Noise.wlt` checks against the engine's own mixture.
- **Syndromes are read as `(1 - ps["Expectation", p])/2`, never `ps["M", q]`,** which
  sidesteps the generator-set-dependence documented in `EngineMeasurementBug.md`. That
  bug is now fixed in the engine (`agExtendToSymplecticBasis` builds genuine
  destabilizers), so the workaround is no longer forced, but reading through
  `"Expectation"` remains a correct and clean choice.

### 3.2 Awkward or missing, worth fixing

- **The flat `{op, q…}` instruction list is a parallel universe to
  `QuantumCircuitOperator`.** Every gadget stores and exposes this list, and **none
  offers a `QuantumCircuitOperator` view or a `["Diagram"]`.** The flat list is the
  right *internal* representation (the frame propagator walks millions of these; a
  `Switch` on a string part is the cheapest dispatch), so the internal choice should
  stay. But there is no bridge out: a user cannot turn a syndrome circuit, a cat gadget,
  or an assembled `FT(C)` into a QF circuit to draw it, compose it, or run it on the
  engine. **This is the highest-value integration addition. Give the shared gadget trait
  a `["QuantumCircuitOperator"]` and a `["Diagram"]`,** built by extending the existing
  `instructionEngineGates` map (which already handles the gate ops; add R -> reset,
  M/MH -> `QuantumMeasurementOperator`). Symmetrically, `QECFaultTolerant` and the
  circuit-consuming entry points should **accept a `QuantumCircuitOperator`** and lower
  it to the flat list, so a user's own circuit composes in rather than having to be
  hand-written as `{{"R",1},{"H",1},…}`.

- **`QECTransversalGate` duplicates the engine's gate matrices.** `$transversalGateMatrix`
  hardcodes X/Y/Z/H/S/V/… and is kept in sync with
  `QuantumOperator[name]["MatrixRepresentation"]` only by a test. **Read the matrices
  from `QuantumOperator` (memoized at first use)** and delete the literal table. The
  test that guards the drift then becomes unnecessary because there is nothing to drift.

- **The phase-carrying Clifford conjugation is a genuine QF capability gap, not just a
  QEC choice.** `transversalAction` recomputes how a Clifford conjugates each Pauli,
  *with the Z4 phase*, by dense matrix and trace overlap, because the engine's frame
  conjugation (`framePropagate`, and `PauliStabilizer`'s own gate updates) drops signs.
  The whole point of the transversal-gate layer is a sign (transversal S is logical
  S†, because `Y^7 = -Ybar`), so the phase cannot be dropped. The engine's
  `StabilizerFrame` already tracks phase-carrying relating Paulis in its `"Paulis"` key,
  so the primitive nearly exists internally. **Longer-term, QF-side: expose a
  phase-carrying "conjugate this Pauli by this Clifford" primitive** from the stabilizer
  engine; then `transversalAction`'s dense-matrix route can be retired. Short-term, keep
  it, but record the dependency so the migration knows what the layer is working around.

- **Noise is simulated by the layer's own frame propagator, in parallel to QF's
  circuit-application path.** This is deliberate and correct for detector-model
  construction (GF(2) frame propagation is the right tool, and it is what makes the exact
  symbolic polynomials possible), so no change is needed. It is worth stating plainly in
  the guide page that the QEC layer runs its own Clifford simulation for noise rather
  than routing through `CliffordChannel`, so a reader does not expect the two to be the
  same code path.

- **Foreign-format returns are fine as they are.** `QECStim` returning a String is the
  right call (it mirrors QASM export); `QECCodeCatalog` returning a `Dataset` is good.
  The only change is the naming and the property-level reachability above.

---

## 4. The redesign

### 4.1 Principles

1. Property dispatch over standalone symbols. Keep the object heads; move utilities and
   constructions onto them.
2. Heads are nouns. Rename the one adjective and the one string-returning "Circuit".
3. One return type per signature. Kill the polymorphic rate; put its decomposition on
   the detector model.
4. Bridge to QF objects at the boundary; keep the fast flat representation inside.
5. Forward-compatible with a bosonic sibling: the shared symbols stay code-type-agnostic;
   the qubit-specific helpers go internal.
6. Update, do not rewrite. Every change below is a rename, a re-home, a guard, or an
   added property. No physics moves.

### 4.2 Keep / merge / internalize / rename

| Current | Action | Becomes |
|---|---|---|
| `QECCode` | **keep**, add guard + fix `"Properties"` + absorb constructions | `QECCode` (qubit stabilizer code) |
| `QECNoiseModel` | keep | `QECNoiseModel` |
| `QECDetectorModel` | keep, add `"LogicalErrorRate"`/`"Acceptance"`/`"Failure"`/`"StimCircuit"` | `QECDetectorModel` |
| `QECSyndromeCircuit` | keep, add `"QuantumCircuitOperator"`/`"Diagram"`/`"Stim"` | `QECSyndromeCircuit` |
| `QECRegister` | keep | `QECRegister` |
| `QECCatState`, `QECPauliMeasurement`, `QECErrorCorrection`, `QECTransversalGate` | keep as heads; share the gadget trait | (unchanged names) |
| `QECFaultTolerant` | **rename** (adjective → noun) | `QECFaultTolerantCircuit` |
| `QECLogicalErrorRate` | keep, **de-polymorphize** the return | `QECLogicalErrorRate` (bare rate; `All` for the report) |
| `QECConcatenate`, `QECRemoveQubit`, `QECPasteCodes` | **merge** onto the head | `QECCode["Concatenate"\|"RemoveQubit"\|"Paste", …]` |
| `QECPauliVector/String/Q/Weight/CommuteQ/Product/Phase` (7) | **merge** into one object | `QECPauli` (+ properties) |
| `QECStimCircuit` | **rename** + also a property | `QECStim` (and `dem["StimCircuit"]`) |
| `QECCodeCatalog` | merge onto head | `QECCode["Catalog"]` |
| `QECClassicalHammingMatrix` | **internalize** | `PackageScope` |
| `QECClearCache` | keep; fix the mechanism (§4.4) | `QECClearCache` |
| `$QECExactEnumerationLimit`, `$QECExactDetectorLimit`, `$QECDecoderSubsetLimit` | keep | (unchanged) |

Public surface: from **28** to about **16** (10 heads + `QECLogicalErrorRate` +
`QECStim` + `QECClearCache` + 3 limits), with `QECPauli` replacing seven and
constructions/catalog folded onto `QECCode`. No object and no capability is lost.

### 4.3 The proposed public API

**`QECCode`: the qubit stabilizer code (the core everything rests on).**

```
QECCode[{gens}] | QECCode[{rows}] | QECCode[ps_PauliStabilizer]     (* construct *)
QECCode[name] | QECCode[name, args]                                 (* named / family *)
QECCode["CSS", hx, hz] | QECCode["Concatenate", outer, inner]
QECCode["RemoveQubit", code, (q)] | QECCode["Paste", c1, c2, s1, s2]  (* constructions, unified *)
QECCode["Catalog"]                                                  (* the Dataset *)
```

Properties unchanged in content, with two fixes: a `QECCodeQ` guard on the dispatch
rule, and a `"Properties"` list derived from the rules (dropping the duplicate
`"Decoder"`). The QF bridges (`"PauliStabilizer"`, `"State"`, `"EncodingCircuit"`) stay.

**`QECPauli`: one object in place of seven functions.**

```
QECPauli["XZZXI"] | QECPauli[row]        (* construct from string or row *)
p["String"] | p["Vector"] | p["Weight"] | p["Phase"] | p["QuditCount"]
p["Times", q]  (or  p ** q)              (* product, Z4 phase carried *)
p["CommutesWith", q]
QECPauliQ[expr]                          (* keep the bare predicate; it is the guard others use *)
```

`QECPauliQ` is the one member worth keeping as a bare function, because it is the
message-free guard the rest of the layer needs (it is the correct fix for the `Quiet`
in §2.3). The other six fold into properties.

**`QECNoiseModel`, unchanged**, with one design note: it carries three levels in one
head via two internal shapes (`"Probabilities"` for code-capacity/phenomenological,
`"Rates"` for circuit). That union is defensible (a code-capacity model is a
phenomenological one with `q=0`, and the same detector machinery consumes both), so keep
it; but if circuit-level rates and per-qubit channels diverge further, splitting
`QECNoiseModel["Circuit", …]` into its own head is the fallback. Keep `"QuantumChannel"`.

**`QECLogicalErrorRate`: one return type.**

```
QECLogicalErrorRate[code, noise, opts]          -> the rate (number or polynomial; Indeterminate if acceptance 0)
QECLogicalErrorRate[code, noise, count, opts]   -> the sampled rate (number)
QECLogicalErrorRate[code, noise, All, opts]     -> <|"Rate","Acceptance","Failure",...|>  (the old Association)
```

and the same decomposition on the detector model: `dem["LogicalErrorRate"]`,
`dem["Acceptance"]`, `dem["Failure"]`. Options unchanged
(`"Decoder"`, `"DecoderReach"`, `"Rounds"`, `"Extraction"`).

**`QECDetectorModel`, unchanged, plus the rate/acceptance/Stim properties** listed in
§4.2, so the memory experiment's conditional outputs live on the object that owns the
detector matrix rather than leaking out of the rate function.

**`QECSyndromeCircuit` and the gadgets** (`QECCatState`, `QECPauliMeasurement`,
`QECErrorCorrection`, `QECTransversalGate`, `QECFaultTolerantCircuit`) keep their
constructors and properties, and all gain, through the shared trait,
`["QuantumCircuitOperator"]` and `["Diagram"]`. `QECFaultTolerantCircuit` and the
syndrome circuit additionally **accept a `QuantumCircuitOperator`** as input.

**`QECStim`** (renamed from `QECStimCircuit`) keeps its signature and `"Observable"`
option and returns a String, and is additionally reachable as `dem["StimCircuit"]`.

**`QECRegister`, `QECClearCache`, the three `$…Limit` globals** are unchanged in name and
signature.

### 4.4 The shared gadget substrate, and the cache

The six circuit-carrying heads should delegate their circuit-like properties to one
internal handler:

```
qecGadgetProp[a_Association, "Instructions"]        := a["Instructions"]
qecGadgetProp[a_Association, "Qubits"]              := a["Qubits"]
qecGadgetProp[a_Association, "Depth"]               := Length[circuitSchedule[a["Instructions"], a["Qubits"]]]
qecGadgetProp[a_Association, "InstructionCount"]    := Length[a["Instructions"]]
qecGadgetProp[a_Association, "GateCounts"]          := Counts[First /@ a["Instructions"]]
qecGadgetProp[a_Association, "Measurements"]        := Count[a["Instructions"], {"M", _}]
qecGadgetProp[a_Association, "Heralds"]             := Count[a["Instructions"], {"MH", _}]
qecGadgetProp[a_Association, "QuantumCircuitOperator"] := gadgetToQuantumCircuit[a]
qecGadgetProp[a_Association, "Diagram"]             := gadgetToQuantumCircuit[a]["Diagram"]
```

Each head keeps its own physics properties (a syndrome circuit's `"MeasurementLabels"`,
a measurement gadget's `"DataWeights"`, a transversal gate's `"LogicalAction"`) and
forwards the shared ones. A single `qecGadgetBox` builds the summary box. Adding a
circuit-like property then happens once.

For the cache, replace the hand-maintained `$QECMemoisedFunctions` list with a single
keyed store: memoize on one internal `Association` (`$qecCache[{function, codeData,
args}]`) rather than on the DownValues of each internal symbol, so `QECClearCache[]`
is `$qecCache = <||>` and no function can escape the reset by being forgotten in a list.
This is the WL "one cache association" pattern and removes the maintenance hazard the
current `Cache.wl` header itself warns about.

### 4.5 Forward compatibility with the bosonic branch

The brief's target is a bosonic code object beside the qubit `QECCode`. The redesign
serves it directly:

- **The shared analysis symbols already have code-type-agnostic names**
  (`QECNoiseModel`, `QECLogicalErrorRate`, `QECDetectorModel`). Keep them generic and
  let them dispatch on the code object's type. A bosonic loss channel is
  `QECNoiseModel["Loss", …]`; the logical-error-rate *question* is the same, even
  though the enumeration under it is not.
- **Internalizing the seven `QECPauli*` functions is what keeps the public layer
  forward-compatible.** They are the only intrinsically qubit-specific symbols in the
  top level. A bosonic code has no Pauli row. Freezing seven qubit-Pauli verbs into the
  public API right before a non-qubit sibling arrives is exactly the thing to avoid;
  folding them into one `QECPauli` object that the qubit `QECCode` uses, and that the
  bosonic code simply does not, is clean.
- **The gadget trait (§4.4) is code-type-agnostic.** "Instructions", "Depth",
  "QuantumCircuitOperator" mean the same thing for any code whose gadgets are circuits.
- Suggested shape when the branch lands: `QECCode` stays the qubit stabilizer code, a
  `QECBosonicCode` (or a bosonic backend selected at construction) sits beside it, and
  both answer `QECLogicalErrorRate[code, noise]`. Nothing in the redesign blocks that;
  the only thing that would have blocked it is qubit-Pauli verbs at the top level, which
  this removes.

### 4.6 Order of work

1. **Cheap and self-contained, do first:** de-duplicate `QECCode["Properties"]`; replace
   the `Quiet` in `Measurement.wl` with a `QECPauliQ` guard; internalize
   `QECClassicalHammingMatrix`. None touches behavior; all three close audit findings
   outright.
2. **The gadget trait and the `QuantumCircuitOperator`/`Diagram` bridge.** Highest user
   value, isolated to the six gadget heads plus one new converter.
3. **De-polymorphize `QECLogicalErrorRate`** and move the acceptance/failure
   decomposition onto `QECDetectorModel`. This changes a return type, so it is the one
   step with a test-migration cost (`Herald.wlt`, `Extraction.wlt`).
4. **Fold the constructions and the catalog onto `QECCode`; collapse the seven Pauli
   functions into `QECPauli`.** Deprecate the old symbols with a one-release alias if
   any notebook depends on them.
5. **Rename `QECFaultTolerant` → `QECFaultTolerantCircuit` and `QECStimCircuit` →
   `QECStim`,** and adopt the QF property-dispatch shape (guard + `undefprop`/`failprop`)
   across the heads.
6. **Read `QECTransversalGate`'s matrices from `QuantumOperator`;** and, on the QF side
   and separately, propose a phase-carrying Clifford-conjugation primitive so the
   dense-matrix `transversalAction` can eventually be retired.

Steps 1-2 are safe to land before the migration into `Kernel/QEC/`. Steps 3-6 are the
ones worth settling before reference pages freeze the names, because every one of them
is a name or a return type a doc page would otherwise pin in place.

---

*The layer is a strong first implementation with the right internal choices: the
engine-compatible row layout, the exact-and-symbolic detector model, the honest
open-assumption accounting, and a functional core with 593 passing cross-checked tests.
The redesign changes none of that. It trims a doubled surface back to its objects,
gives every gadget a way back into a QF circuit, makes one return type per signature,
and takes the qubit-specific helpers out of the top level so the bosonic sibling has
room to stand beside `QECCode` rather than in front of it.*
