# QEC Layer: API Audit and Redesign

*Audit and redesign of the quantum-error-correction layer at
`OngoingProjects/QEC/QECCore/` (Maurice Engelhardt's functional-idiom `Package[]` on
top of the shipped stabilizer engine). This revision goes past the public-surface trim
of the first pass and answers the three questions the review is actually about: are the
function names right (checked against the QEC literature and the reference packages, not
just against Wolfram house style); is the design built on the minimal set of objects that
quantum theory itself uses for error correction; and, given that, what should the
signatures be and what should each function return. The physics under the layer is sound
(593 cross-checked tests pass), so nothing below touches the mathematics. It is about the
shape of the objects and how they compose.*

## 0. What this rests on

- **The full layer, read in place:** the loader and all 21 implementation files, all 20
  test files, and the planning documents, plus a live-kernel dump of the 28 public
  exports and every property list (WL 15.0, paclet 2.1.0). Sanity values agree with the
  code (Steane `{7,1,3}`; bit-flip depolarizing rate `(2p(9−9p+4p²))/9`;
  `QECPauliProduct["X","Z"]` = `{1,1,3}` = −iY; transversal S on Steane = logical `Sdg`).
- **The QEC literature and the reference packages, cross-checked** (details and citations
  in §6). Two threads matter. The **operational-formalism** thread (Knill-Laflamme-Viola;
  Kribs-Laflamme-Poulin operator QEC; Rahn-Doherty-Mabuchi; Cowtan-Burton) says what the
  minimal objects of error correction are. The **software** thread (Gidney's Stim;
  QuantumClifford.jl and its `QECCore.jl` interface package; PyMatching; the Derks-Eisert
  detector-error-model paper; the QUITS simulator) says how a well-built QEC framework is
  layered and what its objects are named. Both were read against the current design.
- **The naming corpus** was taken from the live API of Stim (`stim.Circuit`,
  `stim.DetectorErrorModel`, `stim.PauliString`, `stim.Tableau`, `stim.TableauSimulator`,
  `stim.FlipSimulator`), from QuantumClifford.jl / `QECCore.jl` (`parity_checks`,
  `code_n`/`code_k`/`code_s`, `distance`, `logx_ops`/`logz_ops`, `naive_syndrome_circuit`,
  `naive_encoding_circuit`, `PauliFrame`, `PauliError`, `evaluate_decoder`), and from the
  term-of-art usage in the DEM literature.
- **The `Bosonic-QEC-Plan.md` in the brief does not exist** in the tree or git history;
  the forward-compatibility target is taken from the brief and, below, given a firmer
  footing than "a sibling object" (operator-algebra QEC, §3.5).

The layer grew in five commits, doubling its surface in the last one (`4acb4be7`, the
fault-tolerance layer). That is the right moment to settle names, objects, and return
types, before the migration into `QuantumFramework/Kernel/QEC/` freezes them behind
reference pages.

---

## 1. Inventory: the 28 public symbols

10 property-dispatch objects, 8 free functions, 7 Pauli utilities, 3 global limits.

**Objects** (a head over one `Association`, with a `::noprop` catch-all and a summary box):

| Head | Constructors | Options | Props |
|---|---|---|---|
| `QECCode` | `[{gens}]`, `[{rows}]`, `[name]`, `[name,args]`, `["CSS",hx,hz]`, `[ps_PauliStabilizer]` | none | 34 |
| `QECNoiseModel` | `["Depolarizing"\|"BitFlip"\|"PhaseFlip"\|"BitPhaseFlip",p]`, `[<\|"X"->…\|>]`, `[{pI,pX,pY,pZ}]`, `["Circuit",p\|rates]` | `"MeasurementError"->0` | 13 |
| `QECSyndromeCircuit` | `[code]`, `[code,rounds]` | none | 12 |
| `QECDetectorModel` | `[code,noise]`, `[code,noise,rounds]` | `"Extraction"` | 21 |
| `QECCatState` | `[m]`, `[qubits,check]` | `"Pairs"`, `"Repetitions"` | 13 |
| `QECRegister` | `[code]`, `[code,blocks]` | none | 16 |
| `QECPauliMeasurement` | `[code,P]`, `[code,P,reps]` | `"Pairs"`, `"CatRepetitions"`, `"ErrorCorrection"` | 25 |
| `QECErrorCorrection` | `[code]`, `[code,offset]` | `"Order"` | 21 |
| `QECTransversalGate` | `[code,gate]`, `[code,All]` | none | 12 |
| `QECFaultTolerant` | `[circuit_List,code]` | `"Order"` | 22 |

**Free functions:** `QECLogicalErrorRate[code,noise,(count),opts]` (polymorphic return),
`QECConcatenate`, `QECRemoveQubit`, `QECPasteCodes` (code -> code), `QECStimCircuit`
(returns a String), `QECCodeCatalog` (Dataset), `QECClassicalHammingMatrix` (classical
matrix), `QECClearCache`.

**Pauli utilities (7):** `QECPauliVector`, `QECPauliString`, `QECPauliQ`, `QECPauliWeight`,
`QECPauliCommuteQ`, `QECPauliProduct`, `QECPauliPhase`, over the engine-compatible
`{x1..xn,z1..zn,e}` row (Z4 phase).

**Global limits (3):** `$QECExactEnumerationLimit`, `$QECExactDetectorLimit`,
`$QECDecoderSubsetLimit`.

Two live defects worth carrying forward: `QECCode["Properties"]` lists `"Decoder"` twice,
and `Measurement.wl:212` uses a `Quiet` the package's own house rule forbids (the
message-free `QECPauliQ` is the intended guard). Both are one-line fixes.

---

## 2. Axis 1: are the function names right?

Measured against the field, not just against Wolfram style. The comparison corpus is the
two reference implementations whose designs are cleanest and most used, Stim (the de facto
substrate) and QuantumClifford.jl with its `QECCore.jl` interface package, plus the
term-of-art of the detector-error-model literature.

### 2.1 What the field calls these things

| Concept | Stim | QuantumClifford / QECCore.jl | Literature term of art | Current QEC |
|---|---|---|---|---|
| Pauli operator | `PauliString` (one object) | `PauliOperator` (one object) | Pauli / stabilizer generator | **7 free functions** `QECPauli*` |
| stabilizer code | (not an object; a circuit) | code types `Steane7`, `Shor9`, `Toric`, `Surface` | stabilizer code, `[[n,k,d]]` | `QECCode` |
| check matrix | none | `parity_checks`, `parity_matrix_x/z` | (parity-)check matrix | `code["CheckMatrix"]` |
| code parameters | none | `code_n`, `code_k`, `code_s`, `distance` | `[[n,k,d]]` | `code["Parameters"]`, `["Distance"]` |
| logical operators | none | `logx_ops`, `logz_ops` | logical (X/Z) operators | `code["LogicalX"]`, `["LogicalZ"]` |
| encoding circuit | none | `naive_encoding_circuit` | encoding circuit | `code["EncodingCircuit"]` |
| syndrome-extraction circuit | (in the `Circuit`) | `naive_syndrome_circuit` | syndrome extraction circuit | `QECSyndromeCircuit` |
| noise / error model | (baked into `Circuit`) | `PauliError`, `UnbiasedUncorrelatedNoise` | noise model / error model | `QECNoiseModel` |
| detector error model | **`DetectorErrorModel`** | (via Stim) | **detector error model (DEM)** | `QECDetectorModel` |
| logical error rate | (sampled from the model) | `evaluate_decoder` | logical error / failure rate | `QECLogicalErrorRate` |
| decoder | (external: PyMatching) | `AbstractSyndromeDecoder`, `TableauDecoder` | decoder | `code["Decoder"]` |

The reading is clear and mostly favorable: **the layer's object names are the field's
names.** `QECDetectorModel` is exactly Stim's `DetectorErrorModel`, which the Derks-Eisert
paper calls the standardized interface between a circuit and a decoder; adopting that
noun verbatim is a strength, not an accident. `QECSyndromeCircuit` matches
`naive_syndrome_circuit` (and, correctly, it *is* the non-fault-tolerant bare-ancilla
construction that `naive_` names). `QECNoiseModel`, `QECLogicalErrorRate`, `["CheckMatrix"]`,
`["Parameters"]`, `["Distance"]`, `["LogicalX"]`/`["LogicalZ"]` all sit on standard
terms. The property vocabulary inside the objects (`"UndetectableFaults"`,
`"ObservableRates"`, `"OpenAssumptions"`, `"Overheads"`) is precise and reads well.

### 2.2 The names that the field says to change

- **The 7 `QECPauli*` verbs should be one `QECPauli` object.** This was the first pass's
  recommendation on Wolfram-idiom grounds; the field settles it. Both reference
  implementations model a Pauli as a single object with methods (`PauliString` in Stim,
  `PauliOperator` in QuantumClifford), never as a spray of free functions. `QECPauli`,
  with `p["Weight"]`, `p["Phase"]`, `p["String"]`, `p["Vector"]`, `p["CommutesWith", q]`,
  and product via `p ** q`, is the two-package-precedented shape, and it removes the only
  qubit-specific symbols from the top level (which matters for §3.5).
- **`QECFaultTolerant` is the one adjective head.** No package has a settled noun here
  (Gottesman writes `FT(C)`, "the fault-tolerant simulation of a circuit"). Rename to a
  noun: `QECFaultTolerantCircuit`, parallel to `QuantumCircuitOperator`, or
  `QECFaultTolerantProtocol` for the broader gadget-set object.
- **`QECStimCircuit` returns a String, so "Circuit" over-promises.** The field's verb for
  this is "export to Stim"; Wolfram's own foreign-format head is `QuantumQASM`. Rename to
  `QECStim` and additionally expose it as a property (`dem["StimString"]`), the way Stim
  itself is reached by a method, not a constructor.
- **`QECClassicalHammingMatrix`** is a classical-coding helper in the quantum namespace;
  internalize it (§5).

### 2.3 A naming gap the field exposes

The current layer has no first-class **decoder** object, only `code["Decoder"]`
(a lookup table) and `code["Decoder", noise]` (a coset map). Every serious framework makes
the decoder a first-class, swappable thing: Stim hands its DEM to PyMatching's `Matching`;
QuantumClifford has `AbstractSyndromeDecoder` with `TableauDecoder`, belief-propagation,
and matching implementations behind it; QUITS makes the inner decoder a plug-in. The
literature is explicit that the DEM is the interface precisely so that decoders are
interchangeable behind it. A `QECDecoder` object (even if the only built-ins are the
lookup table and maximum-likelihood, with a documented seam where PyMatching/BP-OSD plug
in) would match the field and is the natural consumer of `QECDetectorModel`. This is an
addition, not a rename, and it belongs with the object-model work of §3.

---

## 3. Axis 2: is the design built on the minimal objects of quantum theory?

This is the load-bearing question. QuantumFramework is organized the way quantum theory
is: a `QuantumState` is a state, a `QuantumOperator` transforms states, a `QuantumChannel`
is a more general transformation, a `QuantumMeasurementOperator` is a measurement, and a
`QuantumCircuitOperator` composes them, so that `op[state]`, `channel[state]`, and
`qc1 /* qc2` all mean what they say. The QEC layer does not join that world. Its objects
are associations with property dispatch; a code does not transform anything, noise is a
separate model rather than a channel you apply, syndrome extraction is a flat instruction
list rather than a measurement, and the correction cycle returns a report rather than a
map. The question the brief asks is whether QEC *can* be expressed in the same object
language. It can, and the reason is that error correction was defined in that language in
the first place.

### 3.1 What quantum theory says the objects are

The operational formulation of quantum error correction, which is the standard one, is
built from exactly the primitives QuantumFramework already ships. A code is an **encoding
isometry** `V : H_L -> H_P` from k logical qubits into n physical ones; its **code space**
is the range of the **projector** `P = V·V†`; physical noise is a **channel** `N` (a CPTP
map); syndrome extraction is a **quantum instrument** `{M_s}`, a measurement that returns
an outcome s and the post-measurement state; and recovery is a **syndrome-indexed family
of channels** `{R_s}`. The whole cycle is a single composition, and the thing it produces
is again a channel, the **effective logical channel**

```
    L  =  D ∘ ( Σ_s  R_s ∘ M_s ) ∘ N ∘ E ,
```

where E is encode (apply V), N is noise, `M_s` is the syndrome branch, `R_s` the recovery
for that syndrome, and D decode (V†). This is Knill-Laflamme-Viola and Rahn-Doherty-Mabuchi
(§6). Two facts from that formulation matter for a software design:

- **The Knill-Laflamme condition** `P E_i† E_j P = λ_ij P` (for noise Kraus operators
  `{E_i}`) is the exact-correctability test, and it is a statement purely about the
  encoder V (through P) and the noise channel N. Correctability is a relation between two
  objects the framework already has.
- **Concatenation is literal composition** of these logical channels: the concatenated
  code's logical channel is the composition of the inner and outer ones. Code constructions
  are not bespoke matrix surgery; they are composition of maps.

Cowtan and Burton make the last point exact for CSS codes: a code *is* an object (a chain
complex over GF(2)), and the constructions the layer already implements are the standard
categorical operations on those objects. Concatenation, the direct sum of two codes, CSS
merging and splitting (lattice surgery), and the LDPC balanced product are all colimits
(pushouts and coequalizers) in the category of codes; a code map is a morphism, and a
morphism corresponds to a physical circuit. In their words the layer's `QECConcatenate`,
`QECPasteCodes`, `QECRemoveQubit`, and CSS constructor are all **morphisms between code
objects**, and morphisms compose.

### 3.2 QuantumFramework already has every one of these objects

| Operational object | QuantumFramework object |
|---|---|
| state `ρ` | `QuantumState` |
| encoding isometry `V : H_L -> H_P` | `QuantumOperator` (rectangular / isometric) |
| code-space projector `P = V·V†` | `QuantumOperator` (a projector), or the subspace itself |
| physical noise `N`, recovery `R_s` | `QuantumChannel` |
| syndrome instrument `{M_s}` | `QuantumMeasurementOperator` (retains outcome + post-state) |
| composition `∘` and the cycle | `QuantumCircuitOperator`, `@`, `/*` |

So the minimal object set for QEC introduces no new primitive type. It is the QF primitive
set, specialized. The layer today touches this in exactly one place, and it is telling:
`code["EncodingCircuit"]` already returns a genuine `QuantumCircuitOperator`, and
`noise["QuantumChannel"]` already returns a `QuantumChannel`. The bridges exist; they are
just not the spine of the design.

### 3.3 What the code object should expose

Make `QECCode` present the operational objects, so the cycle composes with the rest of QF:

```
code["Encoder"]                 -> QuantumOperator     (* the isometry V *)
code["Codespace"] / ["Projector"] -> QuantumOperator   (* P = V·V†, the subspace *)
code["SyndromeMeasurement"]     -> QuantumMeasurementOperator   (* the instrument {M_s} *)
code["Recovery", syndrome]      -> QuantumChannel       (* R_s *)
code["LogicalChannel", noise]   -> QuantumChannel       (* the effective L above *)
```

The first two make a code a subspace you can encode into and project onto; the third makes
syndrome extraction a measurement a user can apply to any `QuantumState`; the last two make
recovery and the whole cycle channels a user can compose, plot, or hand to `QuantumEvolve`.
`code["LogicalChannel", noise]` is the operational payoff: it turns the code-plus-noise into
one QF channel, and `QECLogicalErrorRate` becomes a functional of it (one minus the process
fidelity of L to the identity, or the logical-flip probability). This is the object the
whole layer is implicitly about, and right now it has no name.

### 3.4 The efficiency caveat, stated honestly

There is a real tension, and the design must resolve it rather than ignore it. A
`QuantumChannel` or `QuantumMeasurementOperator` on n qubits is a `4^n`-sized object, and
the entire reason the stabilizer and detector-error-model engines exist is to stay
polynomial. So the operational objects are the **interface and the semantics**, not the
computational representation. The rule is the one QF already lives by: a `QuantumOperator`
may be symbolic or lazy, and `PauliStabilizer` is the efficient backend; likewise a code's
`"Encoder"` is computed and stored as a Clifford tableau or a gate list, materialized as a
dense `QuantumOperator` only when n is small enough to ask for its matrix. The detector
error model is, in fact, already exactly this: an efficient, GF(2) representation of the
composed logical channel L restricted to the syndrome and logical bits (this is the
Derks-Eisert reading of the DEM). So the layer is *already computing* the operational
composition; it is presenting it as a bag of matrices and rates instead of as the channel
it is. Exposing the QF objects is a change of face, not a change of engine.

### 3.5 Why this also fixes the bosonic forward-compatibility

Operator-algebra quantum error correction (Kribs-Laflamme-Poulin; Poulin's stabilizer
formalism for it; Dauphinais-Kribs-Vasmer) is the generalization that replaces the encoding
isometry and the code projector by a **protected subsystem** or, most generally, a
**correctable observable algebra**. It unifies subspace codes (what the layer does today),
subsystem and gauge codes, decoherence-free subspaces, and noiseless subsystems, and it is
the natural frame for the planned bosonic branch, because a bosonic or continuous-variable
code is an algebra of operators, not a set of qubit Paulis. If the shared objects are
`QECCode` presenting an encoder, a codespace, a syndrome measurement, and a logical
channel, then a `QECBosonicCode` sibling answers the same interface with a different
backend, and `QECLogicalErrorRate[code, noise]` reads a channel it does not have to know
the qubit-ness of. The forward-compatibility argument for internalizing the qubit-Pauli
helpers (§2.2) is the shallow version of this; operator-algebra QEC is the deep one.

### 3.6 The verdict on axis 2

Yes, QEC can follow QuantumFramework's object model, and it should, for three reasons that
the literature makes concrete rather than aesthetic. The operational definition of error
correction is already a composition of state, isometry, channel, and instrument, which are
QF's primitives (§3.1, §3.2). Presenting them lets the QEC cycle compose with the rest of
the framework, so a user encodes a `QuantumState`, pushes it through a `QuantumChannel`,
and measures with a `QuantumMeasurementOperator`, instead of learning a separate vocabulary
of associations. And the efficient engine the layer already has (stabilizer tableaux, the
detector error model) is precisely the polynomial realization of that composition, so the
alignment costs nothing at run time (§3.4). The current design has the right *layering*
(it matches Stim/QUITS/Derks-Eisert almost module for module) but the wrong *object
semantics*: it computes channels and measurements and presents them as bags.

---

## 4. Axis 3: signatures, and what each function returns

Once the objects are the operational ones, the signatures and returns mostly write
themselves, and they line up with how the reference packages read.

### 4.1 The principle the packages share

Stim's pipeline is method chaining on objects: `circuit.detector_error_model()`,
`sampler = circuit.compile_detector_sampler()`, `matching.decode(syndrome)`.
QuantumClifford's is `evaluate_decoder(decoder, setup, nsamples)`. Both build a small
number of typed objects and pass them to each other. QuantumFramework's idiom is stronger
than either, because its objects are *applied* and *composed*: `V[state]`, `N[state]`,
`qmo[state]`, `channel1 /* channel2`. So a QEC-in-QF pipeline should read as application and
composition, which is more natural in Wolfram than in Python and is the whole point of
matching the object model.

### 4.2 Constructors take algebra or QF objects; returns are QF objects

- **`QECCode`** constructs from the algebraic data (generators, check matrix, name,
  family) as today, and additionally from native QF objects: `QECCode[ps_PauliStabilizer]`
  already works; add `QECCode[V_QuantumOperator]` (a code *is* its encoding isometry) and
  `QECCode[P_projector]`. Its operational properties return QF objects (§3.3), not
  associations.
- **`QECNoiseModel`** stays the constructor for named/parametric Pauli noise, but the
  function that consumes noise should accept either a `QECNoiseModel` or a bare
  `QuantumChannel`, since noise *is* a channel: `QECLogicalErrorRate[code, channel]` where
  `channel` is either. This is the operational statement that the code does not care how
  the noise was named, only what map it is.
- **`QECLogicalErrorRate`** returns **one type**: a number or a symbolic polynomial (the
  rate), `Indeterminate` when acceptance is zero. The post-selected `Association`
  (`"Rate"`, `"Acceptance"`, `"Failure"`, …) that it returns today only at circuit level is
  the current return-type wart; move that decomposition onto `QECDetectorModel`
  (`dem["LogicalErrorRate"]`, `dem["Acceptance"]`, `dem["Failure"]`), which already owns
  `"DetectorRates"` and `"ObservableRates"`, and offer the full report through the
  `QuantumLinearSolve[m, b, All]` precedent: `QECLogicalErrorRate[code, noise, All]`. The
  bare call then always returns a rate, definable cleanly as a functional of
  `code["LogicalChannel", noise]`.
- **The gadgets** (`QECSyndromeCircuit`, `QECCatState`, `QECPauliMeasurement`,
  `QECErrorCorrection`, `QECTransversalGate`, `QECFaultTolerantCircuit`) keep the fast flat
  instruction list inside, but every one gains `["QuantumCircuitOperator"]` and
  `["Diagram"]` (the highest-value single addition, from the first pass), so the circuits
  return a QF object a user can draw, run, or compose, and the circuit-consuming entry
  points accept a `QuantumCircuitOperator` in.
- **`QECSyndromeCircuit`** additionally returns its measurement as
  `["SyndromeMeasurement"]` -> `QuantumMeasurementOperator`, closing the loop with §3.3.

### 4.3 The signature table, corrected and uniform

| Function | Signature | Returns |
|---|---|---|
| `QECCode` | `[gens\|rows\|name\|ps\|V\|P, args]` | `QECCode` |
| `QECCode["Encoder"\|"Codespace"\|"SyndromeMeasurement"]` | property | `QuantumOperator` / `QuantumMeasurementOperator` |
| `QECCode["Recovery", s]`, `["LogicalChannel", noise]` | property | `QuantumChannel` |
| `QECNoiseModel` | named / rates / `["Circuit", …]` | `QECNoiseModel` (yields a `QuantumChannel`) |
| `QECLogicalErrorRate` | `[code, noise\|channel, (count), opts]`, `[…, All]` | rate (number/polynomial); `All` -> report `Association` |
| `QECDetectorModel` | `[code, noise, (rounds), opts]` | `QECDetectorModel`; `dem["LogicalErrorRate"\|"Acceptance"\|"Failure"]` |
| `QECDecoder` (new) | `[dem]` or `[code, noise]` | `QECDecoder`; `dec["Decode", syndrome]` -> correction |
| `QECStim` (renamed) | `[code, (noise), (rounds), opts]` | `String` |
| gadget heads | `[…]` | object; `["QuantumCircuitOperator"]`, `["Diagram"]` |

Argument order is already uniform and correct in the current layer (`code`, then
`noise`/`channel`, then `count`/`rounds`, then options); the change is in the return types,
not the argument lists.

---

## 5. The redesign, in two tiers

**Tier A, the surface trim** (from the first pass, still valid and now with external
precedent). Fold the four code -> code constructions onto `QECCode` as named constructors,
so building a code is one question with one answer; collapse the seven `QECPauli*` verbs
into one `QECPauli` object (Stim/QuantumClifford precedent); internalize
`QECClassicalHammingMatrix`; move `QECCodeCatalog` onto `QECCode["Catalog"]`; rename
`QECFaultTolerant -> QECFaultTolerantCircuit` and `QECStimCircuit -> QECStim`;
de-duplicate `QECCode["Properties"]`; replace the `Quiet` in `Measurement.wl` with a
`QECPauliQ` guard; share the six gadgets' circuit-property substrate and give it the
`QuantumCircuitOperator`/`Diagram` bridge; de-polymorphize `QECLogicalErrorRate`. This
takes the public surface from 28 symbols to about 16 with no capability lost.

**Tier B, the object-model alignment** (new, from §3-§4, and the more important half).
Make `QECCode` expose the operational objects: `["Encoder"]` and `["Codespace"]` as
`QuantumOperator`s, `["SyndromeMeasurement"]` as a `QuantumMeasurementOperator`,
`["Recovery", s]` and `["LogicalChannel", noise]` as `QuantumChannel`s. Let noise be
consumed as a `QuantumChannel`. Add a first-class `QECDecoder` object as the consumer of
`QECDetectorModel` (§2.3), matching the field's decoder abstraction and marking the
PyMatching/BP-OSD seam. Treat the code -> code constructions as the morphisms they are
(Cowtan-Burton), so they compose. Keep operator-algebra QEC as the frame that makes the
subspace layer, a future subsystem/gauge layer, and the bosonic branch answer one
interface.

### Keep / merge / internalize / rename / add

| Current | Action | Becomes |
|---|---|---|
| `QECCode` | keep; add guard, fix `"Properties"`, absorb constructions, **expose operational objects** | `QECCode` |
| `QECNoiseModel` | keep; **also consumable as a `QuantumChannel`** | `QECNoiseModel` |
| `QECDetectorModel` | keep; add `"LogicalErrorRate"`/`"Acceptance"`/`"Failure"`/`"StimString"` | `QECDetectorModel` |
| `QECSyndromeCircuit` + 5 gadgets | keep; share the trait; add `"QuantumCircuitOperator"`/`"Diagram"`; add `"SyndromeMeasurement"` | (same, `QECFaultTolerant` renamed) |
| `QECFaultTolerant` | rename (adjective -> noun) | `QECFaultTolerantCircuit` |
| `QECLogicalErrorRate` | keep; de-polymorphize; accept a channel | `QECLogicalErrorRate` |
| `QECConcatenate`, `QECRemoveQubit`, `QECPasteCodes` | merge onto the head (they are code morphisms) | `QECCode["Concatenate"\|"RemoveQubit"\|"Paste", …]` |
| `QECPauliVector/String/Q/Weight/CommuteQ/Product/Phase` | merge into one object | `QECPauli` (+ `QECPauliQ` kept as the guard) |
| `QECStimCircuit` | rename + property | `QECStim`, `dem["StimString"]` |
| `QECCodeCatalog` | merge onto head | `QECCode["Catalog"]` |
| `QECClassicalHammingMatrix` | internalize | `PackageScope` |
| `QECClearCache`, 3 `$…Limit` globals | keep | (unchanged) |
| decoder | **add a first-class object** | `QECDecoder` (consumes `QECDetectorModel`) |

### Order of work

1. **Trivial, behavior-preserving, do first:** de-duplicate `QECCode["Properties"]`;
   replace the `Quiet` with a `QECPauliQ` guard; internalize `QECClassicalHammingMatrix`.
2. **The gadget trait and the `QuantumCircuitOperator`/`Diagram` bridge.** Highest
   day-one value; isolated to the six gadget heads plus one converter.
3. **The operational objects on `QECCode`** (`"Encoder"`, `"Codespace"`,
   `"SyndromeMeasurement"`, `"Recovery"`, `"LogicalChannel"`) and noise-as-channel. This is
   Tier B's core and where the design earns its keep.
4. **De-polymorphize `QECLogicalErrorRate`** (move the report onto `QECDetectorModel`), and
   add the `QECDecoder` object.
5. **Fold constructions and the catalog onto `QECCode`; collapse the Pauli verbs into
   `QECPauli`; rename the two heads.** Deprecate old symbols with one-release aliases.
6. **Read `QECTransversalGate`'s matrices from `QuantumOperator`**, and, on the QF side,
   propose a phase-carrying Clifford-conjugation primitive so the dense-matrix
   `transversalAction` can retire (the one genuine QF capability gap the layer works
   around).

Steps 1-2 are safe before the migration into `Kernel/QEC/`. Steps 3-6 change names, return
types, and the object model, so they are the ones to settle before reference pages freeze
them.

---

## 6. References

Operational and categorical formulation (what the minimal objects are):

- E. Knill, R. Laflamme, L. Viola, *Theory of quantum error correction for general noise*,
  Phys. Rev. Lett. 84, 2525 (2000), [doi:10.1103/PhysRevLett.84.2525](https://doi.org/10.1103/PhysRevLett.84.2525).
  The encode/noise/recover composition and the Knill-Laflamme condition.
- B. Rahn, A. C. Doherty, H. Mabuchi, *Exact performance of concatenated quantum codes*,
  Phys. Rev. A 66, 032304 (2002), [doi:10.1103/PhysRevA.66.032304](https://doi.org/10.1103/PhysRevA.66.032304).
  The logical channel, and concatenation as its composition.
- D. Kribs, R. Laflamme, D. Poulin, *Unified and generalized approach to quantum error
  correction*, Phys. Rev. Lett. 94, 180501 (2005), [doi:10.1103/PhysRevLett.94.180501](https://doi.org/10.1103/PhysRevLett.94.180501);
  and *Operator quantum error correction*, Quantum Inf. Comput. 6, 382 (2006),
  [doi:10.1017/CBO9781139034807.008](https://doi.org/10.1017/CBO9781139034807.008).
- D. Poulin, *Stabilizer formalism for operator quantum error correction*, Phys. Rev. Lett.
  95, 230504 (2005), [doi:10.1103/PhysRevLett.95.230504](https://doi.org/10.1103/PhysRevLett.95.230504);
  G. Dauphinais, D. Kribs, M. Vasmer, *Stabilizer formalism for operator algebra quantum
  error correction*, Quantum 8, 1261 (2024), [doi:10.22331/q-2024-02-21-1261](https://doi.org/10.22331/q-2024-02-21-1261).
  The operator-algebra frame for subsystem and bosonic forward-compatibility.
- A. Cowtan, S. Burton, *CSS code surgery as a universal construction*, Quantum 8, 1344
  (2024), [doi:10.22331/q-2024-05-14-1344](https://doi.org/10.22331/q-2024-05-14-1344).
  Codes as objects, constructions as morphisms.
- R. Cleve, D. Gottesman, *Efficient computations of encodings for quantum error
  correction*, Phys. Rev. A 56, 76 (1997), [doi:10.1103/PhysRevA.56.76](https://doi.org/10.1103/PhysRevA.56.76).
  The encoder as an efficiently computable isometry.

Software design and the reference packages (how a QEC framework is layered and named):

- C. Gidney, *Stim: a fast stabilizer circuit simulator*, Quantum 5, 497 (2021),
  [doi:10.22331/q-2021-07-06-497](https://doi.org/10.22331/q-2021-07-06-497),
  [arXiv:2103.02202](https://arxiv.org/abs/2103.02202). `Circuit`, `DetectorErrorModel`,
  `PauliString`, `Tableau`.
- P.-J. H. S. Derks et al. (incl. J. Eisert), *Designing fault-tolerant circuits using
  detector error models*, Quantum 9, 1905 (2025), [doi:10.22331/q-2025-11-06-1905](https://doi.org/10.22331/q-2025-11-06-1905).
  The DEM as the standardized circuit-to-decoder interface, and the three-level circuit
  abstraction.
- M. Kang et al. (incl. K. R. Brown), *QUITS: a modular qLDPC code circuit simulator*,
  Quantum 9, 1931 (2025), [doi:10.22331/q-2025-12-05-1931](https://doi.org/10.22331/q-2025-12-05-1931).
  The decoupled code / circuit / noise / DEM / decoder / rate pipeline.
- O. Higgott, *PyMatching*, ACM Trans. Quantum Comput. 3, 1 (2022),
  [doi:10.1145/3505637](https://doi.org/10.1145/3505637). The decoder that consumes a DEM.
- N. Rengaswamy et al. (incl. H. D. Pfister), *Logical Clifford synthesis for stabilizer
  codes*, IEEE Trans. Quantum Eng. 1, 1 (2020), [doi:10.1109/TQE.2020.3023419](https://doi.org/10.1109/TQE.2020.3023419).
  The transversal/logical-gate layer as synthesis on the code.
- QuantumClifford.jl and its `QECCore.jl` interface package (QuantumSavory):
  `parity_checks`, `code_n`/`code_k`/`code_s`, `distance`, `logx_ops`/`logz_ops`,
  `naive_syndrome_circuit`, `PauliFrame`, `PauliError`, `evaluate_decoder`. The one
  ecosystem that has already factored the code interface into a separate package, which is
  the direction this layer's `Kernel/QEC/` migration points.

Grounding for the code families and the honest scope:

- S. Bravyi et al., bivariate-bicycle qLDPC codes, Nature 627, 778 (2024),
  [doi:10.1038/s41586-024-07107-7](https://doi.org/10.1038/s41586-024-07107-7),
  [arXiv:2308.07915](https://arxiv.org/abs/2308.07915); D. Gottesman, *Stabilizer Codes and
  Quantum Error Correction*, [arXiv:quant-ph/9705052](https://arxiv.org/abs/quant-ph/9705052)
  and the 2026 book the layer is built against; Dennis, Kitaev, Landahl, Preskill,
  *Topological quantum memory*, [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143).

---

*The layer is a strong implementation whose internal choices are right and whose module
layering already matches the best of the field. What it has not done is present itself in
the object language that quantum theory uses for error correction and that
QuantumFramework is built on. A code is an encoding isometry, its noise is a channel, its
syndrome extraction is a measurement, and its correction cycle is a channel: the framework
has all four types already, and the stabilizer and detector-error-model engines are the
efficient way to compute them. The redesign trims the doubled surface to its objects, and
then makes those objects the operational ones, so the QEC cycle composes with the rest of
QuantumFramework instead of standing beside it in a private vocabulary.*
