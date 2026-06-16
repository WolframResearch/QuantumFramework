# The QuantumFramework Stabilizer Subsystem vs the Dedicated Clifford Packages

**What this is.** An evidence-backed comparison of the Wolfram QuantumFramework (QF) stabilizer
subsystem, `Kernel/Stabilizer/` (paclet 2.0.0), against the four packages a practitioner would
otherwise reach for: **Stim**, **QuantumClifford.jl**, **Qiskit** (`Clifford` / `StabilizerState`),
and **Cirq** (`CliffordSimulator`). The emphasis is on *what QF offers that the others do not*, with
the speed trade honestly stated and measured rather than recalled.

**This document is the deep, benchmark-backed companion** to the tutorial
[`OngoingProjects/Stabilizer/Stabilizer-Formalism-By-Computation.md`](../Stabilizer/Stabilizer-Formalism-By-Computation.md)
(which has a short "honest comparison" section). It does not duplicate that section; it expands it
with the full feature matrix, a re-measured cross-package benchmark, per-advantage verified code, and
the complete gap list.

**Provenance.** Every QF claim is read from the live kernel under
`QuantumFramework/Kernel/Stabilizer/` and cited as `file:line`; every benchmark and code snippet in
this document was run and its real output captured. Anchor: repo HEAD `f9dc1cdc` (2026-06-14), which
the QF audit set is synced to. Machine: Apple M2 Pro (Darwin 23.2.0), Wolfram Language 15.0.0, Stim
1.16.0, QuantumClifford.jl 0.11.5 (Julia 1.12), Qiskit 2.4.1, Cirq 1.6.1.

---

## 1. Feature matrix

Legend: **Yes** native; **partial** present but limited or indirect; **No** absent; **bridge** reached
by delegating to an external tool through QF's Python interop. QF cells are verified live against the
2.0.0 working tree (citations in the notes); competitor cells are from the package audits
([`external-packages-audit.md`](../Stabilizer/external-packages-audit.md)) and their documentation.

| Capability | QF `PauliStabilizer` | Stim | QuantumClifford.jl | Qiskit `Clifford`/`StabilizerState` | Cirq `CliffordSimulator` |
|---|:--:|:--:|:--:|:--:|:--:|
| Tableau (Aaronson-Gottesman) simulation | Yes | Yes | Yes | Yes | Yes |
| Destabilizer-augmented tableau (rank tracking) | Yes (always) | Yes (inverse tableau) | Yes (`MixedDestabilizer`) | Yes | Yes |
| Compiled bulk gate fold | Yes (`ApplyCircuit`, compiled-to-C) | Yes (C++) | Yes (Julia) | partial (Python loop) | partial (Python loop) |
| Single-qubit $Z$ measurement | Yes | Yes | Yes | Yes | Yes |
| Arbitrary Pauli-string measurement | Yes (`ps["M","XZZXI"]`) | Yes | Yes | partial | partial |
| Pauli expectation $\langle P\rangle$ | Yes (`ps["Expectation"]`) | Yes | Yes | Yes (`expectation_value`) | partial |
| Bipartite entanglement entropy | Yes (closed form, any $n$) | No (build from tableau) | Yes | No | No |
| Random Clifford sampling | Yes (**Bravyi-Maslov / Mallows**) | Yes (symplectic-uniform) | Yes (symplectic-uniform) | Yes (`random_clifford`) | partial |
| Stabilizer-rank / magic (frames) | Yes (`StabilizerFrame`) | No | partial (`PauliFrame` is for sampling, not rank) | No | No |
| Clifford channel as first-class object | Yes (`CliffordChannel`) | partial (noise channels in circuits) | partial | No | No |
| Graph states + local complementation | Yes (`GraphState`, `LocalComplement`) | partial (Crumble/glue) | Yes (`graphstate`, `Graphs.jl`) | No | No |
| Symbolic / formal-symbol outcomes | **Yes** (`SymbolicMeasure`, SymPhase) | No | No | No | No |
| Exact / rational arithmetic | **Yes** | No (bit / float) | No (bit / float) | No | No |
| Detector sampling | No | **Yes** | Yes (`PauliFrame`) | No | No |
| Detector error model (DEM) export | No | **Yes** | partial | No | No |
| Decoder / PyMatching bridge | No | bridge (sinter/pymatching) | Yes (`QECCore`, decoder exts) | No | No |
| Batched bit-packed Pauli-frame shot sampling | No | **Yes** (`FrameSimulator`) | **Yes** (`pftrajectories`, threaded) | No | No |
| QEC code-family generators | partial (named codes only) | partial (glue) | **Yes** (`QECCore`, ~100 files) | No | No |
| Native qudit ($d>2$) tableau | partial (carries a `d` parameter) | No | No | No | No |
| Conversion into a full object model | **Yes** (State / Operator / Circuit / Channel) | No | partial (to/from Qiskit, Cirq) | within Qiskit | within Cirq |

**Notes and citations (QF cells).**
- Tableau, destabilizer-augmented, $\{2,n,2n\}$ canonical array: `PauliStabilizer.m:24,36`; closed-form
  register constructor `Constructors.m:259`.
- Compiled bulk fold: `Compiled.m:95` (`compiledCliffordFold`, `CompilationTarget :> "C"`/`"WVM"`),
  entry `Compiled.m:224` (`ps["ApplyCircuit"]`) and `PauliStabilizer.m:123` (`PauliStabilizerApply`).
- Pauli-string measurement: `PauliMeasure.m:96` (packed) / `:99` (canonical); expectation
  `InnerProduct.m:156`.
- Entropy (Fattal closed form, $\mathrm{rank}_{\mathbb F_2}-|A|$, $O(n^3)$, no OOM): `Entropy.m:32`.
- Random Clifford (Bravyi-Maslov / Koenig-Smolin Mallows): `RandomClifford.m:53`, sampler
  `RandomClifford.m:14`.
- `StabilizerFrame` (Garcia-Markov Quipu): `StabilizerFrame.m:57`; non-Clifford boundary doubles the
  frame, `GateUpdates.m:182`, `StabilizerFrame.m:109`.
- `CliffordChannel` (Yashin25 Boolean Choi tableau, composition via $\mathbb F_2$ null space):
  `CliffordChannel.m:46,252`.
- `GraphState` / `LocalComplement` (Anders-Briegel): `GraphState.m:48,145`.
- `SymbolicMeasure` / SymPhase (Fang-Ying23): `SymbolicMeasure.m:74`; symbolic-sign tableaux are a
  first-class storage shape, `PauliStabilizer.m:27`.
- Object-model conversion (UpValues so a stabilizer drops into the rest of QF): `Conversions.m:142-149`,
  `HybridInterop.m:187,257`.

The matrix's shape is the thesis: QF is **wider** than the dedicated engines (symbolic, exact, frames,
channels, graph states, full object-model interop), and the dedicated engines are **deeper** on the one
thing they exist for (bulk fault-tolerance Monte Carlo: detectors, error models, decoders,
constant-cost frame sampling).

---

## 2. Performance: re-measured at HEAD, apples to apples

**Task.** An identical deterministic stream of single- and two-qubit Cliffords
(a fixed $H$/$S$/CNOT pattern; $2000$/$10000$/$20000$ gates at $n=100$/$500$/$1000$) starting from
$|0\rangle^{\otimes n}$. The **same op list** drives all three tableau engines: a shared JSON file
read by Stim and QuantumClifford.jl, and the identical pattern in WL. Total wall-clock including the
constructor, minimum of repetitions, cache cleared per repetition. (Qiskit and Cirq are omitted from the
timing table: their stabilizer paths are Python-loop simulators not built for this throughput regime,
and the audits place them with QF's interpreted per-gate route, not the compiled route.)

Scripts: `scripts/qf_final_all.wl`, `scripts/bench_stim_all.py`, `scripts/bench_qc_full.jl`. Streams:
`scripts/stab_ops_{100,500,1000}.json`.

| $n$ / gates | QF `ps["ApplyCircuit"]` | Stim `TableauSimulator` | QuantumClifford.jl `MixedDestabilizer` | QF / Stim | QF / QC.jl |
|---|---:|---:|---:|---:|---:|
| 100 / 2000 | $2.46\,\mathrm{ms}$ | $1.01\,\mathrm{ms}$ | $1.10\,\mathrm{ms}$ | $2.4\times$ | $2.2\times$ |
| 500 / 10000 | $28.2\,\mathrm{ms}$ | $5.31\,\mathrm{ms}$ | $24.4\,\mathrm{ms}$ | $5.3\times$ | $1.2\times$ |
| 1000 / 20000 | $67.7\,\mathrm{ms}$ | $10.9\,\mathrm{ms}$ | $112.9\,\mathrm{ms}$ | $6.2\times$ | $\mathbf{0.6\times}$ |

**Reading it.**
- QF is within a small constant factor of Stim, $2.4\times$ to $6.2\times$ on this run (the
  optimization report's $6\mathrm{f}953\mathrm{ad}5$ measurement was $2.5$ to $5.3\times$; the spread is
  run-to-run noise plus the fact that $n=1000$ rides the high end). **Do not quote the stale "660-3000x"
  numbers**: those predate the June 2026 rework.
- QF **beats QuantumClifford.jl at $n=1000$** ($67.7$ vs $112.9\,\mathrm{ms}$, $1.7\times$ faster) and is
  in the same regime at $n=500$. QC.jl is faster at $n=100$.
- The route is the whole story. The same stream applied **one gate at a time** through the property
  interface (`ps["H",1]["CNOT",1,2]...`) is far slower because each gate pays object wrap/unwrap and
  interpreter dispatch:

| $n$ / gates | QF per-gate fold | QF `ApplyCircuit` | speedup of the bulk route |
|---|---:|---:|---:|
| 100 / 2000 | $73.3\,\mathrm{ms}$ | $2.46\,\mathrm{ms}$ | $30\times$ |
| 500 / 10000 | $428.0\,\mathrm{ms}$ | $28.2\,\mathrm{ms}$ | $15\times$ |
| 1000 / 20000 | $1012.3\,\mathrm{ms}$ | $67.7\,\mathrm{ms}$ | $15\times$ |

  The compiled kernel itself runs at $\sim 0.22$ to $0.27\,\mu\mathrm{s}$/gate, flat in $n$ (Stim's own
  regime); the residual wall-clock is the $\sim 0.7\,\mu\mathrm{s}$/gate encode pass plus the
  pack/unpack boundary (`Compiled.m:180,192`).

**Hot-path operations** (this run, HEAD `f9dc1cdc`; consistent with the bottleneck audit's HEAD
re-verification):

| Operation | $n=100$ | $n=500$ | $n=1000$ | where |
|---|---:|---:|---:|---|
| Register constructor `PauliStabilizer[n]` | $0.13\,\mathrm{ms}$ | $1.46\,\mathrm{ms}$ | $5.35\,\mathrm{ms}$ | `Constructors.m:259` |
| Single-qubit measure `ps["M",1]` | $0.12\,\mathrm{ms}$ | $0.19\,\mathrm{ms}$ | | `Measurement.m:94` |
| Pauli-string measure `ps["M","Z..Z"]` | $4.20\,\mathrm{ms}$ | $51.7\,\mathrm{ms}$ | | `PauliMeasure.m:44` |
| Pauli-form display `ps["Stabilizers"]` | $2.0\,\mathrm{ms}$ | $31.1\,\mathrm{ms}$ | | `Formatting.m:33` |
| Composition `a[b]` | $54.2\,\mathrm{ms}$ | | | `Compose.m:28` |

The remaining gap to Stim is **SIMD lane width** (Stim packs 256-bit vector lanes; QF packs 62-bit
machine words, `Packed.m:53`, `Compiled.m:44`) plus the per-gate encode pass, not an algorithmic deficit.
Composition and AG circuit synthesis remain $O(n^3)$ (`Compose.m:24`, `Conversions.m:91`): correct and
fast enough at these sizes, but not the place to push $n$.

---

## 3. The advantages, made concrete

Every snippet below was run; the output shown is the captured result. Load the working tree first:

```wolfram
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

### 3.1 Symbolic outcomes: SymPhase (`SymbolicMeasure.m:74`)

A measurement outcome can be carried as a *formal symbol* propagated through the tableau, so one circuit
traversal serves arbitrarily many samples and outcome correlations are explicit polynomials. No
floating-point simulator can do this.

```wolfram
bell = PauliStabilizer[{"XX", "ZZ"}];          (* Bell stabilizer group *)
sm   = bell["SymbolicMeasure", 1];             (* measure qubit 1 symbolically *)
sm["Signs"]
sm["SampleOutcomes", 4][[All, "Signs"]]
```
```
{1, 1, 1 - 2 s[1], 1}                            (* the fresh outcome symbol s[1] lives in the sign *)
{{1, 1, -1, 1}, {1, 1, -1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}}   (* sampled by substituting s[1] -> 0/1 *)
```

**Who can do this:** QF only. Stim, QuantumClifford.jl, Qiskit, and Cirq are numeric; an outcome is a
sampled bit, never a symbol you can later substitute or correlate algebraically.

### 3.2 Exact rational arithmetic and bit-reproducible seeded output (`RandomClifford.m:14`)

Coefficients, signs, and the Mallows random-Clifford sampler are exact; the sampler exactifies its
Mallows ratio (`RandomClifford.m:22`, `SetPrecision[..., Infinity]`) so a seeded draw is bit-identical
across runs and there is no floating-point drift.

```wolfram
SeedRandom[42]; a = PauliStabilizer["Random"[6]];
SeedRandom[42]; b = PauliStabilizer["Random"[6]];
a["Tableau"] === b["Tableau"]                   (* True *)
MatchQ[a["Signs"], {(-1 | 1) ..}]               (* True: exact integers, no float sign *)
Head[PauliStabilizer[{"XXX","ZZI","IZZ"}]["Expectation","ZZZ"]]   (* Integer, not Real *)
```
```
True
True
Integer
```

**Who can do this:** QF only at this level. Both reference engines are bit/float; seeded reproducibility
exists but the arithmetic is not exact rational, and an expectation is a machine float. (QF's random
Clifford is also drawn from the mathematically uniform **Mallows** distribution, where Stim and QC.jl
use symplectic-uniform sampling: QF is *ahead* of both on the sampling distribution itself.)

### 3.3 Closed-form entanglement entropy at any $n$ (`Entropy.m:32`)

The Fattal et al. formula $S(A)=\mathrm{rank}_{\mathbb F_2}(\text{generators restricted to } A)-|A|$ is
$O(n^3)$ and independent of $|A|$, so it never materializes a state vector and never runs out of memory.

```wolfram
n = 200;
specs = Join[{"H" -> 1}, Table["CNOT" -> {i, i + 1}, {i, 1, n - 1}]];
ghz = PauliStabilizer[n]["ApplyCircuit", specs];    (* 200-qubit GHZ *)
ghz["Entropy", Range[100]]                            (* half-chain cut *)
```
```
1                                                     (* GHZ has 1 bit of entanglement across any cut; 22.7 ms *)
```

**Who can do this:** QF natively. QuantumClifford.jl can compute stabilizer entanglement too;
Stim, Qiskit, and Cirq leave you to extract it from the tableau yourself.

### 3.4 Stabilizer rank and magic via `StabilizerFrame` (`StabilizerFrame.m:57`)

A non-Clifford gate ($T$, $P[\theta]$) maps a stabilizer state to a superposition of stabilizer states.
The frame represents that superposition exactly; its term count is the stabilizer rank, the resource that
controls magic-state and near-Clifford simulation cost.

```wolfram
plus = PauliStabilizer[1]["H", 1];   (* |+> *)
f    = plus["T", 1];                  (* T|+> as a stabilizer frame *)
f["Length"]                           (* stabilizer rank *)
f["Coefficients"]
```
```
2
{(1 + E^(I Pi/4))/2, (1 - E^(I Pi/4))/2}             (* exact closed-form coefficients *)
```

**Who can do this:** QF only. Stim and Cirq are pure-Clifford; Qiskit's `Clifford` is pure-Clifford;
QuantumClifford.jl's `PauliFrame` is a *sampling* construct, not a stabilizer-rank decomposition object.

### 3.5 `CliffordChannel` as a first-class algebraic object (`CliffordChannel.m:46`)

A Clifford channel (CPTP map sending stabilizer states to stabilizer states) is a Yashin25 Boolean Choi
tableau that composes via $\mathbb F_2$ null-space with AG phase tracking, not a circuit you re-simulate.

```wolfram
cc = CliffordChannel["Identity", 2];
{cc["Rank"], cc["InputQubits"], cc["OutputQubits"]}    (* {4, 2, 2} *)
cc[cc]["Source"]                                        (* composition is itself a CliffordChannel *)
CliffordChannel[PauliStabilizer[{"XX","ZZ"}]]["Rank"]  (* Bell state as a state-prep Choi *)
```
```
{4, 2, 2}
Composition
2
```

**Who can do this:** QF only as a standalone algebraic head. The others model noise as channels *inside a
circuit* to be sampled, not as a composable Choi-tableau object.

### 3.6 Graph states and local complementation (`GraphState.m:48,145`)

```wolfram
gs = GraphState[CycleGraph[5]];
gs["Stabilizers"]                       (* K_i = X_i prod_{j in N(i)} Z_j *)
LocalComplement[gs, 1]["EdgeCount"]     (* toggles edges among vertex 1's neighbors *)
```
```
{XZIIZ, ZXZII, IZXZI, IIZXZ, ZIIZX}
6                                       (* was 5 *)
```

**Who can do this:** QF and QuantumClifford.jl (via `Graphs.jl`); not Stim's core, Qiskit, or Cirq.
(Caveat, stated in §4: QF's `LocalComplement` toggles the graph but does **not** yet track the
vertex-operator Clifford correction, `GraphState.m:162`.)

### 3.7 Conversion into the full QF object model (`Conversions.m:142-149`)

A `PauliStabilizer` becomes a `QuantumState`, `QuantumOperator`, `QuantumCircuitOperator`, or (via
`CliffordChannel`) a channel, so a stabilizer subroutine drops into a larger simulation that also uses
dense states, density matrices, Lindblad evolution, or phase-space methods, without leaving the language.

```wolfram
ghz = PauliStabilizer[{"XXX","ZZI","IZZ"}];
Head[QuantumState[ghz]]                  (* QuantumState *)
Head[ghz["QuantumOperator"]]             (* QuantumOperator *)
Head[ghz["QuantumCircuit"]]              (* QuantumCircuitOperator *)
Normal @ QuantumState[ghz]["StateVector"]
```
```
QuantumState
QuantumOperator
QuantumCircuitOperator
{1/Sqrt[2], 0, 0, 0, 0, 0, 0, 1/Sqrt[2]}   (* exact amplitudes *)
```

**Who can do this:** QF natively, into a single unified object model spanning all of quantum theory.
Stim, QC.jl, Qiskit, and Cirq can interconvert with each other's circuit objects but have no comparable
one-model surface (pure states, density matrices, POVMs, channels, Fock states, Wigner functions, all the
same kind of expression).

### 3.8 Pauli-string measurement and expectation (`PauliMeasure.m`, `InnerProduct.m:156`)

```wolfram
ghz = PauliStabilizer[{"XXX","ZZI","IZZ"}];
{ghz["Expectation","ZZZ"], ghz["Expectation","XXX"], ghz["Expectation","ZII"]}
Keys @ ghz["M","XXX"]                    (* deterministic outcome *)
```
```
{0, 1, 0}                                 (* <ZZZ>=0 (ZZZ not in the group), <XXX>=1, <ZII>=0 *)
{0}                                       (* XXX is a stabilizer: deterministic +1 outcome *)
```

This is the syndrome-extraction primitive; QF has it, and so do Stim and QC.jl. Qiskit and Cirq support it
only partially through their stabilizer-state APIs.

---

## 4. Honest gaps (not soft-pedaled)

These are the places where the dedicated engines win, stated plainly. Several are the entire reason Stim
and QuantumClifford.jl exist.

- **No detector / observable concept.** QF has no notion of a "detector" (a parity of measurement
  outcomes) or a tracked logical observable. Verified absent: no such symbol anywhere in
  `Kernel/Stabilizer/`.
- **No detector error model (DEM) export.** Stim's `DetectorErrorModel` (the compiled error graph that
  decoders consume) has no analogue in QF. `external-packages-audit.md:57` records this gap.
- **No decoder and no PyMatching bridge.** QF does not decode syndromes and ships no minimum-weight
  perfect-matching or BP-OSD binding. QuantumClifford.jl carries `QECCore` plus decoder extensions;
  Stim pairs with sinter/pymatching. QF has neither.
- **No batched bit-packed Pauli-frame sampler with fresh per-shot noise.** This is *the* throughput gap.
  Stim's `FrameSimulator` and QC.jl's `pftrajectories` derive every additional shot at constant cost per
  gate by propagating a Pauli frame diffed against one reference sample (`external-packages-audit.md:354`).
  QF measurement returns a conditional `Association` (`Measurement.m:139`); the user samples by hand, and
  there is no constant-cost frame engine. **Consequence: large-scale QEC Monte Carlo, threshold
  estimation, and $\Lambda$ (logical-error-suppression) campaigns at distance $\gtrsim 7$ are out of
  reach in QF.** `StabilizerFrame` represents *mixtures* and supports `"SampleOutcomes"`, but that is not
  the same primitive.
- **SIMD lane-width gap to Stim.** QF packs 62 tableau bits per machine word (`Packed.m:53`,
  `Compiled.m:44`) to stay inside Int64; Stim uses 256-bit vector lanes. This is most of the residual
  $2.4$ to $6.2\times$ factor in §2 and is a constant, not an algorithmic, deficit.
- **`QuantumCircuitOperator` construction overhead for $10^4$-gate circuits.** Building a circuit *object*
  out of $10^4$ individual gates is slow in the host framework (minutes), independent of the stabilizer
  kernel. Feed `ps["ApplyCircuit", specs]` (or `PauliStabilizerApply`) a gate-spec list directly.
- **The `PauliStabilizer[circuit]` trap is an $O(4^n)$-per-operator path.** `PauliStabilizer[qco]`
  (`Constructors.m:189`) folds an AG *tomographic* decomposition over each operator; it is not the gate
  fold. Measured contrast at $n=12$:

  ```wolfram
  n = 12; specs = Join[{"H" -> 1}, Table["CNOT" -> {i, i+1}, {i, 1, n-1}]];
  First @ AbsoluteTiming[ PauliStabilizer[QuantumCircuitOperator[specs]] ]      (* 6.73 s *)
  First @ AbsoluteTiming[ PauliStabilizer[n]["ApplyCircuit", specs] ]           (* 0.0006 s *)
  ```
  ```
  PauliStabilizer[qc]  n=12: 6732.9 ms
  ps[ApplyCircuit]     n=12: 0.616 ms        (* same state; >10^4x faster *)
  ```
  Use `ps["ApplyCircuit", specs]` or the gate-update chain, never `PauliStabilizer[circuit]`, for
  simulation.
- **Global-phase contract on reconstructed state vectors.** Clifford gate updates deliberately drop the
  global phase (`GateUpdates.m:11-21`): recovering it per gate would cost $O(2^n)$. So a gate-updated
  tableau's `["State"]` equals the true state *up to a global phase*, while the *constructor* path
  (`PauliStabilizer[qs]`/`[qo]`) stores a `"GlobalPhase"` and round-trips exactly:

  ```wolfram
  qs = QuantumState["0"];
  Normal @ PauliStabilizer[qs]["Y", 1]["State"]["StateVector"]                  (* {0, 1} *)
  Normal @ PauliStabilizer[QuantumOperator["Y"][qs]]["State"]["StateVector"]    (* {0, I} *)
  ```
  Both are $Y|0\rangle$ up to the global phase $i$; compare magnitudes or fidelity, not raw amplitudes.
- **Storage-form caveat for equality.** Single gates return a *packed* object, `ApplyCircuit` and
  constructors return a *canonical* one (`Packed.m:80`, `Compiled.m:203`); both are valid everywhere but
  `SameQ` across the two forms fails. Compare `ps1["Tableau"] === ps2["Tableau"]` or `["Stabilizers"]`.
- **`GraphState` and `CliffordChannel` are partial.** `GraphState[ps]` accepts only graph-form
  stabilizers (`GraphState.m:56`, `::nongraph`); `LocalComplement` does not yet track vertex-operator
  Cliffords (`GraphState.m:162`); `CliffordChannel[qc]` handles only deterministic single-Pauli channels
  (`CliffordChannel.m:402`, `::stochastic`).
- **Compiled-kernel packing trap.** A `CompiledFunction` silently falls back to interpreted evaluation
  (~$100\times$ slower, no warning) on unpacked arguments; the kernel forces `Developer`ToPackedArray` at
  every boundary (`Compiled.m:63,192`). Relevant only if you extend the kernel.

---

## 5. Verdict: complementary, not competing

**Use QuantumFramework when** the stabilizer work is *part of a larger physics computation* or needs
something the numeric engines structurally lack:
- exact or symbolic stabilizer analysis (carry free symbols through a circuit; keep signs and
  coefficients rational; bit-reproducible seeded sampling);
- magic / stabilizer-rank accounting via `StabilizerFrame`;
- codes, graph-state surgery, and Clifford-channel algebra as first-class objects;
- a stabilizer subroutine that must hand its state to dense simulation, Lindblad evolution, phase-space
  methods, tomography, or the rest of the Wolfram Language without leaving the object model.

  On raw Clifford throughput QF is now a *real option*, not a toy: $2.4$ to $6.2\times$ Stim and ahead of
  QuantumClifford.jl at $n=1000$, cross-validated against both (`Tests/Stabilizer/CrossPackage_Stim.wlt`,
  `CrossPackage_QuantumClifford.wlt`).

**Use Stim or QuantumClifford.jl when** the task is *high-volume fault-tolerance Monte Carlo and
decoding*:
- distance-$\gtrsim 7$ surface/color-code shot generation (constant-cost bulk Pauli-frame sampling);
- detector error models feeding a decoder;
- threshold and $\Lambda$ campaigns, decoder benchmarking (Stim + sinter/pymatching; QC.jl's `QECCore`
  and decoder extensions).

These are not the same job. The dedicated engines are a few times faster on pure gate throughput and own
the QEC-sampling-and-decoding pipeline outright; QF buys symbolic, exact, frame, channel, graph-state, and
whole-framework interoperability for a small constant in speed. For a research group that already lives in
Wolfram Language, QF removes the round-trip to an external tool for everything *except* large-scale QEC
sampling, where the right move is to call Stim or QuantumClifford.jl and bring the results back.

---

## Reproduce

```
# Shared op streams (already in scripts/; stage to /tmp):
cp "OngoingProjects/Platform Comparison/scripts/"stab_ops_*.json /tmp/

# QF (repo working tree):     wolframscript -file "OngoingProjects/Platform Comparison/scripts/qf_final_all.wl"
# Stim (venv):                /tmp/qcbench/bin/python "OngoingProjects/Platform Comparison/scripts/bench_stim_all.py"
# QuantumClifford.jl:         julia "OngoingProjects/Platform Comparison/scripts/bench_qc_full.jl"
# Advantage + trap snippets:  wolframscript -file "OngoingProjects/Platform Comparison/scripts/qf_stab_verify.wls"
```

Underlying detail: the stabilizer kernel lives in
[`QuantumFramework/Kernel/Stabilizer/`](../../QuantumFramework/Kernel/Stabilizer/) (20 files); the
performance story is in [`QF-Stabilizer-Optimization-Report.md`](QF-Stabilizer-Optimization-Report.md)
and [`QF-Stabilizer-Bottleneck-Audit.md`](QF-Stabilizer-Bottleneck-Audit.md); the cross-package design
analysis is in [`external-packages-audit.md`](../Stabilizer/external-packages-audit.md); the broader
platform context is [`QF-Master-Platform-Comparison.md`](QF-Master-Platform-Comparison.md).

---

## Physicist's summary

The object here is the stabilizer state: instead of $2^n$ amplitudes, you track the $n$-generator Pauli
group that fixes the state, so Clifford evolution is a constant-size rewrite of a binary symplectic
tableau and the whole simulation stays polynomial. The question this document answers is how Wolfram's
implementation of that idea stacks up against the two engines built for it, Stim (C++) and
QuantumClifford.jl (Julia), and against the stabilizer corners of Qiskit and Cirq.

The finding is a clean division of labor. On the one task the dedicated engines optimize, folding a long
Clifford circuit, QF is now within a small constant of Stim (a factor of a few, measured here at $2.4$ to
$6.2\times$ on an identical gate stream) and actually faster than QuantumClifford.jl at a thousand qubits;
the residual gap is hardware vector width (256-bit lanes versus 62-bit machine words), not the algorithm.
What QF adds that none of them have is everything *around* the tableau: measurement outcomes carried as
formal symbols and substituted later (SymPhase), exact rational arithmetic with bit-reproducible seeded
randomness, closed-form entanglement entropy at any size, a stabilizer-frame object that exposes the
stabilizer rank (the cost knob for magic), Clifford channels as composable Choi tableaux, graph states
with local complementation, and, decisively, the ability to hand a stabilizer state straight to a dense
simulator, a Lindblad solver, or a Wigner transform inside one object model. What it lacks is the
fault-tolerance production line: detectors, detector error models, decoders, and the constant-cost
bulk Pauli-frame sampler that makes Stim and QC.jl indispensable for surface-code threshold studies. So
the honest recommendation is complementary use: QF for exact, symbolic, embedded stabilizer analysis;
Stim or QuantumClifford.jl for high-volume QEC Monte Carlo and decoding.

On the code: the subsystem is twenty Wolfram-Language files under
`QuantumFramework/Kernel/Stabilizer/`. Speed comes from two representations, a bit-packed machine-word
tableau (`Packed.m`) and a single compile-to-C gate-fold kernel reached through `ps["ApplyCircuit", ...]`
(`Compiled.m`); correctness is gated by exact-equivalence tests against the canonical array path and by
cross-package suites that replay the same gates through Stim and QuantumClifford.jl. The distinctive heads
(`StabilizerFrame`, `CliffordChannel`, `GraphState`) and the symbolic-measurement path are separate files,
each a short, self-contained extension of the same tableau. Two traps are worth carrying away:
`PauliStabilizer[circuit]` is an $O(4^n)$ tomography path, not the simulator (use `"ApplyCircuit"`), and a
gate-updated tableau reproduces the state only up to a global phase.
