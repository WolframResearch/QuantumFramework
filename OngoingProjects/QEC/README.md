# Quantum error correction

A code layer on top of the framework's stabilizer engine. The engine
(`QuantumFramework/Kernel/Stabilizer/`) holds a *state*: an n-qubit stabilizer state
pushed through Clifford gates in the Heisenberg picture. This layer holds a *code*: a
subspace fixed by a commuting group of Pauli checks, together with everything you would
want to ask about one — its logical operators and distance, its syndromes and decoders,
a circuit that prepares a codeword, a circuit that measures its checks, and whether a
logical qubit survives a given noise model.

Priorities are set by [`QEC-Development-Plan.md`](QEC-Development-Plan.md). Items 1 and 2
of that plan are done; 3 is partly done; 4 and 5 are not started.

**Start here:** [`docs/Tutorials/StabilizerCodes.nb`](docs/Tutorials/StabilizerCodes.nb) —
a tech note in 19 sections with 61 runnable examples, which is the guided tour of
everything below. Its source is the markdown beside it; the notebook is generated.

## What is here

```
QECCore/           the layer: 16 files, ~3,200 lines, 22 public symbols
Tests/             406 tests in 10 MUnit files, plus a local runner
StimCrossCheck/    the external oracle: emit the circuits, compare against Stim
docs/              markdown sources, the build script, and the built notebook
References/        an index of the source material (the PDFs are not in the repo)
EngineMeasurementBug.md   a defect found in the engine, with a minimal reproduction
```

## Running it

From a checkout of this repository, with `$repo` its root:

```wl
PacletDirectoryLoad[FileNameJoin[{$repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

(* the layer is not in the paclet yet, so it is loaded directly *)
Get[FileNameJoin[{$repo, "OngoingProjects", "QEC", "QECCore", "QECCore.wl"}]];
Needs["Wolfram`QuantumFramework`QEC`"];
```

Then, in a later cell — the `Needs` above has to have finished before these names
resolve:

```wl
code = QECCode["SteaneCode"];
code["Parameters"]                                        (* {7, 1, 3} *)
QECLogicalErrorRate[code, QECNoiseModel["Depolarizing", p]]
QECLogicalErrorRate[code, QECNoiseModel["Circuit", 1/1000], "Rounds" -> 1]
```

The tests run in a fresh kernel and print a machine-readable summary:

```
wolframscript -file OngoingProjects/QEC/Tests/RunQECTests.wls
```

The documentation is built from its markdown sources, and each page's example cells can
be evaluated and checked against their recorded outputs first:

```
wolframscript -file OngoingProjects/QEC/docs/validate.wls
wolframscript -file OngoingProjects/QEC/docs/build.wls
```

## What it computes, and what makes it different

Everything a fast simulator does, this does too — but the noise rate may be left
**symbolic**, and then the answer is an exact polynomial rather than an estimate. A
sampler answers "about 0.0022 at p = 0.01"; this answers "10p² − 200p³/9 + 160p⁴/9 −
128p⁵/27, for every p". That is the one move Monte Carlo structurally cannot make, and it
holds at every level of noise, up to and including a circuit whose every gate, reset and
readout can fail.

Some results, all reproduced by the test suite:

| | |
|---|---|
| Exact logical error rate, `[[5,1,3]]`, depolarizing | `10p² − 200p³/9 + 160p⁴/9 − 128p⁵/27` |
| Exact rate, `[[9,1,3]]` Shor, over all 4⁹ errors | `13p² − 340p³/9 + …`, in about a second |
| Maximum likelihood vs a lookup table, Shor, symmetric noise | `13p²` against `31p²` — degeneracy |
| Maximum likelihood vs a lookup table, Steane, biased noise | `651p²/25` against `756p²/25` |
| Exact rate with noisy readout, bit-flip code, one round | `3p² − 2p³ + 6pq − 12p²q + … + q²` |
| Circuit-level rate, bit-flip code, closed form | `32p/15 + 1028p²/225 − …` |

Two facts about fault tolerance fall out of those polynomials that are awkward to state
without them, and both are in the tech note's last two sections. The duality between the
bit-flip and phase-flip codes survives at code capacity and breaks at circuit level
(`32p/15` against `112p/15`), entirely because measuring X-type checks needs eight
Hadamards that measuring Z-type checks does not. And the five-qubit code, distance 3,
stops correcting single faults when its checks are measured with bare ancillas: its rate
goes from order p² to `29p/5`, because 288 circuit faults collapse onto 71 detector
signatures of which 29 carry conflicting logical effects.

## How it is checked

Nothing is trusted because this package computed it. Every fast routine is cross-checked
against a slow, obvious one:

| What | Checked against |
|---|---|
| Pauli algebra, products, phases | Dense 2ⁿ × 2ⁿ matrices built from scratch |
| Row layout | The engine's own tableau, row for row |
| Syndromes | Expectation values on an actual `PauliStabilizer` state |
| Encoding circuits | The state they prepare, check by check |
| Structural results | Gottesman's published values (X̄ = ZIIZX, the [[4,2,2]] generators) |
| CSS and concatenation | The textbook codes they must reproduce, as groups and with signs |
| Noise probabilities | The mixture the engine's own channels return |
| Logical error rate | A closed-form binomial tail, and an independent sampled route |
| The extraction circuit | `code["Syndrome", …]`, and the circuit run as a state on the engine |
| Detectors and faults | Stim, detector by detector, over two million shots |
| Circuit-level rate | The code-capacity polynomial it must reduce to |
| The exported circuit | PyMatching, which decodes it unchanged |

The Stim comparison lives in [`StimCrossCheck/`](StimCrossCheck/README.md) and is the
sharpest of these: five cases across two noise levels, 58 detector firing rates and 5
observable flip rates computed **exactly** here and sampled there, all agreeing within 2.7
standard errors. It needs a Python install, which is why it sits outside the `.wlt` suite.

## Honest limits

- **The extraction circuit is not fault tolerant, and is not claimed to be.** Its ancillas
  are bare rather than cat-state-verified — this is Gottesman §12.1.1, whose title is
  *Non-Fault-Tolerant Measurement of Paulis*. Verified ancillas (Shor, Steane, Knill error
  correction) are roadmap item 5, and the `29p/5` above is what it costs not to have them.
- **Idle noise is not modelled.** Generators are extracted sequentially rather than in
  parallel layers, so the circuit waits more than a real one, not less. Every circuit-level
  number here is therefore optimistic in the direction the circuit is already weakest.
- **The three noise levels are not three estimates of one number.** Their `p` means
  different things, and a phenomenological threshold is not comparable to a circuit-level
  one.
- **Exactness is bounded, in two different ways.** The code-capacity route enumerates 4ⁿ
  errors (`$QECExactEnumerationLimit`); the circuit-level route costs 2^((r+1)m + 2k)
  (`$QECExactDetectorLimit`), which grows in the number of rounds. Past either, the same
  quantity is sampled.
- **Not in the paclet yet.** The files are written to move verbatim into
  `QuantumFramework/Kernel/QEC/`; until they do, there is a transitional loader, the tests
  carry a `Get`, and there are no reference pages behind `?QECCode` or F1.

## What is left

| Plan item | State |
|---|---|
| 1. Rebuild and QA the code core | Done. Remaining: move into the paclet, then the guide page and one reference page per public symbol; harvest the older scattered QEC notebooks |
| 2. Decoder and noise model | Done — three noise levels, the extraction circuit, the detector error model, the memory experiment exactly and by sampling, the Stim bridge |
| 3. Logical error rate and threshold | Rate done at all three levels. `QECThreshold[family, noise]` not written; it needs item 4 to have a family to cross |
| 4. Scalable families | Not started: surface, toric, bivariate-bicycle qLDPC. They slot into `Families.wl` and the rest of the layer needs no change. Matching decoders (PyMatching) plug into the detector model here |
| 5. Fault-tolerant operations, resource counting | Not started |

## Sources

Daniel Gottesman, *Stabilizer Codes and Quantum Error Correction* (Caltech PhD thesis,
[quant-ph/9705052](https://arxiv.org/abs/quant-ph/9705052)) and the 2026 book draft of the
same name; Dennis, Kitaev, Landahl and Preskill,
[quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143) for detectors; Gidney,
[2103.02202](https://arxiv.org/abs/2103.02202) for Stim and the detector error model.
[`References/README.md`](References/README.md) indexes which sections each piece of this
layer rests on, and the tech note's *Where this comes from* section maps every claim to
one of them.
