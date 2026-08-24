# Frontier Snapshot 2026: Reconstruction Briefs

Snapshot date: **2026-08-23**. This file expands Part 16 of `Question-List.md`. It does not crown a
winning platform or code. Each entry identifies a primary result, a bounded reconstruction, and the
claim that the reconstruction is not allowed to make.

## 16.1 Below-threshold superconducting surface-code memory

- **Status:** experimental demonstration.
- **Primary source:** [Quantum error correction below the surface code threshold, Nature (2024)](https://www.nature.com/articles/s41586-024-08449-y).
- **Pinned result:** logical error per cycle decreases through distances 3, 5, and 7; the reported
  distance-5 and distance-7 memories operate beyond break-even, with suppression factor above two
  for a two-unit increase in distance under the experiment's decoder and cycle definition.
- **Reconstruction:** digitize or use source data for logical error per cycle versus distance,
  reproduce the reported suppression factor with uncertainty, and state the exact decoder and
  memory protocol attached to it.
- **Refuting check:** repeat the estimate after changing the fitted distance range and after using
  a different reported decoder. A stable conclusion must survive the stated variations within
  uncertainty.
- **Not established:** a universal surface-code threshold independent of device, circuit, decoder,
  leakage handling, or correlated-error model.

## 16.2 Real-time decoding

- **Status:** experimental demonstration.
- **Primary source:** [Quantum error correction below the surface code threshold, Nature (2024)](https://www.nature.com/articles/s41586-024-08449-y).
- **Pinned result:** a distance-5 memory was operated with an integrated real-time decoder while
  preserving below-threshold performance at a syndrome-cycle time of about $1.1\,\mu\mathrm{s}$.
- **Reconstruction:** compare the decoder-latency distribution, not only its mean, with the rate at
  which syndrome data arrive; determine whether backlog grows with experiment duration.
- **Refuting check:** inspect a high percentile and the latency-versus-round-count scaling. A mean
  below one cycle is insufficient if the tail or backlog diverges.
- **Not established:** that the highest-accuracy offline decoder in the paper is already available
  with the same accuracy and resource use in real time at arbitrary distance.

## 16.3 Colour-code scaling and logic

- **Status:** experimental demonstration.
- **Primary sources:** [Scaling and logic in the colour code on a superconducting quantum processor,
  Nature (2025)](https://www.nature.com/articles/s41586-025-09061-4) and the same-platform surface-code
  comparison reported there.
- **Pinned result:** logical memory errors were suppressed when scaling a colour code from distance
  3 to 5, while transversal logical Clifford gates, magic-state injection, teleportation, and
  lattice-surgery primitives were demonstrated.
- **Reconstruction:** reproduce the reported distance-3-to-5 suppression and construct a two-axis
  comparison: memory suppression on one axis and logical-operation cost or fidelity on the other.
- **Refuting check:** remove any operation that was postselected or not distance-scaled and see
  whether the qualitative comparison changes.
- **Not established:** that either colour or surface codes dominate for every physical error rate,
  gate set, logical workload, or hardware layout.

## 16.4 Neutral-atom fault-tolerant architecture

- **Status:** experimental demonstration plus architectural inference.
- **Primary source:** [A fault-tolerant neutral-atom architecture for universal quantum computation,
  Nature (2025)](https://www.nature.com/articles/s41586-025-09848-5).
- **Pinned result:** arrays of up to 448 atoms were used to study below-threshold surface-code
  correction, fault-tolerant logical operations, teleportation-based synthesis, erasure-aware
  decoding, reset, and reuse.
- **Reconstruction:** make a claim-evidence table in which every advertised architectural
  ingredient is tied to the exact code, distance, rounds, physical qubits, and measured observable
  used to demonstrate it.
- **Refuting check:** mark any edge in the proposed scalable architecture that was inferred from
  separate component demonstrations rather than exercised in one deep error-corrected algorithm.
- **Not established:** execution of an arbitrarily deep useful algorithm at fully scaled logical
  error and overhead.

## 16.5 Experimental bivariate bicycle codes

- **Status:** experimental demonstration.
- **Primary source:** [Demonstration of low-overhead quantum error correction codes, Nature Physics
  (2026)](https://www.nature.com/articles/s41567-025-03157-4).
- **Pinned result:** a 32-qubit superconducting processor implemented repeated nonlocal stabilizer
  extraction for a distance-4 bivariate bicycle code encoding four logical qubits and a distance-3
  punctured code encoding six logical qubits.
- **Reconstruction:** rebuild the two check matrices, verify their ranks and commutation, reproduce
  the syndrome-cycle schedule, and compare measured per-logical-qubit failure with a code-capacity
  baseline.
- **Refuting check:** test whether the measured finite instances show distance scaling or only the
  feasibility of extracting their checks.
- **Not established:** below-threshold qLDPC scaling, a logical gate set, or lower end-to-end
  computational overhead than a surface-code architecture.

## 16.6 Projected low-overhead qLDPC memory

- **Status:** theorem plus circuit-level simulation and architectural proposal.
- **Primary sources:** [High-threshold and low-overhead fault-tolerant quantum memory, Nature
  (2024)](https://www.nature.com/articles/s41586-024-07107-7) and [Asymptotically Good Quantum and
  Locally Testable Classical LDPC Codes](https://arxiv.org/abs/2111.03654).
- **Pinned result:** asymptotically good qLDPC families establish constant rate and linear distance
  in principle; finite bivariate bicycle proposals report competitive circuit-level thresholds and
  substantially reduced encoding overhead under stated connectivity and noise assumptions.
- **Reconstruction:** separate the asymptotic theorem from the finite-code simulation. Reproduce
  one finite code's parameters, check schedule, decoder choice, threshold fit, and overhead ratio.
- **Refuting check:** replace long-range connectivity with a stated routing model and include its
  faults and time cost.
- **Not established:** that asymptotic goodness automatically supplies practical decoding,
  geometrically local checks, low-latency logical gates, or low finite-size overhead.

## 16.7 Concatenated cat-qubit memory

- **Status:** experimental demonstration.
- **Primary source:** [Hardware-efficient quantum error correction via concatenated bosonic qubits,
  Nature (2025)](https://www.nature.com/articles/s41586-025-08642-7).
- **Pinned result:** five biased cat qubits formed an outer distance-5 repetition code; phase-flip
  correction operated below its threshold while increasing cat size suppressed the complementary
  logical bit-flip channel.
- **Reconstruction:** decompose the measured logical error into phase- and bit-flip contributions,
  reproduce their dependence on repetition distance and cat size, and identify the best operating
  point under the paper's cycle definition.
- **Refuting check:** compare distance 5 with distance 3 including every additional fault location;
  do not infer successful scaling merely from the nominal distance.
- **Not established:** full arbitrary-noise protection by the outer repetition code alone or a
  universal fault-tolerant logical gate set.

## 16.8 Beyond-break-even GKP qudits

- **Status:** experimental demonstration.
- **Primary source:** [Quantum error correction of qudits beyond break-even, Nature
  (2025)](https://www.nature.com/articles/s41586-025-08899-y).
- **Pinned result:** GKP-encoded logical qutrit and ququart memories were reported beyond break-even,
  with gains of approximately 1.82 and 1.87 relative to the best corresponding physical memories
  available in the device.
- **Reconstruction:** recover the lifetime or error metric used to define gain, reproduce both gain
  ratios with uncertainty, and identify how reinforcement learning entered protocol optimization.
- **Refuting check:** compare against the best relevant physical qutrit or ququart, not a qubit or an
  average oscillator lifetime chosen after the fact.
- **Not established:** a multi-mode fault-tolerant architecture or that increasing local Hilbert
  space removes the need for outer-code protection.

## 16.9 Erasure-converted logical circuits

- **Status:** experimental demonstration.
- **Primary source:** [Logical qubits with erasure conversion using metastable neutral atoms, Nature
  Physics (2026)](https://www.nature.com/articles/s41567-026-03309-0).
- **Pinned result:** metastable ytterbium qubits supplied located erasure information used in
  adaptive execution, a $[[4,2,2]]$ encoding, and teleportation with conditionally selected
  ancillas.
- **Reconstruction:** write the decoder's information set with and without erasure flags and
  reproduce the logical acceptance or fidelity advantage attributable to the flags.
- **Refuting check:** inject an unlocated Pauli error and confirm that the distance-two code does not
  acquire a correction capability it does not have.
- **Not established:** correction of arbitrary single-qubit Pauli errors by the demonstrated
  $[[4,2,2]]$ protocol.

## 16.10 Learned surface-code decoding

- **Status:** offline decoder demonstration on experimental data and simulated scaling.
- **Primary source:** [Learning high-accuracy error decoding for quantum processors, Nature
  (2024)](https://www.nature.com/articles/s41586-024-08148-8).
- **Pinned result:** AlphaQubit improved logical-error suppression on experimental surface-code
  data relative to the reported matching and tensor-network baselines and was evaluated on
  simulated distances through 11.
- **Reconstruction:** reproduce one fixed-distance test split and one cross-distance generalization
  result with identical event data for all decoders.
- **Refuting check:** exclude training examples from the tested device or distance and report both
  the accuracy change and training-data requirement.
- **Not established:** real-time latency at arbitrary distance, low training cost, or universal
  superiority under noise outside the training distribution.

## 16.11 Trapped-ion qLDPC break-even

- **Status:** experimental demonstration reported in an arXiv preprint.
- **Primary source:** [Breakeven demonstration of quantum low-density parity-check codes
  (2026)](https://arxiv.org/abs/2606.06455).
- **Pinned result:** nine code instances from qLDPC, topological, and concatenated families were run
  on one trapped-ion device; an 18-physical-qubit qLDPC instance encoded four logical qubits and
  achieved a logical error rate reported as up to nine times lower than an earlier superconducting
  implementation of a similar code. Some tested instances had lifetimes comparable to or slightly
  longer than the experiment's physical-qubit baseline.
- **Reconstruction:** rebuild one selected check matrix and verify its parameters, reconstruct its
  memory rounds and lifetime fit, and state the exact physical baseline used for break-even. Compare
  its per-logical error and connectivity assumptions with the superconducting bivariate-bicycle
  experiment in 16.5.
- **Refuting check:** repeat the break-even comparison with the most conservative relevant physical
  lifetime and include state preparation, measurement, reset, cooling, and optical-metastable-ground
  operations in the accounting.
- **Not established:** asymptotic below-threshold qLDPC scaling or universal fault-tolerant
  computation with the demonstrated blocks.

## 16.12 Measurement-free fault-tolerant logic

- **Status:** experimental demonstration reported in an arXiv preprint.
- **Primary source:** [Demonstration of measurement-free universal fault-tolerant quantum
  computation (2025)](https://arxiv.org/abs/2506.22600).
- **Pinned result:** logical state teleportation was demonstrated between two four-qubit
  error-detecting codes without mid-circuit measurement, and a coherently executed universal gate
  set based on state injection was implemented on an eight-qubit code hosting three logical qubits.
  The toolbox was used for an encoded Grover-search demonstration.
- **Reconstruction:** compile the logical circuits, verify their claimed logical actions and fault
  propagation, record every acceptance or postselection condition, and compare with a
  measurement-and-feed-forward implementation under the same physical fault model.
- **Refuting check:** inject one fault at every circuit location and classify it as corrected,
  detected and discarded, or undetected; do not count postselected detection as active correction.
- **Not established:** scalable repeated active QEC without measurement or below-threshold
  performance for arbitrarily deep universal computation.

## 16.13 Logical magic-state distillation

- **Status:** experimental demonstration.
- **Primary source:** [Experimental demonstration of logical magic state distillation, Nature
  (2025)](https://www.nature.com/articles/s41586-025-09367-3).
- **Pinned result:** neutral-atom experiments distilled logical magic states encoded in distance-3
  and distance-5 colour codes and reported higher output logical fidelity than the input logical
  magic states.
- **Reconstruction:** reproduce the protocol's ideal distillation map, acceptance probability, and
  input-to-output fidelity comparison with uncertainty; then account separately for encoding,
  state preparation, logical operations, readout, postselection, and total physical resources.
- **Refuting check:** compare accepted distilled outputs with the best non-distilled logical input
  available under the same resource and postselection accounting, rather than with an easier raw
  baseline.
- **Not established:** a magic-state factory operating at algorithmic target error and throughput,
  or the asymptotic resource scaling of such a factory.

## 16.14 Cross-platform evidence table

- **Status:** synthesis of experimental demonstrations, simulations, and open problems.
- **Primary sources:** the sources in 16.1 through 16.13.
- **Pinned result:** none; the task is methodological.
- **Reconstruction:** build a table with separate columns for physical platform, code, $[[n,k,d]]$,
  corrected noise, rounds, logical observable, decoder, real-time status, distance scaling,
  break-even baseline, demonstrated logical operations, and projected overhead.
- **Refuting check:** remove any scalar "best platform" score and verify that no empty cell is
  silently treated as a zero or a failure.
- **Not established:** a total ordering of platforms whose experiments answer different questions.

## Refresh rule

Revisit this file by **2027-08-23**. A refresh must check retractions or corrections, later versions,
new distance-scaling data, decoder latency, logical operations, and whether proposals have become
experiments. Preserve old snapshots in version control rather than rewriting the historical record.
