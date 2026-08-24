# Frontier Ledger

Snapshot date: **2026-08-23**. The detailed reconstruction contract lives in
`../Frontier-Snapshot-2026.md`; this ledger controls claim status and refresh priority.

| question | topic | status | directly measured or established | main limitation | refresh priority |
|---|---|---|---|---|---|
| 16.1 | superconducting surface-code scaling | experimental demonstration | distance-dependent logical memory error under named decoders | finite distances and device-specific noise | high |
| 16.2 | real-time surface-code decoding | experimental demonstration | decoder latency and logical performance in a distance-5 memory | not the highest-accuracy offline route at arbitrary distance | high |
| 16.3 | colour-code scaling and logic | experimental demonstration | distance-3-to-5 memory suppression and several logical primitives | memory and gate comparisons are workload-dependent | high |
| 16.4 | neutral-atom architecture | experiment plus architectural inference | component demonstrations for correction, logic, erasure use, reset, and reuse | components were not one indefinitely deep useful computation | high |
| 16.5 | experimental bivariate bicycle codes | experimental demonstration | repeated nonlocal checks and finite logical memories | no below-threshold distance scaling | high |
| 16.6 | low-overhead qLDPC projection | theorem, simulation, proposal | asymptotic code parameters and finite circuit-level projections | connectivity, decoder, routing, and gate assumptions | high |
| 16.7 | concatenated cat memory | experimental demonstration | biased bosonic qubits plus outer repetition correction | incomplete protection and limited distance range | medium |
| 16.8 | GKP qudits | experimental demonstration | beyond-break-even qutrit and ququart memories | one oscillator and no full architecture | medium |
| 16.9 | erasure-converted logical circuits | experimental demonstration | adaptive use of located erasures | distance-two code does not correct arbitrary Pauli errors | high |
| 16.10 | learned decoding | offline experiment plus simulation | improved accuracy on specified data and simulated scaling | training cost, distribution shift, and real-time latency | high |
| 16.11 | trapped-ion qLDPC break-even | experimental demonstration in preprint | nine code instances and break-even for selected memories | no below-threshold qLDPC scaling or universal logical gate set | high |
| 16.12 | measurement-free fault-tolerant logic | experimental demonstration in preprint | coherent logical teleportation, universal encoded gates, and encoded Grover search | error-detecting blocks, postselection, and no scalable repeated QEC | high |
| 16.13 | logical magic-state distillation | experimental demonstration | improved output logical fidelity for distance-3 and distance-5 colour codes | no algorithm-scale factory throughput or target-error demonstration | high |
| 16.14 | cross-platform synthesis | open methodological problem | no single measured result | incomparable metrics tempt false rankings | high |

## Status rules

- A theorem is not relabelled as an experiment when a small instance is implemented.
- A simulation is not relabelled as an overhead result until routing, extraction, decoding, and
  logical operations share one stated hardware model.
- `below threshold` requires improvement with increasing distance under a controlled family.
- `beyond break-even` requires a declared physical comparison baseline.
- `fault-tolerant` identifies a property of an operation or architecture under a fault model; it is
  not a synonym for high fidelity.
- A press release, roadmap, or company benchmark may locate a primary source but cannot replace it.
