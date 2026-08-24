# Quantum Error Correction

A prerequisite-ordered collection of pedagogical questions for learning quantum error correction
from classical redundancy and the three-qubit code through fault tolerance, topological codes,
quantum LDPC codes, bosonic encodings, hardware-aware correction, and the experimental frontier.

## Status

- The question collection is drafted: 17 Parts, 188 atomic questions.
- The stable curriculum and the dated frontier are separated.
- The authoring, verification, route, prerequisite, and benchmark contracts are drafted.
- Worked answers have not yet been authored.

## Canonical files

- [`Question-List.md`](Question-List.md): the complete learning sequence.
- [`PIPELINE.md`](PIPELINE.md): binding rules for questions, answers, computation, and review.
- [`Frontier-Snapshot-2026.md`](Frontier-Snapshot-2026.md): dated paper-reconstruction briefs for Part 16.
- [`Conceptual-Oral-Questions.md`](Conceptual-Oral-Questions.md): oral checkpoints that should not be forced into code.
- [`Plans/Prerequisite-Graph.md`](Plans/Prerequisite-Graph.md): dependency structure and entry routes.
- [`Plans/Code-Noise-Decoder-Ledger.md`](Plans/Code-Noise-Decoder-Ledger.md): the carrier and contrast ledger.
- [`Plans/Evidence-Route-Table.md`](Plans/Evidence-Route-Table.md): permitted evidence modes and their traps.
- [`Plans/Pilot-Vertical-Slice.md`](Plans/Pilot-Vertical-Slice.md): the first answer-authoring tranche.
- [`Plans/Frontier-Ledger.md`](Plans/Frontier-Ledger.md): claim status, source, and expiry policy.
- [`Benchmarks/README.md`](Benchmarks/README.md): reproducibility contract for simulations and decoders.
- [`Answers/README.md`](Answers/README.md): answer-file organization and release gate.

## Reader

The stable core assumes finite-dimensional quantum mechanics through density operators, channels,
measurement, tensor products, and elementary quantum circuits. Part 0 supplies the algebraic bridge
needed for readers who have not yet used binary symplectic notation or finite-field linear algebra.

Tags describe the role of a question rather than a university degree:

- `[Bridge]`: prerequisite or transition material.
- `[Core]`: durable knowledge expected of anyone working with QEC.
- `[Advanced]`: graduate or research-facing structure.
- `[Frontier-2026]`: a dated result whose wording must be re-audited as the field changes.
