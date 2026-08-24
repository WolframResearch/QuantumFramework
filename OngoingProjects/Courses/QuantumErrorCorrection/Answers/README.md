# Worked Answers

Worked answers have not yet been authored. They will be added only after the vertical pilot in
`../Plans/Pilot-Vertical-Slice.md` passes the release gate.

## File organization

- `Part-00.md` through `Part-16.md`, one file per Part.
- Question headings copied verbatim from `../Question-List.md`.
- Each answer declares its evidence route from `../Plans/Evidence-Route-Table.md`.
- Simulation artifacts and replay information live under `../Benchmarks/`.
- Frontier answers cite the pinned primary source and snapshot status from
  `../Frontier-Snapshot-2026.md`.

## Per-answer rhythm

1. Physical question.
2. Defining objects and conventions.
3. Derivation or literal computation.
4. Adversarial check.
5. Failure edge.
6. Physical meaning.

## Release gate

An answer is complete only when its equations and examples have been independently reproduced, all
literal code runs, the central check can reject a wrong answer, and no fatal or misleading physics
defect remains. A missing answer is preferable to an unverified one.

