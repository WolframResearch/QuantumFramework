# Independent WL quality review — round 1

Verdict: PASS on all seven axes. No artifact edits. `OPEN ISSUES: 0`.

- **A — Documented symbols: PASS.** Live `Information` and `Options` agree
  with the package APIs and return structures. Installed WL 15 documentation
  supports conditional `ForAll`, real polynomial `Resolve` and ties-to-even
  `Round`. The pinned Kickmix source supports the exact modulo additions,
  little-endian registers, lookup cleanup and phase arguments. No new QF APIs
  are used by this optimization.
- **B — Idiomatic WL: PASS.** Row mapping, `MapIndexed`, associations and
  lexical `With` express the transformations. There are no unnecessary
  imperative loops. Row mapping for column-dependent precision is justified.
- **C — Terse reasoning: PASS.** Quantization, signed error, interval accuracy
  and capacity form a coherent argument. The README states the error signs,
  positivity required before squaring, monotone accumulators and scaling.
- **D — Invariants: PASS.** Independent kernel probes reproduce `e0-t*e1`,
  the quantified envelope, normalization identity, unit phase modulus and
  monotone partial-accumulator theorem. All selected interval/capacity
  certificates reproduce. Compiled tests verify restoration, cleanup, phases
  and register dimensions. Hamiltonian Hermiticity, energy conservation and
  eigen-relations are inapplicable to this diagonal phase oracle.
- **E — Nontrivial examples: PASS.** Unequal precisions, two dependent stages,
  recycled coefficient storage, 192 intervals, every valid input and rejected
  precision choices extend the zero-error base case. All nine Python tests
  pass, including 262,144 official small-domain shots.
- **F — Symbolic-first: PASS.** Rational coefficients, symbolic errors,
  quantified inequalities, arbitrary row-count bounds and exact scaling
  precede numerical screening. The general row bound and precision limit
  reproduce. Fixed interval count and source follow the frozen design.
- **G — Regimes: PASS.** Every row of the general/special/limit/numerical/
  failure table has executable evidence. Independent exhaustive screening
  reproduces zero phase failures and zero capacity violations. A separate
  70-digit Decimal check at 387 selected points agrees with the integer
  reference, including its rounding allowance.

The reviewer reproduced all 21 WL checks and all eight certificate sets:
four accepted and four rejected. It independently inspected the terminal
official report: PASS on 200,540,928 shots, actual cost 527,287.988387595,
142 qubits, 3,460 Toffolis, one T gate and 37 fine rotations.

A rebuild of the current default generator reproduced circuit SHA-256
`3e8da3474cb20e7db77bb51a86514c3ef422c9a5673b0955a8106d76969396b4`.
The generator, dependency, saved binary, verifier and all twelve bound
evidence hashes matched. This resolved a provisional concern that the
generator docstring claimed full-run evidence before that run finished.
The final artifact has the required matching evidence.

OPEN ISSUES: 0
