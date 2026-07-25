# Shared Brief: Plan Outlines for Parts 4 through 23

Pinned 2026-07-23 after the one-time adversarial-draft interrogation; format revised the same day
on direct feedback (see Revision log). This brief governs every `Part-NN-Plan.md` in this folder.
Parts do not re-negotiate it; a change here propagates to all parts.

PIPELINE.md (one level up) binds every file, including the plans themselves: no em dashes, all math
valid TeX with braced subscripts, no process narration, pure Wolfram Language with $\hbar=m=1$.

## Brief

- **Goal.** Each `Part-NN-Plan.md` is an outline from which the Part's solutions can be authored
  without re-deriving any decision: reading an entry shows how the code will be done.
- **Reader.** Equally (a) a cold authoring session with zero memory of the planning conversation,
  so every example, route, and check is named explicitly, and (b) Mads reviewing; the entry must
  read as a suggested approach, not as bookkeeping.
- **Genre.** One Markdown plan per part at `Plans/Part-NN-Plan.md` (zero-padded), outline level
  only; plus one consolidated `All-Parts-Plan.md` (a straight concatenation, built at aggregation).
- **Scope.** In: Parts 4 through 23, all 130 questions. Out: solution prose, kernel probing inside
  parts (routes cite `Route-Table.md`), any change to Parts 1 through 3.
- **Examples.** Pinned centrally in `Example-Ledger.md`, one primary example per question,
  collision-checked against the PIPELINE Appendix A reserved states and across all 20 plans.
  Escape valve: an example the kernel refuses at authoring time may be replaced within its stated
  class; the swap is recorded in the plan and in the ledger, never silently.
- **Named-trivial policy.** When the question itself names the trivial object (the free Gaussian of
  5.3, the Gaussian ground state of 4.5, the Gaussian trial of 7.3), that object stays the carrier
  and the entry adds one discriminating contrast that exposes what the trivial case hides.
  Everywhere else the PIPELINE rule stands: never a Gaussian as a shortcut, never the trivial base case.
- **Routes.** Every DE-solving entry cites its class in `Route-Table.md` and inherits that class's
  binding traps in its prose. Class C0 (no differential equation) names its WL machinery directly;
  any symbol the planner is not certain exists carries the tag `(verify at authoring)`.
- **Verification.** Every entry states a check that could refute the result and that reuses the
  definitions (never a re-typed literal), woven into the approach prose. The five regimes (general
  symbolic; exactly solvable anchor; asymptotic or limiting; numerical reference; failure or edge)
  are a thinking checklist for the paragraph, not a table artifact: cover the ones that are
  meaningful for the question.

## Entry format (prose approach)

Each entry is the question as a heading, then ONE flowing paragraph (roughly 4 to 8 lines) that
says what to do: the pinned example with its parameters (and the contrast where the ledger names
one), the strategy with the named WL functions, the class traps where they bite, and the
verification, ending with the closing move (a limit, an experiment, or a contrast). No labeled
sublines, no tables. The format exemplar:

```
### 6.2 [BSc] How do I compute the tunneling transmission $T(E)$ through a rectangular barrier?

Set up plane waves in the three regions of a barrier of height $V_0$ and width $L$, match $\psi$
and $\psi'$ at both walls, and Solve the linear system for the transmission amplitude; $T=|t|^2$
comes out in closed form via ComplexExpand and FullSimplify under $k,\kappa,L>0$ (region by
region, because DSolve refuses the piecewise potential silently). Verify $R+T=1$ from the
probability currents, recover the over-barrier resonances by the continuation $\kappa\to iq$
(with $T=1$ exactly at $\sin qL=0$), and cross-check the whole thing against a transfer-matrix
product. Close with the opaque limit $T\sim e^{-2\kappa L}$.
```

## Charge and sign conventions (family-wide)

Minimal coupling is written for a unit positive charge, $i\partial_t\psi=[\tfrac12(\hat p-A)^2+\phi]\psi$,
so the Pauli term that follows from $\tfrac12[\vec\sigma\cdot(\hat p-A)]^2$ is $-\tfrac{g}{2}\vec B\cdot\vec S$,
negative. Any entry writing a magnetic Hamiltonian uses this sign; an entry that wants the opposite
charge says so and carries the consequence through its stated spectrum and degeneracy labels.

## Cross-references between entries

A plan entry may name a sibling entry to divide labour ("4.4 owns the truncated matrices") or to
place a result in the Part. That is planning coordination and is allowed here. It is not prose to
carry into the authored answer: PIPELINE section 5 requires each finished answer to define the state
or operator it uses and to be readable alone, so an authored answer never points at another answer.
Treat every such pointer in these plans as a note to the author, to be spent and dropped.

## Per-plan structure (lean)

1. Header: part title, question count, class census with a Route-Table cite (two lines).
2. Common ground: one short paragraph of the foundations the part rests on, each with its defining
   expression (this seeds the Part opening required by PIPELINE section 5).
3. The per-question entries, in order. Nothing after them.

## Critique protocol (per part)

An outline-level `/wl-quality` critic pass checks: every entry names its pinned example (matching
the ledger), cites or clearly follows its class route with the binding traps inherited, and states
a refuting check that reuses definitions; the meaningful regimes are covered in the prose; no
invented WL symbols (uncertain ones tagged); plan prose conforms to PIPELINE (no em dashes, valid
TeX, no process narration). Maximum 2 rounds per part. Findings still open after round 2 are
reported verbatim in the aggregate index, not silently dropped and not fixed by a third lap.

Back edge: a critic finding that refutes a route goes to `Route-Table.md` as a class-level revision
(one per class, kernel-evidenced), then propagates to every affected part. No re-routing inside a
single part.

## Ownership

Part drafters write only their own `Part-NN-Plan.md`. The coordinator alone writes the shared files
(this brief, `Example-Ledger.md`, `Route-Table.md`, `All-Parts-Plan.md`, the aggregate index) and
commits, straight to main, user identity only.

## Objections preempted (from the pressure-test)

1. A pinned example may fail the kernel: covered by the escape valve above (swap within class, recorded).
2. Parallel drafting versus a global no-repeat ledger: examples are pinned centrally before fan-out;
   drafters inherit, never invent.
3. Verification talk can go vacuous at outline level: every entry must name the concrete check, and
   the critic grades the check's ability to refute, not its presence.
4. C0 has no probed gates: uncertain symbols carry `(verify at authoring)` and the critic hunts
   invented ones.

## Revision log

- 2026-07-23: format revised on direct feedback after the shared files shipped: entries are
  question-heading plus one prose approach paragraph (labeled Example/Route/Verify sublines
  dropped); the per-file regimes table and ledger-additions sections are dropped (regimes stay as
  a prose checklist; ledger bookkeeping lives in `Example-Ledger.md` alone); a consolidated
  `All-Parts-Plan.md` is added at aggregation.
