# PIPELINE: Authoring DNA for "Quantum in Continuous Space" Worked Solutions

This file is the standard every answer in the solution set is held to, and the shape each answer
takes. It was distilled from the finished Part 1 (nine answers, refined through direct feedback and
two adversarial review loops) and is meant to govern Parts 2 through 23 without further re-derivation.

The source of truth is the literate-Markdown file (`Solutions-*.md`). The notebook (`.nb`) is a build
output produced by md2nb with `"Evaluate" -> False`, so the reader runs the cells. Conventions are
pure Wolfram Language (no QuantumFramework) and natural units $\hbar=m=1$.

Throughout this file: no em dashes; every formula is valid TeX; subscripts and superscripts braced.

---

## 0. The governing principle

**Earn every claim with a visible computation, using the actual defining object of the problem, in
the smallest honest step, on a nontrivial and diverse example, inside a narrative that teaches by
computing.**

Every rule below is a corollary of this one. The recurring failure mode is settling for what
"technically answers the question" or the path of least resistance (an inherited file, a hand-typed
consequence, a repeated state) instead of this deeper standard. Guard against it by gating each
answer on the deeper standard, not the minimal pass.

When the reader asks "why did you do X" or "why not Y", the expected response is to **opine honestly
on the root cause first, then resolve**, never to defend.

---

## 1. Physics content and honesty

- **Pure Wolfram Language, no QuantumFramework; $\hbar=m=1$** throughout.
- **Scope matches the question's declared scope.** Part 1 is one-dimensional because the questions
  are; do not reach into higher dimensions, spin, or many-body content that later Parts own. Read the
  Part title and the question text, and answer exactly that.
- **Never assert physics.** Every physical claim about a state (that it *is* the bound state of some
  potential, that its energy is $E$, that it is normalized, that it solves the evolution equation) is
  either derived from its defining differential equation or verified in a cell. Part 1 examples:
  1.4 earns "this is the hydrogen $1s$ reduced radial" from the $l=0$ Coulomb radial equation; 1.6
  earns "solves the free time-dependent Schrodinger equation" from a zero residual.
- **The defining object must enter a computation and do real work.** If a physical mechanism is the
  point of the question, compute *with* it, do not smuggle in its consequence. In 1.1 the contact
  potential $V=-g\,\delta(x)$ enters as its own energy $\langle V\rangle=-g\int\delta(x)|\psi|^2\,dx$
  via an actual `Integrate`; the decay rate and binding energy then fall out of minimizing
  $\langle H\rangle(\kappa)$, rather than a jump condition typed in by hand.
- **Cash the abstraction out into the lab where you can.** End an answer by saying what an experiment
  would see (1.9: "a position measurement sees nothing; a momentum measurement sees the entire
  boost"). Name physical scales (Bohr radius, beat frequency, recoil form factor).

---

## 2. The explicit example (the state)

- **Never a Gaussian as a shortcut. Never the trivial base case**: not the ground state alone, not a
  single stationary level with zero current, not $\theta=0$, not a fixed quantum number standing in
  for the general one.
- **Every state is nontrivial and drawn from a recognizable physical field**, and **diverse across
  answers**: no state reused as a crutch. Track what has been used (Appendix A) and pick something new.
- **Choose the state to expose structure the trivial case cannot**: a node, a complex overlap, a
  cusp, a heavy $1/x$ tail, non-normalizability, a genuinely time-dependent current. Prefer a state
  that is a *harder, more honest test* than the smooth default (1.7 uses a kinked tent function
  precisely because "a piecewise-linear trial state is a more honest grid test than a smooth bell
  curve").
- **Ladder recognizable to consequential.** When the richest example needs heavy caveats (the
  non-normalizable Airy beam in 1.6), lead with a simple recognizable normalizable one (the sloshing
  two-level well) and let the exotic case escalate from it. Give both when both teach.
- **Generalize rather than assert a general law.** If a formula holds for all $n$, widen from the
  single case to the family and read the pattern off (1.5 goes from the displaced $n=1$ state to the
  whole Hermite ladder $\{\tfrac12,\tfrac32,\tfrac52,\tfrac72\}$).
- **Discriminate, do not merely illustrate.** Contrast the chosen state with its neighbors so the
  example distinguishes the concept: 1.8 sets the chirp's position-dependent velocity against a plane
  wave's uniform flow; 1.3 sets the sinc's slow $1/x$ tail against exponential decay.

Aim, per answer, to touch these five regimes (from the quality loop):

| regime | what it means here |
|---|---|
| general symbolic | parameters ($\kappa,a,L,q,\dots$) kept symbolic, closed forms |
| exactly solvable | a special case with a known answer to anchor against |
| asymptotic or limiting | a `Limit` or `Series` (weak coupling, large kick, wide well) |
| numerical reference | an independent `NIntegrate` or grid spot-check |
| failure or edge | the boundary where it breaks (non-normalizable, point measure, wall) |

Not every answer fills every row, but a shallow answer usually fills only the first.

---

## 3. Symbolic-first, and cross-checked

- **Full closed form wherever the kernel allows.** Fall to numeric only after checking that symbolic
  genuinely fails, and say so. Keep every parameter symbolic; pin a number only for a plot.
- **Cross-check every nontrivial closed form independently, and the check must be able to refute it.**
  A numeric spot-check that **reuses the definitions** (never a re-typed literal, per the quality
  loop), and/or a limit or asymptotic. In Part 1, 1.2/1.3/1.4/1.5 each pair a closed form with an
  independent `NIntegrate` at one point; 1.3 adds boundary limits, 1.4 adds a large-$q$ `Series`.
- **Two representations must agree.** Symbolic versus numeric; a quantity computed two ways; a general
  proof beside a concrete instance (1.6 proves continuity for arbitrary $f,g$, then for the Airy beam).
- **Probe the kernel, do not recall it** (the `pde-route` discipline). Test a symbol's behavior in a
  live kernel before relying on it. Hard-won this session, all reusable:
  - `D[Abs[x], {x, 2}]` yields **no** `DiracDelta` in current WL (the `Abs` to delta chain is closed);
    put the delta in through an `Integrate`, not a derivative.
  - `NIntegrate` silently misses a Dirac point measure (returns `0.`); `Integrate` is exact. So a
    numeric cross-check is the wrong tool for a delta and the right tool everywhere else.
  - `ComplexExpand` is needed to reduce `Abs[complex]` or `Abs[sum]^2` to real trigonometric form (the
    density of a superposition, 1.6).
  - The norm of a time-dependent state needs `t \[Element] Reals` in the assumptions, or `Abs` keeps
    an `Im[t]` factor.
  - `Minimize` returns a `Piecewise` on the sign of a parameter; for a clean variational condition use
    `Solve[D[E, \[Kappa]] == 0, \[Kappa]]` and confirm the second derivative is positive.
  - `Series[f, {q, Infinity, n}]` for large-argument falloff; `Limit` for boundary behavior.
  - Always give `Integrate`/`Simplify` the positivity and reality assumptions ($\kappa>0$,
    `Element[{x,t}, Reals]`) that produce a clean closed form.

---

## 4. Cell rhythm (learning by computing)

- **One displayed result per cell.** The cell culminates in a single unsuppressed output (its last
  line, no `;`). Supporting definitions and intermediate assignments above it end in `;`. A result you
  will reuse is captured as `name = expr` with **no** trailing `;`, so it both displays and binds. Two
  independently interesting results are two cells. A numeric-comparison cell may bundle setup lines and
  end in a single displayed pair (1.7).
- **A bridge sentence before every cell, ending in a colon, naming the object symbolically**: "Compute
  $\langle V\rangle$:", "Verify continuity $\partial_t\rho+\partial_x j=0$:". No orphan cells.
- **State the closed-form formula in the prose before the code**, then compute it. The reader should
  know what to expect before the cell runs.
- **Drop the trailing `;` when the outcome is interesting** so the value shows. Keep `;` only on pure
  definitions and intermediate assignments.
- **After the cell, interpret the meaning woven into a sentence**, not the raw output token, and let it
  raise the next step. A bare verification ($0$, $1$, `True`) may stand silently; anything else earns
  a sentence about what it means physically.
- **Reuse names across cells** (`ovl`, `dens`, `j`, `T`, `Vexp`); never re-type a wavefunction literal.
- **Idiomatic functional WL, no loops.** Vectorize (1.7 sums a Born probability with
  `Total[Abs[\[Psi][grid]]^2] dx` over `Subdivide`, relying on listability), use `Table` for a family,
  replacement rules for a general proof. No `Do`, `For`, `While`, `AppendTo`.

---

## 5. Prose voice and structure

- **No process or workflow narration in teaching text.** No "we derived rather than assumed", no
  contrast with a rejected draft, no meta about the method used to get here.
- **No out-of-context labels.** Do not call a calculus object a "soliton". Every physical name must fit
  the setting of the question.
- **Self-contained.** Each answer defines the state or operator it uses, computes from that definition,
  and can be read or run on its own. No cross-references to other answers, no reliance on a symbol set
  elsewhere.
- **No em dashes anywhere; all math is valid TeX.** Use colons, commas, parentheses, periods,
  semicolons. Unicode math and raw `x^2`-in-prose are not allowed; wrap in `$...$`.
- **Discuss a result as a sentence with its closed form inline.** Reserve display math ($$...$$) for
  matrices or multi-line derivations, always bridged and interpreted, never a bare `$$label=formula$$`.
- **The section heading is the question verbatim, with its `[BSc]`/`[MSc]` tag.** The tag sets depth:
  BSc answers are foundational and clean; MSc answers go further (operator domains, ordering, deeper
  limits).
- **Part openings give concise, precise common ground.** One paragraph after the Part heading that
  represents the core foundations the Part's questions rest on, each with its defining expression
  (for Part 1: inner product, Born rule, expectation value, current and continuity, phase), very
  concise, not a table of contents.

### 5a. Prose sentences are claims, and are held to the same standard as cells

The earn-it rule of §1 applies to *framing* sentences too, not only to claims a cell can check. A
framing sentence is the most load-bearing thing in an answer, because it is what a reader carries
away after forgetting the algebra. Two tests, both cheap:

- **The stranger test.** Every sentence must survive being read by someone who never saw a previous
  draft and does not know the author exists. "It is not chosen by hand", "rather than assert that",
  "where the delta earns its keep", "not a convenience but a necessity" all fail instantly. A
  correction received must change the physics silently, never appear in the text as self-defense.
- **The motivation-first test.** Before any computation, the reader must already know what it is for
  and why *this* form of it. A cell that arrives before its motivation is a defect even when its
  output is correct. State the strategy ("build $\langle H\rangle(\kappa)$ and minimize"), then compute
  the pieces.

Conflations that repeatedly produced false sentences, all of them worth checking by name:

- the **state** (a ray, fixed only up to a global phase) vs an observable's **spectrum** (the possible
  results) vs the **Born probability** of those results vs the **state update** on measurement vs the
  **ensemble** reading (one run returns one number; a density is what many runs reproduce). Never call
  $|\psi|^2$ "the observable content": it is the outcome distribution of one observable, and the phase
  it discards carries the momentum, the current, and the energy.
- **Hermitian vs self-adjoint.** Hermiticity is necessary, not sufficient; an observable must be
  self-adjoint, and the difference is a statement about domains.
- **Binding energy vs eigenvalue.** The eigenvalue is $-g^2/2$; the binding energy is $+g^2/2$.
- **Representation bias.** $|\tilde\psi(p)|^2$ is as much "the physics" as $|\psi(x)|^2$; say so or
  do not privilege position silently.
- **Unstated hypotheses.** Continuity needs $V$ real; $L^2$ membership does not imply pointwise decay;
  an operator identity needs a domain on which both orderings are defined; a named effect ("form
  factor", "Wannier function of a flat band", "Breit-Wigner in position") must be the right object,
  not a near neighbour.

---

## 6. Routing a differential equation (pde-route)

For any question that solves, routes, or chooses a method for a differential equation, eigenproblem,
or evolution, follow the `pde-route` procedure rather than reaching for a remembered solver:

1. **Classify** (ODE/PDE/DAE; order; linear; elliptic/parabolic/hyperbolic; does superposition hold).
2. **Enumerate** all admissible routes (exact `DSolve`; symbolic eigenproblem `DEigensystem`;
   transform; numeric eigenproblem `NDEigensystem`; numeric evolution `NDSolve`; assemble-then-solve).
3. **Gate by probing**, not recall. Known traps: `NDEigensystem` returns smallest in magnitude, not
   most negative (wrong physics silently for a negative spectrum, fix with a spectral shift);
   `DEigensystem` needs homogeneous boundary conditions and the domain form at $\pm\infty$;
   `NIntegrate` misses point measures.
4. **State the criterion before scoring**: exactness, then falsifiability, then independence, then
   inspectability, then cost.
5. **Pick a primary route plus an independent cross-check** (different machinery, so shared errors
   cannot hide).
6. **Verify** with a refinement sweep and a benchmark. A limit every candidate satisfies is not a
   check. Compare $|\psi|^2$, absolute overlaps, or projectors, never raw eigenvector components
   (the phase is arbitrary).

---

## 7. Self-containment mechanics

Answers run in one kernel in document order for verification, yet each must also run alone. Both
work if symbols are managed carefully:

- **Each answer redefines its own state.** A same-arity redefinition (`\[Psi][x_] := ...`) cleanly
  overwrites the previous one, so `\[Psi]` can be the state in every answer.
- **A different arity does not overwrite, it coexists and can pollute.** If a symbol is used with one
  arity in some answers, do not introduce it with a different arity elsewhere: the well eigenfunctions
  in 1.6 are named `\[Chi][n_, x_]` precisely because `\[Phi]` is a single-argument overlap partner in
  1.2 and 1.4.
- **Do not assign the reserved free symbols** used as generic unknowns (for example a free ODE unknown,
  or generic test functions `f,g` in a proof). Keep them symbolic.
- **Verify the whole document runs in order.** Extract the code cells and run them end to end; the run
  must reach the end with no error before an answer set is declared done.

---

## 8. Verification and build workflow

- **Do not self-certify.** After writing or changing an answer, run the two adversarial loops:
  `/wl-verify` for correctness (a fresh kernel reproduces every claim and worked example) and
  `/wl-quality` for depth and idiom (a fresh critic grades doc-grounding, functional style, physics
  invariants, non-triviality, symbolic-first, and the five regimes). Apply the **deeper** fix, not the
  cheap one: extend to the nontrivial case, redo the derivation symbolically, or cover the missing
  regime. Earn `OPEN ISSUES: 0` from a genuine critic pass before recording it.
- **Run a third loop, aimed at the prose, not the cells.** Both loops above grade code, so both will
  pass a document whose code is right and whose physics sentences are wrong. Part 1 passed a
  three-round `/wl-quality` and `/wl-verify` with zero issues while carrying roughly ten conceptual
  physics errors: a false rigid-transport claim, a survival probability described as an ejection
  probability, two different potentials given the same name, an incorrect self-adjointness argument,
  a continuity proof silently assuming $V=0$ and real, a mis-attributed form factor, a wrong Wannier
  band, and a binding energy with the wrong sign. Point a fresh critic at the **framing and
  foundational prose** with an explicit mandate: claims that are false, overreaching,
  representation-biased, resting on unstated hypotheses, or contradicted by the document's own later
  content. Have it cross-read the parts against each other, since the sharpest disproof of a Part 1
  sentence was sitting in Part 1 itself.
- **Run the extract-and-run harness** on the whole document (all cells, document order) and confirm it
  finishes cleanly. Any changed answer invalidates a prior quality pass on that Part.
- **Build the notebook with md2nb**, `"Evaluate" -> False`, from the Markdown source. Rebuild after any
  content change so the `.nb` stays in sync. Do not hand-edit the `.nb`.

---

## 9. Per-answer template

```
### N.k [BSc|MSc] <the question, verbatim>

<Concept in 2 to 4 sentences: define the quantity and its meaning. Name the chosen state, its
physical field, and one sentence on why it is nontrivial and what structure it exposes.>

<Bridge sentence naming the object symbolically, ending in a colon:>

```wl
<definition ending in ;>
<one computation, displayed (no ; if it is the interesting result), or captured as name = expr>
```

<Interpret the result as physics, in a sentence. If a physical claim was made, the cell just earned
it. Forward-reference the next step.>

<... repeat: bridge, cell, interpret ...>

<Close by cashing the result into a limit, an experiment, or a contrast with a neighbor.>
```

---

## 10. Pre-ship checklist

- [ ] State is nontrivial, field-grounded, and not previously used (Appendix A updated).
- [ ] Every physical claim about the state is derived from its equation or verified in a cell.
- [ ] The defining object of the problem actually appears in a computation.
- [ ] Closed form wherever possible; each nontrivial one cross-checked by an independent route that
      reuses the definitions.
- [ ] At least a limit or an asymptotic where one is meaningful; the five regimes considered.
- [ ] One displayed result per cell; a bridge sentence with a colon before each; formula in prose
      first; `;` dropped on interesting outputs; output interpreted, not echoed.
- [ ] Functional WL, no loops; assumptions supplied for clean closed forms.
- [ ] Self-contained: own definitions, no leaks, non-colliding names.
- [ ] No process narration, no out-of-context labels, no em dashes, valid TeX.
- [ ] Whole document runs in order (clean exit); `/wl-verify` and `/wl-quality` run; notebook rebuilt.

---

## Appendix A. Used-states ledger (extend, never repeat)

Track the primary state or system of each answer so later answers do not reuse it. A transformation
question (a boost, a Fourier transform) may take a known packet as its *object*; that is not a reuse
of the *concept*.

Part 1:

| answer | state / system | field |
|---|---|---|
| 1.1 | Morse oscillator ground state $z^{\lambda-1/2}e^{-z/2}$, $z=2\lambda e^{-ax}$ | molecular vibration |
| 1.2 | Lorentzian $1/(x^2+a^2)$; plane wave | Breit-Wigner resonance; continuum |
| 1.3 | sinc $\sin x/(x\sqrt\pi)$ | flat-band Wannier / single-slit diffraction |
| 1.4 | hydrogen $1s$ reduced radial $2\kappa^{3/2}x\,e^{-\kappa x}$ | atomic |
| 1.5 | displaced Hermite $n=1$ | harmonic oscillator first excited |
| 1.6 | infinite-well two-level superposition; Berry-Balazs Airy beam | particle in a box; accelerating beam |
| 1.7 | tent / Rayleigh-Ritz trial state | finite-element building block |
| 1.8 | chirped sech $\operatorname{sech}(x)\,e^{i\beta x^2}$ | position-dependent phase |
| 1.9 | boosted Hermite $n=1$ | Galilean boost of a known packet |

Part 2 (only answers where the state is the object; the operator-identity answers 2.1, 2.2, 2.4, 2.5,
2.8, 2.9 use generic or convenient test functions, exempt):

| answer | state / system | field |
|---|---|---|
| 2.3 | boosted first-excited $\operatorname{sech}(x)\tanh(x)\,e^{ikx}$ | Poschl-Teller well excited state (node) |
| 2.6 | infinite-well ground state $\sqrt{2/L}\sin(\pi x/L)$ | particle in a box |
| 2.7 | $\operatorname{sech}(x)/\sqrt2$ vs Gaussian | minimum-uncertainty contrast |

Part 3 and beyond: append rows here as answers are written.

---

## Appendix B. Reference exemplar

The finished Part 1 in `Solutions-Parts-1-3.md` is the worked reference for every rule above. When a
rule is ambiguous for a new answer, match the corresponding Part 1 answer: its cell rhythm, its
earn-the-claim structure, its cross-checks, and its closing move.
