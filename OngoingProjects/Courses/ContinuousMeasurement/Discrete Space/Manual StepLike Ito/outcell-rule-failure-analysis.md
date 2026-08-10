# Why the Out-cell rule was violated and the gate missed it

Case study: `continuous-measurement-trajectories-V2.md` in this folder. The rule at
issue: *the Out cell holds only the value; a string-keyed Association is allowed only
when every key is a bare canonical name/symbol (or a bare parameter value) so the keys
read as a legend; a key that describes the value's method/condition/meaning, or an
Association that staples two independent facts together, is banned.* Stated in both
`~/.claude/CLAUDE.md` and the `learning-by-computing` skill; the `/essay-verify` loop
ran three adversarial rounds and did not catch the surviving violations.

This report separates the authoring failure from the verification failure, enumerates
every Out-cell Association in V2 with a corrected form, and proposes safeguards. A
working deterministic linter (`outcell_lint.py`, in this folder) is included; it flags
all three violations and zero of the compliant cells, and is clean on the exemplar.

---

## 1. One-paragraph verdict

**Authoring.** The three convergence-study cells return `<|"errors" -> <table>,
"orders" -> <table>|>`: a dict of dicts, keyed by role-words, stapling an error table
to a table derived from it. That is the generic "return a labeled dict of results"
engineering shape, not a cell regenerated against the exemplar (whose Out associations
are always flat legends keyed by a bare symbol, ket, or stat name). **Verification.**
The rule the loop enforces is *mechanical* (decidable from the cell's shape), but it was
handed to an adversarial *same-model* reviewer as prose to re-judge. On that axis the
reviewer shares the author's aesthetic prior ("a labeled dict reads as a legend"), so a
fresh context stripped the author's private rationalization but reproduced the shared
one. The wrong tool was aimed at the rule: a mechanical rule needs a mechanical gate,
and the adversarial loop should be reserved for correctness and argument soundness,
where a fresh context genuinely differs.

---

## 2. Every Out-cell Association in V2, classified

"Out cell" = the expression a cell displays (its last statement, not terminated by `;`
and not a `:=` definition). A `<|...|>` returned from inside a function body, or an
intermediate assignment ending in `;`, is a *data carrier*, not an Out cell, and the
rule does not govern it. That distinction is load-bearing and is where a blunt
"grep every `<|`" reading goes wrong.

| Line | Cell / role | Displayed? | Keys | Verdict |
|---|---|---|---|---|
| 205 | inside `simulate := Module[... <|...|>]` | no (return value) | Times, States, Record, Noise | compliant (data carrier; canonical field names) |
| 218 | inside `infer := <|...|>` | no (return value) | Times, States, Record | compliant (data carrier) |
| 583 | Part VII, `AssociationThread[{1.0,0.5,0.0}, ...]` | **yes** | `1.0, 0.5, 0.0` | **compliant** (bare parameter-value keys; one swept quantity) |
| 639 | inside `observedOrders := AssociationThread[...]` | no (helper body) | inherited step sizes | compliant (data carrier) |
| 643 | `deterministicErrors = AssociationThread[...];` | no (`;`-suppressed) | step sizes | compliant (data carrier; bare param keys) |
| 649 | Part VIII, `<|"errors"->…, "orders"->…|>` | **yes** | `errors, orders` | **VIOLATION** |
| 673 | inside `coupledFamily := … AssociationThread[…]` | no (helper body) | step sizes | compliant (data carrier) |
| 683 | inside `coupledFamily` Module, `<|"States"->…,"Noise"->…|>` | no (return value) | States, Noise | compliant (data carrier) |
| 688 | inside `adjacentErrors := AssociationThread[…]` | no (helper body) | step sizes | compliant (data carrier) |
| 696 | inside `rootMeanSquare := AssociationThread[…]` | no (helper body) | inherited keys | compliant (data carrier) |
| 706 | Appendix pathwise, `<|"errors"->…, "orders"->…|>` | **yes** | `errors, orders` | **VIOLATION** |
| 730 | Appendix innovations, `<|"mean"->…, "variance"->…, "lag-one correlation"->…|>` | **yes** | `mean, variance, lag-one correlation` | **compliant** |
| 744 | Appendix non-commuting, `<|"errors"->…, "orders"->…|>` | **yes** | `errors, orders` | **VIOLATION** |

**Three violations: lines 649, 706, 744.** All three share one shape,
`<|"errors" -> X, "orders" -> observedOrders[X]|>`, and fail on two independent grounds:

- **Staple of a quantity with a transform of itself.** `"orders"`'s value is
  `observedOrders[X]`, computed *from* `"errors"`'s value `X`. Two dependent derived
  objects packed into one Out cell to save a cell. This is precisely the banned
  "staples two independent facts together."
- **Dict of dicts.** Each value (`deterministicErrors`, `maximumErrors`, `ncErrors`,
  and their `observedOrders[…]`) is itself an Association. The outer keys are section
  headers on two sub-tables, i.e. annotations, not a legend over parallel scalar values.
  The tell that the Out cell is being used as a labeled report.

The keys `"errors"`/`"orders"` are role-words. They are not descriptive-token keys of
the vivid `"grid-substituted"` kind, which is exactly why a keyword-surface check waves
them through.

**Two compliant cells a blunt reading would wrongly flag** (this is where I disagree
with "every keyed dict is suspect"):

- **730** `<|"mean" -> Mean[innovations], "variance" -> Variance[innovations],
  "lag-one correlation" -> Correlation[Most[innovations], Rest[innovations]]|>`. Three
  *parallel scalar* reductions of *one* object (the innovation sequence), each key the
  bare canonical name of its slot's value. This is the `{"Mean","Variance"}` legend the
  rule explicitly allows, extended by one more standard statistic. None of the values is
  an association; none is a transform of another; together they answer one question (are
  the innovations standard-normal and white). Compliant. The exemplar uses this exact
  form (`introduction-to-QIS-revised.md:2883,2891`).
- **583** `AssociationThread[{1.0, 0.5, 0.0}, <Frobenius distances>]`. Keys are the bare
  parameter values (the efficiencies), values one derived quantity swept over them. The
  parametric legend. Compliant.

The distinction that separates the three violations from 730 and 583 is **structural,
not lexical**: a violation is a dict whose values are themselves tables, or where one
value is computed from a sibling. A compliant legend maps bare keys to *parallel
scalars* of one object, or a bare parameter value to its one result. A linter can decide
this without judging whether a word "sounds canonical" (see §6).

---

## 3. Corrected form

The payoff of each convergence cell is the observed *order*; the error table is
supporting scaffolding. Give each its own cell under its own bridge, each Out a flat
legend keyed by the bare parameter value. For the Part VIII deterministic cell (V2
lines 636 to 652), replace the single stapled cell with two:

Bridge: *"...halve the step repeatedly and read the trace distance from the
master-equation answer at each step size:"*

```wl
deterministicSteps = {0.08, 0.04, 0.02, 0.01};
deterministicErrors = AssociationThread[deterministicSteps,
   traceDistance[
      Last@simulate[initialState, hamiltonian, measuredOperators, {0.}, {},
         #, finalTime]["States"],
      masterSolution[finalTime]] & /@ deterministicSteps];
deterministicErrors
```

*As one can see, the error falls in a regular way as the step shrinks.*

Bridge: *"Estimate the local order from the base-2 log of the ratio of successive
errors:"*

```wl
observedOrders[deterministicErrors]
```

*The estimated order sits near one: halving the step roughly halves the error, the mark
of a first-order scheme.*

Each Out cell is now a bare legend `<|0.08 -> …, 0.04 -> …, …|>` keyed by the step size
(a bare parameter value), and the two distinct quantities live in two cells. The
appendix cells at 706 and 744 take the identical fix (`maximumErrors` then
`observedOrders[maximumErrors]`; `ncErrors` then `observedOrders[ncErrors]`), with the
"order below one" reading kept as prose beneath the second cell. A terser single-cell
alternative is to display only `observedOrders[…]` (the payoff) and describe the error
decay in the bridge; the two-cell form is more faithful because it lets the reader see
the errors shrink. Cells 730 and 583 need no change.

---

## 4. Root cause: the six preliminary points, tested

**(1) Write-from-prior. HOLDS.** The three cells are built by helper functions
(`observedOrders`, `rootMeanSquare`) that return Associations, then wrapped in a two-key
report dict. That is the "assemble a labeled result bundle" reflex. The exemplar has no
dict-of-dicts Out cell anywhere (checked: 20 associations, all flat legends keyed by
symbol/ket/stat-name). So the shape did not come from regenerating against the exemplar.

**(2) Propagation. HOLDS, with a refinement.** The `<|"errors"->X,"orders"->…|>` template
is reused verbatim three times. But it did *not* infect 730, which was written to a
different, compliant mold in the same appendix. So propagation is real but *template-local*:
the author can produce compliant Out cells and did, a few lines away. The failure is not a
global inability, it is a single bad template copied without re-examining its class.

**(3) Patch-the-instance. HOLDS (structurally); the specific old token is brief-sourced.**
The current keys `"errors"`/`"orders"` are exactly what a token-only relabel of a
previously-descriptive key leaves behind: the banned *structure* (dict-of-dicts staple)
survives untouched. I could not find `"rms maximum trace distance"` in git history (V2 was
committed once, `38cc9ea1`), so the precise round-3 rename is known only from the session
brief, not independently confirmable. What *is* confirmable: V1 (the pre-rewrite sibling)
is full of descriptive-token keys, including the literal banned example
`"largest purity error at eta = 1"` (V1:816), plus `"recorded channels per step"`,
`"state dimensions"`. V2 removed those but kept the structural staple. The cleanup tracked
the rule's *examples*, not its intent.

**(4) Shared blind spot. HOLDS, and is the deepest point.** Rounds 1 to 2 reportedly
rationalized "stays inside the `{Mean,Variance}` idiom." That rationalization is available
*because a legitimate instance of the idiom exists in the same essay* (730). The allowed
exception is the handle: any labeled dict can be waved through as "just the legend idiom."
A fresh context removes the author's ego but not this analogical over-extension, because it
is the same model with the same aesthetic prior. This is the strongest argument for taking
the Out-cell decision *away* from model judgment and giving it to a deterministic gate.

**(5) Checklist-surface-matching. PARTLY HOLDS; the framing is too kind to the prompt.**
The claim is that the verifier prompt is a keyword filter. It is not *only* that: check 4
in `essay-verifier.md` does contain the correct clause "or an Association that staples two
independent facts together." So the rule was present. But that clause is abstract and
carries no example, while the four vivid examples right before it are all
descriptive-token keys (`"grid-substituted"`, `"purity at eta = 1"`). The agent anchored on
the concrete examples and left the abstract clause inert. The fix is therefore *not* "add
the staples clause" (it is there) but "make it operational": give the structural test and
force a per-cell answer to the deep question. Diagnosing this as "the prompt was just a
keyword list" would misdirect the fix.

**(6) Cost-as-false-confidence. HOLDS directionally; the mechanism is a missing gate, not
a crowded-out one.** Three expensive rounds did manufacture "surely it is clean." But no
cheap deterministic check was ever *in* the pipeline to be displaced; the expense filled a
vacuum where a ten-line linter should have run first. The lesson is ordering: run the
mechanical gate before the adversarial loop, so the loop never has to spend judgment on the
mechanical class at all.

**Additions the six miss:**

**(7) Tool-to-rule mismatch (subsumes 4, 5, 6).** The Out-cell rule is decidable from
cell shape. Enforcing a decidable rule by adversarial same-model review is the category
error. Mechanical rules get mechanical gates; adversarial review is for the undecidable
axes (is the argument sound, is the number real, is the coinage a real synonym).

**(8) The exception is the attack surface.** "Banned except when the keys read as a
legend" has no decision procedure attached, so the exception swallows the rule by analogy.
Any "X banned except Y" needs Y mechanically decidable, or it launders violations. Here Y
should be stated structurally: *values are parallel scalars of one object, or a parameter
value maps to its one result; never a table, never a transform of a sibling.*

**(9) No anti-cosmetic-patch step.** The loop re-spawns a fresh verifier on the revised
file but never asks "was the prior finding resolved structurally, or only relabeled?" A
rename that preserves the banned structure passes, because the new verifier sees only the
new token. The re-spawn must carry the prior finding and demand structural resolution.

---

## 5. Safeguards, by leverage

1. **Deterministic pre-check, run before the adversarial loop (highest leverage).** A
   linter that enumerates every *displayed* Out-cell Association and hard-fails the
   mechanical classes: (A) a value that is a subterm of a sibling value (staple of a
   quantity and its transform), (B) a value that is itself an Association/table (dict of
   dicts), (C) a key embedding a condition or relation (`"purity at eta = 1"`). Built and
   validated: `outcell_lint.py` (§6). This alone catches the exact class that survived
   three rounds, in well under a second, with zero exemplar false positives. Wire it as
   step 0 of `/essay-verify` and as a fast local check the author runs while writing.

2. **Re-scope the adversarial verifier off the mechanical axis and onto intent.** Once the
   linter owns staples and condition-keys, the verifier's Out-cell check becomes: *for
   every displayed Out-cell Association, quote the rule, then answer one question in
   writing: is this cell a legend (bare keys over parallel scalars of one object, or a
   parameter value over its one result), or a labeled report? If you cannot say "legend"
   in one sentence without hedging, it is a report.* Require the answer per cell, not a
   grep summary. Forbid the phrase "reads as a legend" as a conclusion without naming the
   one object whose parallel values the keys label.

3. **Prompt the loop against cosmetic patches.** When re-spawning after a fix, pass the
   prior finding verbatim and instruct: *verify the previous violation was resolved
   structurally (the shape changed), not by renaming a key or token. A relabel that keeps
   the structure is an open finding.*

4. **Write the rule's exception as a decision procedure, not a vibe** (edits in §7). In
   both CLAUDE.md and the skill, replace "keys read as a legend" with the structural test,
   so author and verifier apply the same decidable criterion.

5. **A standing note on when the adversarial loop is the wrong tool.** Mechanical/structural
   rules (Out-cell shape, em dashes, fence structure, `Print` bans) get deterministic gates.
   The same-model adversarial loop earns its cost only where a fresh context genuinely
   judges differently: correctness (does it run, is the number real), argument soundness,
   and semantic coinage (is this noun a real synonym for a named quantity). Spending three
   rounds on a shared-prior aesthetic axis buys false confidence, not coverage.

---

## 6. The linter: reach, limits, results

`outcell_lint.py` (this folder) is pure static analysis, no Wolfram kernel. It extracts
` ```wl ` fences, finds each cell's *displayed* Out expression (skips `;`-terminated
statements and `:=` definitions, so function-internal data carriers like `simulate`'s
return are correctly ignored), and checks any Out-cell `<|...|>` or `AssociationThread`.

**Hard rules (deterministic, zero false positives on the exemplar):**
- `STAPLE-TRANSFORM`: one top-level value is a subterm of another (transform of a sibling).
- `NESTED-DASHBOARD`: a top-level value is an Association/table (literal, or a symbol/helper
  known to return one, tracked from assignments in the file).
- `CONDITION-IN-KEY`: a string key embeds a relational operator or a condition clause
  (`at`/`when`/`given`/`for` + parameter).

**Soft advisory (heuristic, routes attention, never gates):**
- `KEY-LEXICON`: a string key that is not obviously a bare canonical token (symbol, ket,
  stat name, short token). Flags for human/verifier confirmation.

**Results:**

| File | Hard | Soft | Meaning |
|---|---|---|---|
| exemplar `introduction-to-QIS-revised.md` | **0** | 0 | the compliant reference passes cleanly |
| `continuous-measurement-trajectories-V2.md` | **7** | 1 | all on 649/706/744; soft note on "lag-one correlation" |
| `continuous-measurement-trajectories-V1.md` (control) | 3 | 13 | condition-keys caught hard; flat descriptive staples soft |

The 7 hard hits on V2 localize to exactly the three violating cells and nothing else;
730 and 583 are not hard-flagged. The soft note on 730's `"lag-one correlation"` is the
honest boundary: a two-word canonical name trips the advisory channel and a human confirms
it in one look.

**Honest limits.** The hard gate is *sound* (no exemplar false positives) but not
*complete*. A flat Out cell stapling several *independent* facts with plausibly-canonical
short keys (`<|"step sizes"->…, "record dimensions"->…, "state dimensions"->…|>`) is not
mechanically separable from a legitimate multi-statistic legend
(`<|"mean"->…, "variance"->…|>`) without semantics: both are flat dicts of scalars with
short keys. That residue is exactly what the re-scoped verifier (safeguard 2) must judge.
The division of labor: the linter kills the decidable classes cheaply and certainly; the
verifier judges the semantic residue and is told not to re-litigate what the linter owns.

---

## 7. Proposed edits (for approval; not applied)

These touch the user's global config and skills. Proposed, not made.

### 7a. `~/.claude/CLAUDE.md`, "Computational Essays" section

Append a decision procedure to the Out-cell bullet, so the exception is decidable:

> **Out-cell Association, decide structurally.** A displayed Out-cell Association is a
> *legend* only if either every key is a bare canonical name/symbol over a *parallel
> scalar* of one object (`{"Mean","Variance"}` over `{Mean[x],Variance[x]}`), or every
> key is a bare parameter value over its one derived result
> (`AssociationThread[{1.0,0.5,0.0}, distances]`). It is a *labeled report* (banned) if
> any value is itself an Association/table (dict of dicts), or one value is computed from
> a sibling value (a quantity beside its own transform), or a key embeds a condition
> (`"purity at eta = 1"`). A report is two cells, each a legend under its own bridge. Run
> the deterministic Out-cell check before `/essay-verify`; it owns this class.

And in the enforcement line, add: *before the adversarial loop, run the deterministic
Out-cell linter; the loop is for the semantic residue, not the mechanical class.*

### 7b. `learning-by-computing/SKILL.md`, "Code DNA", the Out-cell bullet

Replace the current bullet's exception wording with the structural test above (same text),
so the skill and CLAUDE.md state one decidable criterion rather than "keys read as a
legend," which is the phrase the reviewer used to wave the violations through.

### 7c. `essay-verify/SKILL.md`, Procedure

Add **step 0** before round 1:

> 0. **Deterministic pre-check.** Run `outcell_lint.py` on the changed essay. Fix every
>    HARD finding before spawning any verifier; a mechanical violation must never consume
>    an adversarial round. Carry the SOFT findings into the verifier prompt as cells to
>    judge.

And amend step 4 (re-spawn):

> When re-spawning after a fix, include the prior round's findings verbatim and instruct
> the new verifier to confirm each was resolved *structurally* (the cell's shape changed),
> not by renaming a key. A relabel that preserves the structure is an open finding.

### 7d. `~/.claude/agents/essay-verifier.md`, check 4

Replace the keyword-leaning check 4 with an intent-first, per-cell one:

> 4. **The Out cell is a legend, not a report.** For *every displayed* Out-cell
>    Association (ignore `<|...|>` returned from inside a function or suppressed by `;`),
>    quote the rule, then answer in one written sentence: *what single object's parallel
>    values (or what parameter's one result) do these keys label?* If you cannot, it is a
>    labeled report and a finding. Mechanical staples (a value that is itself a table; a
>    value computed from a sibling; a key embedding a condition) are the linter's job and
>    should already be gone; if you see one, the pre-check was skipped, so report it. Your
>    judgment is for the residue the linter cannot decide: a flat dict of short-keyed
>    scalars that are *not* parallel reductions of one object (`"step sizes"`,
>    `"record dimensions"`, `"state dimensions"` are three different objects, not a
>    legend). Do not conclude "reads as a legend" without naming the one object.

---

*Linter: `outcell_lint.py` in this folder. Report and linter written on `main`; the essay
itself is unmodified.*
