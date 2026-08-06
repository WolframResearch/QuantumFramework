# The Question-Chain Method (portable prompt)

You are given a research topic, an informal research question, or a draft paper. Your task is to compile it into a **formal question chain**: a single self-contained document, structured like a Lean/theorem-prover file, that captures the topic's minimal set of definitions and open problems so completely that a competent reader (or another AI) given the document and nothing else can state and begin solving the problems.

## Why this method exists (the goals)

1. **Kill vague questions.** Informal question lists fail two tests: (a) handed to a stranger with no context, the questions are not actionable ("what survives that choice?", "what is the answer good for?" are pronouns wearing question marks); (b) the questions do not force the object they are about: hidden choices live in surrounding prose. The method fixes both by the theorem-prover discipline: every reference is a *named, previously defined item*, never a pronoun, so questions can be short AND self-contained at the same time.
2. **Make the backbone load-bearing.** Every section of downstream work must answer exactly one chain item; orphan content signals a missing question, and a question no section answers signals dead weight.
3. **Make derivability checkable.** The central problem must be solvable *from the document alone*, by mechanical computation, with no undeclared inputs. This is the document's acceptance test, and it is executed, not asserted.
4. **Expose choices instead of hiding them.** Any convention the answer depends on (a chosen witness, gauge, representation, normalization) must appear as an explicit named item, with a separate problem auditing what depends on it.

## The artifact: structure

Produce one document with these layers, in this order:

- **Preamble**: the discipline in three sentences; audit status (rounds run, findings, repairs; never claim more than was verified).
- **Primitives**: the ambient standard mathematics, declared once (sets, spaces, operations, special constants, with explicit formulas for anything conventional, e.g. named matrices). Primitives are ambient: never listed in dependency lines.
- **The chain**: numbered items D1, D2, ... interleaved with problems P1, P2, ...
  - **Definition (Dn)**: introduces one object, using only primitives and earlier items. A definition may assert a membership; that assertion is a proof obligation discharged trivially or in its lemma block.
  - **Lemma (na, nb, ...)** attached to its item: a true, provable statement (stated, not proven). Tag each lemma's role: *certification* (justifies a name or locates an object), *interpretation* (gives meaning; consumed by no problem's computation), or *substrate* (supplies coordinates/types the computations use).
  - **Problem (Pn)**: a well-formed task over earlier items only, with a stated deliverable format ("give in closed form", "give necessary and sufficient conditions on X", never a bare "characterize"). Include a solvability note when the problem is a finite computation from the file.
  - **Reading** line under every item: one to three sentences of interpretive prose for the target expert audience, explicitly outside the formal chain. Readings must *unpack*, never merely *name* (explaining a definition by the jargon it formalizes is the anti-pattern). Readings never carry information a problem's solution needs.
- **Meta-problems (M1, M2, ...)**: the non-formalizable assessment duties (what is it good for; where has this been studied), explicitly quarantined outside the formal chain.
- **Minimality record**: per item, why it survives the deletion test (remove it: what breaks?) and the collapse test (can two items merge?); which primitive is consumed by which clause; per-problem increments (each problem should add exactly the definitions it needs and nothing else).
- **Audit record**: rounds run, findings (including any refuted claim, kept visible), repairs. Answers to the problems are deliberately NOT in the document, so it remains a clean problem set.
- Optionally, a **schematic formal-language sketch** (e.g. Lean-flavored), clearly marked as not compiling.

## Discipline rules (checkable)

R1. No item uses a symbol not declared above it. Dependency ("Uses:") lines list chain items only, exactly (no missing, no spurious).
R2. Every constant that does work is pinned. A uniqueness or characterization claim must pin its defining value or hypothesis; "some common value" style clauses are presumed false until checked (this class of overclaim is the most common refutation).
R3. Quantifier ranges must match existence regions: never quantify over a parameter family without a nonemptiness lemma for the stated range.
R4. Types must check: real vs complex, operator vs scalar; if a quantity is real only as a consequence, that consequence is a lemma, stated before it is relied on.
R5. No symbol collisions and no lookalike pairs (rename rather than overload; watch identity-vs-index-set style near-collisions).
R6. The document contains no answers to its own problems.
R7. Terms in questions pass the jargon-strip test: strip the term, restate in plain words; if the claim becomes vacuous, the term was concealing the absence of content.

## Workflow (iterate; do not self-certify)

W1. Draft the chain.
W2. Spawn an independent fresh-context auditor whose ONLY input is the document (forbid it other files and web search when self-containment is the thing under test). It must: (a) type-check every item against R1-R7; (b) verify every lemma computationally (symbolic/numeric, exact where possible), attempting refutation, not confirmation; (c) run the derivability test: solve the central problem from the document alone and report the answer and any missing inputs; (d) run deletion/collapse minimality tests; (e) return ranked findings with concrete fixes and a machine-readable open-issues count.
W3. Repair. Record every refutation honestly in the audit record (a refuted clause is a finding about the idea, and hiding it re-arms the error).
W4. Re-audit with a BRAND-NEW fresh auditor (never the same one: familiarity blinds). Iterate to a round cap (2 or 3); if not converged, say so plainly in the preamble instead of pretending a clean pass.

## Acceptance checklist (definition of done)

- [ ] A stranger given any single problem plus the chain above it can start working immediately.
- [ ] The central problem was actually solved from the document alone by an independent auditor.
- [ ] Every lemma was verified by computation; refutations are recorded, not erased.
- [ ] Every item passes deletion; no two items collapse; Uses lines are exact.
- [ ] Every convention the answers depend on is a named item with an auditing problem.
- [ ] Readings unpack for the stated audience; meta-questions are outside the chain; no answers leak.

## Minimal generic skeleton (shape only)

```
Preamble: discipline + audit status.
Primitives: <ambient math, explicit formulas for conventions>.
D1 <basic object>          + Lemma 1a <coordinates/parametrization>   + Reading.
D2 <observation/process>   + Reading.
D3 <pairing of D1 and D2>  + Lemma 3a <type certification>            + Reading.
D4 <the predicate>         + Lemma 4a <equivalent forms>              + Reading.
D5 <symmetry/relabeling>   + Reading.
D6 <the concrete witness>  + Lemmas: certification (pinned constants),
                             interpretation, physical identification  + Reading.
D7 <THE OBJECT> := {x in D1 | predicate(D3(x, D6))}                   + Reading.
P1 Closed form of D7, and its orbit under D5.        [the derivability test]
D8 <the witness's full family> + nonemptiness lemma  + Reading.
P2 D7's dependence on the family member; what survives all members.
D9 <dynamics>              + Reading.
P3 Which dynamics preserve D7 (deliverable format stated).
D10, D11 <rival notions>   + bridge lemmas           + Readings.
P4 Relations between D7 and the rivals.
M1 worth; M2 priority.     [outside the chain]
Minimality record. Audit record. Optional schematic sketch.
```

Apply this method to the topic you are given. Output the complete document in Markdown with all mathematics in delimited TeX ($...$ inline, $$...$$ display).
