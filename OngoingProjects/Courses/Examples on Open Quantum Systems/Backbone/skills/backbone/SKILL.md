---
name: backbone
description: Write or revise a backbone, the code-free statement layer of a physics essay: a one-page dependency graph of setup assumptions (R), choices (C), the equation of motion as a derived node (S0) or posited contract (E), derived statements (S), a headline with a two-stage audit line (H), and recorded non-assumptions (N), with weakened-assumption statements and an abstract read off the graph. Use whenever Mads asks for a backbone, a statement layer, an assumption graph, a minimal-assumption or foundations-first treatment, the DNA method, an axiomatic skeleton, or "what does this essay really assume", for an essay or a research piece, in any repository; covers both extracting a backbone from an existing essay and writing one from scratch for a topic with no written case. Do not fire for ordinary computational essays without such a request (learning-by-computing owns those); this treatment is opt-in.
---

# Backbone: the statement layer of a physics essay

A backbone is the upstream artifact an essay is later built on: a single page holding the graph of assumptions, choices, and derived statements, each given in words and formulas only, nothing executable. Reading it alone conveys the core idea; the essay (prose, code, figures) is downstream and supplies the computational checks. The failure mode this method exists to kill is writing from a generic prior: plausible, fluent, false statements. Everything below is arranged against that.

## Read these before writing anything

Generate against the exemplars, never from memory of the rules. The canonical files (this machine):

- **The rules, authoritative (nine rules, v3)**: `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/Backbone/backbone-notation-key.md`. On any conflict with this SKILL.md, the key wins.
- **Exemplar pages**: `backbone-cavity-cat.md` and `backbone-rapid-purification.md` in that same Backbone folder; their drawn graphs are `backbone-cavity-cat-graph.svg` and `backbone-rapid-purification-graph.svg`.
- **Downstream companion** (where verifications live): `two-examples-from-minimal-assumptions.md`, one folder up, beside the catalog. It predates the key and keeps its own local labels; the key's closing note gives the mapping.
- Background on why the method is shaped this way (optional): `/Users/mohammadb/Documents/GitHub/TheoremProving/lean-vs-wolfram-theorem-proving.md`.

Read the key and at least one exemplar page in full before drafting. If the paths are missing (another machine), ask Mads for the current location rather than reconstructing from memory.

## The two layers

Every page separates model construction from analysis. In the **construction layer**, physical roots R (facts about the arrangement; existence claims and closure claims are different kinds and say so) are translated into the equation of motion S0 by cited standard theory, with each root visible in the equation's shape. In the **analysis layer**, every consequence consumes the equation as a contract, so the analysis survives any construction that yields the same equation. A page may instead posit the contract alone, labeled E, with a precise scope statement and construction declared out of scope, but only when every planned weakening is expressible at the contract level; the moment one weakening reaches beneath the equation (a lost record share, an opened channel, a warmed bath), the construction layer is owed. The weakening set dictates the layer, never taste.

## The shape of a page

Layer order, each node with a bold header naming its label, kind, and parents ("from:"), then words and formulas, then an italic tail:

1. Preamble: what the page is, the conventions, the label key, where the checks live.
2. **The setup, and the choice**: roots R1..Rn (one line of scope each), choices C1.. (meanings we fixed; another choice was possible and is compared).
3. **The equation**: S0 (from roots; cited) or E (posited; scope statement).
4. **Consequences**: statements S1.., each naming its parents in the header.
5. **The headline** H: the core claim; the gloss states the mechanism, not just the formula; ends with the two-stage audit line (the contract and analysis-layer premises, then the roots realizing the contract or E's scope statement, with any conjecture in the closure surfaced by name).
6. **When an assumption is weakened**: primed statements, each doubling as an independence proof; recorded non-assumptions N with the reason and the boundary of the deletion.
7. **The abstract, read off the graph**: the node glosses linearized in layer order; nothing added.

Budget: 6 to 12 nodes, one page. Past that it has become the essay again.

## The tails: provenance, then verification

Every statement's italic tail answers two separate questions.

**Provenance** (where the statement came from), one of three:
- *Derived in place*: the algebra is in the node, checkable by eye.
- *Cited*: a pinned source (book with chapter, or paper with journal or arXiv id). Never a naked "standard".
- *Conjectured*: stated and believed, unverified in origin. Propagates to every descendant regardless of verification (a numerically checked conjecture is still a conjecture) and surfaces in the headline's audit line.

**Verification** (how it was tested), per clause:
- *Exact decision*: a decision procedure or canonical form returned True on exact input (Reduce, Resolve, FullSimplify, canonical forms), with the exact call recorded. A verdict, not a certificate.
- *Analytic cross-check*: an independent closed form or derivation compared.
- *Numerical check*: computation downstream. Single-point checks say so, naming what was varied and what was held fixed ("at a single efficiency, other values untested").
- *None*.

The tail names the clause scope (whole node, or which clause got which check) and what the checking route shares with the finding route; route independence is the method's substitute for a proof checker, since the structural claims physics makes have none. The tail is the compressed form only: the full manifest of each check (the exact test, what varied, what was held fixed, corner cases, clause scope) lives beside the check in the downstream companion.

## The two workflows

**Extraction (an essay exists).** Compress it: find the headline, read the assumptions off what the existing derivations and computations actually consumed, type the choices honestly, wire the weakening nodes to the essay's existing turning-off computations, and mark provenance and verification from where each claim's truth already lives.

**From scratch (no written case).** The order matters, because it is what prevents assumption padding and prior-recall:

1. State the candidate headline in one sentence.
2. Attempt the derivation by hand, in node-sized steps. The assumption list is read off what the derivation consumed; never write the list first. An assumption nothing consumed gets cut; a step that fails exposes a missing one.
3. Define the words before the list uses them; nested terms (thermal inside Gaussian) get a separating example; ambiguous adverbs ("faster", "purifies") get their inequivalent meanings split.
4. Any candidate statement inside a decidable fragment (an algebraic identity, a real inequality over finitely many parameters) must earn its *exact decision* verdict via wolframscript on exact input before it becomes a node. Wrong candidates die at the door.
5. What is neither derivable in place, nor pinnable, nor decidable gets provenance *conjectured*, visibly; a scratch computation may still give it a numerical check, and the pair (conjectured, numerically checked) is an honest and expected state for a young node. Scratch computations in a throwaway kernel during drafting are encouraged (normal careful work, not a gated loop); the artifact stays code-free.
6. The essay later written against the backbone upgrades provenance where derivations land, supplies the full manifests, and may send corrections upward; the backbone is a versioned hypothesis, not scripture.

## The linter, and the graph

Before delivering any page, run the deterministic checker and fix what it finds:

```
python3 "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/Backbone/backbone_lint.py" PAGE.md
```

It derives the graph from the node headers and diffs the authored metadata against it: root consumption, header-versus-tail consistency, audit-line closure and conjecture surfacing, tail grammar (the provenance and verification kinds are fixed typed vocabularies the linter parses, not free prose), symbol-definition closure (a symbol defined with an equals sign in one node and used in another demands an edge path), and, where the page declares its companion (`<!-- companion: path -->` in the preamble), a reachable manifest anchor (`<!-- manifest: S1 -->` in the companion) for every tail claiming a verification. Its `--status` flag prints each node's provenance, verification kinds, conjecture taint, and manifest state. Errors must be fixed; warnings must be fixed or explicitly justified to Mads; a deviation from a key rule may be declared in the preamble ("Declared deviation (rule N): reason") and then surfaces as a notice instead of an error, loud but legal. The reason this exists: from-lists written by the author drift from the derivation the author wrote, and the three historical defects that prompted the linter were all of that kind.

Graphs are **generated, never hand-laid**:

```
python3 .../backbone_lint.py PAGE.md --render PAGE-graph.svg
```

The generator owns the layout and the typography, which retires both the page-versus-drawing drift class and the SVG formula-mangling class. Each node's box text comes from a display caption written in the page as an invisible comment inside the node block: `<!-- box: first line | second line -->`, a curated summary of the node (its key inequality, rate, or verdict), not the gloss's opening words; use `_{...}` and `^{...}` for letter indices (the generator emits correct tspans), Greek and digit sub/superscripts as bare glyphs, about 40 characters per line. A node without a caption falls back to its gloss's opening words, which is legible but far less useful, so caption every node. After rendering, open the SVG in the browser pane and screenshot it once to confirm the render is sane. If you must hand-edit an SVG anyway (legacy graphs only): no LaTeX markup or combining accents in SVG text; Greek letters, ⊗, ⟨⟩, ≥, ′ and digit sub/superscripts (², ⁻¹, ₂) are safe bare glyphs; modifier letters (ᵐ, ₜ) are banned; letter indices need tspan pairs with dy negated afterward.

## Guardrails

- Never the em dash character, anywhere.
- No coined nouns; standard names or literal descriptive phrases.
- Prose glosses interpret the mechanism; no computed numbers hard-coded in text (fixed inputs you chose are fine).
- Self-audit before delivering: every root consumed by a non-weakening statement; roots pairwise separable by corner models; every tail carrying exactly one provenance and an honest per-clause verification; the audit line's closure correct by tracing the from-lists and naming the layer boundary.
- Adversarial review gates (inquiry-critic and kin) are ask-first, always: propose, never launch.
