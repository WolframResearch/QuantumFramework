---
name: innovation
description: Farm candidate claims for novelty testing on a physics topic by procedure rather than by asking the model to be creative: enumerate variations with predictions stated in advance, arbitrage inequivalent meanings of load-bearing words, transplant structures from distant formal systems through an explicit dictionary, mutate the problem's representation and test the mutation by compression, and map the literature onto a backbone's dependency graph to kill rediscoveries and surface unmapped cells; every candidate faces the strongest affordable falsifier, and Mads receives the survivors, the instructive deaths, and the tradeoffs. Use whenever Mads asks for novel angles, new questions, innovation, unasked questions, research directions, "what could be new here", "what has nobody done", or candidate claims for an essay, backbone, or research topic. Operates on a backbone page (the `backbone` skill's artifact); if none exists, build one first. Do not use for ordinary drafting, verification of known content, or literature summary.
---

# Innovation: farming candidate claims by procedure

The premise, stated honestly: a language model generates from the distribution of what has been written, so its best spontaneous ideas are the field's existing ideas, and asking it to "be creative" produces novelty-shaped text whose fluency masks that it is known, wrong, or vague. Scientific novelty has many sources this skill does not reach (new instruments, new data, impossibility results, technological change); this workflow concentrates on the sources a model can systematically exploit: mechanical variation with priced expectations, collisions whose failed predictions are kept as first-class objects, structure transplanted across distant fields, representation changes tested by compression, and retrieval used honestly. The model's role is volume and breadth; taste stays with Mads, and the loop is engineered to maximize well-tested candidates per hour of his attention, not to simulate his judgment.

One distinction governs everything here: this skill produces **candidates for novelty testing**, never certified novelty. A claim has three independent statuses, and collapsing them is the failure mode this revision exists to prevent:

```
Validity:      untested | refuted | supported | derived
Prior art:     unchecked | known (with citation) | no close prior found (protocol recorded)
Significance:  trivial | incremental | substantial | undecided
```

"Refuted" and "known" are different deaths: one means the physical reasoning failed, the other that it succeeded and arrived somewhere already occupied. Both stay on the ledger; a rediscovery with a citation is evidence about the search's calibration, not an embarrassment to hide.

## Prerequisite: an address space

This skill operates on a backbone page (see the `backbone` skill; rules at `/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/Examples on Open Quantum Systems/Backbone/backbone-notation-key.md`). The backbone gives the topic a finite address space of nodes and edges, which turns "find an unasked question" into an algorithm, and its known limitation is faced by Move 5 below: local coordinates cannot reveal what lies outside their own representation. If no backbone exists, build one first, or ask Mads whether a rough one should be drafted as scaffolding.

## The candidate record

Fields at birth, for every candidate:

```
Claim:                     one falsifiable sentence; if it cannot fail, it is not a candidate
Move:                      which move produced it
Scope:                     regime, state family, system class it speaks about
Expectation before testing: what the first collision is predicted to show, on record BEFORE anything runs
Reason it might hold:      the physical or structural argument, one or two sentences
Cheapest falsifier:        the least costly test that could kill it, with a cost estimate
Validity / Prior art / Significance:  the three statuses above
```

Fields added when a candidate survives its collisions and a prior-art pass:

```
Closest prior result:      the nearest published thing, cited
Exact delta:               precisely what this claim adds beyond it
Why the delta matters:     conceptual or practical significance, stated clinically
Novelty-search protocol:   databases, queries, terminology variants, citation trails, date; and a confidence
Remaining falsifier:       the strongest test not yet run
Next test:                 what would be run next, at what cost
Decision:                  retain | reject | rediscovery | defer   (Mads's call)
```

The expectation field is the surprise detector: the model under-recognizes its own surprises, so the prediction must be on record for a miss to be mechanical to spot, and a corrected wrong expectation is historically where the best statements have come from. The ledger is the evidence record and may be as detailed as it needs to be; compression is the backbone's discipline, not the ledger's.

## The five moves

**Move 1: variation sweep with priced expectations.** Enumerate what the backbone makes mechanical: each root weakened in each physically meaningful direction, each pair's corner models, regimes and crossings the page does not price (long times, strong coupling, low efficiency, and their intersections, where idealized analyses stop because their idealization does). Expectation first, then the test.

**Move 2: literature mapping on the graph.** Runs AFTER the generative moves and the cheap kills (see the staging below), in a dual role: attach every published question found to the backbone node or edge it interrogates, killing rediscoveries with their citations, and surface the cells left over. A cell with nothing attached is **unmapped**, never "unaddressed": the search may have missed terminology, the result may live as a lemma inside a paper about something else, the databases have edges, or the cell is a corollary nobody thought worth stating. The honest product is "no close prior found under the recorded protocol", with the protocol written down. Retrieval hierarchy: the Undermind connector where available (call get_orientation first; deep searches take minutes and are worth it); otherwise arXiv, Crossref, and Semantic Scholar with terminology variants; then citation chaining from the closest papers found; and when coverage is still thin, a recorded limitation instead of a confident absence.

**Move 3: meaning-split arbitrage.** The backbone's definitional discipline splits load-bearing words into inequivalent meanings (deadline versus race, pathwise versus mean versus passage). For each split, determine which meaning the published claims silently assume; the other meaning's version of the claim is a candidate, and this move-type has produced real papers.

**Move 4: transplant through a structural dictionary.** Name the two or three genuinely adjacent formal systems and push a central invariant or discipline through an explicit dictionary. Every transplant fills, before it becomes a candidate:

```
Source objects:            Source transformations:      Source invariant or theorem:
Dictionary into target:    Assumptions preserved:       Assumptions lost:
New testable consequence:
```

A dictionary that maps only vocabulary and yields no falsifiable consequence dies immediately; that filter is what separates transplant from metaphor.

**Move 5: mutate the address space.** The complement of a graph finds only holes in the graph; this move searches for a better graph. Candidate mutations: replace the state variable or the objective; move from ensemble dynamics to trajectories, or from time evolution to a stopping-time question; hunt for a monotone or an invariant; quotient out degrees of freedom the headline never uses; pose the dual optimization; exchange operational and geometric descriptions; ask which apparently distinct nodes collapse under a new quantity. A representation candidate has a collision test like any other, and it is **compression**: it survives only if it re-derives the same headline from fewer roots, merges nodes the old coordinates kept apart, or turns a numerically checked tail into a derived-in-place one. Compression is measurable on the graph; a representation that does not compress is a change of notation.

## Staging and collisions

1. Run Moves 1, 3, 4, 5 before any literature contact, so generation is not anchored to what retrieval surfaces.
2. Apply the cheap kills to everything: dimensional analysis, limiting cases, algebraic checks; claims inside decidable fragments get exact decisions through wolframscript on exact input, call recorded. Mads's standing computation rules apply (they load with every session): wolframscript -file, no Quiet, no Print, throwaway kernels, code discarded.
3. Run Move 2 once, in its dual role.
4. Compute or derive the survivors, at whatever depth each claim warrants; estimate cost per candidate rather than by test type, since a five-minute search can save an hour of simulation and a one-line counterexample can save both.
5. Deep literature and citation-chain searches on serious survivors only; fill their survivor fields.
6. Independent verification of finalists (a route sharing nothing with the route that found them).

The requirement is coverage, not mortality: every candidate must face the strongest affordable falsifier relevant to its claim. A pass where everything survives is a diagnostic to investigate, not a quota to correct; and manufacturing harshness (weak candidates bred to die, minor qualifications booked as deaths) is the mirror-image failure and equally disqualifying.

## Selection and the interface back to the backbone

Present Mads the survivors and the instructive deaths in one table, with the tradeoffs visible per candidate: expected correctness, novelty confidence, significance, testing cost, scope, and relevance to the current essay; headline-impact is one column, not the ranking, since a claim can leave the headline intact and still matter more than one that overturns it. Mads decides; the skill never launches adversarial review loops on its own (inquiry-critic is the later gate for a chosen candidate, ask-first, always). A retained survivor enters the backbone as a node with provenance *conjectured* and whatever verification its collisions earned; the ledger is saved beside the backbone as `<backbone-name>-candidates.md` so later passes inherit the kill list and the search protocols.

## Boundary, stated plainly

This skill does not promise representation-level innovation; Move 5 attempts it, with compression as the bar. Any candidate of that kind receives the strongest correctness, prior-art, and significance scrutiny available, and the verdict comes from that scrutiny, not from the fact that a model proposed it: source-based disbelief would violate the collision philosophy this skill runs on. What can be promised is the reliable product: unpriced corners, unasked-so-far-as-the-recorded-search-shows questions, meaning arbitrage, disciplined transplants, and an honest ledger of what died and why.

## Guardrails

- Never the em dash character; no coined nouns; claims in standard vocabulary or literal descriptive phrases.
- No candidate without its birth fields; no collision without its prior expectation; no silent smoothing of a miss.
- Every field clinical, the three status axes included; excitement adjectives ("groundbreaking", "remarkably") are banned everywhere in the ledger, and enthusiasm is expressed by what a claim survived, cited.
