# Backbone Notation Key

The reference card for reading a backbone page and its dependency graph. A backbone is the statement layer of an essay: a graph of assumptions, choices, and derived statements, given in words and formulas only, with the essay (prose, code, figures) built on it downstream.

## Node kinds

| Label | Name | What it is | Example |
|---|---|---|---|
| R &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#E7EEF9" stroke="#5B7DB1" stroke-width="1.5"/></svg> | root | A setup assumption: a physical fact posited about the arrangement, with its scope in one line. | R2 (cat): the environment is at zero temperature, it absorbs and never returns. |
| C &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#FDF2DA" stroke="#C09A4A" stroke-width="1.5"/></svg> | choice | A meaning we fixed rather than a fact; another choice was possible, and the page compares them. | C1 (cat): coherence means the cross-term weight ratio. |
| S0 &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#E9F4EA" stroke="#5E9C6B" stroke-width="1.5"/><rect x="4" y="4" width="22" height="10" rx="2" fill="none" stroke="#5E9C6B" stroke-width="0.8"/></svg> | equation | The equation of motion, derived from the roots by standard theory and cited; never itself a root. Each root should be visible in the equation's shape. | S0 (cat): dρ/dt = γ D[a]ρ. |
| S &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#F4F4F4" stroke="#909090" stroke-width="1.5"/></svg> | statement | A derived claim; its header names its parents. | S1 (cat): the dyadic law, the overlap raised to the lost fraction. |
| H &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="2" y="2" width="26" height="14" rx="4" fill="#FFFFFF" stroke="#333333" stroke-width="2.5"/></svg> | headline | The terminal statement, the essay's core claim; its gloss states the mechanism, not just the formula. | H (cat): Γ_dec/γ = d²/2, interference dies before the light dims. |
| N &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#FAFAFA" stroke="#999999" stroke-width="1.5" stroke-dasharray="4 3"/></svg> | recorded non-assumption | A premise deliberately absent, recorded with the reason it is not needed and the boundary of the deletion. | N (cat): a free drive commutes with the loss; nothing rests on stillness. |
| ′ (prime) | weakened variant | The same statement under one weakened root. | S1′ (cat): coherence under a bath with memory: flat onset, revival, d² survives. |

## Edge and border conventions

| Mark | Meaning |
|---|---|
| solid arrow &nbsp;<svg width="58" height="14" style="vertical-align:-2px"><line x1="2" y1="7" x2="44" y2="7" stroke="#666666" stroke-width="1.8"/><path d="M 44 2.5 L 53 7 L 44 11.5 z" fill="#666666"/></svg> | Derives from: the target consumes the source. |
| dashed arrow &nbsp;<svg width="58" height="14" style="vertical-align:-2px"><line x1="2" y1="7" x2="44" y2="7" stroke="#A6699C" stroke-width="1.8" stroke-dasharray="6 4"/><path d="M 44 2.5 L 53 7 L 44 11.5 z" fill="#A6699C"/></svg> | Weakens: the target states what holds when the source root is weakened. |
| double border &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#E9F4EA" stroke="#5E9C6B" stroke-width="1.5"/><rect x="4" y="4" width="22" height="10" rx="2" fill="none" stroke="#5E9C6B" stroke-width="0.8"/></svg> | Cited: the node's derivation is imported from standard theory, not reproved on the page. |
| dotted border, disconnected &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#FAFAFA" stroke="#999999" stroke-width="1.5" stroke-dasharray="4 3"/></svg> | Recorded non-assumption: nothing rests on it; the detachment is the content. |
| italic tail of a node | Where the truth comes from: *derived in place* (the algebra is in the node), *cited* (a pinned source), *decided* (a decision-procedure verdict on exact input, call recorded), *checked* (verified downstream along a route independent of the one that found it), or *conjectured* (unverified; propagates to every descendant and surfaces in the headline's audit line). |

## Rules of the graph

1. The equation of motion is never a root. It enters as S0, translated from the setup roots by standard theory and cited, because the physically meaningful weakenings are acts on the setup, and only the setup layer makes their effect on the equation canonical.
2. Every root must be consumed by at least one statement that is not its own weakening. A root that fails this test is not an assumption of the theory but only a direction along which the model can be varied.
3. Roots must be pairwise separable by corner models: for independent roots A and B there is a model with A and not B, and one with B and not A. Every weakening node doubles as one of these independence proofs.
4. Every statement carries its truth source (the italic tail), one of five marks. *Derived in place*: the algebra is in the node, checkable by eye. *Cited*: a pinned source. *Decided*: a decision procedure or canonical form returned True on exact input, with the exact call recorded; a verdict, not a certificate. *Checked*: verified computationally downstream, naming what the checking route shares with the finding route; route independence is this method's substitute for a proof checker, since the structural claims physics essays make have none. *Conjectured*: stated and believed, unverified, and never wearing one of the other marks.
5. The page is the source of truth; the drawn graph is a rendering of its from-lists and should eventually be generated from the page mechanically, so the two cannot drift apart.
6. The abstract is the graph linearized: read the node glosses in layer order, from the setup through the headline to the weakenings, and the dependency notation drops away.
7. The headline ends with an audit line: the complete transitive set of roots and choices it rests on, read off the from-lists, with any conjecture anywhere in that closure surfaced by name. Conjecture status propagates to every descendant; a clean audit line and one carrying a conjecture describe different objects, and the difference must be visible at the top.
8. While drafting, any candidate statement that falls inside a decidable fragment (an algebraic identity, a real inequality over finitely many parameters) must earn its *decided* verdict before it becomes a node, so wrong candidates die at the door instead of entering the graph.

## Current instances

| Page | Graph | Checked downstream in |
|---|---|---|
| backbone-cavity-cat.md | backbone-cavity-cat-graph.svg | two-examples-from-minimal-assumptions.md |
| backbone-rapid-purification.md | backbone-rapid-purification-graph.svg | two-examples-from-minimal-assumptions.md |
