# Backbone Notation Key

The reference card for reading a backbone page and its dependency graph. A backbone is the statement layer of an essay: a graph of assumptions, choices, and derived statements, given in words and formulas only, with the essay (prose, code, figures) built on it downstream.

## Node kinds

| Label | Name | What it is | Example |
|---|---|---|---|
| R &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#E7EEF9" stroke="#5B7DB1" stroke-width="1.5"/></svg> | root | A setup assumption: a physical fact posited about the arrangement, with its scope in one line. | R2 (cat): the environment is at zero temperature, it absorbs and never returns. |
| C &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#FDF2DA" stroke="#C09A4A" stroke-width="1.5"/></svg> | choice | A meaning we fixed rather than a fact; another choice was possible, and the page compares them. | C1 (cat): coherence means the cross-term weight ratio. |
| S0 &nbsp;<svg width="30" height="18" style="vertical-align:-4px"><rect x="1" y="1" width="28" height="16" rx="4" fill="#E9F4EA" stroke="#5E9C6B" stroke-width="1.5"/><rect x="4" y="4" width="22" height="10" rx="2" fill="none" stroke="#5E9C6B" stroke-width="0.8"/></svg> | equation | The equation of motion, derived from the roots by standard theory and cited; never itself a root. Each root should be visible in the equation's shape. | S0 (cat): dρ/dt = γ D[a]ρ. |
| E | contract | The same equation node when posited without its construction: an equation-level axiom with a precise scope statement, legitimate only when no planned weakening reaches beneath it (rule 1). | No instance yet; both current pages expose their construction, so their equation node is S0. |
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
| italic tail of a node | Provenance, then verification. Provenance: *derived in place*, *cited* (pinned), or *conjectured* (propagates to every descendant and surfaces in the audit line). Verification, per clause: *exact decision* (call recorded), *analytic cross-check*, *numerical check* (single-point checks say so), or *none*; the tail names the clause scope and what the checking route shares with the finding route. Full test manifests live beside the checks in the companion. |

## Rules of the graph

1. The page has two layers. In the construction layer, physical roots R are translated into the equation of motion S0 by cited standard theory. In the analysis layer, every consequence consumes the equation as a contract, and therefore survives any construction that yields the same equation. The equation is never a root by default; a page may posit the contract alone (label the node E, with a precise scope statement, construction declared out of scope) exactly when every planned weakening is expressible at the contract level. The moment one weakening reaches beneath the equation (a lost record share, an opened channel, a warmed bath), the construction layer is owed: the weakening set dictates the layer, never taste.
2. Every root must be consumed by at least one statement that is not its own weakening. A root that fails this test is not an assumption of the theory but only a direction along which the model can be varied.
3. Roots must be pairwise separable by corner models: for independent roots A and B there is a model with A and not B, and one with B and not A. Every weakening node doubles as one of these independence proofs.
4. Every statement's tail states its provenance, one of three: *derived in place* (the algebra is in the node, checkable by eye), *cited* (a pinned source, never a naked "standard"), or *conjectured* (stated and believed, unverified in origin). Provenance *conjectured* propagates to every descendant regardless of verification, because a numerically checked conjecture is still a conjecture.
5. The same tail states its verification, per clause: *exact decision* (a decision procedure or canonical form returned True on exact input, call recorded), *analytic cross-check* (an independent closed form or derivation compared), *numerical check* (computation downstream; single-point checks say so, naming what was varied and what was held fixed), or *none*. The tail names the clause scope (whole node, or which clause) and what the checking route shares with the finding route; route independence is the method's substitute for a proof checker, since the structural claims physics makes have none. The full manifest of each check (the exact test, what varied, what was held fixed, corner cases, clause scope) lives beside the check in the downstream companion; the page carries only the compressed form. This applies forward; earlier companions are retrofitted when next touched.
6. The page is the source of truth; the drawn graph is a rendering of its from-lists and should eventually be generated from the page mechanically, so the two cannot drift apart.
7. The abstract is the graph linearized: read the node glosses in layer order, from the setup through the headline to the weakenings, and the dependency notation drops away.
8. The headline ends with a two-stage audit line: first the contract and the analysis-layer premises it rests on, then the roots realizing the contract by the cited translation (or the scope statement of a posited E), with any conjecture anywhere in the closure surfaced by name. A clean audit line and one carrying a conjecture describe different objects, and the difference must be visible at the top.
9. While drafting, any candidate statement that falls inside a decidable fragment (an algebraic identity, a real inequality over finitely many parameters) must earn its *exact decision* verdict before it becomes a node, so wrong candidates die at the door instead of entering the graph.

## Current instances

All backbone pages and graphs live in this Backbone folder; companions live beside the essays they check, one level up unless a fuller path is given.

| Page | Graph | Checked downstream in |
|---|---|---|
| backbone-cavity-cat.md | backbone-cavity-cat-graph.svg | ../two-examples-from-minimal-assumptions.md |
| backbone-rapid-purification.md | backbone-rapid-purification-graph.svg | ../two-examples-from-minimal-assumptions.md |
| backbone-feedback-cooling.md | (none yet) | ../open-system-simulation-catalog.md (Feedback Cooling section) |
| backbone-catalog-spine.md | backbone-catalog-spine-graph.svg | ../open-system-simulation-catalog.md (toolkit and all twenty examples) |
| backbone-entanglement-sudden-death.md | backbone-entanglement-sudden-death-graph.svg | backbone-entanglement-sudden-death-checks.wl (exact-decision battery; no companion essay yet) |
| ../../ContinuousMeasurement/Continuous Space/Caldeira-Leggett/backbone-gaussian-bath-exact-tcl2.md | backbone-gaussian-bath-exact-tcl2-graph.svg (same folder as its page) | ../../ContinuousMeasurement/Continuous Space/Caldeira-Leggett/gaussian-bath-exact-tcl2-conditions.md (the source draft) |

Tooling: `backbone_lint.py`, in this folder, derives each page's graph from its node headers and checks the mechanical rules (root consumption, header-versus-tail consistency, audit-line closure and conjecture surfacing, tail grammar, symbol-definition closure as warnings), with declared deviations downgraded to notices; its `--render` flag regenerates the page's graph SVG from the derived structure, so graphs are generated, never hand-laid, and cannot drift from their pages. Run it before delivering a page; the fixtures in `lint-fixtures/` reproduce the historical defects it must always catch.

A note on labels across layers: the companion essay predates this key and keeps its own local labels (A1 through A5 with C1 in its part one, B1 through B4 with C2 and C3 in its part two); this key governs backbone pages only. Where the mapping is not obvious: the essay's A1 (no memory) is the cat page's R3 and its A2 the page's R2; the essay's B2 (efficiency) is the purification page's R2 and its B3 the page's R3; the essay's numerical assumptions (A5, the cutoff; B4, the step) have no backbone counterpart on purpose, being downstream-layer concerns.
