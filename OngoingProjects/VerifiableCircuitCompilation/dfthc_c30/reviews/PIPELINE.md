# Review pipeline

Deliverable class: code, with no notebook build. The wl-ship graph completed
two correctness visits (one round each) and two quality rounds, with one
quality-to-correctness back edge. Every review used a fresh independent agent.
The first quality round found redundant MapThread around Listable BitXor.
The direct array operation replaced it; the revised source returned to
correctness and then to a fresh quality review.

Correctness: OPEN ISSUES: 0. The independent reviewer rebuilt identical
circuits, factors and transformed words, reran the full official verifier,
reproduced all WL and dirty-target probes, and added boundary cases through
4096 rows. See correctness-1.md and correctness-v2-1.md; the final snapshot
is bound in correctness-v2-1-sha256.json.

Quality: OPEN ISSUES: 0. The independent critic confirmed all axes A–G,
including documented symbols, functional WL, physics invariants, nontrivial
structure, exact symbolic arguments and required regimes. See quality-2.md.
The historical quality-1.md retains its original one-issue verdict.

The final evidence manifest binds the public sources, results and reviews.
The preceding project audit was pushed as df53dba; the C30 loader stage
remains local. No factors were trained and no Kaggle submission was made.
