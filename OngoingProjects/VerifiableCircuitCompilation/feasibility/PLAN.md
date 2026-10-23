# Feasibility graph

Scope: compare three architectures before selecting a production implementation.
The incumbent and its full-verifier results are immutable inputs.

| Node | Artifact | Exit signal |
|---|---|---|
| Shared audit | documentation ledger, stage counts, exact contract | AUDITED |
| A: seed + Newton | dyadic linear/quadratic seeds, continuum certificates, resource probes | CERTIFIED_ESTIMATE or REJECT |
| B: direct phase approximation | dyadic quadratic/cubic polynomials, certificates, resource probes | CERTIFIED_ESTIMATE or REJECT |
| C: phase propagation | exact parity identities, retained partial-product phases, compiled resource probes | EXACT_BLOCK or REJECT |
| Correctness | independent kernel reproduction | OPEN ISSUES: n |
| Quality | independent seven-axis review | OPEN ISSUES: n |
| Decision | evidence-labelled comparison | RECOMMENDATION or INCONCLUSIVE |

Seeds: experiment.py:94 varies only five reference widths; verify-qf.wls:33
checks the Newton identity, without deriving an alternative architecture.
The supplied generator already implements truncated partial products and HMR uncomputation.

Shared audit fans out to A/B/C. Failed inequalities return to approximation;
excessive cost returns to architecture. A/B have two specified architecture probes
each, rather than an unbounded search. Each symbolic query has a timeout and an
unresolved query never counts as a certificate. Correctness precedes quality;
executable quality changes return to correctness. Each review uses at most three
fresh rounds; at most three correctness visits. No new full benchmark submission
is chosen or represented as verified by this feasibility exercise.

Terminal class: code plus research report; no generated notebook is required.
