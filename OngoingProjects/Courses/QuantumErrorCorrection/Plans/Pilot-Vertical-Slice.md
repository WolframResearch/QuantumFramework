# Pilot Vertical Slice

These fifteen answers are authored before the full collection. They test whether the question
contract, tooling, notation, and verification gates work across the whole causal chain. They remain
in prerequisite order even though they sample distant Parts.

| order | question | purpose | minimum acceptance evidence |
|---:|---:|---|---|
| 1 | 0.3 | establish symplectic commutation | exhaustive agreement with direct matrices for all two-qubit Paulis |
| 2 | 1.1 | establish classical failure accounting | exact polynomial and Monte Carlo agreement within uncertainty |
| 3 | 2.2 | establish syndrome without logical-information leakage | symbolic $\alpha,\beta$ calculation and full single-$X$ table |
| 4 | 2.5 | establish error discretization | arbitrary one-qubit operator expanded and recovered by linearity |
| 5 | 3.4 | expose coherent-versus-stochastic noise | matched average infidelity with different repeated/logical scaling |
| 6 | 4.1 | establish general correctability | independent derivation of Knill-Laflamme plus a failed counterexample |
| 7 | 5.4 | establish syndrome as coset information | complete small-code table with logical-equivalence classification |
| 8 | 7.2 | establish CSS construction | Steane ranks, commutation, logical Paulis, and $[[7,1,3]]$ verified |
| 9 | 8.2 | establish circuit-level fault propagation | every single fault propagated through two gate orderings |
| 10 | 10.5 | establish honest statistical scaling | synthetic known-threshold dataset recovered with uncertainty and finite-size drift shown |
| 11 | 11.6 | establish a topological patch from algebra | check matrix, qubit count, logical strings, and distance for two sizes |
| 12 | 12.3 | establish decoder implementation | matching result checked against exhaustive logical-class likelihood at small distance |
| 13 | 13.8 | establish finite qLDPC construction | literal matrices, commutation, ranks, and independently verified parameters |
| 14 | 14.6 | establish finite-energy bosonic honesty | ideal limit plus a finite-squeezing logical-error calculation |
| 15 | 16.1 | establish frontier reconstruction | source-data fit, uncertainty, decoder label, and explicit nonclaim |

## Pilot release gate

- All literal code runs from a clean kernel or environment.
- Every question uses its declared evidence route and records tool/library versions.
- No two answers introduce conflicting Pauli, syndrome, or rate conventions.
- Every answer includes a computation or argument capable of refuting its central claim.
- The fifteen answers receive an independent physics pass before Part-wide authoring begins.
- Any question that cannot fit the answer rhythm without becoming a chapter is split in the main
  list before the stable core is frozen.
