# Evidence Route Table

Every worked answer declares one primary route. Secondary routes are cross-checks and must not share
the same hidden assumption as the primary route.

| class | primary evidence | required check | silent trap | useful independent cross-check |
|---|---|---|---|---|
| Q0 | exact $\mathrm{GF}(2)$, Pauli, stabilizer, or symplectic algebra | ranks, commutators, logical anticommutation, equivalence modulo stabilizer or gauge | doing binary algebra over the reals; losing Pauli phase where phase matters | direct matrix representation for a small instance |
| Q1 | exact finite-dimensional channel or recovery calculation | CPTP conditions and action on a logical-reference entangled state | checking only $\lvert 0_{L}\rangle$ and $\lvert 1_{L}\rangle$; replacing coherent noise by a Pauli channel | Choi matrix or independent Kraus representation |
| Q2 | literal circuit propagation and exhaustive bounded faults | enumerate every claimed fault location and classify residual data and logical error | verifying the abstract code but not the extraction circuit; omitting idle, reset, measurement, or leakage faults | stabilizer propagation plus direct small-Hilbert-space simulation |
| Q3 | decoder built for a fixed syndrome model | exact lookup or brute-force maximum likelihood on a small case | treating the lowest-weight physical error as the most likely logical class; using mismatched weights | exhaustive coset summation on a reduced instance |
| Q4 | Monte Carlo, finite-size scaling, latency, or resource model | stored configuration, uncertainty, convergence, distance sweep, and failure-edge search | common random data generated differently for each decoder; declaring a crossing from too few sizes; hiding routing or factory time | exact enumeration at small distance or a second simulator |
| Q5 | proof or operational conceptual argument | explicit hypotheses, conclusion, and counterexample when a hypothesis is removed | citing a theorem name as an explanation; confusing asymptotic with finite-size statements | finite example that realizes the theorem and one that violates a dropped premise |
| Q6 | bounded reconstruction of a primary result | source data or digitization, declared metric, uncertainty, status label, and limitation | repeating an abstract or press claim; combining separately demonstrated components into one claimed experiment | independent reanalysis or an alternative decoder on the released data |

## Route assignment by Part

| Parts | dominant routes | reason |
|---|---|---|
| 0-1 | Q0, Q1, Q3 | algebra and classical inference are exact at the pinned sizes |
| 2-4 | Q1, Q5 | first codes and general recovery must connect channels to operational arguments |
| 5-7 | Q0, Q1 | stabilizer structure and finite codes admit exact checks |
| 8-9 | Q2 with Q0/Q1 cross-checks | fault propagation is a property of literal circuits |
| 10 | Q4 and Q5 | thresholds join asymptotic reasoning to finite statistical evidence |
| 11 | Q0, Q2, Q3 | topology, extraction, and syndrome graphs must agree |
| 12 | Q3 and Q4 | decoding is inference plus benchmark discipline |
| 13 | Q0, Q3, Q4, Q5 | finite constructions must remain separate from asymptotic theorems |
| 14 | Q1 and Q4 | bosonic encodings join exact channel structure to finite-energy numerics |
| 15 | Q2-Q4 | hardware adaptation couples circuits, information, and inference |
| 16 | Q6 | the task is reconstruction, not literature narration |

## Cross-check independence test

Before calling two routes independent, list their shared inputs. They are not independent when they
use the same truncated error set, calibrated noise model, decoder weights, finite box, routing
assumption, or source digitization. Agreement can validate implementation consistency while leaving
the common premise untested; the prose must say so.
