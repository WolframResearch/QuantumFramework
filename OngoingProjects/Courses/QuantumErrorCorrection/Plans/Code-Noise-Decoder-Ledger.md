# Code-Noise-Decoder Ledger

This seed ledger pins the primary carrier and discriminating contrast for each Part before worked
answers are authored. Question-level answer plans must extend it; they may not silently swap a code,
noise model, decoder, or success observable because a computation becomes inconvenient.

| Part | primary carrier | noise or fault model | inference or recovery | observable | discriminating contrast |
|---:|---|---|---|---|---|
| 0 | one- and two-qubit Paulis and channels | one explicit Kraus channel | exact algebra | commutators, Choi spectrum, trace preservation | real arithmetic versus $\mathrm{GF}(2)$ |
| 1 | three-bit repetition and a small linear code | binary symmetric, nonuniform flip, erasure | majority, maximum likelihood, exact BP | block failure probability | flip versus located erasure |
| 2 | three-qubit repetition and Shor code | one Pauli and one coherent rotation | ideal syndrome plus recovery/frame | logical channel | stochastic Pauli versus coherent rotation |
| 3 | one qubit and one extraction circuit | Pauli, dephasing, damping, leakage, bias, correlated burst | none or fixed diagnostic | channel metrics and logical error per cycle | equal average infidelity, different worst-case or logical behavior |
| 4 | repetition code plus one erasure example | a stated correctable or approximate operator span | explicit, Petz-type, and optimized recovery | Knill-Laflamme matrix, entanglement fidelity, coherent information | exact versus approximate correction; capacity versus memory claims |
| 5 | small qubit and prime-power-qudit stabilizer codes | bounded generalized-Pauli set plus one non-Pauli error | syndrome table and normalizer | code dimension and logical class | binary tableau versus trace-symplectic representation |
| 6 | Steane, five-qubit, degenerate, and local-code examples | bounded Pauli and erasure sets | exhaustive distance search and cleaning argument | $[[n,k,d]]$, $d_{X}$, $d_{Z}$, locality trade-off | abstract distance versus circuit-level or geometric constraints |
| 7 | Shor, Steane, five-qubit, Bacon-Shor, damping code | common circuit-level model where compared | code-specific recovery | logical error and extraction cost | CSS versus non-CSS versus subsystem |
| 8 | weight-four check, Steane block, and a single-shot example | every single circuit fault, selected pairs, noisy checks | flag/cat/encoded or redundant syndrome | residual data class and malignant fault sets | repeated extraction versus single-shot; gadget test versus exRec proof |
| 9 | Steane block and surface-code patches | single faults and circuit-level magic-state noise | teleportation, cultivation, distillation, frame tracking | logical process, acceptance, space-time volume | transversal/locality-preserving versus surgery or resource-state routes |
| 10 | concatenated toy model plus synthetic code family | code-capacity, phenomenological, circuit-level, rare correlated events | fixed decoder plus declared rare-event sampler | threshold, pseudo-threshold, $\Lambda$, tail probability, overhead | suppression versus break-even; direct versus biased sampling |
| 11 | toric, rotated planar, colour, Floquet, and cluster-state instances | Pauli data, preparation, and measurement faults | syndrome or measurement-history graph | logical operators across space, time, and measurement schedule | static stabilizer versus dynamically generated code |
| 12 | distance-3/5 surface instances and one CSS qLDPC instance | identical stored event streams | lookup, matching, union-find, BP family, tensor-network, learned model | logical error, approximation error, latency, memory | exact logical-coset probability versus fast approximation |
| 13 | verified hypergraph-product, expander, and bivariate-bicycle instances | code-capacity and a declared extraction or logical-operation circuit | exact small decoder, small-set-flip, and BP-based decoder | parameters, logical error, connectivity and logical-gate cost | asymptotic family versus usable finite instance |
| 14 | cat, binomial, and finite-energy GKP encodings | loss, dephasing, displacement, and ancilla faults | parity/modular syndrome, autonomous stabilization, fault-enumerated gadget | logical channel per time, cycle, or gate | passive code-space protection versus fault-tolerant operation |
| 15 | photonic fusion, XZZX/repetition, and $[[4,2,2]]$ erasure examples | fusion failure, photon loss, biased Pauli, leakage, flagged erasure | information-aware decoder | logical gain from located side information | flags retained versus discarded; located versus residual unlocated faults |
| 16 | released data and circuits from primary sources | exactly the source's stated model or experiment | source decoder plus one bounded alternative where possible | source-specific measured quantity | directly measured result versus extrapolation |

## Collision rules

1. A code may recur when the question changes a load-bearing element: noise, syndrome circuit,
   decoder, logical operation, or performance observable.
2. Repeating the same single $X$ error on the same code is not a new example.
3. A named canonical code remains when it is the subject; add a contrast that exposes what the
   canonical example hides.
4. Comparisons use a single stored event stream or a reproducibly coupled sampling procedure.
5. Finite qLDPC parameters are computed from the literal matrices; family labels are not trusted.
6. Frontier data keep the paper's units and denominator until an explicitly documented conversion.

## Question-level extension fields

Each future answer plan adds one row with:

`question | code instance | physical parameters | noise | information | decoder/recovery | output |
refuting check | seed/data provenance | permitted replacement class`
