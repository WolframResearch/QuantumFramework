# Prerequisite Graph

The course is a directed acyclic learning graph, not merely a numbered list. The stable spine is

$$
0 \rightarrow 1 \rightarrow 2 \rightarrow 3 \rightarrow 4 \rightarrow 5 \rightarrow 6
\rightarrow 7 \rightarrow 8 \rightarrow 9 \rightarrow 10,
$$

with stabilizer methods from Part 5 also feeding the topological and qLDPC branches, and channel
methods from Parts 3 and 4 feeding the bosonic branch.

## Part-level dependencies

| Part | requires | new capability unlocked |
|---:|---|---|
| 0 | finite-dimensional quantum mechanics | binary Pauli algebra, channels, explicit measurement records |
| 1 | 0.1 and elementary probability | syndrome decoding as inference over cosets |
| 2 | 0 and 1 | quantum redundancy, error discretization, first logical channel |
| 3 | 0.5 and 2 | honest physical and circuit-level noise models; operational performance metrics |
| 4 | 2 and 3 | exact and approximate correctability, optimized recovery, channel-rate distinctions |
| 5 | 0.1-0.4 and 2 | qubit and prime-power-qudit stabilizers, syndromes, normalizers, logical Paulis |
| 6 | 4 and 5 | parameters, distance, erasure capability, general and locality-constrained bounds |
| 7 | 3-6 | canonical CSS, non-CSS, subsystem, and channel-adapted codes |
| 8 | 3, 5, and 7 | repeated and single-shot extraction; gadget-level and exRec fault proofs |
| 9 | 5, 7, and 8 | protected logical gates, locality restrictions, and universality resources |
| 10 | 3, 8, and 9 | threshold language, finite-size and rare-event estimation, resource accounting |
| 11 | 5, 6, and 8 | static and Floquet topological codes, cluster states, homology, lattice surgery |
| 12 | 3, 5, 10, and 11 | exact logical-coset and approximate decoders, uncertainty, latency, fair benchmarks |
| 13 | 5-7 and 12 | qLDPC constructions, expander decoding, finite instances, connectivity, logical operations |
| 14 | 3 and 4 plus oscillator quantum mechanics | cat, binomial, GKP, autonomous QEC, fault-tolerant bosonic gadgets |
| 15 | 3, 8, 11, 12, and 14 | erasure and fusion architectures, hardware-noise shaping, code-noise-decoder co-design |
| 16 | the branch relevant to each source | bounded reconstruction of current results, including qLDPC break-even and logical resources |

## Branch structure

### Stabilizer-to-topological branch

`5 -> 6 -> 8 -> 11 -> 12`

The learner must be able to distinguish stabilizers, logical operators, and fault propagation before
surface-code pictures become meaningful. Decoder algorithms arrive after the syndrome graph is
derived, not before.

### Stabilizer-to-qLDPC branch

`5 -> 6 -> 7 -> 12 -> 13`

qLDPC begins only after CSS commutation, code parameters, degeneracy, and belief-propagation traps
are available. Asymptotic parameters never substitute for a verified finite instance.

### Channel-to-bosonic branch

`3 -> 4 -> 14 -> 15`

The reader needs Kraus channels and exact versus approximate correctability before finite-energy GKP
states, biased cat noise, or autonomous stabilization can be discussed honestly.

### Architecture branch

`8 -> 9 -> 10 -> 12 -> 15 -> 16`

This branch distinguishes an abstract code from a physical implementation: extraction, gates,
decoder latency, routing, reset, fusion, and resource overhead all enter before an architecture
claim can be audited. Fusion-based routes also require the cluster-state material in Part 11.

## Entry routes

- A reader who completed `QuantumInDiscreteSpace` may begin with Part 0 as a diagnostic and skip only
  questions whose outputs they can reproduce.
- A classical coding reader may begin with Part 1 but still needs Part 0.2-0.6 before Part 2.
- A stabilizer expert may begin with Part 8 only after passing the Part 3 noise and Part 6 distance
  diagnostics.
- A bosonic-QEC reader may take `0 -> 2 -> 3 -> 4 -> 14`, but Part 15 and the frontier still require
  decoder and fault-tolerance vocabulary from Parts 8, 10, and 12.
- No reader begins with Part 16. Frontier reconstruction without the relevant branch produces paper
  summaries, not understanding.

## Back edges forbidden

- Do not use the surface-code threshold to define a threshold in Part 10.
- Do not use a decoder library before the learner has constructed one exact lookup decoder.
- Do not invoke qLDPC asymptotics to justify a finite-code overhead claim.
- Do not explain general correctability by quoting stabilizer syndromes; Knill-Laflamme is broader.
- Do not call a flag circuit fault-tolerant before enumerating the fault set and residual logical class.
