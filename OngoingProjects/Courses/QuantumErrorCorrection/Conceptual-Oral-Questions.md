# Conceptual and Oral Checkpoints

These checkpoints test whether the learner can reason about QEC without hiding behind a simulator.
They complement the computational bank; they are not answerless versions of coding exercises.

## Protection and information

1. Why does QEC require redundancy but not copies of an unknown state?
2. What information does an ideal syndrome reveal, and what logical information must it not reveal?
3. Why does correcting $X$, $Y$, and $Z$ suffice for an arbitrary one-qubit error but not make
   arbitrary physical noise equivalent to a stochastic Pauli channel?
4. Explain correctability as absence of logical information from the environment.
5. Why is known-location erasure correction stronger than unknown-location error correction at the
   same code distance?
6. Give a case where a recovery succeeds on $\lvert 0_{L}\rangle$ and $\lvert 1_{L}\rangle$ but still
   fails on a superposition or on entanglement with a reference.

## Codes and logical equivalence

7. Why does a syndrome label an error coset rather than an error?
8. When are two physical Pauli errors logically equivalent?
9. How can a degenerate code assign the same syndrome to distinct physical errors without ambiguity?
10. Why is the smallest-weight syndrome representative not always the most likely logical class?
11. What changes when a stabilizer code is replaced by a subsystem code?
12. Why do excellent asymptotic rate and distance say little by themselves about one finite code's
    syndrome circuit and decoder?

## Fault tolerance and scale

13. Why is correcting the final data state insufficient to certify a syndrome-extraction circuit as
    fault-tolerant?
14. Distinguish code distance from circuit-level distance using a hook error.
15. Why can a transversal gate be fault-tolerant but not belong to a universal transversal set?
16. What does Eastin-Knill forbid, and what does it leave open?
17. Distinguish threshold, pseudo-threshold, error suppression with distance, and break-even.
18. Why can a logical-error floor remain even far below the threshold predicted by an independent
    stochastic model?

## Decoders, hardware, and evidence

19. When does a decoder need the whole syndrome history rather than the latest syndrome?
20. Why is decoder latency part of the quantum architecture rather than a classical afterthought?
21. What information is gained by an erasure flag, and what physical error can still remain hidden?
22. Why can a code optimized for biased noise lose its advantage under bias-destroying gates?
23. What distinguishes an experimental demonstration from an architectural projection assembled
    from several separate demonstrations?
24. Why is there no defensible one-number ranking of current QEC platforms?
25. Why can two physical channels with the same average gate infidelity have sharply different
    worst-case and logical effects, and why is average fidelity alone insufficient for a threshold
    claim?
26. Why do surface-code-like schemes repeat noisy syndrome measurements, whereas a single-shot code
    can use one sufficiently redundant noisy syndrome record?
27. What remains protected when a Floquet code has no single fixed stabilizer group throughout its
    measurement cycle?
28. In a photonic fusion architecture, which failures become located erasures and which residual
    faults must still be treated as unlocated Pauli errors?
29. Why do high rate and distance not make a qLDPC code computationally useful without a practical
    decoder, syndrome circuit, and fault-tolerant logical operations?
30. Why are positive coherent information, a below-threshold memory, and experimental break-even
    logically different claims?

## Oral-answer standard

A satisfactory answer names the relevant object and assumption, gives a concrete example, and
identifies a failure edge. Fluency with a slogan or theorem name is not sufficient.
