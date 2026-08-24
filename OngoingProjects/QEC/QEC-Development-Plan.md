# QEC Development Plan

Quantum error correction protects one *logical* qubit by spreading it across many *physical* qubits, in a subspace chosen so that the errors we care about push the state out of that subspace in a way we can detect by measuring a fixed set of parity checks, and then undo. Everything below is organized around that one idea: build the code, read out the error's fingerprint, guess the error, undo it, and ask whether the logical qubit actually survived the noise.

---

## 1. Where we are now

There are two layers, and they are in very different shape.

### The stabilizer engine (shipped in the paclet: solid)

The foundation in `Kernel/Stabilizer/` is good and trustworthy. In physics terms it lets you:

- hold an n-qubit stabilizer state compactly (the state is stored as the Pauli operators that fix it, not as 2^n amplitudes);
- push it through Clifford gates in the Heisenberg picture (each gate conjugates the stabilizers), with a fast compiled path;
- **measure** any qubit or check operator, getting a definite outcome when the state already fixes that operator and a random 50/50 outcome otherwise; this projection is exactly the syndrome-extraction step QEC runs on;
- read expectation values, inner products, and entanglement entropy;
- work with the textbook codes as states (5-qubit, Steane/7, Shor/9, in their |0_L> and |1_L> versions), with graph states and Clifford channels alongside.

It is cross-checked against Stim and QuantumClifford.jl, so the numbers can be trusted. This is a sound base to build QEC on.

### The code layer (new, in `OngoingProjects/QEC/`: a rough first pass)

On top of the engine, Maurice's package treats a set of stabilizer generators as a *code* and computes its structure: how many logical qubits it holds and how far apart its codewords are (the distance) with a concrete witness error; the logical operators; a circuit that prepares a codeword; the CSS, concatenation, and qubit-surgery constructions; and, for a single-qubit error, the lookup from syndrome to correction.

Honest status, and it matters for planning: **the mathematics is right on the small cases (distances, logical operators, and encoders all check out against an independent re-derivation, and the package's own 179 tests pass), but this is a pre-QA prototype and needs heavy refactoring and revision before it can be the foundation.** Specifically:

- it is written in an imperative style (explicit loops, building lists in place) rather than the functional idiom the rest of the paclet uses;
- Pauli operators are carried as strings and the algebra is done with dense mod-2 matrices, which is fine for a 9-qubit code but does not carry to the code sizes we will actually want;
- distance is found by brute-force search over error weights (unavoidable in general, since computing exact distance is NP-hard, but the search is not organized for larger codes);
- several exported names are too generic and will collide (`Syndrome`, `CodeDistance`, `RemoveQubit`, `LogicalOperators`, `PasteCodes`);
- a couple of the messages are wrong or misleading, and it only corrects a single error;
- and it is not part of the paclet, has no documentation, and its tests are not in the standard test harness.

So: trust the numbers on small codes, but treat the code layer as something to rebuild cleanly, not to ship as it stands. Separately, older QEC material (5/7/9-code notebooks, magic-state distillation, the `Courses/QuantumErrorCorrection` set) exists but is scattered and unreviewed; harvest the good parts into whatever we build rather than leaving them loose.

The honest one-line summary: we can represent and measure codes, and analyze the small ones exactly, but we cannot yet answer the question QEC is actually about, which is whether a code protects a qubit under realistic noise.

---

## 2. What to develop

Ordered the way the physics builds up. Each item says what the physics is, then the concrete function that would deliver it and how it works.

### (a) A code is a subspace that hides a logical qubit. Make that the clean core.

*Physics.* A stabilizer code is a subspace fixed by a commuting group of Pauli checks; the logical information lives in the Pauli operators that commute with every check but are not themselves in the check group. We already compute all of this, but across loose functions. We want one trustworthy object that *is* the code and carries its checks, its logical operators, its code space, and its parameters [[n, k, d]].

*Deliver.* A rebuilt code object, `QECCode[...]`, constructible from stabilizer generators, from a parity-check matrix, or by name, exposing `["Distance"]`, `["LogicalOperators"]`, `["Parameters"]`, and an encoding circuit that prepares a codeword. This is the clean-up of the current package into a core the rest can rely on.

### (b) An error leaves a fingerprint; the decoder guesses the error from it.

*Physics.* A physical error either commutes or anticommutes with each check. The pattern of anticommutations is the *syndrome*, and it is the only thing we are allowed to learn: measuring the checks reveals the syndrome without disturbing the logical qubit. Correcting is then an inference problem: given the syndrome, find the most likely physical error that could have produced it. Our current lookup only inverts single-qubit errors.

*Deliver.* A decoder that covers every error the code guarantees to fix (up to weight ⌊(d-1)/2⌋) for small codes by most-likely inference, and, for the surface-code family where errors form chains on a lattice, hands the syndrome to minimum-weight matching (through PyMatching). `DecodeSyndrome[code, syndrome]` returns the inferred correction; physically it is choosing the shortest set of flips consistent with the fingerprint.

### (c) Noise is probabilistic. Give it a physics object.

*Physics.* Real errors are not chosen adversarially; each qubit independently suffers X, Y, or Z with some probability p (depolarizing, bit-flip, or phase-flip), and the check measurements themselves misreport with some probability q. A code only means something relative to such a noise model. The engine has Clifford *channels*, but not "each qubit flips with probability p" as something you can sample.

*Deliver.* A noise object you attach to a code and draw errors from, at two honesty levels: the simplest, where each qubit flips with probability p between rounds (code-capacity noise); and the realistic one, where the syndrome-extraction gates and measurements are themselves noisy (circuit-level noise). `NoiseModel[...]` plus a sampler that returns a random error. For the circuit-level case, emit the syndrome-extraction circuit and let Stim sample it quickly.

### (d) The real question: run the cycle under noise and see if the logical qubit survives.

*Physics.* This is the whole point of a code. Prepare a codeword, hit it with a sampled error, measure the checks to get the syndrome (with measurement noise included), decode, apply the correction, and ask whether the *logical* operator got flipped. Repeat many times, and the fraction of logical flips is the **logical error rate**, the code's true figure of merit.

*Deliver.* `LogicalErrorRate[code, noise]` returning that probability, built on a `CorrectionCycle` upgraded to run under a noise model over one or many rounds. For small codes we can also sum the exact error probabilities in closed form, which is the exact check on the sampled number.

### (e) The threshold: when adding qubits starts to help.

*Physics.* The deepest fact in QEC. Below a critical physical error rate p_th, a larger code (bigger distance d) has an exponentially *smaller* logical error rate, so adding qubits suppresses errors; above p_th, adding qubits makes things worse. Measuring p_th, and showing you are below it, is the headline result of every serious QEC experiment (Google's Willow demonstrated a surface code operating below threshold).

*Deliver.* Run (d) across a family of code sizes and physical error rates and find the crossing point. `Threshold[family, noise]`, together with the characteristic figure: logical error rate versus p, one curve per distance, all crossing at p_th.

### (f) The codes that actually scale: the surface code and qLDPC.

*Physics.* Our named codes are the textbook [[3,1,1]], [[5,1,3]], [[7,1,3]], [[9,1,3]], all fixed and tiny. Fault tolerance is built on the surface code, a two-dimensional lattice of qubits whose distance grows with the lattice size and whose checks are local and easy to measure, and, more recently, on quantum-LDPC codes such as the bivariate-bicycle "Gross" code [[144, 12, 12]] that pack many logical qubits into far fewer physical ones. We have none of these.

*Deliver.* Constructors that build these families at a chosen size and return the same code object, so distance, decoding, and thresholds all work on them unchanged: `SurfaceCode[d]`, `ToricCode[...]`, `BivariateBicycleCode[...]`.

### (g) Computing on the protected qubit without breaking it (fault tolerance).

*Physics.* Storing a qubit safely is only half the job. You must also apply logical gates and measurements so that a single fault cannot grow into an error larger than the code can fix. That is transversal gates, lattice surgery (merging and splitting surface-code patches to entangle logical qubits), and magic-state distillation (manufacturing the one resource transversal gates cannot supply). And you must be able to count the cost: how many physical qubits and how much time a given logical computation needs.

*Deliver (longer term).* Logical-operation constructors, and a resource estimate (physical qubit count, time, magic-state budget) for a target logical circuit.

### The one thing that is genuinely ours, running through all of the above

Everything we compute exactly on a small or symbolic code, we check against Stim, the community's fast simulator. The one move Stim structurally cannot make is to carry a *symbolic* parameter and certify a whole family at once: a code parametrized by its size n, or a check left at a symbolic angle. That exactness and symbolic generality, not speed and not scale, is what Wolfram adds. Stim simulates fast and large; we derive and certify exactly and symbolically. The two are complementary, and we should lean on that difference rather than fight it.

---

## 3. What we should not build

Not a fast large-scale simulator (Stim owns that), not a production or real-time decoder (PyMatching and Riverlane own that), not lattice-surgery compilation at industrial scale (Qiskit and tket own that), not GPU decoding. Connect to those; do not rebuild them.

---

## 4. Order of work

1. Rebuild and QA the code core (a), and fold in the good scattered material. Nothing else is reliable until this is clean.
2. Add the decoder (b) and the noise object (c), so a code can be *tested* under noise, not just described.
3. Add the logical error rate (d) and the threshold (e): the physics figures of merit, and the first results worth showing anyone.
4. Add the scalable families (f).
5. Fault-tolerant operations and resource counting (g), longer term.

Check every exact result against Stim as you go, and put the weight on the exactness-and-symbolic angle wherever it is genuinely ours.

---

## References

- Gottesman, *Stabilizer Codes and Quantum Error Correction* (thesis), and Nielsen & Chuang chapter 10: the stabilizer formalism and the code constructions the engine and the package rest on.
- Gidney, *Stim: a fast stabilizer circuit simulator*, [arXiv:2103.02202](https://arxiv.org/abs/2103.02202); Higgott & Gidney, sparse-blossom matching (PyMatching).
- Google Quantum AI, below-threshold surface code (Willow), [Nature 2024](https://www.nature.com/articles/s41586-024-08449-y).
- Bravyi et al., IBM bivariate-bicycle qLDPC codes / the Gross code [[144,12,12]], [arXiv:2308.07915](https://arxiv.org/abs/2308.07915).
- Minimum-distance of stabilizer codes is NP-hard: [arXiv:2203.04262](https://arxiv.org/abs/2203.04262).
</content>
