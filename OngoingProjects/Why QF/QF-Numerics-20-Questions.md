# Twenty quantum questions to develop with QF

Rebuilt: 2026-10-03. Status updated: 2026-10-05. This collection contains **twenty proposed computational studies for review**. Each specifies a physical question, a bounded investigation, an answer, and substantive QF work. N17 now has a [worked draft](N17-Atomic-Revival-and-Wigner-Negativity.md); the other cards remain proposals in this collection.

**Selection context.** The cards were chosen for scientific payoff, numerical substance, substantive QF contribution, breadth across fields and capabilities, and manageable scope for readers with master's-level quantum mechanics. Develop a selected question through physical setup, visible computation and interpretation, with one core action per code cell. The linked papers support the models and mechanisms; they do not by themselves establish QF feasibility, completed results or an advantage over Qiskit. The collection is a candidate pool, not an implementation priority list.

Existing IDs are preserved for continuing subjects. N01, N07, N12, N13 and N16 are retired from this collection; N21–N25 are new. N19 has a different, operational switch question. Numbering below counts the twenty cards; it is not a ranking or an implementation order.

| No. | Card | Question |
|---:|---|---|
| 1 | N21 | How does subsystem size limit local equilibration? |
| 2 | N23 | Can correcting the exchange frequency rescue a spin-only model? |
| 3 | N10 | How much dephasing helps an excitation reach its destination? |
| 4 | N14 | Which pulse overlap gives the best transfer through a lossy level? |
| 5 | N11 | How does reaction lifetime limit directional contrast? |
| 6 | N09 | Can a convincing spectral envelope hide wrong vibrational intensities? |
| 7 | N17 | How much cavity Wigner negativity can the atomic population record guarantee? |
| 8 | N08 | How many spatial modes retain the information about source separation? |
| 9 | N02 | Should an illumination probe be maximally entangled? |
| 10 | N22 | Does entanglement in a reconstructed state establish entanglement in the experiment? |
| 11 | N03 | When does imperfect purification improve pair quality, and at what yield? |
| 12 | N04 | When should the same error syndrome trigger a different correction? |
| 13 | N06 | How many phase-estimation bits distinguish two known energies? |
| 14 | N05 | Does the best simulation depth depend on what we measure? |
| 15 | N25 | Does maximizing the average cut maximize the chance of an optimal cut? |
| 16 | N24 | Can resetting a qubit erase its past? |
| 17 | N15 | Which trapped states survive a change in interface shape? |
| 18 | N18 | What heralding cost buys a better optical signal? |
| 19 | N19 | How much control coherence does a switch need to preserve entanglement? |
| 20 | N20 | How much extractable energy does a measurement record make accessible? |

## 1. N21: How does subsystem size limit local equilibration?

A closed quantum system remains pure under unitary evolution, yet measurements on a small part can become nearly stationary. Do larger blocks remain as close to their time-averaged states as smaller ones over the same observation window?

Start an eight-spin chain in an alternating up/down product state. Use nearest-neighbor Ising interactions with both longitudinal and transverse fields, avoiding the integrable transverse-field limit. Compare blocks of one through four adjacent spins with their reductions of the global infinite-time average, the **diagonal ensemble**. Keep the Hamiltonian and observation window fixed.

**What answers it.** Reduced-state trace distance over time, its late-window fluctuations, and block entropy. These distinguish a mixed subsystem from one that actually stays near its time average; persistent fluctuations and recurrences belong in the answer.

**QF work.** Build labeled spin interactions, evolve the pure joint state, and obtain block states and entropy through subsystem operations.

**Learning payoff.** Determine the limits of local equilibration without mistaking it for global mixing or thermalization. [Physical basis: thermalization through entanglement](https://doi.org/10.1126/science.aaf6725); the spin-chain comparison is an adaptation.

## 2. N23: Can correcting the exchange frequency rescue a spin-only model?

Two opposite-spin fermions in neighboring wells can exchange places by tunneling through configurations with both particles in one well. At strong repulsion this is approximated by a two-spin Heisenberg interaction. As tunneling grows, does that approximation merely get the oscillation frequency wrong, or does it omit essential charge motion?

Use a symmetric, half-filled Hubbard dimer, initially one opposite spin per site. Vary the repulsion-to-hopping ratio. Compare exact evolution with two spin models: the usual exchange coupling $4J^2/U$, and the exact low-energy singlet–triplet splitting. Here $J$ is hopping energy and $U$ is on-site repulsion.

**What answers it.** Spin-exchange probability, double-occupancy probability, and error against each embedded spin-model state over the first exchange cycles. Retain population outside the spin subspace.

**QF work.** Compose explicitly encoded fermionic operators, evolve the occupation state, and use a rectangular isometry to compare the four-dimensional spin space with the sixteen-dimensional occupation space. Encode the fermionic anticommutation signs in the hopping operators.

**Learning payoff.** Separate a repairable frequency error from failure of the retained state space. [Model and exact limits](https://arxiv.org/abs/1502.02194).

## 3. N10: How much dephasing helps an excitation reach its destination?

In a symmetric transport network, part of an excitation can occupy a superposition that never reaches the exit. Dephasing can release that trapped component, but strong dephasing also slows transfer. Which rate gives the best chance of collection before recombination removes the excitation?

Use three equal-energy sites: sites 1 and 2 couple equally to site 3, which feeds an irreversible trap. Start on site 1. Include recombination into a separate absorbing state. Hold couplings, trap rate, recombination rate and collection time fixed; vary only local site dephasing.

**What answers it.** Collected probability at the fixed deadline, together with trapped, surviving and lost populations. Representative time histories should explain the low- and high-dephasing limits.

**QF work.** Combine coherent site operators and physical jumps in one five-state model; obtain capture and loss directly from quantum evolution and measurements.

**Learning payoff.** Identify when destroying coherence improves transport and when it prevents arrival. [Mechanism](https://doi.org/10.1088/1367-2630/11/3/033003); this three-site network is a minimal adaptation.

## 4. N14: Which pulse overlap gives the best transfer through a lossy level?

Two long-lived levels communicate through an intermediate level that can decay. Driving the destination coupling first can guide the state through a superposition with little intermediate occupation. Too little pulse overlap, or evolution that is too rapid, can defeat that protection.

Start in the initial level of a resonant three-level system, with decay from the intermediate level into a fourth, absorbing state. Compare both orders of the same Gaussian pulses. Fix their width, peak coupling, decay rate and a common observation interval that contains every pulse pair; vary their separation.

**What answers it.** Final destination and lost populations versus signed pulse delay, supported by intermediate-level histories. The best transfer settings must be distinguished from settings that merely avoid loss by doing little.

**QF work.** Assemble the time-dependent transition operators and loss jump, evolve directly, and measure populations on the resulting quantum states.

**Learning payoff.** Understand how finite-time dynamics limits dark-state transfer. [Pulse mechanism and experimental context](https://doi.org/10.1038/ncomms10628); the absorbing-loss model is a generic adaptation.

## 5. N11: How does reaction lifetime limit directional contrast?

A chemical reaction can depend on whether two electron spins are in their singlet or triplet sector. An anisotropic coupling to a nuclear spin makes their conversion sensitive to field direction. If recombination happens too early, how much directional information reaches the reaction products?

Use two electron spins and one unpolarized spin-half nucleus, initially with the electrons in a singlet. Adopt the axial hyperfine reference-and-probe model: only one electron couples to the nucleus. Hold the field magnitude and hyperfine coupling fixed; vary field direction and the common singlet/triplet lifetime.

**What answers it.** The singlet probability integrated against the exponential reaction-time distribution, and its maximum-minus-minimum angular contrast versus lifetime. Compare reaction and spin-precession timescales to explain the result.

**QF work.** Compose electron–nuclear interactions on labeled subsystems, evolve their joint state, and measure the electronic singlet sector. WL supplies the lifetime integral.

**Learning payoff.** Connect coherent spin dynamics to a chemical yield, without attributing the model to a particular organism. [Model and physical scales](https://doi.org/10.1103/PhysRevA.85.040304).

## 6. N09: Can a convincing spectral envelope hide wrong vibrational intensities?

An electronic transition changes a molecule's equilibrium geometry, vibrational frequencies and normal-mode directions. Neglecting mode mixing may distort individual transition strengths even when broadening makes the overall spectrum look convincing. How does spectral broadening hide those intensity errors?

Use the published two-mode harmonic model for the $\mathrm{SO}_2\rightarrow\mathrm{SO}_2^+$ transition, initially in the vibrational ground state. Compare its full coordinate transformation with an independent-mode approximation. Keep frequencies and the displacement in the stated coordinate convention unchanged; remove only the mode rotation. Apply the same broadening to both spectra.

**What answers it.** Resolved line-strength differences and the discrepancy between broadened spectra as resolution changes. Highlight errors in appreciable lines, rather than large relative errors in negligible ones.

**QF work.** Compose the correctly ordered squeezing, mode mixing and displacement operators; joint Fock probabilities give line strengths. Control the finite Fock cutoff.

**Learning payoff.** Test whether agreement with a smooth measured envelope validates its microscopic assignments. [Published model and parameters](https://doi.org/10.1039/C7SC04602B).

## 7. N17: How much cavity Wigner negativity can the atomic population record guarantee?

**Status:** A [worked computational draft](N17-Atomic-Revival-and-Wigner-Negativity.md) is available in this folder. It records the implemented model, results and numerical checks.

An atom exchanging excitations with a cavity shows oscillations that collapse and revive. The atom's population record is easy to measure; the field's Wigner function is not. Given only that record, how much negative volume of the unconditional cavity field can be guaranteed?

Start an excited two-level atom with a coherent cavity field of mean photon number four, under resonant Jaynes-Cummings dynamics with cavity loss. Compare it with a phase-randomized coherent field that has the same photon-number distribution.

**What answers it.** Excitation-number conservation makes the record depend only on the photon-number distribution, so the two fields give the same record at every loss rate. The phase-averaged field's negative volume is a lower bound for every field with that distribution (convexity), and the phase-randomized preparation attains it. The record determines the distribution once it extends past the collapse. Compare the guaranteed and the true negative volume over time and loss, and identify when an atomic superposition makes the record phase sensitive.

**QF work.** Evolve the joint qubit-Fock state for both preparations, reduce to atom and field, and evaluate the field's Wigner function through its angular harmonics, cross-checked against QF's Wigner representation.

**Learning payoff.** Establish exactly what an accessible subsystem signal can and cannot certify about another subsystem. [Atom-field mechanism](https://doi.org/10.1103/PhysRevLett.65.3385).

## 8. N08: How many spatial modes retain the information about source separation?

Two equally bright incoherent sources can produce an almost featureless merged image. A detector that separates spatial wave modes measures different information from one that records photon positions. How much detector complexity is needed to retain that information?

Assume a Gaussian point-spread function and a known centroid. Compare position detection with resolving the first $K$ Hermite–Gaussian modes and collecting all remaining modes into one residual outcome. Use the same detected-photon budget and vary source separation relative to image width.

**What answers it.** Separation Fisher information for both detectors and the smallest $K$ retaining 95% of the complete mode measurement's information. Evaluate the approach to coincidence from positive separation; do not substitute a singular zero-probability outcome into a regular error formula.

**QF work.** Form the incoherent mixture of displaced single-photon spatial states and apply alternative complete measurements to that state. WL supplies mode overlaps, position integrals and information calculations.

**Learning payoff.** Turn an abstract resolution advantage into a finite measurement-design decision. Spatial modes are not photon-number states. [Imaging model and mode measurements](https://doi.org/10.1103/PhysRevX.6.031033).

## 9. N02: Should an illumination probe be maximally entangled?

A detector click may come from background rather than a reflecting target. If the background occupies one signal mode more often than another, an unentangled probe can favor the quieter mode. A retained partner offers correlations, but maximum entanglement need not balance these two advantages best.

Work within detected single-photon events in two distinguishable modes. A missing target returns the known background; a present target returns the prepared signal with fixed effective weight $\eta$, and background otherwise. At equal target priors, vary the background imbalance. Compare the best product probe, a Bell probe, and a probe with adjustable Schmidt weights in the background eigenmodes: $\sqrt{\lambda}|00\rangle+\sqrt{1-\lambda}|11\rangle$. Assume ideal partner storage.

**What answers it.** The best receiver and decision error for each preparation, together with the preferred Schmidt weights. All strategies send one signal photon; the entangled strategies additionally retain a partner.

**QF work.** Apply the channel only to the signal, retain joint correlations, and construct the minimum-error measurement from the received-state difference.

**Learning payoff.** Design the entanglement for the background instead of maximizing it indiscriminately. [Model and optimal probe](https://doi.org/10.1038/s41534-020-00303-z).

## 10. N22: Does entanglement in a reconstructed state establish entanglement in the experiment?

Tomography returns a physical density operator, but finite counts can move its estimated entanglement across the separable boundary. How often would interpreting that estimate as a certificate lead to the wrong conclusion, and what can the same counts establish with stated confidence?

Generate repeated synthetic experiments for two-qubit Bell–white-noise mixtures on both sides of the boundary. Measure all nine local Pauli settings with equal counts per setting. Compare the negativity of QF's maximum-likelihood estimate with a prespecified Bell-fidelity witness evaluated directly from those counts, using a one-sided 95% confidence bound. Assume independent preparations and ideal calibrated measurements.

**What answers it.** False-positive and detection frequencies versus sample size for the naive point-estimate rule and the confidence-based rule. Certification requires a lower confidence bound on Bell fidelity above $1/2$; an entangled point estimate alone does not suffice. Inconclusive data are an allowed outcome.

**QF work.** Simulate quantum measurement counts, reconstruct physical states, and evaluate their entanglement; WL supplies the explicit statistical bound and repeated-experiment summary.

**Learning payoff.** Separate state estimation from evidence for a physical claim. MLE remains a useful estimator; its point estimate is not a confidence statement. [Statistical issue and direct bounds](https://doi.org/10.1103/PhysRevLett.114.080403).

## 11. N03: When does imperfect purification improve pair quality, and at what yield?

Purification consumes two noisy entangled pairs and keeps one only for selected measurement outcomes. Better accepted pairs can be valuable, but the local gates can spoil the improvement and failed attempts consume resources.

Start with two identical Bell-diagonal pairs: overlap $F$ with the desired Bell state and equal populations in the other three. Apply one bilateral-CNOT round, measuring the sacrificial pair and accepting matching outcomes. After each local CNOT, apply independent depolarization of equal strength to its two qubits; keep readout ideal. Vary input overlap and gate error.

**What answers it.** Accepted-pair overlap alongside acceptance probability and surviving pairs per input pair. Identify where quality improves and which quality–yield combinations are attainable. Do not assume the accepted state retains the input's symmetric Bell populations.

**QF work.** Compose the four-qubit circuit and gate channels, carry accepted branches with their weights, and reduce to the surviving pair.

**Learning payoff.** Judge a purification step by both its delivered quality and its pair consumption. [Recurrence protocol](https://doi.org/10.1103/PhysRevA.54.3824); imperfect local gates are the specified adaptation.

## 12. N04: When should the same error syndrome trigger a different correction?

A three-qubit phase-flip code cannot distinguish a single phase flip from the complementary two-qubit flip: they have the same syndrome. The right correction therefore depends on which physical error is more likely. Can a decoder calibrated only to individual error rates make the encoded qubit worse?

Mix independent phase flips of marginal probability $p$ with a burst channel that applies one uniformly selected two-qubit phase flip with total probability $3p/2$. Both components give every qubit the same marginal error $p$. Vary their mixture weight, with $0<p\leq1/4$. Keep encoding, syndrome readout and correction ideal.

**What answers it.** Conditional error probabilities for each syndrome, the boundary at which its preferred recovery changes, and decoded entanglement fidelity for conventional and correlation-aware decoders versus an unencoded qubit.

**QF work.** Preserve an untouched reference, encode, apply the two noise models, recover full syndrome branches, and decode before evaluating the reference–logical state.

**Learning payoff.** Show why equal physical error rates do not imply the same recovery rule. [Correlated-noise context](https://arxiv.org/abs/quant-ph/9702001); the matched-marginal burst model is defined here.

## 13. N06: How many phase-estimation bits distinguish two known energies?

A phase-estimation register produces a distribution of binary readouts for each energy. Two energies with the same separation can be easier or harder to distinguish depending on their alignment with that binary grid. What register size makes the label sufficiently reliable?

Use two known eigenphases separated by $1/32$ of a cycle, with equal prior weights. For a fixed base evolution time $t_0$, the energy gap obeys $\Delta E\,t_0/(2\pi\hbar)=1/32$. Vary their known common offset while keeping both phases in one alias-free interval. Prepare each eigenstate separately to obtain conditional readout distributions; their equal mixture describes an input with equal eigenstate weights.

**What answers it.** The minimum Bayes assignment error for each register size and offset, the smallest size reaching 5% error, and the longest and total controlled evolution required in units of $t_0$. Distinguish best tested results from analytical sufficient bounds.

**QF work.** Keep target and readout registers explicit, compose controlled propagator powers and inverse Fourier transformation, then measure the readout register.

**Learning payoff.** Translate energy resolution into an operational error and coherent-time cost, rather than counting visibly separated peaks. [Phase-estimation derivation](https://doi.org/10.1098/rspa.1998.0164).

## 14. N05: Does the best simulation depth depend on what we measure?

Increasing the number of product-formula steps can reduce the error from splitting noncommuting evolution, but require more imperfect gates. The best compromise may differ when we need one local observable rather than the entire final state.

Use a four-spin chain with nearest-neighbor $\hat X_j\hat X_{j+1}$ interactions and a transverse field proportional to $\hat Z_j$, initially in an alternating computational-basis state. Fix couplings and evolution time. Compare first-order product formulas with different step counts, applying independent local depolarization of equal strength to both participating spins after every two-spin factor.

**What answers it.** Error in the first spin's polarization $\langle\hat Z_1\rangle$ and whole-state trace distance from exact ideal evolution, versus step count at fixed gate-error strengths. Report the best tested count under each criterion and explain discrepancies using the resulting states.

**QF work.** Construct the Hamiltonian and circuit factors from the same labeled operators, compose noise, and measure both requested outputs without changing representations.

**Learning payoff.** Choose simulation effort for the physical quantity required. An accidental cancellation in one observable is not whole-state accuracy. [Noisy product-formula context](https://arxiv.org/abs/2504.10247); this small first-order comparison is an adaptation.

## 15. N25: Does maximizing the average cut maximize the chance of an optimal cut?

A variational circuit for graph partitioning produces a distribution of cuts. Improving its average cut does not necessarily maximize the chance of obtaining an optimal cut in a limited number of measurements. Which objective better serves that goal on the same circuit?

Use depth-one QAOA on a six-vertex ring with two extra edges, $(1,3)$ and $(2,5)$; all eight edges have unit weight. Prepare six $|+\rangle$ qubits and apply one cost layer and one mixer layer. Compare angle choices maximizing mean cut with those maximizing the average of the best 25% of the probability distribution, the upper-tail CVaR. Use exact probabilities to isolate objective choice from sampling and optimizer noise.

**What answers it.** Full cut distributions, optimal-cut probability and repetitions needed for 95% success. Compare against uniform sampling and the best tested success probability; report CVaR ties rather than choosing a favorable one silently.

**QF work.** Build the parameterized cost/mixer circuit and obtain its quantum output distribution. WL evaluates cut scores and searches the two-angle domain.

**Learning payoff.** Distinguish a convenient variational objective from the success criterion of the algorithm. [Objective design](https://doi.org/10.22331/q-2020-04-20-256); this graph is an illustrative instance, not a quantum-speedup claim.

## 16. N24: Can resetting a qubit erase its past?

If a qubit is reset to the same state, should its later behavior be independent of what happened before? That expectation can fail when the environment retains information and interacts with it again, even while an uninterrupted evolution shows no revival of distinguishability.

Prepare a system qubit in either $|0\rangle$ or $|1\rangle$, with the same maximally mixed environment qubit. Apply a partial swap, discard and reprepare only the system in $|0\rangle$, then apply a second identical partial swap with the retained environment. Compare with replacing the environment by a fresh one. Use $\hat U_\theta=\cos\theta\,\hat I-i\sin\theta\,\widehat{\mathrm{SWAP}}$, with $0<\theta\leq\pi/4$, so uninterrupted distinguishability decreases throughout both collisions.

**What answers it.** The difference between the two final system states and a measurement that detects it, alongside uninterrupted trace distance. The reset state and future operations must be identical for both histories.

**QF work.** Apply a subsystem replacement channel to the correlated joint state, retain the environment, and compose the later interaction and measurement.

**Learning payoff.** Detect memory through an intervention, rather than identifying it only with revivals in a decay curve. [Operational criterion and partial-swap construction](https://doi.org/10.1103/PhysRevLett.120.040405).

## 17. N15: Which trapped states survive a change in interface shape?

A quantum walker can remain near an interface between different topological phases. A smooth interface can also bind additional states whose existence depends on its shape. How much observed trapping comes from the protected state, and how much from these extra states?

Use a finite one-dimensional split-step walk with a two-state coin. Hold its two bulk walks and the initial localized state fixed. Change the interface from sharp to smooth by varying the second coin rotation in position, while the first remains uniform so the relevant chiral symmetry is preserved.

**What answers it.** Localized eigenphases, their initial-state overlaps and boundary-residence probabilities versus interface width. Separate modes pinned at zero or $\pi$ from additional paired modes. Keep remote-boundary effects outside the observation window.

**QF work.** Compose conditional shifts and coin rotations on the coin–position space, inspect the step operator, and obtain position probabilities by subsystem reduction.

**Learning payoff.** Distinguish protection of a state's existence from protection of an observed trapped population. [Model and ordinary interface modes](https://arxiv.org/abs/1003.1729).

## 18. N18: What heralding cost buys a better optical signal?

Loss weakens an optical coherent signal. A quantum-scissors device can conditionally amplify its vacuum and one-photon components, but it truncates higher components and rejects many attempts. What improvement is physically attainable, and at what success rate?

Fix the original brightness and pure-loss transmission, leaving the signal phase unknown. Use one ancillary photon, vacuum, a tunable beam splitter, a balanced beam splitter and ideal photon counting. Accept the two one-photon herald patterns with their corresponding phase corrections; vary only the tunable splitter.

**What answers it.** Fidelity to the original signal versus total heralding probability, with direct transmission as the baseline. Derive the accepted instrument from the optical network before reducing it to a two-dimensional output. Establish phase covariance rather than assuming the unknown phase is available for correction.

**QF work.** Carry optical operations and weighted measurement branches into a rectangular map from the input Fock space to the vacuum–one-photon space; embed that output for comparison with the target.

**Learning payoff.** Distinguish a physical conditional improvement from an arbitrarily rescaled nonunitary matrix. [Single-photon scissors and success weights](https://doi.org/10.1103/PhysRevA.102.063715).

## 19. N19: How much control coherence does a switch need to preserve entanglement?

Applying complete dephasing in the $X$ and $Z$ bases successively erases all information about the input state. A coherent control of their order can make one control outcome identify a correctable branch. How much coherence must that control retain for the heralded output to preserve entanglement with a reference?

Send half a Bell pair through complete $X$- and $Z$-dephasing channels. Construct their switch from the four pairs of scaled Kraus operators. Initialize the order-control qubit in $|+\rangle$; vary its remaining off-diagonal coherence factor $v$ from zero to one. Measure it in the $X$ basis, accept $|-\rangle$ and apply $\hat Y$ to the transmitted qubit. Compare with either sequential order and with an incoherent order control.

**What answers it.** Heralding probability, conditional Bell fidelity and negativity versus control coherence, including the entanglement-loss boundary. More heralds need not mean a better link.

**QF work.** Combine named operator switches into an explicit Kraus channel, retain the untouched reference, and evaluate weighted measurement branches.

**Learning payoff.** Turn coherent order into an operational link-quality question. The comparison concerns this heralded protocol, not universal capacity superiority. [Channel construction and ideal branch](https://arxiv.org/abs/1809.06655).

## 20. N20: How much extractable energy does a measurement record make accessible?

Alice measures one of two interacting thermal spins and sends the outcome to Bob, who controls the other. Her measurement can inject energy. To isolate the value of the record, Bob must be compared with someone acting on the same disturbed ensemble without knowing the outcome.

Use equal local fields proportional to $\hat Z_A$ and $\hat Z_B$, and a coupling proportional to $\hat X_A\hat X_B$. Alice measures $\hat X_A$; Bob chooses a rotation generated by $\hat Y_B$. At fixed coupling, vary temperature and optimize the same rotation family with and without the record. Neglect free evolution during the measurement-and-control sequence.

**What answers it.** Bob's average extracted energy in both cases, Alice's injected energy and the initial thermal entanglement. Determine whether the record remains useful after the initial state becomes separable. Evaluate energy using the full interacting Hamiltonian.

**QF work.** Construct the Gibbs state, carry weighted measurement branches into conditional local rotations, and evaluate energy on the resulting quantum states.

**Learning payoff.** Separate usable information, measurement disturbance and entanglement in one energy account. [Two-spin protocol](https://doi.org/10.1088/1751-8113/46/45/455304); the matched post-measurement baseline is our stated adaptation.

## Undermind research references

Saved discovery search used to develop this collection. All 201 titles/relevance summaries were inspected during selection; selected papers received targeted or full-text examination. This does not mean that all 201 papers were read in full or that the proposed QF studies were executed.

- [Compact computational questions across quantum science](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3?path=/Compact%20computational%20questions%20across%20quantum%20science): 201 candidate-paper records spanning compact studies across quantum science.

Workspace: [Roadmap for Wolfram quantum](https://app.undermind.ai/projects/bcedc042-7a0d-4f3d-b4ca-4e1c138571a3); workspace ID: `bcedc042-7a0d-4f3d-b4ca-4e1c138571a3`. Links checked on 2026-10-04. For reuse, inspect the saved search goal and results before launching a new search, follow the primary-source links on each card, and check the current QF implementation when developing a study.
