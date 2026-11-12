# Talk script: Continuous Quantum Measurement in the Wolfram Language

For each slide: the model, the question, how to read the figure, and the physics it shows.

## How the Wolfram Language handles a watched quantum system

Each slide asks one physical question and answers it twice, in closed form and by sampling quantum trajectories, and checks one against the other. Five ideas recur. The stochastic master equation is written once as an `ItoProcess`, which samples trajectories and applies Itô's lemma to any function of the state. Its average over records is a Lindblad master equation, built and evolved in QuantumFramework. Exact results come from the Liouvillian: its spectrum gives decay rates and frequencies, its discriminant locates exceptional points, and a linear solve gives the noise spectrum. A filter is Bayes' rule folded over the measurement record, checked against the full conditional state driven by the same record. Where possible, the equations come from operator algebra instead of matrices: the canonical commutation relation for a particle, the Pauli relations for a code.

## 1. A driven qubit under continuous measurement

**Model.** A Rabi-driven qubit under continuous σ_z measurement with strength k and efficiency η; single quantum trajectories against the Lindblad average.

**Question.** Does measurement stop the Rabi oscillation?

**Figure.** ⟨σ_z⟩ for three trajectories and the Lindblad average in red dashed, at η of 0.3 and 1; the controls set k and the detuning Δ, both in units of Ω.

**Physics.** The drive rotates the state and the measurement backaction pushes it toward an eigenstate of σ_z. For weak measurement each trajectory keeps oscillating, and the average decays only by dephasing between trajectories. For strong measurement each trajectory is pinned near a pole and jumps rarely: the quantum Zeno effect, which the average shows only as slow decay. With η below one the conditional state stays mixed and never reaches the poles. Detuning reduces the Rabi amplitude, so under weak measurement the record moves the state slowly between the poles.

## 2. The quantum Zeno crossover as an exceptional point

**Model.** The Liouvillian of the same measured qubit, whose eigenvalues give every decay rate and oscillation frequency of the averaged dynamics.

**Question.** What happens at the crossover from Rabi oscillation to Zeno freezing?

**Figure.** Top: the logarithm of the discriminant D of the characteristic polynomial, for complex k; the dots are exceptional points. Bottom: the slowest decay rate and the oscillation frequency for real k; dashed lines mark the exceptional points. Columns: four detunings.

**Physics.** On resonance the qubit is a damped harmonic oscillator, with k as the damping. Below the exceptional point it is underdamped; above it, it is overdamped and the slowest rate falls as k grows, which is the Zeno effect. The exceptional point, where k equals Ω, is critical damping: the frequency vanishes and the relaxation is fastest. Detuning splits it into two exceptional points, with no oscillation between them. At a detuning of one over two root two they merge and then move into the complex k plane, out of reach of any real measurement strength; the frequency never vanishes and the peak in the rate becomes a smooth maximum.

## 3. A charge qubit read by a point contact: when the Rabi peaks merge

**Model.** A double-dot charge qubit monitored by a quantum point contact; the current noise spectrum S(ω) in closed form, from the same Liouvillian.

**Question.** Do the Rabi peaks of S(ω) merge at the exceptional point?

**Figure.** S(ω) for several k, for the largest k with two peaks, and at each exceptional point; panels are level offsets ε; gray lines mark the bare Rabi frequency.

**Physics.** Coherent tunnelling makes the current fluctuate at the Rabi frequency, and measurement-induced dephasing broadens the peaks. Two broadened peaks merge when their width exceeds their splitting, which on resonance happens before the exceptional point. At an offset ε equal to Ω the peaks merge although there is no exceptional point at all. Peak merging in S(ω) is a line-shape effect and does not locate the exceptional point.

## 4. Dispersive readout pointer states

**Model.** A transmon coupled dispersively to a readout resonator with dispersive shift χ and linewidth κ, driven by a readout pulse (TransmonLab). The interactive circuit shows the Hamiltonian at three levels of approximation.

**Question.** How do the pointer states for the two qubit states separate as the ratio of κ to χ changes?

**Figure.** Left: the intracavity coherent states for the two qubit states in the IQ plane, with their steady states. Right: the separation of the two output fields. The control sets the ratio of κ to χ.

**Physics.** The qubit-dependent shift detunes the resonator in opposite directions, so the two coherent states acquire opposite phases. A narrow resonator separates them strongly inside but leaks slowly; a wide one leaks fast but separates them weakly. The output separation therefore settles fastest at large κ and is largest at intermediate κ, near the known optimum where κ is twice χ.

## 5. Born's rule as a martingale

**Model.** A qubit under perfect σ_z measurement, starting in a pure state with ⟨σ_z⟩ equal to z₀; 400 trajectories per start.

**Question.** With what probability does a trajectory end at the upper pole?

**Figure.** Left: the fraction of trajectories ending at the upper pole against Born's prediction, with error bars. Right: four trajectories and the mean of 400. Controls: z₀ and drive.

**Physics.** Without a drive, z has zero drift, so it is a martingale: its mean stays at z₀. Every trajectory still ends at a pole. A bounded martingale that ends at the poles must split between them in proportion to z₀, which is exactly Born's rule, so the rule follows from the dynamics. A drive adds drift, the mean leaves z₀, and the trajectories keep jumping.

## 6. Quantum Bayesian update

**Model.** A sequence of weak measurements alternating between σ_z and a tilted Pauli operator, each applied as a Kraus operator, Korotkov's quantum Bayesian update.

**Question.** How does each readout change the mixedness, the determinant of ρ?

**Figure.** Left: the factor by which each readout changes det ρ, against the inverse square of w, the probability of that readout; markers are tilt angles. Right: the mixedness over 20 readouts for low and high detector noise D.

**Physics.** Every Kraus operator has unit determinant, so each readout changes det ρ by a factor fixed by w alone, and all points lie on the diagonal. Likely readouts purify the state, unlikely ones mix it, and on average it purifies. At low D the purity grows by orders of magnitude in 20 readouts; at high D each readout carries little information and purification is partial.

## 7. Measurement purifies a qubit

**Model.** The Bloch-vector length r of a measured qubit, derived with Itô's lemma in spherical coordinates and sampled.

**Question.** Does continuous measurement purify the state, and what does an efficiency below one cost?

**Figure.** The stochastic differential equation for r, then r for three trajectories at drive strengths Ω of 0, 1 and 5 in units of k; the control sets η; the dashed line marks the square root of η.

**Physics.** Neither Ω nor Δ appears in the equation for r: a Hamiltonian only rotates the Bloch vector. The drift of r is positive inside the dashed line, where information gain beats dephasing. At η of one the line sits at the surface, the noise vanishes there, and every trajectory reaches the surface and stays pure. With η below one the drift at the surface is negative except at the poles, so a trajectory is pure only at a pole. Without a drive it reaches one; a strong drive keeps it away and it stays mixed.

## 8. A monitored qubit stays in the Bloch ball

**Model.** A driven, relaxing qubit read out along σ_z, with measurement-induced dephasing Γ_d and readout noise split into an informational part Γ_CI and a backaction part Γ_BA; their sum relative to Γ_d sets η.

**Question.** Does every trajectory stay a valid density matrix, inside the Bloch ball?

**Figure.** The noise of r² on the sphere and the condition for outward drift; then r for η of 0.5, 1 and 3, with the region outside the ball shaded.

**Physics.** On the sphere the noise of r² vanishes, so only the drift decides whether a trajectory leaves. The drift points outward somewhere exactly when the rate of information gain exceeds the total dephasing rate. A physical detector has η at most one, so its trajectories stay in the ball; η of three, an unphysical detector, takes them out. Positivity of the conditional state is the statement that η cannot exceed one.

## 9. When the readout record fixes the population

**Model.** A qubit with its integrated readout Q, integrated in Stratonovich form along one smooth noise path.

**Question.** Does the integrated readout Q alone determine z on every trajectory?

**Figure.** The Itô equation for the difference between arctanh z and the scaled Q; then z (solid) against the value predicted from Q (dashed), with and without drive Ω and relaxation γ₁.

**Physics.** arctanh z is half the log-likelihood ratio of the two states, and Bayes' rule adds the integrated record to it. The equation has zero noise for all rates and zero drift without drive and relaxation, so then the record alone determines z and the curves coincide. Drive and relaxation change z without entering the record, and the curves separate.

## 10. A trapped particle watched continuously

**Model.** A harmonic oscillator under continuous position measurement. The conditional state stays Gaussian, and its covariances obey a deterministic Riccati equation, solved in closed form.

**Question.** Can measurement squeeze Σ_x below the zero-point variance?

**Figure.** Σ_x in units of the zero-point variance for three measurement strengths μ, with the closed-form steady state dashed; panels η of 1, 1/2 and 1/5; each title gives the squeezing condition.

**Physics.** Position measurement reduces Σ_x through information gain, while backaction heats the momentum, and the trap rotates Σ_p into Σ_x. At η of one the conditional state is squeezed for every strength. With η below one the heating is unchanged but the information gain is reduced, so squeezing needs μ above a threshold that grows as η falls.

## 11. A watched Gaussian: the center wanders, the widths do not

**Model.** The five Gaussian moments ⟨x⟩, ⟨p⟩, Σ_x, C and Σ_p of the conditional state, derived from the canonical commutation relation with NonCommutativeAlgebra and sampled.

**Question.** Do the covariances fluctuate along a trajectory like the means?

**Table.** Drift and noise of each moment: the noise of all three covariances is zero, and every steady state has purity equal to the square root of η.

**Figure.** ⟨x⟩ for four trajectories (top) and their purity (bottom), at η of 1/5 and 1.

**Physics.** For a Gaussian state the conditional variance after a measurement depends on the measurement precision, not on the outcome. The means therefore diffuse and the trajectories fan out, while the covariances follow one deterministic curve and the purity curves coincide. The steady-state purity is set by the detector efficiency alone.

## 12. The quantum Kalman filter

**Model.** The conditional wavefunction on a 256-point grid under position measurement, and a Kalman filter that propagates five Gaussian moments, both driven by the same record.

**Question.** Can the Kalman filter track ⟨x⟩ of the full conditional state?

**Figure.** ⟨x⟩ from the wavefunction (solid) and from the filter (dashed); rows: strength μ; columns: quartic anharmonicity ε.

**Physics.** In a harmonic trap the conditional state is exactly Gaussian, so the five moments are a sufficient statistic and the filter is exact. Anharmonicity makes the state non-Gaussian and the filter's model wrong. Under weak measurement the estimate drifts; under strong measurement each innovation corrects it, and the error shows only at large ε.

## 13. How much cavity Wigner negativity the atomic record certifies

**Model.** An excited atom coupled to a lossy, detuned cavity mode (Jaynes-Cummings); only the atomic inversion ⟨σ_z⟩ is recorded. A coherent field is compared with its phase-randomized version, which has the same Fock-state distribution p_n.

**Question.** How much Wigner negativity of the field does the atomic record certify?

**Figure.** Top: ⟨σ_z⟩ and the negative volume of the Wigner function for the coherent field (solid) and its phase-randomized version (dashed); controls: κ and Δ in units of g. Bottom: Wigner functions of both fields at four times.

**Physics.** The Jaynes-Cummings interaction conserves the excitation number, so the inversion depends only on p_n, not on the coherences between Fock states. The two records are therefore identical. Phase averaging can only remove Wigner negativity, so every field with the same p_n has at least the negative volume of the phase-randomized field. That lower bound is what the record certifies. The Wigner snapshots show it is small: a thin negative ring, against the strong interference fringes of the coherent field.

## 14. Two qubits entangled by a measurement alone

**Model.** Two non-interacting qubits, initially in a product of plus states, read by a half-parity meter: its mean output distinguishes the two even-parity states but not the two odd-parity states. Each readout is applied with Korotkov's Bayesian update.

**Question.** Can the measurement alone entangle the qubits?

**Figure.** The concurrence of five records (lines), the closed form evaluated at each record's running mean (markers), and its ceiling (dashed); panels: dephasing rate γ of the odd-parity coherence.

**Physics.** A record whose running mean stays near zero projects onto the odd-parity subspace without distinguishing its two states, which leaves a Bell state. Records whose running mean drifts to either side project onto an even-parity product state and lose the entanglement. The meter also dephases the odd-parity coherence, so entanglement requires the measurement rate to exceed γ; once γ wins, no record entangles the pair.

## 15. Continuous syndrome tracking

**Model.** The three-qubit bit-flip code with stabilizers Z₁Z₂ and Z₂Z₃ measured continuously. The full 8×8 conditional state is compared with a filter for the four syndrome probabilities, both driven by the same records.

**Question.** Are four syndrome probabilities enough to track the syndrome?

**Figure.** The four syndrome probabilities from the filter; rows: efficiency η; columns: flip rate γ in units of κ. Each title gives the largest deviation from the full state at two time steps.

**Physics.** The records depend only on the syndrome, and bit flips only move probability between syndromes, so the four probabilities obey a closed equation. The deviation halves with the time step, so it is integration error. With high η and rare flips the filter jumps cleanly between syndromes; with low η or frequent flips it rarely settles.

## 16. The noise term in the syndrome equation of a watched code

**Model.** The equation for the four syndrome probabilities p_j, derived in the Pauli algebra from a general state with 64 coefficients and no matrices, then sampled.

**Question.** What term does the measurement record add to the equation for p_j?

**Table.** One row per syndrome: the flip term, the stabilizer backaction, which is zero, and the record terms multiplying each Wiener increment.

**Figure.** The code-space probability for three trajectories, the mean of 1000, and the flips-only solution (dashed); rows η, columns γ in units of κ.

**Physics.** The stabilizers commute with the syndrome projectors, so measuring them does not change the syndrome probabilities on average. The record term moves probability toward the syndromes the record favours, in proportion to p_j, and has zero mean. Single trajectories therefore jump between zero and one, sharper at higher η, while the ensemble mean follows the flips-only solution.

## 17. Markovian feedback that protects the code space

**Model.** The bit-flip code with each flip channel monitored by homodyne detection and the current fed back as a Hamiltonian, the Wiseman-Milburn feedback master equation.

**Question.** How well does the feedback stabilize an encoded state, and what does an efficiency below one cost?

**Figure.** The two identities the design rests on; then the code-space population (top) and fidelity (bottom) from the logical zero and logical plus states, for several η and without feedback (dashed).

**Physics.** The feedback is chosen so that the effective jump operator annihilates the code space: a flip is undone before it can act. At η of one every encoded state is a fixed point. With η below one the feedback also returns detector noise to the qubits: the plateau drops with η, and at low η the feedback does worse than none. From the logical zero the fidelity falls below the code-space population, because the noisy feedback causes a logical bit flip that the stabilizers cannot detect.

## Concluding remarks

A single trajectory shows what the ensemble average hides: oscillation, Zeno jumps, collapse and purification. Every measurement both gains information and causes backaction, and the efficiency η sets the balance; it decides positivity, the steady-state purity, the squeezing threshold and the cost of feedback. When the conditional state stays in a small family, a few numbers updated by Bayes' rule track it exactly. The Zeno crossover is an exceptional point of the Liouvillian. And the record is a resource: it creates entanglement and, fed back, protects encoded states.
