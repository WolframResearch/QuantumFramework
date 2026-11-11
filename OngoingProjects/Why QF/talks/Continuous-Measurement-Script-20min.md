# Talk script (20 min): Continuous Quantum Measurement in the Wolfram Language

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

## 5. Born's rule as a martingale

**Model.** A qubit under perfect σ_z measurement, starting in a pure state with ⟨σ_z⟩ equal to z₀; 400 trajectories per start.

**Question.** With what probability does a trajectory end at the upper pole?

**Figure.** Left: the fraction of trajectories ending at the upper pole against Born's prediction, with error bars. Right: four trajectories and the mean of 400. Controls: z₀ and drive.

**Physics.** Without a drive, z has zero drift, so it is a martingale: its mean stays at z₀. Every trajectory still ends at a pole. A bounded martingale that ends at the poles must split between them in proportion to z₀, which is exactly Born's rule, so the rule follows from the dynamics. A drive adds drift, the mean leaves z₀, and the trajectories keep jumping.

## 8. A monitored qubit stays in the Bloch ball

**Model.** A driven, relaxing qubit read out along σ_z, with measurement-induced dephasing Γ_d and readout noise split into an informational part Γ_CI and a backaction part Γ_BA; their sum relative to Γ_d sets η.

**Question.** Does every trajectory stay a valid density matrix, inside the Bloch ball?

**Figure.** The noise of r² on the sphere and the condition for outward drift; then r for η of 0.5, 1 and 3, with the region outside the ball shaded.

**Physics.** On the sphere the noise of r² vanishes, so only the drift decides whether a trajectory leaves. The drift points outward somewhere exactly when the rate of information gain exceeds the total dephasing rate. A physical detector has η at most one, so its trajectories stay in the ball; η of three, an unphysical detector, takes them out. Positivity of the conditional state is the statement that η cannot exceed one.

## 11. A watched Gaussian: the center wanders, the widths do not

**Model.** The five Gaussian moments ⟨x⟩, ⟨p⟩, Σ_x, C and Σ_p of the conditional state, derived from the canonical commutation relation with NonCommutativeAlgebra and sampled.

**Question.** Do the covariances fluctuate along a trajectory like the means?

**Table.** Drift and noise of each moment: the noise of all three covariances is zero, and every steady state has purity equal to the square root of η.

**Figure.** ⟨x⟩ for four trajectories (top) and their purity (bottom), at η of 1/5 and 1.

**Physics.** For a Gaussian state the conditional variance after a measurement depends on the measurement precision, not on the outcome. The means therefore diffuse and the trajectories fan out, while the covariances follow one deterministic curve and the purity curves coincide. The steady-state purity is set by the detector efficiency alone.

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

## 17. Markovian feedback that protects the code space

**Model.** The bit-flip code with each flip channel monitored by homodyne detection and the current fed back as a Hamiltonian, the Wiseman-Milburn feedback master equation.

**Question.** How well does the feedback stabilize an encoded state, and what does an efficiency below one cost?

**Figure.** The two identities the design rests on; then the code-space population (top) and fidelity (bottom) from the logical zero and logical plus states, for several η and without feedback (dashed).

**Physics.** The feedback is chosen so that the effective jump operator annihilates the code space: a flip is undone before it can act. At η of one every encoded state is a fixed point. With η below one the feedback also returns detector noise to the qubits: the plateau drops with η, and at low η the feedback does worse than none. From the logical zero the fidelity falls below the code-space population, because the noisy feedback causes a logical bit flip that the stabilizers cannot detect.
