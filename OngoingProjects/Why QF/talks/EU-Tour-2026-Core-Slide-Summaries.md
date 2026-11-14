# EU Tour 2026 core deck: slide summaries

### Correlated decay: bright, dim and dark states

1. **Core message.** When two atoms share the same light field, how fast one decays depends on the state of the pair. The symmetric shared excitation radiates faster than a lone atom (superradiance), the antisymmetric one slower (subradiance), and at perfect correlation the antisymmetric state does not radiate at all: it is a dark state.

2. **Physical setup.** Two identical two-level emitters with decay rate γ couple to common field modes. The correlation r between 0 and 1 is how much of their emission goes into shared modes. In free space r depends on separation and dipole orientation and approaches 1 when the emitters are much closer than a wavelength (the Dicke limit). Realizations include trapped-ion pairs at adjustable distance, atoms in a shared cavity or in front of a mirror, and superconducting qubits on one waveguide, where a spacing of a whole number of wavelengths gives r close to 1. The coherent dipole-dipole exchange only shifts the energies of |T⟩ and |S⟩, so it does not change the excitation number and is dropped.

3. **What the computation gives.** QuantumEvolve solves the master equation exactly for symbolic r, and the kernel returns closed forms for all four initial states (t in units of 1/γ): |T⟩ decays as e^{-(1+r)t} and |S⟩ as e^{-(1-r)t}; |eg⟩, half singlet and half triplet, decays as [e^{-(1-r)t} + e^{-(1+r)t}]/2; |ee⟩ per excitation decays as [(1+r)² e^{-(1+r)t} + (1-r)² e^{-(1-r)t} − 4r² e^{-2t}] / (2(1−r²)). At r = 0 all four curves fall on the lone-atom curve e^{-t}, a check that independent atoms behave as single atoms. At r = 1, |eg⟩ levels off at 1/2: the triplet half radiates and the singlet half is trapped forever, and that flat line is the dark state. |ee⟩ per excitation becomes (1+t)e^{-2t}, a two-step cascade |ee⟩ → |T⟩ → |gg⟩ with each step at the doubled rate 2γ. The |ee⟩ curve starts with the lone-atom slope and only pulls ahead later: with both atoms excited there is no shared phase yet, so superradiance from |ee⟩ builds up as the decay creates the correlation. The shading shows which side of the lone-atom curve each state lies on.

4. **Immediate applications.** At r = 1 a qubit stored as α|gg⟩ + β|S⟩ is shielded from shared emission: both states are dark, so the pair forms the simplest decoherence-free subspace, the principle behind atom-array and waveguide quantum memories. Related uses that need further modelling: superradiant states for faster, brighter emission (readout, narrow-line superradiant lasers), and neighbouring qubits that share a bath, whose correlated errors matter for error-correction decoders that assume independent errors.

### When a noisy circuit matches the master equation

1. **Core message.** The usual way device simulators model an idle period, one gate followed by an amplitude-damping and a phase-damping channel, is exact only when the Hamiltonian commutes with the noise. With no drive (or a drive whose effect the readout cannot see) one gate plus two channels reproduces the Lindblad equation for every parameter value. With a drive off resonance the circuit must be cut into short slices, and how short depends on the drive, the detuning and the dephasing together.

2. **Physical setup.** A Ramsey experiment on one qubit: Hadamard, a free wait of length t, Hadamard, read P(e). During the wait the qubit is detuned by δ from the reference frequency, relaxes at rate γ₁ (T₁ decay), dephases at rate γ_φ, and may feel a residual drive Ω (a leftover tone, crosstalk from a neighbouring line, or an imperfectly switched-off pulse). This is the standard T₂* calibration on superconducting qubits, trapped ions and spin qubits.

3. **What the computation gives.** The master equation is solved exactly with QuantumEvolve and compared with the circuit (𝒫_μ ∘ 𝒜_λ ∘ 𝒰_{t/m})^m for m = 1 and m = 8, with λ = 1 − e^{−γ₁t/m} and μ = 1 − e^{−4γ_φ t/m} chosen so each channel matches its Lindblad term exactly over one slice. The circuit diagrams show one slice and eight slices.
   - Without drive, dots sit on the line in every panel: the rotation about z, the decay and the dephasing all commute, so one gate and two channels are exact. The fringe oscillates about 1/2 at frequency δ and its envelope decays at 1/T₂ = γ₁/2 + 2γ_φ.
   - With drive Ω and detuning δ ≠ 0 the one-step circuit misses badly (in the shown panel it settles near 1/2 while the master equation settles near 0.6 to 0.7, because the driven steady state is not at the equator); eight slices come close but not all the way.
   - At δ = 0 the routes still give the same fringe even with the drive on. The reason: the second Hadamard turns the x component of the Bloch vector into the measured population, and a drive about x leaves x untouched, so x decays at the same rate in both routes.
   - Dephasing shrinks the swing of the fringe and the one-step gap, but turning it on widens the gap left by eight slices.

4. **Immediate applications.** Deciding when a gate-plus-channel noise model of a simulator is trustworthy for idling qubits that still feel a drive (crosstalk, spectator tones, always-on couplings), and how finely to slice the idle time when it is not. It also tells an experimentalist fitting Ramsey data which rates the fringe can reveal: without drive only 1/T₂ and δ, not γ₁ separately.

### Fastest relaxation at the Zeno exceptional point

1. **Core message.** A driven qubit under continuous σ_z measurement relaxes fastest at one particular measurement strength, the exceptional point where the Liouvillian's two Rabi eigenvalues and their eigenvectors coincide. Weaker measurement leaves damped Rabi oscillation; stronger measurement freezes the qubit (Zeno effect). Detuning splits, moves and finally removes this point from the real axis: for |Δ/Ω| > 1/(2√2) no measurement strength reaches it.

2. **Physical setup.** A qubit driven at Rabi frequency Ω, detuned by Δ, while a detector measures σ_z with strength k, averaged over the detector record so the state obeys a Lindblad equation with dephasing operator σ_z. In practice: a transmon read out dispersively through a cavity while a Rabi tone plays, with k set by the readout power; also trapped ions or atoms with a weak fluorescence probe. Superconducting-qubit experiments have used exactly this setup to reach Liouvillian exceptional points.

3. **What the computation gives.** The Liouvillian comes straight from QuantumOperator["Liouvillian"], and Reduce on the discriminant of its characteristic polynomial (steady-state root removed) returns the exact condition for an exceptional point at some k > 0: −1/(2√2) ≤ Δ/Ω ≤ 1/(2√2).
   - On resonance the Rabi pair is λ = −k ± √(k² − Ω²). The relaxation rate (bottom row, blue) rises linearly to a sharp cusp at k = Ω, where the oscillation frequency (orange) drops to zero; beyond it the slow rate falls toward Ω²/(2k), the Zeno slowdown. The top row shows the single zero of the discriminant at k = Ω on the real axis.
   - At Δ/Ω = 1/4 the zero splits into two real ones: the oscillation stops between them and restarts beyond the second, now at a frequency set by the detuning.
   - At Δ/Ω = 1/(2√2) the two zeros merge, and the relaxation peak is at its sharpest and highest of the four panels.
   - At Δ/Ω = 1 the zeros leave the real axis as a complex-conjugate pair, the oscillation never stops, and the cusp becomes a smooth hump.

4. **Immediate applications.** Choosing the measurement strength that makes the Liouvillian gap, the slowest decay rate, as large as possible, so the qubit approaches its steady state fastest at long times. That steady state is the maximally mixed state, so this is not reset or cooling. The optimum is a sharp exceptional point only for small detuning; at larger detuning a smooth optimum remains. Exceptional points are also proposed for sensing and for chiral state transfer, and this slide gives the exact parameter window in which a single measured qubit has one.

### Entanglement sudden death

1. **Core message.** Local decay can destroy entanglement completely at a finite time, even though every single-qubit population decays only exponentially. Whether that happens depends on the initial state, not just on the noise: it needs both qubits to decay and more weight on |ee⟩ than on |gg⟩.

2. **Physical setup.** Two qubits, A and B, prepared in cos θ |gg⟩ + sin θ |ee⟩, each coupled to its own independent zero-temperature bath, so each loses energy by spontaneous emission at its own rate γ_A, γ_B with no correlated decay. Realized with pairs of photons passing through controllable amplitude-damping setups, with trapped ions or atoms, and with superconducting qubits, each relaxing through its own T₁.

3. **What the computation gives.** QuantumEntanglementMonotone gives the concurrence of the evolved state for general θ, rates and t:
   C(t) = 2 e^{−(γ_A+γ_B)t/2} sin θ · max(0, cos θ − sin θ √((1−e^{−γ_A t})(1−e^{−γ_B t}))).
   For 0 ≤ θ ≤ π/2 it reaches zero, and stays there, when (1−e^{−γ_A t})(1−e^{−γ_B t}) = cot²θ.
   - In the panel γ_B = 0 every curve only fades: the product is zero, so it never reaches cot²θ.
   - In the panels γ_B/γ_A = 3/10 and 1, the curves θ = π/3 and 5π/12 hit zero at the marked dots, earlier for larger θ and earlier for faster B. The curves θ = π/8 and π/4 (dashed) fade without dying.
   The mechanism: decay feeds population into |eg⟩ and |ge⟩ incoherently. When |ee⟩ starts large, that mixed population soon outweighs the remaining |gg⟩ ↔ |ee⟩ coherence, and the state becomes separable.

4. **Immediate applications.** Entanglement distribution and storage (quantum repeaters, networked memories): a protocol that uses the pair, such as distillation or teleportation, has to finish before the death time, not merely within a few T₁. Within this family of states, the same formula says that θ ≤ π/4, with |gg⟩ dominant, makes the pair fade gradually instead of dying.

### Watched qubits pin and jump

1. **Core message.** A single run of a continuously measured, driven qubit looks nothing like the average. On resonance and with a perfect detector, under weak watching each run keeps Rabi-oscillating, and the average decays only because the runs drift out of phase. Under strong watching each run sits at a pole and makes rare, sudden jumps to the other one (quantum Zeno effect with quantum jumps); the average shows only a slow decay.

2. **Physical setup.** The same driven qubit as the exceptional-point slide (Rabi Ω, detuning Δ, σ_z measured at strength k), but now each run is conditioned on its own detector record, with detection efficiency η. This is weak continuous measurement as done with a transmon read out through a cavity and a near quantum-limited amplifier; the record is the amplified output signal.

3. **What the computation gives.** The Itô stochastic master equation is written as Bloch-vector equations and sampled with ItoProcess / RandomFunction; three runs are plotted against the Lindblad average from QuantumEvolve.
   - At k/Ω = 0.1 and η = 1 each run swings between the poles with a full ±1 amplitude, while the red average is a damped oscillation: it decays only because the runs lose phase with each other.
   - At η = 0.3 the runs no longer reach the poles. A detector that misses part of the signal leaves each run mixed.
   - On resonance, as k grows past Ω the oscillation stops. At k/Ω = 10 the runs stay near z = 1 with brief spikes and an occasional full jump to the other pole, and stronger watching makes the jumps rarer.
   - Detuning shrinks the Rabi swing to a fast wiggle near the starting pole. Under weak watching the record then carries the runs slowly between the poles; under strong watching they look as they do on resonance.

4. **Immediate applications.** Real-time tracking of a qubit for measurement-based feedback (stabilizing a Rabi oscillation or a state), catching and reversing a quantum jump mid-flight, and Zeno-type protection of a state by frequent measurement. The efficiency η sets how faithfully the record tracks the state.

### A monitored qubit stays in the Bloch ball: Itô's lemma on the sphere

1. **Core message.** A stochastic master equation for a monitored qubit keeps every run a valid density matrix exactly when the detector learns no more than it disturbs: η Γ_d ≤ Γ_d + γ_φ. For any real detector (η ≤ 1) this holds automatically, whatever the drive and relaxation. A model with η Γ_d > Γ_d + γ_φ pushes the Bloch vector outside the ball.

2. **Physical setup.** A Rabi-driven qubit with relaxation γ₁ toward z = 1 and extra dephasing γ_φ, read out continuously along σ_z, as in dispersive readout of a transmon with a phase-sensitive amplifier. The detector dephases the qubit at rate Γ_d. Its noise enters in two ways: an informational part Γ_CI that updates z, and a back-action part Γ_BA that kicks the phase. The efficiency is η = (Γ_CI + Γ_BA)/(2Γ_d).

3. **What the computation gives.** ItoProcess supplies the Itô derivative of r² = x² + y² + z², and evaluated on the sphere r = 1 it gives two exact results:
   - The noise term of r² is zero on the sphere: a measurement kick only slides a pure state along the surface, so only the drift can push it in or out.
   - Reduce finds that the drift points outward somewhere exactly when η Γ_d > Γ_d + γ_φ. Neither γ₁ nor Ω appears.
   In the three panels:
   - At η = 0.5 the runs wander inside the ball.
   - At η = 1 they get close to r = 1 (the record nearly purifies the state) but stay inside.
   - At η = 3 they leave the ball within about one unit of time.
   The mechanism: the detector's dephasing pulls inward at rate Γ_d (plus γ_φ), its information pushes toward purity at a rate set by ηΓ_d, and the sign of the difference decides.

4. **Immediate applications.** A consistency check on fitted detector models. If an estimate of Γ_CI, Γ_BA and Γ_d from readout data gives η > 1, the detector model is unphysical, and once η Γ_d > Γ_d + γ_φ its trajectories leave the Bloch ball. The same check applies to any quantum-trajectory or filtering code used for feedback, which must never produce r > 1. Extra dephasing that carries no information only widens the safety margin.

### Randomized benchmarking: same decay, different floors

1. **Core message.** Randomized benchmarking reports one number, the decay p of survival with sequence length, and that number sees only the average gate fidelity of the noise. Two very different errors with the same fidelity, incoherent amplitude damping and a coherent over-rotation, give the same decay. What tells them apart is the floor the survival settles to and, above all, the spread from one random sequence to the next.

2. **Physical setup.** Standard single-qubit randomized benchmarking. Start in |0⟩, apply m Cliffords drawn at random from the 24-element group, then the one Clifford that undoes them, and measure. Every gate carries the same error. The error is either energy relaxation during the gate (amplitude damping with probability γ, a T₁ effect) or a miscalibrated pulse amplitude (a small extra X rotation by ε), with ε chosen so both have the same average fidelity F. This is routine on superconducting, trapped-ion and spin-qubit devices.

3. **What the computation gives.** Each noise is reduced to its action on the Bloch vector, r ↦ M r + c, read off QuantumChannel and QuantumOperator acting on Bloch states. Averaging over the Clifford rotations of the cube leaves a shrink factor p = Tr M/3 = 2F − 1, which gives the exact survival P(m) = (1 + c_z + M_zz p^m)/2. Equal fidelity means cos ε = √(1−γ) − γ/2. The simulation runs 200 random sequences at each length and compares them with this formula.
   - The two mean curves (cyan for damping, magenta for over-rotation) each lie on their dashed exact curve. The two curves differ, but they decay at the same rate p and settle to different floors.
   - Damping settles at the floor (1+γ)/2, slightly above 1/2, because it keeps pushing toward |0⟩. The over-rotation settles at 1/2.
   - The error bars give the spread between single sequences: small for damping, very large for the over-rotation.
   - On the sphere and in the side cuts, damping shrinks the runs into the ball and they stay clustered near the axis, drifting toward a point above the centre. The over-rotation keeps every run pure, on the surface, and the random gates scatter them over the whole sphere, so only their average reaches the centre.

4. **Immediate applications.** Reading randomized-benchmarking data correctly. The reported error rate cannot say whether the error is coherent, but coherent errors build up much faster in long circuits and in error correction, so the difference matters. A large variance between sequences is a cheap extra sign of a coherent error; a floor of 1/2 alone is not, since depolarizing noise gives it too. If the variance points to one, the likely next step is recalibrating the pulse, not improving T₁.

### Dephasing: GHZ gains nothing, squeezing a bounded factor

1. **Core message.** Without noise, an entangled GHZ probe of n qubits measures a frequency n times better (in Fisher information) than n independent qubits: this is the Heisenberg limit. With independent dephasing on each qubit, that advantage disappears completely once the Ramsey time is chosen well. A spin-squeezed probe still helps, but only by a fixed factor of at most e ≈ 2.7, however many qubits are used.

2. **Physical setup.** This is Ramsey interferometry, the measurement behind atomic clocks and atomic or spin magnetometers. n two-level systems (atoms, ions, NV centres) pick up a phase ωt and each loses coherence independently at rate γ (pure dephasing, T₂ noise). A fixed total time T is split into T/t repeated runs of length t. The precision is set by the quantum Fisher information F_Q per unit time, F_Q/t, and the run length t is free to optimize.

3. **What the computation gives.** QuantumEvolve solves one dephasing qubit symbolically, and the quantum Fisher information follows from its eigen-decomposition. A GHZ state behaves as one qubit whose phase and dephasing rate are both multiplied by n, so symbolic n comes for free:
   - F_Q(GHZ) = n²t² e^{-4nγt}
   - F_Q(product) = n t² e^{-4γt}

   Both reach the same best value of F_Q/t, n/(4eγ). GHZ reaches it at t = 1/(4nγ) and the product state at t = 1/(4γ). The plot shows information per qubit per unit time:
   - Every GHZ curve peaks at the same height as the product curve, shifted to runs n times shorter. The phase grows n times faster but the coherence dies n times faster.
   - The dashed squeezed-state curves (one-axis twisting, optimal squeezing, a lower bound on F_Q) rise above the product peak, more so at n = 100 than at n = 10.
   - All the squeezed curves stay under the dotted bound t e^{-4γt}/(1 − e^{-4γt}), which holds for any state at large n. That bound approaches 1/4 as t goes to 0, exactly e times the product peak 1/(4e). So e is the largest gain in Fisher information any entangled state can give under this noise, a factor √e in standard deviation.

4. **Immediate applications.** This sets the realistic target for quantum-enhanced clocks and magnetometers whose noise is local dephasing:
   - Use a moderately squeezed state with a well-chosen Ramsey time, which can gain up to a factor e.
   - Do not build large GHZ states for this, since they cost a great deal and gain nothing.
   - GHZ states pay off only when the run length is capped by something other than the qubits' own dephasing (for example a noisy local oscillator), or when the noise is not independent.

   The same few lines give the optimal run length for any n directly.

### Teleportation with the classical channel as record wires

1. **Core message.** Teleportation moves an unknown qubit state from Alice to Bob using one shared Bell pair and two classical bits. Both of Bob's corrections are needed to transmit an arbitrary unknown state. Each missing correction erases a specific part of the state's Bloch vector, so the failure has a recognizable signature.

2. **Physical setup.** There are three qubits: Alice's input, plus a Bell pair split between Alice and Bob. Alice performs a Bell-basis measurement (CNOT, then H, then measure both qubits). Bob applies X if m₂ = 1 and Z if m₁ = 1. In hardware the corrections are classical feed-forward: mid-circuit measurement, then a gate chosen by the result. This is done on trapped ions, superconducting processors and photonic links. Here the two outcomes are kept as record qubits, and the corrections become gates controlled by those records. By the deferred-measurement principle this gives the same statistics as real-time feed-forward.

3. **What the computation gives.** Bob's Bloch vector is computed in closed form for a general input (β, φ) in each of four cases:
   - **Both corrections:** Bob's Bloch vector equals Alice's exactly. On the sphere every blue dot sits on its orange input.
   - **X only:** the Z that should have been applied is applied at random, which flips the x and y components. Only the z component survives, so the dots fall onto the vertical axis.
   - **Z only:** a random X flips y and z, so only the x component survives and the dots fall onto the x axis.
   - **No corrections:** Bob holds the maximally mixed state, a dot at the centre of the sphere, for every input.

   Without the classical bits, Bob's state carries no information about the input. This is why teleportation cannot send information faster than light.

4. **Immediate applications.**
   - **Testing feed-forward in hardware.** The pattern of surviving Bloch components shows directly which classical correction is missing.
   - **Building blocks.** The same circuit with record-controlled gates is the core of gate teleportation, of entanglement swapping in quantum repeaters, and of measurement-based computation. In all of these, the routing of the corrections is the part most easily got wrong.

### Noise caps Grover's speedup

1. **Core message.** Grover's search finds one marked item among N in about (π/4)√N oracle calls, against about N for guessing. A small phase-flip error on every qubit at every iteration erodes this. The cost stays near its noiseless value only while a whole search collects well under one error, so the error per qubit per iteration that keeps the full quadratic speedup is about 1/(n√N), which shrinks as the search grows. Beyond that scale the advantage shrinks gradually; it does not vanish at once.

2. **Physical setup.** Grover search runs on n = log₂N qubits, with one oracle accepting x* = 5. Each qubit is dephased before every iteration: Z is applied with probability p. Dephasing is the dominant error in many qubit platforms. A run of k iterations is checked with one extra oracle call and repeated until it succeeds. The cost is the expected number of oracle calls, (k+1)/P_k, minimized over k and compared with N for guessing.

3. **What the computation gives.** QF's Grover iteration is used as a matrix. Phase flips on every qubit damp each density-matrix element ρ_xy by (1−2p)^{d(x,y)}, where d(x,y) is the Hamming distance between x and y. This gives the exact noisy evolution.
   - **Left plot (N = 256).** Without noise, the success probability oscillates and reaches 1 at about 12 iterations. With noise each peak is lower and arrives earlier. At p = 0.03 the best success is only about 0.23, so the search must be repeated several times.
   - **Right plot.** The expected cost stays flat near the noiseless value until p reaches about 1/(n√N): about 0.02 for N = 64 and about 0.008 for N = 256. Above that it climbs steeply toward the classical cost. It crosses the dashed guessing line only at heavy noise, p around 0.2, and only because every run must make at least one iteration; running none is guessing itself, which costs N.

4. **Immediate applications.** This is a dephasing budget for amplitude amplification, the subroutine used in Grover search, amplitude estimation and many quantum speedups, with the oracle and diffuser themselves taken as perfect. The dephasing that keeps the full quadratic speedup falls like 1/(n√N), which grows harder to meet as the problem grows. So keeping the full speedup at scale needs error-corrected logical qubits. The same model gives, for a given hardware error rate, the largest search that still beats classical search, and the best number of iterations to run before restarting.

### Shor's algorithm factors 15 and 21

1. **Core message.** Shor's algorithm reduces factoring N to finding the order r of a number a modulo N: the smallest r with a^r ≡ 1 (mod N). Phase estimation produces one probability distribution whose peaks sit at k/r. Continued fractions read r off those peaks, and gcd(a^{r/2} ± 1, N) gives the factors, except for even r with a^{r/2} ≡ −1 or for odd r.

2. **Physical setup.** There is a work register of 4 qubits (for 15) or 5 qubits (for 21), holding |1⟩. A genuine modular multiplier U_a|y⟩ = |ay mod N⟩ is built as a permutation matrix. A register of t = 8 or 9 control qubits runs phase estimation (controlled powers of U_a, then the inverse quantum Fourier transform). Factoring 15 and 21 is the standard small-scale demonstration on NMR, photonic, trapped-ion and superconducting hardware, though many such demonstrations simplify the multiplier using knowledge of the answer.

3. **What the computation gives.** The circuit diagram shows the order-finding circuit, drawn with 4 of the 8 control qubits. The six panels show the exact output distributions:
   - **N = 15**, where 256/r is an integer, so every peak is exact:
     - a = 2 has r = 4, with four equal peaks at 0, 1/4, 1/2, 3/4 and factors 3 and 5.
     - a = 11 has r = 2, with factors 5 and 3.
     - a = 14 ≡ −1 has r = 2, but the gcds are 1 and 15, so it fails and another a must be chosen.
   - **N = 21**, where 512/r is not an integer for a = 2 and a = 4:
     - a = 2 has r = 6, and the peaks at k/6 spread slightly over neighbouring outcomes. Continued fractions of the outcomes above the cutoff still return denominators whose least common multiple is 6, giving factors 7 and 3.
     - a = 4 has odd order 3, so no factors and a restart.
     - a = 13 has r = 2, and 512/2 is an integer, so its peaks are exact, giving factors 3 and 7.
   - Peaks at 0 and at k/r with k sharing a factor with r give only a divisor of r (open circles). This is why the order is taken as the least common multiple over several outcomes.

4. **Immediate applications.**
   - **A reference for hardware demonstrations.** These exact distributions are what an honest order-finding experiment on 15 or 21 should produce with a full modular multiplier. Measured histograms can be compared against them; whether a demonstration built the answer into its circuit can be judged only from the circuit itself.
   - **Checking the classical post-processing.** The same code tests the continued-fraction and least-common-multiple steps, and the choice of the number of control qubits t, before scaling up.

### QAOA on MaxCut: the exact landscape and a certified optimum

1. **Core message.** For small graphs, the QAOA expected cut is an exact trigonometric function of the angles. Its one-layer global maximum can therefore be found and certified exactly, not just estimated by an optimizer. How many layers reach the true maximum cut depends on how much of the graph each layer can see: rings need depth that grows with their length.

2. **Physical setup.** This is MaxCut on four graphs: the triangle, K₄, and rings of 5 and 6 vertices, with one qubit per vertex. Each QAOA layer applies ZZ rotations by angle γ on every edge (the cost) and X rotations by θ on every qubit (the mixer), starting from |+⟩^n. The score is the expected number of cut edges. QAOA is run on superconducting processors and, for related graph problems, on Rydberg-atom arrays.

3. **What the computation gives.** E(γ, θ) is derived in closed form, and the density plot shows the full one-layer landscape as a fraction of the maximum cut. The one-layer maxima, certified by Maximize, are:
   - **triangle:** 2, which is the maximum cut;
   - **K₄:** 3.70, an algebraic number (root of the polynomial in the title), below the maximum cut of 4;
   - **both rings:** exactly 3/4 of their edges (3.75 of 5 and 4.5 of 6).

   With more layers, local searches from random starting angles find:
   - **K₄ at p = 2:** 4, the maximum cut;
   - **ring of 5 at p = 2:** 4, the maximum cut;
   - **ring of 6:** 5 at p = 2 (that is 5/6), and the full cut of 6 at p = 3.

   This matches the rule that a p-layer circuit sees only 2p+1 consecutive edges of a ring. Until it sees the whole ring it cuts (2p+1)/(2p+2) of the edges. A search that reaches the maximum cut to numerical precision has found the best score, because no state can cut more.

4. **Immediate applications.**
   - **Exact benchmarks for QAOA.** Certified optimal angles and values let one test whether classical angle optimizers get stuck in local maxima, and show the ideal value that noisy hardware runs should reach.
   - **Depth planning.** The ring result gives a concrete rule: on rings, QAOA needs depth comparable to the ring's length before its expected cut reaches the maximum cut. This guides the choice of circuit depth on real devices.

### An infinite Ising chain from 2d − 1 sites: the light cone of a Trotter circuit

1. **Core message.** In a circuit of nearest-neighbour gates, a local observable after d Trotter steps depends only on a finite window of sites, its light cone. For the transverse-field Ising chain started from all spins up, the window is 2d − 1 sites wide, because the first ZZ layer only adds a phase. So the magnetization of a site in an infinitely long chain, after d Trotter steps, is computed exactly on 2d − 1 qubits.

2. **Physical setup.** This is a quench of the transverse-field Ising chain, H = J Σ ZᵢZᵢ₊₁ + h Σ Xᵢ, from the product state |0…0⟩. It is approximated by first-order Trotter steps of size Jδt = 0.2: a layer of ZZ rotations on neighbouring pairs, then X rotations. This is the standard digital quantum simulation run on superconducting and trapped-ion processors, including large kicked-Ising experiments whose results are checked with this kind of light-cone reduction.

3. **What the computation gives.** ⟨Z_c⟩ is computed on windows of m = 1, 3, …, 2d + 3 sites. In the plot, every window with m ≥ 2d − 1 (blue) gives the same value, the infinite-chain value. Smaller windows (red) cut the light cone and give a wrong value: for d = 2 and h/J = 1, about 0.74 against about 0.70 for one site.
   - The circuit drawing colours each gate: purple gates affect Z_c, grey gates cancel against their inverses, and orange gates of the first step only multiply |0…0⟩ by a phase.
   - Keeping only the purple gates (the star) gives the same value with far fewer gates: 6 against 26 for d = 2.
   - Stepping through d and h/J, the magnetization falls faster as the transverse field grows, and at h/J = 2 it changes sign within the times shown.

4. **Immediate applications.**
   - **Fewer qubits for the infinite chain.** Simulating the thermodynamic limit of a quench needs only about 2d qubits for d Trotter steps, and the circuit can be trimmed to its light cone, which cuts both qubit count and gate count on hardware.
   - **Exact references for large experiments.** Large digital-simulation experiments can be checked site by site against this exact light-cone value, for example when assessing error mitigation.

### A $1/\sqrt{x}$ phase circuit: its error, rounding included, proved for every input

1. **Core message.** A quantum arithmetic circuit can be given an error bound that is a proof, not a sampling estimate. The bound covers every input value and every way the fixed-point registers can round, and it is decided exactly over the real numbers.

2. **Physical / experimental setup.** The circuit is a phase oracle, |x⟩ → e^{2πi/√x}|x⟩, of the kind used in fault-tolerant algorithms. The computation runs in several stages:
   - **Normalization.** Leading-one detection and a barrel shift write x = 4^k m with 1 ≤ m < 4.
   - **Table lookup.** A QROM table look-up loads the three coefficients of a quadratic for the piece of [1, 4) that m falls in. There are 3·2^r pieces, and every coefficient is a w-bit number.
   - **Evaluation.** Horner's rule evaluates the quadratic with fixed-point multiply-adds that round after each step.
   - **Phase and cleanup.** A shift by 2^{-k} restores the scale, a phase-gradient step writes e^{2πi y}, and every work register is uncomputed.

   The diagram on the slide is this block structure, not a gate-level circuit.

3. **What the computation gives.** The top-left panel shows that the six quadratics (r = 1) lie on top of 1/√m. The bottom-left panel shows the problem: at the samples alone, the relative error is known only at those points and could hide a larger value in between.

   On each interval the claim becomes a polynomial inequality in three real variables: the position t and the two rounding errors δ₁ and δ₂. Resolve decides it exactly, with rational numbers and no rounding in the proof. The six-panel grid (r = 1, 2 and w = 8, 12, 16) shows the worst case over rounding as a solid line under each interval's proved bound η (dashed). Every interval in every panel is proved.

   The panels show which error dominates as the word length grows:
   - **At 8 bits** coefficient and step rounding set the error. It jumps from piece to piece and barely improves when the pieces are halved.
   - **At 12 bits** rounding and interpolation error are comparable on the wide pieces, and rounding takes over again on the narrow ones.
   - **At 16 bits** the smooth interpolation ripple sets the error where it is largest, on the wide pieces and at small m. Where it has shrunk, on the narrow pieces and at larger m, rounding sets the error again.

   The model treats t as exact and the shift by 2^{-k} as free of error.

4. **Immediate applications.**
   - Fixed-point 1/√x and 1/x phase oracles are a large part of the cost of fault-tolerant chemistry and physics simulation, for example the Coulomb potential in grid-based (first-quantized) simulation, and of linear-system algorithms.
   - A proved bound per (r, w) lets a resource estimate choose the smallest register width and table size that meets a target accuracy, with no safety factor taken on trust.
   - The same Resolve step certifies the rounding of any piecewise-polynomial block whose error condition is a polynomial inequality, as it is for 1/√x and 1/x. For arcsine (amplitude loading) or exponentials, the gap between the polynomial and the function first needs a separately proved remainder bound.

### A thousand-qubit GHZ state: entanglement entropy as a rank over F₂

1. **Core message.** For stabilizer states, entanglement entropy is a rank over the two-element field. It can be computed exactly for a thousand qubits, with no state vector. The GHZ and cluster states also show that entanglement can depend on the shape of a cut, not only on its size.

2. **Physical / experimental setup.** The GHZ state comes from H followed by a chain of CNOTs. The 1D cluster state comes from H on every qubit followed by CZ on each neighbouring pair. Both are the standard benchmark entangled states on trapped-ion and superconducting processors. Cluster states are also the resource for measurement-based quantum computing, including in photonic platforms.

   Each state is stored as an n×2n binary tableau, one row per stabilizer generator. The entropy of a set A is S(A) = rank G_A − |A|, where G_A keeps the columns of A. A Hadamard on qubits 2 to n turns GHZ into the star graph state, and the cluster state is the chain graph state. In both, S(A) is the rank over F₂ of the block of the adjacency matrix linking A to the rest, whose ones are the edges that cross the cut.

3. **What the computation gives.** The interactive n = 20 panel shows the cut edges in red and S(A) read off the tableau:
   - **Star (GHZ).** Every cut edge ends on the centre, so S(A) = 1 however many edges are cut.
   - **Chain (cluster).** A block is crossed by one edge, so S(A) = 1. A scattered set is crossed by many separate edges and has a large S(A), though not always one bit per cut edge (the example has 8 edges cut and S = 6).

   At n = 1000 the log plot confirms the pattern. GHZ carries exactly one bit across every cut, block or random, of every size: its entanglement is one shared parity. The cluster state carries one bit across a block boundary, but a random set carries hundreds of bits, peaking near k = n/2.

4. **Immediate applications.**
   - Checking the entanglement structure of large Clifford circuits and graph-state resources exactly, at sizes where state-vector methods are impossible.
   - Choosing qubit partitions for distributed quantum computing and for cutting circuits between modules, where the shape of the cut sets the communication cost.
   - Predicting how much entanglement a tensor-network or MPS simulation must carry for a given qubit ordering.

### Repetition code: the error rate of the decoded bit is a polynomial in the flip probability

1. **Core message.** A longer repetition code helps only when each qubit flips with probability below one half, and the decoded error rate is an exact polynomial in that probability. Here it is computed from n + 1 stabilizer tableaux instead of 2^n error patterns.

2. **Physical / experimental setup.** One bit is stored as |0…0⟩ in an odd number n of qubits. Each qubit goes through an independent bit-flip channel with probability p, and the bit is decoded by majority vote. Because the vote depends only on how many qubits flipped, one tableau per flip count is enough.

   This is the memory experiment run as a first error-correction benchmark on superconducting and trapped-ion processors. It is also the outer code of cat-qubit architectures, where the hardware suppresses one error type and a repetition code handles the other.

3. **What the computation gives.** P_bit(n) is computed exactly up to n = 101 and equals the regularized incomplete beta function I_p((n+1)/2, (n+1)/2).
   - **Linear plot.** All curves cross at p = 1/2 and sharpen toward a step there as n grows. Reduce proves that the longer code is better exactly for 0 < p < 1/2.
   - **Log-log plot.** At small p each code is a straight line with slope (n + 1)/2, the smallest number of flips that fools the vote. So P_bit ≈ C(n, (n+1)/2) p^{(n+1)/2}, and n = 3 already has no linear term.
   - **Closed form.** The difference between neighbouring codes factors as P_bit(n+2) − P_bit(n) = C(2m+1, m)[p(1−p)]^{m+1}(2p − 1), whose sign is that of 2p − 1. This is the simplest exact threshold theorem.

4. **Immediate applications.**
   - The exact small-p exponent and prefactor are the ideal reference for repetition-code memory experiments: the suppression per added pair of qubits, and the crossing at p = 1/2 for independent flips read out perfectly.
   - In biased-noise hardware (cat qubits, where the repetition code corrects the remaining phase flips, and dephasing-dominated spin qubits) this polynomial is the logical error rate of the protected error type per round; with a per-round noise model it sets how many qubits a target lifetime needs.
   - The idea of one tableau per error class is the starting point for computing exact rates of larger stabilizer codes.

### Maximum-likelihood decoding: exact logical error rates of three codes

1. **Core message.** For these three codes, a decoder that picks the most probable class of errors beats one that picks the lightest error only when the noise is biased: under depolarizing noise the two give the same exact polynomial for the five-qubit, Steane and Shor codes. Under Z-biased noise they differ, and maximum likelihood lowers the leading coefficient.

2. **Physical / experimental setup.** One logical qubit is stored in the [[5,1,3]], [[7,1,3]] Steane or [[9,1,3]] Shor code. Each physical qubit goes through an independent Pauli channel, and the checks are read without error, the code-capacity model. Two noise models are used:
   - depolarizing, with p_X = p_Y = p_Z = p/3;
   - Z-biased, with p_X = p_Y = p/10 and p_Z = p, typical of dephasing-dominated qubits.

   The check matrices (X part in blue, Z part in orange) show where degeneracy comes from. The Shor code's weight-two Z checks make many light errors equivalent.

3. **What the computation gives.** Summing over all 4^n errors gives P_L as an exact polynomial in p for each decoder (the table). The maximum-likelihood decoder fixes its choice per syndrome at small p, and that choice stays the most probable up to p = 0.2.
   - **Depolarizing noise.** The two decoders agree term by term, with leading coefficients 10p², 13p² and 49/3 p² for the five-qubit, Shor and Steane codes. For these codes the class the lightest error points to is also the most probable one.
   - **Z-biased Steane code.** Maximum likelihood gives 651/25 p², against 273/10 p² for minimum weight. Some two-error syndromes can come from an XZ pair or from an XY pair, ten times rarer, in a different class. Counting weight alone, the minimum-weight decoder sometimes undoes the rare pair.
   - **Top row.** P_L/p² shows that the gap opens at p = 0 and that the maximum-likelihood curve never lies above the other.
   - **Bottom row.** P_L is set against the unencoded error probability p_X + p_Y + p_Z, and the shading marks where encoding helps. The pseudothresholds under depolarizing noise are 0.14, 0.1 and 0.081 for the five-qubit, Shor and Steane codes, in the order of their p² coefficients. Under bias, maximum likelihood moves the pseudothreshold up slightly.

4. **Immediate applications.**
   - For hardware with biased noise, the exact polynomials say how much is gained by switching from a lookup or matching decoder to a maximum-likelihood one (a tensor-network or table decoder) before building it.
   - Exact P_L(p) is a test oracle for Monte Carlo decoder benchmarks at small p, where sampling is too slow to resolve the leading coefficient.
   - The comparison of pseudothresholds guides the choice of code for small early demonstrations.

### A surface code measured by noisy gates: the exact logical error rate

1. **Core message.** For a distance-3 surface code with faulty CNOTs and faulty readout, the failure probability of one round (the logical qubit flipped, or a detector pattern the decoder refuses) is an exact polynomial in p_C and p_M. It has no linear term, so every single fault is corrected. A decoder that also knows pairs of faults lowers the rate several times and moves the crossing with ε_L = p, in this model, to larger p.

2. **Physical / experimental setup.** The code has nine data qubits on a 3×3 grid and eight ancillas, one per check: four X-type and four Z-type checks of weight two or four. This is the 17-qubit layout used in surface-code experiments on superconducting processors.

   One round of check measurement uses Hadamards around the X checks and CNOTs in a fixed order. The noise is circuit-level:
   - each CNOT is followed by one of the 15 two-qubit Paulis with probability p_C/15;
   - each readout flips with probability p_M;
   - single-qubit gates and idling are noiseless.

   A perfect final measurement follows the round, so there are 16 detectors plus the two logical bits.

3. **What the computation gives.**
   - **Fault table.** A stabilizer simulation of X and Z inserted after each CNOT gives every single fault's detector pattern and logical flips. The table shows how many of the 368 single faults fire 0 to 4 detectors and how many flip the logical qubit. Faults that fire no detector never flip it. Together with the missing linear term in ε_L, which rules out two single faults with one pattern and different logical effects, this is the circuit-level statement that this CNOT order keeps distance 3.
   - **Exact rate.** Faults are independent and their effects add bit by bit, so a Walsh-Hadamard transform over the 2^18 detector-and-logical patterns gives ε_L exactly. Each CNOT contributes a factor 1 or 1 − 16p_C/15, and each readout 1 or 1 − 2p_M.
   - **Single-fault decoder.** ε_L ≈ (2213/15) p_C p_M + (48326/225) p_C² + 28 p_M². The 28 p_M² is just "two of eight readings wrong", since this decoder refuses any pattern it does not know.
   - **Pair decoder.** ε_L ≈ p_C p_M/5 + (2438/45) p_C². It has no p_M² term, because two wrong readings are always recognized. Readout errors alone first defeat it at three.
   - **Log-log plot.** The pair decoder sits several times lower, and its crossing with ε_L = p moves from about 0.0027 to about 0.015.

4. **Immediate applications.**
   - Exact circuit-level rates are a reference for validating fast samplers and decoders (Stim plus matching or lookup) at small p, where Monte Carlo cannot resolve the leading coefficients.
   - The fault table shows directly whether a chosen CNOT schedule produces hook errors that lower the distance. That is how gate order is chosen for surface-code hardware.
   - The polynomial says how much readout fidelity compared with CNOT fidelity matters, and how much a better decoder buys, for a 17-qubit device before running it.

### Continuous syndrome tracking: four numbers follow the records

1. **Core message.** When the stabilizers of the bit-flip code are measured weakly and continuously, four probabilities (one per syndrome) are enough to track the syndrome from the two noisy records. They agree with the full 8×8 conditional state up to the integration error.

2. **Physical / experimental setup.** Three qubits carry the bit-flip code, and each flips at rate γ. The parities Z₁Z₂ and Z₂Z₃ are monitored continuously at rate κ with detector efficiency η. Each record is the parity value plus white noise.

   This is continuous parity measurement as done with dispersive readout of superconducting qubits through a shared resonator. The full model is the stochastic master equation for the conditional state. The filter is the classical Wonham filter for a hidden four-state Markov chain. Flips relax the sector probabilities toward 1/4 at rate 4γ (every nontrivial syndrome mode decays at that rate), and each new pair of readings multiplies each probability by its Gaussian likelihood.

3. **What the computation gives.**
   - **Left plot.** The true syndrome history is shaded and the line is the probability the filter gives it, with the fraction of time the filter's best guess is right:
     - With a good detector and rare flips, the filter jumps to the new syndrome soon after each flip and is almost always right.
     - With a poor detector and frequent flips it rarely settles before the next flip.
   - **Right plot.** The largest gap between the filter and the sector populations of the full 8×8 state, driven by the same record (a positivity-preserving Rouchon-Ralph step), falls in proportion to the time step. So the difference is integration error, and the four numbers follow the sector populations exactly.

   They close on themselves because both stabilizers are diagonal in the sectors and each flip only moves probability between sectors. For the same reason, records simulated from a classical flip history have the same statistics as the quantum model.

4. **Immediate applications.**
   - Real-time syndrome tracking for continuous error correction. Four numbers per step fit easily on an FPGA controller, unlike the full conditional state. They give the current syndrome only: errors with the same syndrome can differ by a logical flip, so deciding the recovery needs the syndrome together with the logical value, eight numbers in this code.
   - The detection delay and accuracy as functions of η and γ/κ show how measurement strength and detector efficiency set how quickly and how reliably a flip is seen.
   - The same reduction applies to any code whose monitored checks commute and whose errors map syndromes to syndromes, including the larger codes that feedback protocols target.

### Markovian feedback that protects the code space, and what the detector costs

1. **Core message.** Feeding a measured current straight back to the qubits can freeze an encoded state exactly, but only if the detector catches every photon. Any inefficiency feeds detector noise back onto the qubits, and at low efficiency the feedback does more harm than no feedback at all.

2. **Physical setup.** Three qubits hold the three-qubit bit-flip code, |0_L⟩ = |000⟩ and |1_L⟩ = |111⟩, and each qubit flips at rate γ. Each flip channel c_i = √γ X_i is engineered so that its output is monitored by homodyne detection of efficiency η; this is a bit-flip channel, not ordinary spontaneous emission, whose jump operator lowers the qubit. Each detector current drives a correcting Hamiltonian I_i F_i at once, with no filtering and no delay (Wiseman-Milburn direct feedback). F_i is built from the code projector so that it pushes a flipped state back into the code space. The closest laboratory setting is continuous error correction on superconducting qubits read out through quantum-limited amplifiers, though there the measurement is of parities rather than of the flip channels.

3. **What the computation gives.** The plot title shows that the two operator identities, c_i F_i + F_i c_i = 0 and (c_i − i F_i) Π_code = 0, hold to exactly zero. Together they make every code state a fixed point when η = 1: the blue lines stay at 1 for both starting states. For η < 1 the extra term (1−η)/η D[F_i] acts like a noise source, and its strength grows without bound as η drops. The code-space population then settles on a plateau that falls with η. At η = 1/10 and 1/50 the population drops much faster than with no feedback, though it ends near the same value. Without feedback the three qubits flip independently until each is equally likely up or down, so the code-space population tends to 2/8 = 1/4. The bottom row separates two kinds of failure:
   - From |0_L⟩ the fidelity keeps sinking below the population. The noisy loop also turns |000⟩ into |111⟩, a logical error the code cannot see.
   - From |+_L⟩ fidelity and population coincide, because |+_L⟩ is unchanged by that logical flip.

   Without feedback the |+_L⟩ fidelity also tends to 1/4 rather than 1/8, because bit flips leave XXX unchanged.

4. **Immediate applications.** The model sets the detector efficiency a continuous-feedback error-correction scheme needs before it beats doing nothing, which is a design target for measurement chains on superconducting and photonic qubits. It also shows that one should check both code-space population and logical fidelity. A loop can keep the state inside the code space while steadily scrambling the logical information, and population alone does not reveal that.

### Compiling a rotation from H and T gates: each T gate buys a fixed factor in precision

1. **Core message.** In a fault-tolerant computer an arbitrary single-qubit rotation must be built from Clifford gates and T gates. The precision you get is set by the number of T gates. Each T gate doubles the number of distinct circuits, and counting them sets the best possible rate, an error falling as 2^{−2T/3}, about one decimal digit for every five T gates. The enumeration shows the circuits reach that rate.

2. **Physical setup.** This is the compilation step of any error-corrected processor (surface-code superconducting, trapped-ion or neutral-atom machines). Clifford gates H and S are cheap there. Each T gate consumes a magic state that must be distilled, at a cost of many physical qubits and many rounds of error correction. A target rotation U must be approximated by a word in H, S and T. The error used is ε = 1 − |Tr W†U|/2, which ignores global phase and grows as the square of the leftover rotation angle.

3. **What the computation gives.** Every Matsumoto-Amano normal form up to T count 16 is enumerated, 196,606 forms times 24 Cliffords, 4,718,544 circuits in all. For each of 200 Haar-random targets the best error at each T count is kept.
   - **Single targets (grey):** each curve falls in irregular steps. A target can be far from every circuit with one T count and close to one with the next.
   - **Across targets:** the median (blue) follows the orange counting line 2^{−2T/3} over more than three decades, and the 10% to 90% band is a few T gates wide.

   The slope comes from geometry. Single-qubit rotations up to phase form a three-dimensional space, so 2^T points can come no closer than about 2^{−T/3} to a typical target, and ε grows as the square of that distance. Reaching this floor also needs the points to spread evenly over the rotations, which the median curve confirms. Setting a factor of 10 equal to 2^{2T/3} gives T = (3/2) log₂10 ≈ 5. This matches the known optimal scaling of approximate synthesis methods: about 3 log₂(1/δ) T gates for an operator-norm error δ, with ε ∝ δ².

4. **Immediate applications.** The result gives a T-count budget per rotation, which feeds directly into resource estimates for fault-tolerant algorithms full of small rotations: phase estimation for chemistry, Trotterized simulation, and quantum Fourier transforms. Combined with the cost of one distilled magic state, it turns a target accuracy into qubits and runtime. The spread across targets also shows how much is gained by choosing rotation angles, or the size of the leftover error, wisely.

### A beam splitter with symbolic photon numbers

1. **Core message.** A lossless beam splitter is a rotation of a spin. N photons in two modes form a spin N/2 (the Schwinger representation), and the splitter rotates it by 2θ. Every photon-counting probability is therefore a squared Wigner D-function, obtained once for all photon numbers.

2. **Physical setup.** Two optical modes, n₁ and n₂ photons in, meet at a beam splitter with mixing angle θ (θ = π/4 is the 50:50 splitter) and phase φ. Photon-number-resolving detectors count m photons in one output. The other output must then hold n₁ + n₂ − m, since a lossless splitter conserves photon number. This is the basic element of every linear-optics experiment, from bulk optics to integrated photonic chips.

3. **What the computation gives.** QuantumInterferometer returns the amplitude symbolically in n₁, n₂, m, θ and φ. Its square is |d^{(n₁+n₂)/2}_{(n₁−n₂)/2, m−(n₁+n₂)/2}(2θ)|². The phase φ only multiplies the D-function by a phase, so it drops out of every counting probability. The three panels evaluate this one formula:
   - **|1,1⟩:** at θ = π/4 the probability of m = 1 (one photon in each output) goes to zero. This is the Hong-Ou-Mandel effect; the photons leave together, as |2,0⟩ or |0,2⟩, each with probability 1/2.
   - **|2,2⟩:** at the balanced point every odd m is zero, and the probabilities are 3/8, 1/4 and 3/8 for m = 0, 2 and 4.
   - **|1,3⟩:** at the balanced point m = 2 is zero, so the splitter never leaves two photons in each output. This zero comes from a node of the D-function, not from a symmetry of the input.

4. **Immediate applications.** The closed form gives output statistics for any Fock input in one formula. Uses include:
   - designing and calibrating photon-number-resolving interferometers;
   - predicting multiphoton interference in boson-sampling and photonic-chip experiments;
   - Holland-Burnett style phase sensing with twin Fock states |N,N⟩, where the zeros at odd m are the output signature of the balanced splitter; sub-shot-noise phase estimation also needs a phase-encoding and readout protocol.

### Hong-Ou-Mandel with partially distinguishable photons

1. **Core message.** Two-photon interference measures how identical the photons are. The coincidence rate is a weighted mix of the curve for identical photons and the curve for distinguishable ones. The weight is the squared overlap of their internal states, here their polarizations.

2. **Physical setup.** Two single photons enter the two ports of a beam splitter (angle θ, phase φ). One is horizontally polarized; a waveplate turns the other's polarization by χ. The detectors count photons but not polarization. Each port therefore carries two modes, H and V, and the one splitter acts on the H pair and the V pair alike. This is the standard Hong-Ou-Mandel test used to qualify single-photon sources such as quantum dots, down-conversion sources and trapped-atom emitters.

3. **What the computation gives.** The four-mode network gives the coincidence probability in closed form:
   - P_coinc = cos²2θ cos²χ + (cos⁴θ + sin⁴θ) sin²χ, independent of φ.
   - **χ = 0 (identical photons):** the curve reaches zero at the balanced splitter, a perfect Hong-Ou-Mandel dip.
   - **χ = π/2 (orthogonal photons):** the curve is the classical one, with value 1/2 at the balanced splitter.
   - **In between:** at θ = π/4, P_coinc = sin²χ/2, so the dip visibility 1 − 2P_coinc equals cos²χ, the squared overlap of the two polarization states.

   The plotted family shows the dip filling in smoothly as χ grows, while the end points θ = 0 and π/2 stay at 1 for every χ.

4. **Immediate applications.** The depth of the dip is the standard figure of merit for photon indistinguishability. It is what one quotes for a single-photon source before using it in boson sampling, photonic quantum computing or entanglement swapping. The same weighted-sum structure covers any internal degree of freedom (timing, frequency or spatial mode mismatch), so for ideal single photons the formula converts a measured visibility directly into an overlap (for mixed internal states, into Tr(ρ₁ρ₂)).

### The nonlinear sign gate of Knill, Laflamme and Milburn

1. **Core message.** Linear optics alone cannot make photons interact. Adding one helper photon and keeping only runs with the right detector clicks produces an effective nonlinearity: it flips the sign of the two-photon amplitude, which no phase shifter can do. The gate succeeds with a probability that does not depend on the input.

2. **Physical setup.** A signal mode carries up to two photons, c₀|0⟩ + c₁|1⟩ + c₂|2⟩. Two auxiliary modes carry one photon and the vacuum. A three-mode network of beam splitters and a phase shifter mixes them; this is the KLM transfer matrix from the slide. A run is kept only when the auxiliary detectors count one photon and zero photons. This is the building block of measurement-induced gates in photonic quantum computing on bulk or integrated optics.

3. **What the computation gives.** QuantumInterferometer decomposes W into three beam splitters and one phase shifter (the diagram) and gives the exact amplitudes ⟨m,1,0|U|m,1,0⟩.
   - **m = 0, 1, 2:** the amplitudes are 1/2, 1/2 and −1/2, so the kept signal is (c₀|0⟩ + c₁|1⟩ − c₂|2⟩)/2.
   - **What that means:** the sign of c₂ is flipped and the overall factor 1/2 is the same for every input. The gate therefore works with success probability (1/2)² = 1/4, whatever the input state.
   - **m = 3 to 5:** the amplitudes are unequal, (4√2 − 5)/2, (16√2 − 23)/2 and (56√2 − 79)/2, so the gate is valid only on inputs with at most two signal photons.

4. **Immediate applications.** Two such gates, one in each arm of an interferometer, make a controlled-Z between two photonic qubits with probability 1/16. This was the original route to universal photonic quantum computing. The exact amplitudes give the success probability and show what leaks in with higher photon numbers, which matters for multiphoton sources and for comparing with newer schemes based on fusion of small entangled states.

### A Schrödinger-cat state in phase space

1. **Core message.** A superposition of two coherent states shows interference fringes in phase space that dip below zero, which no classical state can do. The value of πW at the origin is exactly the photon-number parity, so it is fixed in closed form by the relative phase and the overlap of the two coherent states.

2. **Physical setup.** A single field mode, such as a microwave cavity mode or an optical mode, is prepared in N(|α⟩ + e^{iφ}|−α⟩). Here α is the real amplitude and φ the relative phase. In superconducting circuits and in cavity QED these cats are made routinely. Their Wigner function is measured by displacing the field and then measuring the photon-number parity with a coupled qubit or atom.

3. **What the computation gives.** QF builds the cat on 40 Fock levels and computes πW over phase space for α from 0.5 to 3 and φ = 0, π/2 and π. All panels share one color scale.
   - **Shape:** the two coherent lobes sit at x = ±√2 α. Between them, fringes along p oscillate with a spacing that shrinks as α grows, and they dip below zero; the odd cat is negative at the origin for every α, and at α = 1 the fringes are already negative while the lobes still overlap.
   - **The origin:** each panel's printed πW(0,0) matches the closed form ⟨Π⟩ = (e^{−2α²} + cos φ)/(1 + e^{−2α²} cos φ).
   - **φ = 0 and φ = π:** these are the even and odd cats, parity eigenstates with πW(0,0) = +1 and −1 for every α.
   - **φ = π/2:** the parity is e^{−2α²}, which goes to zero quickly with α. At α = 2 it is about 3×10⁻⁴, and the panel shows 0.000.
   - **Small α:** the cats look like the vacuum (even) or a single photon (odd).

4. **Immediate applications.** Cat states are the logical states of bosonic cat qubits. In the four-component cat code both codewords share one parity, so a lost photon flips parity and parity measurement is the error syndrome; in a two-component cat, parity is the logical value itself, and a lost photon is a logical error. The identity πW(0,0) = ⟨Π⟩ is the principle behind direct Wigner tomography in cavities. The closed-form parity gives the expected signal at any α and φ, and the fringe spacing shows how quickly larger cats become sensitive to photon loss and dephasing.

### GHZ and W: arbitrarily close, never connected

1. **Core message.** GHZ and W are the two inequivalent kinds of genuine three-qubit entanglement. No local operation, even one allowed to fail with some probability, turns one into the other, yet W can be approached as closely as you like from inside the GHZ class. W lies on the boundary of the GHZ class.

2. **Physical setup.** Three qubits are held by three separate parties. They could be photons, trapped ions or superconducting qubits, all platforms that prepare GHZ and W states routinely. Each party may apply any invertible local operation: a unitary or a local filter (a weak measurement whose successful outcome is kept). This is stochastic local operations and classical communication (SLOCC). The quantity that tells the classes apart is the Cayley hyperdeterminant of the 2×2×2 amplitude array. Its normalized form is the three-tangle τ₃, the part of the entanglement that no pair of qubits carries.

3. **What the computation gives.** All twelve entries of the three local maps are kept as symbols.
   - After the maps, the hyperdeterminant of GHZ is exactly (det A det B det C)²/4, never zero. The hyperdeterminant of W stays identically zero.
   - Since a local map only rescales the hyperdeterminant, the zero of W and the nonzero value of GHZ can never be exchanged. This proves the two classes are separate.
   - The left plot repeats the test on 300 random integer maps: every GHZ point lies on the diagonal and every W point sits at zero.
   - For the family made by the map A_ε on every qubit, the kernel gives τ₃ = 4ε²/(ε⁴+3ε²+3)² and infidelity with W equal to ε²(ε²+3)/(ε⁴+3ε²+3). Both vanish as ε² (as 4ε²/9 and ε²) when the determinant ε of each map goes to zero, which the right plot shows as two parallel lines of slope 2.

   So W is reached only in the limit where the map stops being invertible.

4. **Immediate applications.** For pure three-qubit states with genuine three-party entanglement, the three-tangle tells GHZ from W, a question pairwise measures cannot answer. It is also zero for product and biseparable states, so it is not a classifier on its own, and mixed lab states need more than this test. The result also warns that a prepared state can be arbitrarily close to W in fidelity while belonging to the GHZ class, so fidelity alone does not certify the class. Outside physics, the same fact (W has tensor rank 3 but can be approximated by rank-2 tensors) is the border-rank phenomenon used in fast matrix-multiplication algorithms.

### The AKLT ring: from a matrix-product state to an exact ground state

1. **Core message.** A spin-1 ring whose state is fixed by one small 2×2×3 tensor is an exact ground state of the AKLT Hamiltonian. Its spin correlations decay exponentially, falling by a factor of 3 per site. Both facts follow from contracting small closed tensor networks, without ever writing the 3ⁿ amplitudes.

2. **Physical setup.** A ring of n spin-1 sites is coupled to its neighbours by the AKLT Hamiltonian: a sum of projectors onto total spin 2 for each neighbouring pair, which equals a Heisenberg term plus a biquadratic term. Each spin 1 is built from two spin-1/2's, and neighbouring spin-1/2's form singlets. The bond dimension 2 of the matrix-product state is exactly that pairing. This is the textbook model of the Haldane phase of spin-1 chains. That phase is seen in quasi-one-dimensional nickel-chain magnets, and quantum simulators aim to prepare states in it.

3. **What the computation gives.** The energy network contracts to exactly 0 for n = 6, 12 and 24. Since the Hamiltonian is a sum of projectors, energy zero means an exact ground state.
   - The correlation ⟨S^z_1 S^z_{1+r}⟩ has magnitude (4/3)3^{-r} with alternating sign, a correlation length of 1/ln 3 ≈ 0.91 sites. The plot shows this as the dotted straight line.
   - The contracted points sit exactly on the ring formula. Near r = n/2 the ring formula bends up, because the correlation running the other way round the ring adds to it on an even ring. On a larger ring the bend moves out.
   - The network drawings show what is contracted: the ring of A tensors joined to its mirror ring, with the projector or the two S^z operators inserted in red.

   The cost grows linearly with n, so hundreds of sites would stay exact.

4. **Immediate applications.** The AKLT state is the reference state for symmetry-protected topological order: an exactly solvable ground state with a gap, and spin-1/2 edge modes on an open chain. It is used to check quantum simulators and numerical methods. AKLT-type states are also resources for measurement-based quantum computation, where single-site measurements drive the computation. The method shown (exact expectation values by closed-network contraction) is how matrix-product-state results are checked against exact answers.

### The entanglement profile of an Ising ring, and the bond dimension it needs

1. **Core message.** Away from the critical point, the entanglement of a block saturates (an area law). The bond dimension needed then stays flat over the sizes shown. At the critical point, the entropy grows logarithmically with the central charge c = 1/2 of the Ising critical point, and the needed bond dimension grows with the system.

2. **Physical setup.** A ring of n spins-1/2 sits in the transverse-field Ising model: ferromagnetic ZZ coupling J and transverse field h. There are three regimes: the ordered phase (h/J = 1/2), the quantum critical point (h/J = 1) and the paramagnet (h/J = 2). Experimental versions include Rydberg-atom arrays, trapped-ion chains and the quasi-one-dimensional magnet CoNb₂O₆. The ground state is found by exact diagonalization and written as an untruncated matrix-product state. At each cut, the Schmidt values give the entanglement entropy S(ℓ) and the bond dimension χ needed to keep the discarded weight below 10⁻⁶.

3. **What the computation gives.**
   - **Entropy for n = 16.** The ordered phase has S ≈ ln 2 at every cut: the finite ring's symmetric ground state is close to a cat state of all-up and all-down, which carries one bit. The paramagnet has S ≈ 0.17, nearly flat. At criticality S(ℓ) follows the Calabrese-Cardy profile (c/3) ln[(n/π) sin(πℓ/n)] + s₁, and the fit gives c = 0.51, close to the exact 1/2.
   - **Bond drawing.** The bond thicknesses show the same picture. In the gapped phases χ stops growing a few sites in from the ends (6 in the middle). At criticality it keeps growing toward the middle, up to 12.
   - **Size scan.** For rings of 8 to 18 sites, χ at the middle cut stays flat at 6 in both gapped phases. At criticality it rises steadily, from 8 to 13.

   The growth of χ at criticality is read off the computed Schmidt values, not derived from the entropy.

4. **Immediate applications.** Before simulating a device or a material with tensor networks, this tells how large the bond dimension must be. It is cheap away from criticality and grows with size near a phase transition. Measuring c from an entropy profile identifies the universality class of a transition, a test now possible in quantum simulators that measure entanglement entropy. The same numbers tell whether a matrix-product-state simulation of a given circuit or chain is under-resolved.

### Leakage of a fast π pulse on a transmon

1. **Core message.** A fast π pulse on a transmon leaves error in two places: population leaked to |2⟩ and above, and population left behind in |0⟩. DRAG (a second quadrature that follows the slope of the pulse) suppresses leakage by orders of magnitude. The error left over then mostly reflects a shift of the qubit frequency, which a small detuning or a different DRAG strength can correct.

2. **Physical setup.** The device is a standard transmon with E_J/E_C = 50 and E_C/h = 250 MHz. The TransmonLab spectrum gives ω₀₁/2π = 4.74 GHz, ω₁₂/2π = 4.45 GHz and anharmonicity α/2π = −0.287 GHz. The circuit drawing shows the drive line, the junction with its shunt capacitor, and the coupling to a readout resonator. The pulse model is:
   - the lowest five levels of an anharmonic oscillator;
   - a raised-cosine π pulse of length τ, with |α|τ = 5, 10 and 20, about 2.8, 5.5 and 11 ns;
   - a drive detuning δ = 0 or ±0.05|α| (about ±14 MHz);
   - DRAG strength λ = 0, 1/2 and 1.

   The model is solved in the rotating frame, after the rotating-wave approximation.

3. **What the computation gives.** The plots show leakage P(n ≥ 2) (solid) and P(0) (dashed) during the pulse.
   - **Leakage.** Leakage always rises through |1⟩, because the drive only couples neighbouring levels. Only the shortest pulse without DRAG keeps leaking to the end. Otherwise leakage peaks mid-pulse and partly returns.
   - **Effect of DRAG.** At |α|τ = 10 and δ = 0, the final leakage drops from 9.9% (λ = 0) to 2.7% (λ = 1/2) to 0.12% (λ = 1). At |α|τ = 20 with λ = 1 it is below 10⁻⁶.
   - **Population left in |0⟩.** This is the remaining error. At |α|τ = 10 and λ = 1 it is about 8%, against 0.3% at λ = 1/2. The |2⟩ level and the DRAG term shift the qubit transition in opposite directions, and near λ = 1/2 the two shifts nearly cancel.
   - **Detuning.** A detuning removes what is left: a negative one helps without DRAG, a positive one helps at λ = 1.

4. **Immediate applications.** This is the design problem for every single-qubit gate on transmon processors (IBM, Google, Rigetti and others). It shows how to choose pulse length, DRAG strength and drive frequency together, and it explains why calibrations tune both λ and a frequency offset. Leakage matters most for error correction, because a leaked qubit is not corrected by the code and spreads errors. A model like this sets how short gates can be before leakage dominates.

### Entanglement held by a pump on an IBM processor: predicted, then measured

1. **Core message.** A third qubit that is reset without being read can drive two qubits toward an entangled cluster state. On ibm_fez the pair ends well above the separable bound F = 1/2 after three pumped cycles. A noise model built only from the processor's published calibration predicts the result before the run, and gets the right ordering of the three stabilizers. It is optimistic, so the run loses something the calibration-based model does not capture.

2. **Physical setup.** The hardware is three neighbouring superconducting qubits on IBM's Heron processor ibm_fez (qubits 123, 122, 121), with the ancilla in the middle; this line had the highest predicted fidelity among the ten best-calibrated lines. Each half-cycle copies one stabilizer (X₁Z₂ or Z₁X₂) onto the ancilla, uses it to flip the pair back into the +1 branch when needed, and resets the ancilla unread. That reset pumps entropy out of the pair, the same principle as reservoir engineering or autonomous error correction. After three cycles the pair is rotated into one of three bases and read out. There are three jobs of 4096 shots each, and the fidelity with the cluster state is F = (1 + ⟨X₁Z₂⟩ + ⟨Z₁X₂⟩ + ⟨Y₁Y₂⟩)/4. The model compiles the circuit to IBM's native gates (RZ, √X, CZ, reset) and inserts calibrated noise channels: depolarizing noise after each gate, amplitude and phase damping from T₁ and T₂ during idle times (the pair idles through every ancilla reset), a reset error set by the readout error, and readout errors.

3. **What the computation gives.** The circuit drawings show one pump cycle, the three submitted jobs, and the noisy native-gate cycle as QF simulates it.
   - **Parities.** The bar chart compares three sets of values. The exact cluster state gives 1 for all three. The model predicted roughly 0.93, 0.96 and 0.92 for X₁Z₂, Z₁X₂ and Y₁Y₂. The device gave roughly 0.86, 0.93 and 0.85.
   - **Fidelity.** Measured F = 0.910 ± 0.003, against a prediction near 0.95. Both are far above 1/2, so the entanglement is certified.
   - **What the model gets right.** Z₁X₂, pumped last in each cycle, holds best, and the other two lose more.
   - **What it misses.** The gap between prediction and measurement is many times the shot noise, so it is real error the model misses. Model approximations and calibration drift can contribute; possible physical sources are crosstalk and ZZ coupling during the long ancilla resets, leakage, and reset imperfections beyond the readout figure.

4. **Immediate applications.** This is a working example of measurement-free entanglement pumping (showing stabilization would need a comparison run without the pump), a building block for autonomous error correction on hardware where mid-circuit readout is slow. The prediction-before-measurement procedure is a device-validation tool. The size and pattern of the gap between calibration model and hardware show how much error the calibration-based model misses and on which stabilizers; locating it in particular gates or idle periods needs further experiments. This tells where error budgets for error-correction experiments need work.

### The first Bianchi identity from the Riemann tensor's Young diagram

1. **Core message.** All of the Riemann tensor's symmetries, the first Bianchi identity included, are one statement: R lies in the image of the Young projector of shape {2,2} for the tableau with rows {1,3}, {2,4}. That projector enforces antisymmetry in each pair, symmetry under pair exchange and Bianchi, and nothing else.

2. **Physical setup.** There is no lab experiment here. The object is the curvature tensor of a d-dimensional spacetime or Riemannian manifold, as used in general relativity. The four index slots of a tensor are permuted by S₄. The Young projector for the tableau with rows {1,3}, {2,4} (symmetrize the rows, antisymmetrize the columns) picks out one irreducible representation of GL(d). The test works because in d = 4 the 24 permutations act faithfully on four-index tensors (Schur-Weyl duality), so an identity among these maps that holds there holds in every dimension.

3. **What the computation gives.** The projector and each relation are built as exact 256×256 matrices.
   - **The four Riemann relations.** Each relation matrix times the projector is exactly zero, so every tensor in the image obeys all four Riemann relations. A control relation that is not a Riemann symmetry gives 2520 nonzero entries, so the test can fail.
   - **Counting components.** For symbolic d, the projector's rank s_{2,2}(d) = d²(d²−1)/12 equals the count left by the pair symmetries and Bianchi, N(N+1)/2 − C(d,4) with N = d(d−1)/2. For example, it is 20 in four dimensions.
   - **The plot.** Direct counts for d = 2 to 6 (dots) sit on the closed forms (lines). The gap between "pair symmetries only" and the Riemann count is C(d,4): Bianchi removes nothing for d ≤ 3, one component in four dimensions, and grows quickly after that.

4. **Immediate applications.** Computer algebra for gravity and field theory uses this result to reduce curvature expressions to canonical form: projecting onto the Young shape replaces a list of hand-coded symmetry rules. Component counts like 20 in four dimensions (10 Weyl plus 10 Ricci) fix how many independent quantities a curvature measurement or a numerical-relativity code must carry. The same Schur-Weyl machinery is used in quantum information for the Schur transform, symmetric subspaces and permutation-invariant states, so the construction carries over directly to multi-qubit problems.

