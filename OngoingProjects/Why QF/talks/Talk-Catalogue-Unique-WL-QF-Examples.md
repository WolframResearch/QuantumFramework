# What WL and QF let a quantum researcher do: a catalogue for the talk

Built from a full read of 180 notebook sources (every `.md` with a built `.nb` beside it) in continuous-measurement, open-quantum-systems, QIS book, QUALink, TransmonLab, TN courses, QuantumFramework/OngoingProjects and TensorNetworks, plus two kernel probes on 2026-10-28. "Unique" means the usual stack (QuTiP, Qiskit, Cirq, Strawberry Fields, scqubits, ITensor, SymPy, SciPy) cannot do it or can only sample it numerically.

**One fact governs preparation.** md2nb sources carry no Out cells, so most results below exist only as claims in prose. The "Status" column says which ones have recorded output; every pick must be run on QF HEAD before the talk.

Status key: **K** = rerun in a kernel this session; **O** = output recorded in the source; **P** = claimed in prose only; **F** = fix needed first.

## 1. A symbolic parameter survives the whole quantum computation (QF)

| Example | What comes out | Source | Status |
|---|---|---|---|
| QAOA MaxCut landscape | closed form in (γ, θ), then `Maximize` gives the exact algebraic optimum | QF-Showcase Claim 1; QF-Webinar | K (triangle, 2.3 s); K4 is P, fix `QuantumOperator["ZZ"/2, #]` |
| Heisenberg pair, entanglement entropy vs t | binary entropy of cos²2Jt from `QuantumEvolve` + `QuantumPartialTrace` | QF-Showcase; Webinar | K (0.17 s) |
| Qubit Lindblad for every ρ₀ and both rates | full solution, then `Reduce` proves T2 ≤ 2T1 | QF-Showcase | P |
| Quantum Fisher information from `FubiniStudyMetricTensor` | n for product probes, n² for GHZ, exact | Courses/QuantumInDiscreteSpace/Answers/Part-17 §17.3 | P |
| Coherent information of amplitude damping | h₂((1−γ)p) − h₂(γp), optimum log₂(4/3) at p = 4/9 for γ = 1/4 (hand-checked) | Part-14 §14.4 | P |
| Data-processing inequality across the named channel library | `Reduce` finds no violating p, one line per channel | Part-14 §14.5 | P, may be slow (2322-leaf expression) |
| SU(2) double cover from a controlled-RX circuit | sin²(α/4), period 4π, from `mkCirc[α][]["Probabilities"]` | QIS book ch20 | P; hardware fringe recorded (§9) |
| Rabi and Landau-Zener propagators | closed form; LZ in ₁F₁/parabolic-cylinder functions with Stueckelberg oscillations | QF-Showcase | P |

## 2. Linear and Gaussian optics with symbolic photon numbers (QF `QuantumInterferometer`)

QuantumInterferometer is the newest and most active part of QF (20 kernel commits, Oct 17-27, ending at ed916377 with real lossy probabilities). It is the strongest single section for a physics audience: every item below runs in under a second, and the first five were rerun on HEAD today. Sources: Quantum Interferometer/QI-symbolic-showcase (symbolic reference), interferometer-showcase (14 experiments), Boson Sampling/why-boson-sampling-is-hard and boson-sampling.

| Example | What comes out | Status |
|---|---|---|
| Beam splitter, arbitrary n₁, n₂ → t₁, t₂ | `KroneckerDelta[n1+n2, t1+t2] WignerD[{(n1+n2)/2, (n1−n2)/2, (2t1−n1−n2)/2}, φ, 2θ, −φ]`: the splitter is an SU(2) rotation of Schwinger spin (n₁+n₂)/2 | K (0.06 s) |
| Hong-Ou-Mandel for any n | `WignerD[{n, 0, k−n}, 0, π/2, 0]`; the n = 1 coincidence is exactly 0 | K |
| Loss with symbolic η | (1−η)^{n−k} η^k Binomial[n, k] for symbolic n and k | K |
| Lossy HOM, both arms at η | probabilities fully symbolic in η: coincidence stays 0, bunched outputs η²/2 each, total 1 | K (new at HEAD) |
| Two squeezers into a balanced splitter | coincidence amplitude cancels for equal squeezing, peaks for opposite (two-mode squeezed vacuum) | K |
| Tritter (3-mode Fourier) on \|111⟩ | exact rationals: every 2+1 pattern is 0, \|111⟩ survives at 1/3, \|300⟩-type at 2/9 | K |
| Cascaded splitters | inactive sum with `Piecewise` photon-flow conditions, activated and matched to the concrete query | O |
| Two-mode squeezed vacuum | P(n, n) = tanh^{2n} r / cosh² r, normalizes for symbolic n; r → 0 is vacuum | O |
| KLM nonlinear sign gate with symbolic input c₀, c₁, c₂ | heralded output {c₀, c₁, −c₂}/2, success 1/4 for every input; `"ParameterSpec"` survives | O |
| Fock-space bridge | `qi["CircuitOperator", d]` turns the interferometer into a Fock-space `QuantumCircuitOperator` | O |
| Bosons, distinguishable particles, fermions | permanent vs determinant on one HOM splitter; cost of one amplitude set by the occupation pattern | O; benchmark runs 148 s at 28 photons, precompute |
| Also numeric | Bloch-Messiah, GBS form, Reck/Clements meshes, waveguide quantum walks, samples | O |

Known limits to say out loud: loss combined with squeezing and lossy sampling are not supported; one general odd-output WignerD amplitude does not simplify; `"Probability"` on truncated Fock states renormalizes and hides lost norm.

| Related Fock-space item | What comes out | Source | Status |
|---|---|---|---|
| Coherent state under the oscillator, α, ω, t symbolic | `uosc @ CoherentState[] == E^(-I ω t/2) CoherentState[][α E^(-I ω t)]` is True | QF-Showcase; QO book ch2 | P (closed form at fixed cutoff) |
| HOM with polarization overlap χ | cos⁴θ + sin⁴θ − ½ sin²2θ cos²χ, χ → 0 and π/2 limits | boson-sampling Part I | O |

## 2b. Open systems as formulas (QF `QuantumEvolve` with symbolic rates)

Source: Why QF/showcase/QF-OpenSystems-Showcase (found by the keyword search; prose-only results, physics hand-checked).

| Example | What comes out | Status |
|---|---|---|
| Entanglement sudden death | exact death time t* = (1/γ) ln(1/(1 − cot θ)) and its existence condition θ > π/4, from `Reduce` on the symbolic solution; matches `QuantumEntanglementMonotone` | P |
| Discrete-Wigner master equation for a qutrit | the Lindblad flow becomes a real Markov generator on 9 phase-space points; the strange state's negativity dies exactly at 2 ln 2/γ | P, F (drop the "no classical sampler" claim; mana citation is Veitch et al. 2014) |
| Generator, propagator, channel, circuit as one algebra | semigroup law in symbolic t₁, t₂; master equation proved equal to `QuantumChannel["AmplitudeDamping"[1 − e^{−γt}]]` on a symbolic ρ₀; Ramsey fringe ½(1 + e^{−(γ+4γφ)t/2} cos δt) | P |
| Superradiance and a decoherence-free subspace | collective decay at 1 ± r; the singlet is exactly dark at r = 1 | P (symbolic r fails in DSolve) |
| Coherent error corrected by record-controlled gates | e^{−iθX} digitized by syndrome measurement, corrected with no classical channel; syndrome weights cos²θ, sin²θ | numeric fidelity O; symbolic fidelity not yet run (Why QF/QF-Showcase-Claim6-Feedforward-Implementation) |

## 3. Theorems by quantifier elimination (WL `Resolve`/`Reduce` on physics)

| Example | What is proved | Source | Status |
|---|---|---|---|
| No two-qubit bound entanglement | identity 16 det ρ^{T_B} = Π(2λₖ − Σλ) in all 16 Pauli coefficients, plus `Resolve[ForAll …]` | Part-11 §11.12b | P |
| Robertson-Schrödinger for every qubit state and observable pair | slack = \|n̂₁×n̂₂\|²(1−\|r\|²); `Reduce` finds no violation | Part-06 §6.1 | P, 7-variable Reduce untested |
| Helstrom optimality | `Reduce[ForAll[t, …]]` | Part-12 §12.7 | P |
| Continuously measured trapped particle | conditional width below zero-point for every strength; purity √η; threshold μ* = 2√(1−η)/η^{3/2} | continuous-measurement/Continuous Space/measurement-in-a-potential, Part IV | P |
| Monitored transmon stays in the Bloch ball | no outward radial drift for any η ∈ [0,1] and nonnegative rates | Discrete Space/NDSolve vs Ito/NDSolve-WhenEvent-Reproduction-Report; Regularized-SDE-NDSolve-Simulation | P |

## 4. Stochastic calculus as symbolic objects (WL `ItoProcess`, `StratonovichProcess`)

| Example | What comes out | Source | Status |
|---|---|---|---|
| Ito's lemma on the collapse SME into spherical coordinates | exact dr equation ∝ (r²−1): pure stays pure, mixed purifies; then sampled with `RandomFunction` | QIS book ch06; intro-to-QIS-revised Part IV | P |
| Stratonovich → Ito for a noisy Hamiltonian | the −2γ dephasing drift appears symbolically (hand-checked) | QIS book ch05 | P |
| Laser phase diffusion | mean field e^{−κt}cos ωt exactly, Lorentzian by `FourierTransform`, non-stationary two-time correlation | ch05; intro-to-QIS | P |
| Quantum Bayesian (Korotkov) Kraus algebra at symbolic readings | non-commuting kicks differ by 2i sinh sinh σ_y; step is unimodular so purity is exact for every reading; POVM resolved by `Integrate`; half-parity meter concurrence ceiling tanh(δv²t/2s) | Discrete Space/Korotkov Bayesian/Korotkov-Quantum-Bayesian §10-20 | P |
| Quantum Kalman-Bucy filter as an `ItoProcess` | the purity invariant is kept without renormalization | Continuous Space/simulating-the-measurement-record | P |

## 5. Operator algebra and exact propagators

| Example | What comes out | Source | Status |
|---|---|---|---|
| Gaussian pulse, area theorem | `DSolveValue` on a 2×2 matrix unknown gives the Erf propagator, proved equal to exp(−iθσx/2) | QIS book ch05 | P |
| Ehrenfest and Mandelstam-Tamm for an arbitrary H(t) | identities with undetermined functions, no state solved for | Part-07 §7.5b, §7.7 | P |
| Exact Jaynes-Cummings propagator | interaction-picture `MatrixExp`, g and t symbolic, cutoff 50; dressed-state splitting 2g√(n+1) | QO book ch4 | O |
| Syndrome populations of continuous QEC | closed drift γ(1−4pⱼ) from `SolveAlways` over a general ρ | open-quantum-systems/syndrome-populations-from-the-algebra | P (plain WL) |
| Conditional-variance noise cancels under Gaussian closure | `NonCommutativeAlgebra` with [x,p] = iħ plus ordered Wick | bad essays/simulating-continuous-position-measurement | P, F (lives in "bad essays") |

## 6. Symbolic tensors (WL + TensorNetworks)

| Example | What comes out | Source | Status |
|---|---|---|---|
| Symmetry proofs at symbolic dimension n | antisymmetric·symmetric = 0, Ricci symmetric from Riemann symmetries, curl∘grad = 0, no components | TN courses/wolfram_doc_pages/ArrayDot_extended §6.1 | O, F (symbol reuse; fresh kernel) |
| Riemann identities derived, not assumed | `YoungProject` onto {2,2}; right-annihilator search in C[S₄] finds first Bianchi | TensorNetworks/…/SymmetrySubcontextTutorial; Pillar0 Lesson 2 | P |
| Keyl-Werner spectrum estimation | Schur-sampling probabilities via isotypic projectors (hand-checked) | TensorNetworks/…/young_tableaux | P |

## 7. Stabilizers

| Example | What comes out | Source | Status |
|---|---|---|---|
| Measurement outcomes kept as formal bits | `PauliStabilizer[…]["SymbolicMeasure", …]`: random outcomes are symbols, determined ones Boolean functions of them; Stim called in the same notebook as cross-check | QEC/verification-artifact/exact-stabilizer-codes-cross-checked | P (needs Python stim) |
| Logical error rate as a polynomial in a QF channel parameter | `QECLogicalErrorRate[QECCode["Repetition", 3], QuantumChannel["BitFlip"[p]]]` = 3p² − 2p³ | QEC/docs/Tutorials/StabilizerCodes | P, QECCore unshipped |

## 8. Hardware from the same notebook (QUALink, Arbel)

| Example | What it shows | Source | Status |
|---|---|---|---|
| Double-cover fringe on Arbel | measured Pₑ over 0-4π matches sin²(α/4); compiler emits `x $0` at 2π and an empty program for +I | QUALink/ClaudeNotes/…/SU2DoubleCover | O (recorded 2026-06-25) |
| Bell run against a depolarized null | χ² test finds a directional error channel | QUALink/ClaudeNotes/why symbolic/symbolic-wl-for-hardware | O |
| Pulse program as a symbolic object | a DRAG coefficient swept by ReplaceAll; the device IR queried and rewritten with `Cases`/`MapAt` | same | O |
| CZ derived symbolically, then calibrated on chip | {11, 02} block extracted from a symbolic two-transmon Duffing Hamiltonian | QUALink/…/Tutorials/ArbelNativeGatesAtPulseLevel | O for the algebra; hardware unrun (about 1700 jobs) |
| Active reset with feedback written as WL `If` | real-time conditional on an OPX1000 | QIS book ch18 | P (eval: false) |
| QF circuit → QASM → Arbel counts | full chain from a QF object | QUALink/…/Tutorials/QFCircuitsOnArbel | P, F (dephasing prose is wrong physics) |

## 9. Physics answers worth naming
- **Wigner negativity the atomic record can certify** (N17): only the phase-averaged field's negative volume; most of the negativity near half the revival is invisible to the record. QF, O.
- **Zeno as an exceptional point**: the symbolic Liouvillian spectrum {0, −2k, −k ± √(k²−Ω²)} (open-system-simulation-catalog, Zeno entry). Plain WL, P.

## Suggested running order (45 min, 12 items)
1. QAOA landscape and `Maximize` (§1). 2. Two-spin entropy, then differentiate it (§1). 3. Beam splitter for symbolic n, Wigner D, HOM for any n (§2). 4. Loss law, then lossy HOM symbolic in η (§2). 5. Entanglement sudden death time and its θ > π/4 condition (§2b). 6. QFI n vs n² (§1). 7. No two-qubit bound entanglement (§3). 8. Ito's lemma on the collapse SME (§4). 9. Korotkov Kraus algebra at symbolic readings (§4). 10. Squeezing below zero-point for every strength (§3). 11. Formal measurement bits in the stabilizer engine (§7). 12. Double cover: closed form, compiler, then the recorded Arbel fringe (§1, §8).

## Keywords for the search
The search ran over all 1,824 `.md`, `.wl`, `.wls`, `.m` files in the eight folders. Each keyword marks a capability the usual stack lacks; the count is files hit.

| Capability | Keywords | Files |
|---|---|---|
| Theorems by quantifier elimination | `Resolve[ForAll`, `Reduce[ForAll`, `ForAll[`, `SolveAlways[` | 23 |
| Symbolic stochastic calculus | `ItoProcess[`, `StratonovichProcess[`, `ItoProcessChangeVariables` | 16 |
| Non-commutative and bosonic algebra | `NonCommutativeAlgebra`, `BosonicNormalOrder`, `ExpBosonOrder`, `BosonicBCH`, `ZassenhausTerms`, `MagnusTerms`, `DysonTerms` | 22 |
| Symbolic optics and phase space | `QuantumInterferometer[`, `WignerD[`, `WignerFunction[`, `HusimiQFunction[`, `SOrderedFunction[`, `Hafnian`, `Permanent[` | 16 |
| Metrology | `FubiniStudyMetricTensor`, `QuantumFisher`, `Heisenberg limit`, `Cramér-Rao` | 12 |
| Exact matrix propagators | `DSolveValue[… Matrices[`, `∈ Matrices` | 10 |
| Symbolic tensors | `TensorReduce[`, `YoungProject[`, `ArraySymbol[`, `SymbolicIdentityArray` | 64 |
| Stabilizer engine | `PauliStabilizer[`, `SymbolicMeasure` | 87 |
| Control theory on quantum moments | `RiccatiSolve[`, `LyapunovSolve[` | 4 |
| Exact optimization and limits | `Maximize[`, `Minimize[`, `Asymptotic[`, `FindSequenceFunction[` | 42 |
| Geometry | `Berry phase`, `geometric phase`, `exceptional point`, `double cover` | 44 |
| Entanglement | `QuantumEntanglementMonotone`, `Concurrence`, `LogNegativity`, `bound entangle` | 29 |
| Hardware | `QUAExecute`, `QUAConnect`, `IBMJobSubmit`, `QolabExperiment`, `Arbel` | 314 |

The search surfaced six files outside the 180 notebook sources that add to the catalogue: QF-OpenSystems-Showcase (§2b), QF-Showcase-Claim6-Feedforward-Implementation (§2b), exact-relations-from-monitored-sdes (heating law and moment-closure test from `NonCommutativeAlgebra`; fix its (λη)^{−1/4} vs λ^{−1/4}η^{−3/8} contradiction), the TensorNetworks symmetry demos (Schur-Weyl blocks of the Heisenberg chain; demo3's Keyl-Werner plan uses a pure product state and is wrong), Arbel-Geometric-Phase-Investigation (live hardware fringe, but a software phase gate, not a geometric phase), and VerifiableCircuitCompilation (192 interval proofs by `Resolve[ForAll …]` for a 1/√x circuit; the most unique WL capability found, but compilation arithmetic rather than physics).

The rest of the tables is the handout.

## Fix before showing
- Every P item: run on QF HEAD and record output and time; several run long (DPI `Reduce`, 7-variable Robertson `Reduce`, K4 FullSimplify).
- intro-to-QIS Lindblad identity binds `ρ` but uses `r` (cell will not return True).
- ch20 / SU2DoubleCover: "visibility ≈ Bell fidelity²" contradicts the recorded run (0.983² vs V = 0.73).
- QF defects to avoid on stage (Part-11, Part-14): symbolic `"EntanglementEntropy"` silently returns 0 bits; `QuantumDistance` "RelativeEntropy" gives S(ρ) for pure σ and can swap arguments (not reproduced on QF 2.1.3, 2026-11-10: S(I/2‖|0⟩) = ∞, S(|0⟩‖I/2) = 1 bit, S(|+⟩‖|0⟩) = ∞, S(ψ‖ψ) = 0, all correct; recheck before citing); `QuantumDistance[a, b, "Fidelity"]` is the distance 1 − √F (|0⟩,|+⟩ gives 1 − 1/√2), by design; for the fidelity use `QuantumSimilarity[a, b, "Fidelity"]`, which returns √F.

## Leave out
TransmonLab tutorials and CPB reference pages (their printed numbers contradict the reference outputs: charge dispersion, α at EJ/EC = 20, g, a ground energy below −EJ). DRAG (all numeric). The named-circuit tour (the Fujitsu talk). bell-unentangled-photons (circular Schmidt check). Quantum wake (unnormalized, no analysis). The QF community Ito stepper (η and V·V terms wrong). Archived drafts. Plain-WL open-systems catalog entries other than Zeno, unless the point is WL rather than QF.
