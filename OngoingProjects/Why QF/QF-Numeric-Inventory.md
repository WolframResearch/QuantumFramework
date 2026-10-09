# QF numeric capability inventory

What Wolfram QuantumFramework can do numerically today, measured rather than claimed, and what stands in the way of each problem class ranked in the companion survey, [Quantum software numerics priorities](reports/Quantum%20software%20numerics%20priorities.md).

**Environment.** QF 2.1.1 working tree, Wolfram Language 15.0, 4 October 2026 (the P0 defects re-measured after their fixes on 5 October), 48 GB macOS machine. Swap was 97% full and other kernels were running during part of the session, so absolute timings carry roughly ±2× noise; ratios measured inside one kernel are reliable. Every number below comes from a script in [numeric-inventory-probes/](numeric-inventory-probes/), runnable from that folder against the repository's paclet. The per-process probes' result lines are in `results.log` (kernel message dumps removed); `batchA`, `batchB2`, `batchB3`, `batchD` and `t0` print to the terminal. Run probes one per kernel with `./runprobes.sh probe.wls`, which kills a kernel that outlives its time limit: `TimeConstrained` and `MemoryConstrained` did not stop several of these computations.

## Readiness against the survey's ranking

| # | Problem class | Verified in QF | What blocks a study | Readiness |
|---|---|---|---|---|
| 1 | Driven-dissipative Kerr: bistability, Liouvillian gap | Lindblad evolution converges in the cutoff; exact steady state from the Liouvillian null space, which also exposes multiple steady states; 28-digit evolution | A solve at N = 50 takes 3.6 s, 17 to 140 times the published numbers (see Dynamics) | **Ready** to N of about 60 |
| 2 | Transmon beyond two levels: leakage, DRAG | 4-level transmon, Gaussian pi pulse; the DRAG coefficient that minimizes leakage is found numerically in 1.2 s | | **Ready** |
| 3 | Small codes beyond Pauli noise | Density matrices, channels and qudits exist | The default circuit method takes 0.57 s at 12 qubits, 0.76 s at 16 and 1.5 s at 20 with TensorNetworks 1.1.0; no trajectory engine found | Repetition codes only |
| 4 | Cavity-QED spectra | Damped Jaynes-Cummings with the jump operator on the cavity alone matches the exact one-excitation solution | No correlation or spectrum routine | Dynamics ready; spectra are a gap |
| 5 | Exponentially small splittings and gaps | 28 correct digits in time evolution | Not yet run on the double-well benchmark | Promising, verify next |
| 6 | Rydberg arrays | Not tested | Same ODE and circuit performance limits as 1 and 3 | Untested |
| 7 | Quench dynamics, light cones, Trotter error | Full-chain light-cone contraction exact to 16 qubits | No automatic light-cone pruning: 24 qubits did not finish | Small chains only |
| 8 | Lossy phase estimation, Fock-space photonics | NOON quantum Fisher information matches n^2 eta^n to 6 digits up to n = 10 in 0.4 s; exact cat-state Wigner values | No bosonic loss channel | **Ready** |
| 9 | Trajectories, continuous measurement | Not found | No trajectory solver | Gap |
| 10 | Lindblad vs Bloch-Redfield vs HEOM | Lindblad and Kossakowski forms | No Bloch-Redfield or HEOM | Gap |
| 11 | Characterization and tomography | 9 Pauli settings x 1000 shots reconstructed to fidelity 0.9998 in 0.6 s; Bayesian estimator present | A fidelity call failed at 10,000 shots | **Ready** |
| 12 | Molecular dissociation | Not tested | | Validation only |

Variational studies are slow, consistent with the survey's advice to concede them: one energy of a 6-qubit, 12-parameter ansatz costs 0.35 s, and `FindMinimum` took 263 s.

## What to build first

1. **Transmon leakage and DRAG (class 2).** Ready now, no maintained incumbent since Qiskit Dynamics was archived, and the answer is decided numerically: the optimal DRAG coefficient is 0.825, between the two textbook values 0.5 and 1. Extend with T1 and T2 during the pulse, which is cheap at 4 levels, and an idle-gap sequence run with explicit step control.
2. **Lossy phase estimation (class 8).** The match leg already reproduces n^2 eta^n. The extend leg is to optimize the input state under loss and compare it with the coherent-plus-squeezed bound sqrt((1 - eta)/(eta N)) of Demkowicz-Dobrzanski et al., Fig. 12. The optimal state is found numerically, not derived.
3. **Exponentially small gaps (class 5).** Run QuantumToolbox.jl's double-well benchmark. If QF's arbitrary precision holds there, the study is a match plus a verification win, because that incumbent documents its own failure at N = 150.
4. **Tomography with credible regions (class 11).** How many shots does it take to certify entanglement at 95% posterior probability? QF's Metropolis-Hastings estimator answers that directly.

Kerr (class 1) is the survey's top problem and now runs at N up to about 60 in seconds per solve; cavity-QED spectra (class 4) still need a correlation routine.

## Defects and gaps, by priority

### P0: block a ranked study

Resolved on 5 October 2026; what each now does:

1. **Lindblad performance.** Above `$QuantumEvolveSparseThreshold` = 1,024 unknowns a numeric `QuantumEvolve` keeps its operators sparse, splitting a time-dependent one into numeric arrays with time-dependent coefficients; below it NDSolve integrates the dense equations faster. The Kerr oscillator end to end takes 0.42 s at N = 30, 0.88 s at 35, 1.41 s at 40, 3.6 s at 50 and 7.8 s at 60, agreeing with the matrix exponential of the Liouvillian to 7e-7 at N = 50. The incumbents' published N = 50 times are 0.026 to 0.21 s, on different parameters (QuantumToolbox.jl paper, arXiv:2504.21440). What remains is NDSolve's per-evaluation overhead: at N = 50 an explicit Runge-Kutta solve takes 1,171 steps and 11,707 right-hand-side evaluations at about 180 µs each, against 47 µs for the sparse matrix-vector product alone.
2. **Jump operators on a subsystem** act as the identity on the Hamiltonian's other qudits. Damped Jaynes-Cummings with `AnnihilationOperator[12, {2}]` as the jump operator matches the one-excitation solution e^(-kappa t/2) (cos w t + kappa/(4 w) sin w t)^2, w = sqrt(g^2 - kappa^2/16), to 1e-6, and writing each jump of a two-mode loss problem on its own mode gives the same density matrix as writing it on the full space.
3. **A matrix on an order** spans as many qudits as the order names: `QuantumOperator[IdentityMatrix[9], {2}]` is one 9-dimensional qudit, a dimension splits evenly when it can and takes its smallest divisors first when it cannot (12 on {1, 2} gives 2 x 6), and a matrix too small for its order is broadcast over it, as a named operator is. The two-mode NOON loss problem at Fock cutoff 9 that segfaulted runs in 0.08 s, with the |8,0><0,8| coherence matching e^(-n t)/2 to 5e-8.

### P1: degrade a study or its credibility

4. **Short pulses**, the most reported failure across QuTiP users (issue #2814), resolved on 5 October 2026. A numeric `QuantumEvolve` reads landmarks off the symbolic Hamiltonian and restarts its integration at each one: every discontinuity, and two and four widths either side of the centre of every Gaussian or exponential profile, which covers square, Gaussian, sech and tanh pulses. A pi pulse at t = 50 in [0, 100], with no options set:

   | Width | Square `Piecewise` | Square `UnitStep` | Gaussian | sech |
   |---|---|---|---|---|
   | 1.0 | captured | captured | captured (skipped before) | captured |
   | 0.1 | captured | captured | captured (skipped before) | captured |
   | 0.01 | captured (skipped before) | captured | captured (skipped before) | captured |
   | 0.001 | captured | captured | captured | captured |

   Each solve takes about 0.02 s; a train of 50 Gaussian pi/2 pulses of width 0.04 lands on the exact population to 2e-10 in 0.07 s. A smooth drive has no landmarks and runs as before. A pulse given as data (an `InterpolatingFunction`) is not read, so `MaxStepSize` remains the control there.
5. **Dense circuits at 20 qubits**, resolved on 8 October 2026. A 4-step Trotter quench applied with QF's default tensor-network method ran past 4 GB at 20 qubits, because TensorNetworks 1.0.11 contracts the whole network as one product. TensorNetworks 1.1.0 contracts along a greedy path, and QF puts diagonal gates on their wires' shared indices, so the default method now takes 0.57 s, 0.76 s and 1.54 s at 12, 16 and 20 qubits, against 1.8 s, 2.9 s and past 4 GB before (see [Circuits before and after](#circuits-before-and-after-the-work-of-6-to-8-october-2026)). The `"Schrodinger"` method takes 12.4 s and 73.3 s at 12 and 16 qubits and took 1,656 s at 20. Results agree throughout.
6. **`"MergeInterpolatingFunctions" -> False`** now starts the evolution at the start of the time range, resolved on 5 October 2026: a Rabi flop over [-1, 1] gives cos^2 2 = 0.1732.
7. **Circuit parameters do not bind by `ReplaceAll`.** `ansatz /. Thread[pars -> values]` silently leaves the circuit symbolic, and every later step computes symbolic expressions. Parameters must be declared (`"Parameters" -> pars`) and bound by calling the circuit, `ansatz @@ values`.
8. **Named operators with structure go through general algorithms.** QF stores every operator as a `SparseArray`, so a named operator whose structure fixes an answer still gets it from a general algorithm. At their default arguments, 23 of the 108 named operators are diagonal (Z, the phase gates, S, T, CZ, CPHASE, `"Z"[d]`, `"Diagonal"`), 18 are permutations (X, `"X"[d]`, SUM, SWAP, CNOT, Toffoli, Fredkin), 6 are permutations with phases (Y, CY), and the Fourier operator is the Fourier matrix. Of the 43 named circuits, PhaseOracle is diagonal, BooleanOracle and SimonOracle are permutations, and the Fourier circuit is the Fourier matrix, as is the Fourier basis. The exponential of a diagonal operator, and any function of one, already reads the diagonal (`191f1ac3`). These do not, measured on 6 October 2026 on main at `f46ec628` (probe `structured.wls`):

   | Named operator | Operation | QF | From the structure |
   |---|---|---|---|
   | `"Diagonal"[h]`, 10 qubits | `"Eigenvalues"` | 0.32 s, printing `Eigensystem::arh` | the diagonal, 0.1 ms |
   | `"Z"[128]` | `"Eigenvalues"` | 2.4 s, printing `Eigensystem::arhm` | the diagonal |
   | `"X"[64]` | `"Eigenvalues"` | 2.0 s | the 64th roots of unity, from its one cycle |
   | `"X"[64] @ "Z"[64]` | `"Eigenvalues"` | 1.6 s | the 64th roots of the product of its phases, from its one cycle |
   | Fourier circuit as an operator, 4 qubits | `"Eigenvalues"` | 8.6 s | F^4 = 1 |
   | `QuantumBasis["Fourier"[16]]` | change of basis | did not finish in 30 s | F^-1 = F^dagger, 0.5 ms |
   | `"Fourier"[14]` circuit | applied to a state | 4.1 s | `Fourier` of the amplitudes, 0.5 ms, agreeing to 1.5e-15 |
   | `"PhaseOracle"`, random function of 10 variables | build, then apply | 3.0 s + 11.6 s, 226 gates | the sign vector of its truth table, 0.5 ms, agreeing exactly |

   The change of basis is resolved on 6 October 2026 (`2ca84551`): an exactly unitary basis matrix inverts by its conjugate transpose, so a state goes into the Fourier basis in 0.01 s at d = 8, 0.02 s at 16, 0.12 s at 32 and 0.98 s at 64, where d = 16 to 64 did not finish before (the probe's second run). Machine, symbolic and non-unitary basis matrices keep the general inverse.

   Measuring a diagonal observable is slow too, 18.3 s for 64 outcomes on 6 qubits, but its eigensystem takes 0.01 s of that. The time goes into applying a measurement with many outcomes, which a computational-basis measurement does in 0.27 s, so the structure alone does not remove it. The task briefs in [Structured Arrays/tasks](../Structured%20Arrays/tasks/README.md) cover the diagonal and permutation spectra (briefs 2 and 4), the Fourier basis (brief 1) and the QFT on a state (brief 5); the permutations with phases and the phase oracle are not briefed.

### P2: friction or unsupported claims

9. No named bosonic loss channel; loss has to be written as Lindblad.
10. Light-cone locality is not automatic. The showcase draft's claim that the same contraction gives the same number "for a site deep inside a chain of fifty thousand" holds only because the draft builds the 9-qubit cone by hand.
11. `WhenEvent` works through `"AdditionalEquations"`, but only by referring to the undocumented internal state symbol `\[FormalS]`.
12. The QuEST backend is unavailable without a runtime download from qtechtheory.org.
13. `QF-Stabilizer-vs-Packages.md`, the source cited for the showcase's Stim comparison, is not in the repository.
14. Stabilizer per-gate cost is flat only up to about 3,000 qubits (see Circuits); the half-chain entropy dominates beyond that.
15. An evolved state read before its time is bound (`psi["StateVector"]` rather than `psi[t]`) failed with `Interpolation::inddp` when NDSolve repeats a grid point, which it does at a discontinuity of the Hamiltonian, because the Wolfram/Arrays lazy container rebuilt its interpolations over the whole grid. Resolved in Arrays 1.4.1, which QF now requires: the read matches `psi[t]` on both sides of a pulse edge.

## Measurements

### Dynamics

| Probe | Result |
|---|---|
| Transmon pi pulse, 4 levels, alpha/2pi = -200 MHz, Gaussian, sigma = T/4, no DRAG | leakage 0.475 at T = 4 ns, 0.0188 at 8 ns, 5.8e-4 at 16 ns; 0.07 to 0.3 s per solve |
| DRAG scan at T = 8 ns, lambda = 0 to 1.5 in steps of 0.25 | 0.0188, 0.00891, 0.00304, 6.98e-4, 1.15e-3, 3.56e-3, 7.07e-3 |
| `FindMinimum` over lambda | lambda* = 0.8254, leakage 5.78e-4, 1.2 s |
| Rabi flop with `WorkingPrecision -> 40` | 28 correct digits of cos^2 1; 7 at default settings |
| `WhenEvent` stopping when P0 first drops below 1/2 | integration ends at 0.7853981668; pi/4 = 0.7853981634 |
| Steady state of the Kerr oscillator at N = 20 from the Liouvillian null space | 0.84 s, null space of dimension 1 |
| Kerr oscillator, `QuantumEvolve` end to end, t in [0, 10] | 0.06 s at N = 10, 0.10 s at 20, 0.42 s at 30, 0.88 s at 35, 1.41 s at 40, 3.6 s at 50, 7.8 s at 60; mean photon number converged to 2.0732 by N = 20 |
| Damped Jaynes-Cummings, qubit x cavity of 12, g = 1, kappa = 0.05, jump on the cavity alone | population of the excited state with an empty cavity 0.2960, 0.1557, 0.3981 at t = 1, 2, 4, each within 1e-6 of the exact one-excitation solution |

### Bosonic

| Probe | Result |
|---|---|
| Mach-Zehnder at Fock cutoff 30 by direct composition | 1.8 s |
| Even cat, alpha = 2, Fock 40: W(0) | 0.6366198 = 2/pi, in 0.007 s; odd cat gives -2/pi |
| 21 x 21 Wigner grid at Fock 40 | 3.35 s; minimum -0.392 from the interference fringes |
| Lossy NOON at eta = 0.9, quantum Fisher information vs n^2 eta^n | n = 2: 3.2400002 vs 3.24; 4: 10.497601 vs 10.4976; 6: 19.131879 vs 19.131876; 8: 27.549908 vs 27.549901; 10: 34.867857 vs 34.867844; each under 0.5 s |

### Circuits

| Probe | Result |
|---|---|
| Stabilizer, 20,000 random H, CNOT and S gates | `ApplyCircuit` 1.32 s at 1,000 qubits (66 µs per gate), 0.85 s at 3,000 (42 µs), 5.19 s at 10,000 (259 µs) |
| Stabilizer half-chain entropy after those gates | 0.79 s, 4.2 s and 45.1 s |
| Light-cone magnetization, full chain, depth 4, greedy path | 2.15 s at 9 qubits, 3.32 s at 16, both -0.1711558 as with the hand-built cone; 24 qubits did not finish in 10 min at 6 to 10 GB |
| Dense `"Schrodinger"` engine, same circuit | 13.9 s at 12 qubits, 72.9 s at 16, 1,656 s at 20 |
| Default tensor-network method, same circuit | 0.57 s at 12 qubits, 0.76 s at 16 and 1.54 s at 20 with QF `abe5a3d7` and TensorNetworks 1.1.0; 1.8 s, 2.9 s and past 4 GB with QF `0bd75989` and TensorNetworks 1.0.11 (see the next section) |

### Circuits before and after the work of 6 to 8 October 2026

Re-measured on 9 October 2026 with probe `before-after.wls`, run by `before-after.sh`; the raw results are in `before-after.tsv` and the table comes from `before-after-table.py`. Each cell is the median of three runs, each in a fresh kernel that times only the workload, with the three environments alternating. The machine was an Apple M3 Max with 48 GB, under a load average of 4.2 to 11.7 (median 7.0) from other kernels. One repeat ran during a load spike, which the medians discard; 39 of the 47 timed cells have all three runs within 15% of their median.

- **Before**: QF `0bd75989`, WolframResearch/main before the work, with the TensorNetworks 1.0.11 and Arrays 1.4.1 it required.
- **QF before, new TN and Arrays**: the same QF with TensorNetworks 1.1.0 and Arrays 1.4.2, which separates what those two releases contribute.
- **After**: QF `abe5a3d7` with TensorNetworks 1.1.0 and Arrays 1.4.2.

| Workload | Before | QF before, new TN and Arrays | After | Before / after | Same result |
|---|---|---|---|---|---|
| QFT, 14 qubits: build (s) | 1.77 | 1.77 | 0.332 | 5.3x |  |
| Phase oracle, 10 variables: build (s) | 3.31 | 3.35 | 0.809 | 4.1x |  |
| Controlled phase gate of the QFT: build (ms per gate) | 16.1 | 16.9 | 2.37 | 6.8x |  |
| QFT, 14 qubits, applied to the zero state with qc[] (s, first / second) | 2.3 / 0.373 | 2.2 / 0.275 | 0.662 / 0.278 | 3.5x / 1.3x | yes |
| QFT, 14 qubits, applied to a random state (s, first / second) | 2.46 / 0.563 | 2.31 / 0.249 | 0.548 / 0.198 | 4.5x / 2.8x | yes |
| Phase oracle, 10 variables, applied to the zero state (s, first / second) | 12.9 / 5.01 | 8.78 / 0.945 | 2.01 / 0.956 | 6.4x / 5.2x | yes |
| Phase oracle, 10 variables, applied to the uniform state (s, first / second) | 12.8 / 4.99 | 8.73 / 0.943 | 1.64 / 0.569 | 7.8x / 8.8x | yes |
| Trotter quench, 4 steps, 12 qubits, default method (s) | 1.8 | 1.64 | 0.573 | 3.1x | yes |
| Trotter quench, 4 steps, 16 qubits, default method (s) | 2.9 | 1.95 | 0.763 | 3.8x | yes |
| Trotter quench, 4 steps, 20 qubits, default method (s) | out of 4 GB | 5.56 | 1.54 | before failed | yes |
| Change of basis into QuantumBasis["Fourier"[16]] (s) | over 60 s | over 60 s | 0.028 | before failed | yes |

The builds gain from QF alone: the middle column matches the first. QF builds gates in their final form (`af52593b`), and states and operators read their basis's data and shape off the basis (`22c88aaf`, `531cb2e6`, `659478eb`).

The applies gain from both:
- The TensorNetworks and Arrays releases bring the phase oracle's second applies from 5.0 s to 0.94 s, and make the 20-qubit quench fit in 4 GB.
- QF's apply-path changes (`f11c72f9`, `3900fe9e`, `bfeec465`) bring the first applies down: 2.2 s to 0.66 s for the QFT, 8.8 s to 2.0 s for the oracle.
- Diagonal gates on their wires' shared indices (`2d73fb17`, `cf951c97`) bring the quench down at every size (1.64 s to 0.57 s at 12 qubits, 5.56 s to 1.54 s at 20) and the oracle on the uniform state from 0.94 s to 0.57 s.

The change of basis is `c9100bb4`: an exactly unitary basis matrix inverts by its conjugate transpose. Every environment gives the same state: the amplitudes agree to 1e-10 and <Z> to 1e-12.

### Measurement, estimation, optimization

| Probe | Result |
|---|---|
| 9 two-qubit Pauli settings x 1,000 shots, then `QuantumStateEstimate` | 0.51 s to simulate, 0.11 s to estimate; maximum-likelihood infidelity 2.1e-4 (QF's `"Fidelity"` distance is 1 - F); shot sampling is multinomial, so its cost does not grow with the shot count |
| 6-site transverse-field Ising VQE, 12-parameter hardware-efficient ansatz | 0.35 s per energy (2.45 s on the `"Schrodinger"` engine); `FindMinimum` 263 s to -7.2342 against the exact -7.2962, a 0.85% gap set by the shallow ansatz |
| Present but not exercised | Bayesian estimate (Metropolis-Hastings with a Bures prior), quantum natural gradient, stochastic parameter-shift gradients, `QuantumAdiabaticEvolve`, QSVT, `QuantumLinearSolve`; `FubiniStudyMetricTensor` returned an unexpected shape and is unverified |
