# QF numeric capability inventory

What Wolfram QuantumFramework can do numerically today, measured rather than claimed, and what stands in the way of each problem class ranked in the companion survey, [Quantum software numerics priorities](reports/Quantum%20software%20numerics%20priorities.md).

**Environment.** QF 2.1.1 working tree, Wolfram Language 15.0, 4 October 2026 (the P0 defects re-measured after their fixes on 5 October), 48 GB macOS machine. Swap was 97% full and other kernels were running during part of the session, so absolute timings carry roughly ±2× noise; ratios measured inside one kernel are reliable. Every number below comes from a script in [numeric-inventory-probes/](numeric-inventory-probes/), runnable from that folder against the repository's paclet. The per-process probes' result lines are in `results.log` (kernel message dumps removed); `batchA`, `batchB2`, `batchB3`, `batchD` and `t0` print to the terminal. Run probes one per kernel with `./runprobes.sh probe.wls`, which kills a kernel that outlives its time limit: `TimeConstrained` and `MemoryConstrained` did not stop several of these computations.

## Readiness against the survey's ranking

| # | Problem class | Verified in QF | What blocks a study | Readiness |
|---|---|---|---|---|
| 1 | Driven-dissipative Kerr: bistability, Liouvillian gap | Lindblad evolution converges in the cutoff; exact steady state from the Liouvillian null space, which also exposes multiple steady states; 28-digit evolution | A solve at N = 50 takes 3.6 s, 17 to 140 times the published numbers (see Dynamics) | **Ready** to N of about 60 |
| 2 | Transmon beyond two levels: leakage, DRAG | 4-level transmon, Gaussian pi pulse; the DRAG coefficient that minimizes leakage is found numerically in 1.2 s | | **Ready** |
| 3 | Small codes beyond Pauli noise | Density matrices, channels and qudits exist | Dense circuit engine takes 14 s at 12 qubits and 28 min at 20; no trajectory engine found | Repetition codes only |
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
5. **Dense circuit engine.** A 4-step Trotter quench takes 13.9 s at 12 qubits, 72.9 s at 16 and 1,656 s at 20, about 10.6 s per gate at 2^20 amplitudes. Results are correct throughout.
6. **`"MergeInterpolatingFunctions" -> False`** now starts the evolution at the start of the time range, resolved on 5 October 2026: a Rabi flop over [-1, 1] gives cos^2 2 = 0.1732.
7. **Circuit parameters do not bind by `ReplaceAll`.** `ansatz /. Thread[pars -> values]` silently leaves the circuit symbolic, and every later step computes symbolic expressions. Parameters must be declared (`"Parameters" -> pars`) and bound by calling the circuit, `ansatz @@ values`.

### P2: friction or unsupported claims

8. No named bosonic loss channel; loss has to be written as Lindblad.
9. Light-cone locality is not automatic. The showcase draft's claim that the same contraction gives the same number "for a site deep inside a chain of fifty thousand" holds only because the draft builds the 9-qubit cone by hand.
10. `WhenEvent` works through `"AdditionalEquations"`, but only by referring to the undocumented internal state symbol `\[FormalS]`.
11. The QuEST backend is unavailable without a runtime download from qtechtheory.org.
12. `QF-Stabilizer-vs-Packages.md`, the source cited for the showcase's Stim comparison, is not in the repository.
13. Stabilizer per-gate cost is flat only up to about 3,000 qubits (see Circuits); the half-chain entropy dominates beyond that.
14. An evolved state read before its time is bound (`psi["StateVector"]` rather than `psi[t]`) failed with `Interpolation::inddp` when NDSolve repeats a grid point, which it does at a discontinuity of the Hamiltonian, because the Wolfram/Arrays lazy container rebuilt its interpolations over the whole grid. Resolved in Arrays 1.4.1, which QF now requires: the read matches `psi[t]` on both sides of a pulse edge.

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

### Measurement, estimation, optimization

| Probe | Result |
|---|---|
| 9 two-qubit Pauli settings x 1,000 shots, then `QuantumStateEstimate` | 0.51 s to simulate, 0.11 s to estimate; maximum-likelihood infidelity 2.1e-4 (QF's `"Fidelity"` distance is 1 - F); shot sampling is multinomial, so its cost does not grow with the shot count |
| 6-site transverse-field Ising VQE, 12-parameter hardware-efficient ansatz | 0.35 s per energy (2.45 s on the `"Schrodinger"` engine); `FindMinimum` 263 s to -7.2342 against the exact -7.2962, a 0.85% gap set by the shallow ansatz |
| Present but not exercised | Bayesian estimate (Metropolis-Hastings with a Bures prior), quantum natural gradient, stochastic parameter-shift gradients, `QuantumAdiabaticEvolve`, QSVT, `QuantumLinearSolve`; `FubiniStudyMetricTensor` returned an unexpected shape and is unverified |
