# Open-quantum-system and quantum-dynamics libraries: numerical capabilities, published benchmarks and documented limits (QuTiP 5 family, QuantumOptics.jl, QuantumToolbox.jl, Dynamiqs)

Scope note: this file describes what the other libraries do and how well. It does not evaluate QF. Sources were read on 2026-10-04: the three library papers as full-text PDFs, the repositories and docs through the GitHub API, and the benchmark scripts and result JSONs. "Verified" means read in the primary source. Items under **Inferences** are my own reasoning. Library versions on 2026-10-04 come from each project's GitHub "latest release": QuTiP **v5.3.1** (2026-08-04), qutip-qtrl **v0.2.0** (2026-06-29), qutip-qoc **v0.2.0** (2026-03-24), qutip-qip **v0.4.2** (2026-06-23), qutip-jax **v0.1.1** (2025-05-29), QuantumToolbox.jl **v0.50.0** (2026-10-02), HierarchicalEOM.jl **v2.16.2** (2026-09-24), QuantumOptics.jl **v1.2.10** (2026-09-11), Dynamiqs **v0.3.6** (2026-06-08). Sources: [qutip releases](https://github.com/qutip/qutip/releases), [qutip-qtrl](https://github.com/qutip/qutip-qtrl/releases), [qutip-qoc](https://github.com/qutip/qutip-qoc/releases), [qutip-qip](https://github.com/qutip/qutip-qip/releases), [qutip-jax](https://github.com/qutip/qutip-jax/releases), [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl/releases), [HierarchicalEOM.jl](https://github.com/qutip/HierarchicalEOM.jl/releases), [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl/releases), [dynamiqs](https://github.com/dynamiqs/dynamiqs/releases). GitHub stars on 2026-10-04, a rough measure of adoption: qutip 2081, QuantumOptics.jl 624, dynamiqs 325, QuantumToolbox.jl 172, HierarchicalEOM.jl 49 ([qutip](https://github.com/qutip/qutip), [QuantumOptics.jl](https://github.com/qojulia/QuantumOptics.jl), [dynamiqs](https://github.com/dynamiqs/dynamiqs), [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl), [HierarchicalEOM.jl](https://github.com/qutip/HierarchicalEOM.jl)).

Primary papers:
- **QuTiP 5**: Lambert et al., "QuTiP 5: The Quantum Toolbox in Python", arXiv:2412.04705 (v1 2024-12-06, v2 2025-10-01), *Physics Reports* 1153 (18 Jan 2026) 1–62, DOI 10.1016/j.physrep.2025.10.001 ([arXiv](https://arxiv.org/abs/2412.04705), [DOI](https://doi.org/10.1016/j.physrep.2025.10.001)).
- **QuantumToolbox.jl**: Mercurio, Huang, Cai, Chen, Savona, Nori, *Quantum* 9, 1866 (2025), arXiv:2504.21440 ([arXiv](https://arxiv.org/abs/2504.21440), [DOI](https://doi.org/10.22331/q-2025-09-29-1866)).
- **QuantumOptics.jl**: Krämer, Plankensteiner, Ostermann, Ritsch, *Comput. Phys. Commun.* 227, 109–116 (2018), arXiv:1707.01060. This paper is old: it describes v0.4.1 ([arXiv](https://arxiv.org/abs/1707.01060), [DOI](https://doi.org/10.1016/j.cpc.2018.02.004)).
- **Dynamiqs**: there is no published paper. The README cites Guilmin, Bocquet, Genois, Weiss, Gautier (2025) as "in preparation" ([README](https://github.com/dynamiqs/dynamiqs)).

---

## Q1. Which numerical problem classes do the official galleries, documentation and papers feature, and how well does each library serve them?

### Takeaway
QuTiP has by far the widest set of problem classes:
- Lindblad, Monte Carlo trajectories, Bloch–Redfield, Floquet, HEOM, non-Markovian Monte Carlo and stochastic master equations
- PIQS for permutation-invariant spin ensembles
- optimal control and pulse-level circuits

Its gallery is built on the cavity-QED canon: damped Jaynes–Cummings and vacuum Rabi, Kerr, cat states, resonance fluorescence, correlation functions, optomechanics and superradiance.

QuantumToolbox.jl mirrors QuTiP's scope in Julia, gets HEOM through HierarchicalEOM.jl, and adds:
- adaptive Fock-truncation solvers (DFD/DSF)
- low-rank evolution
- parameter-sweep maps
- arbitrary precision

QuantumOptics.jl is broad and mature (semiclassical, particle/FFT and many-body models, stochastic, Bloch–Redfield, correlations) but has no HEOM.

Dynamiqs is deliberately narrow. It covers the Schrödinger and Lindblad equations, Floquet, propagators and jump/diffusive stochastic unravellings, and is optimised for GPUs, batching and gradients. Its API reference lists no steady-state, correlation/spectrum, Bloch–Redfield or HEOM solver.

### Cited Findings

**QuTiP 5 (core) solver set and showcased problems**
- The QuTiP 5 paper's table of contents lists these solver sections:
  - `sesolve`/`mesolve`, including time-dependent systems, and JAX/GPU via Diffrax
  - `steadystate`
  - `mcsolve` (quantum trajectories)
  - `nm_mcsolve` (Monte Carlo for non-Markovian baths)
  - `brmesolve` (Bloch–Redfield)
  - Floquet methods
  - `smesolve` (stochastic master equation)
  - `HEOMSolver`
  - excitation-number-restricted (ENR) states
  - JAX automatic differentiation
  - MPI
  - QuTiP-QOC (GRAPE, CRAB, GOAT, JOPT)
  - QuTiP-QIP

  — [QuTiP 5 paper](https://arxiv.org/abs/2412.04705)
- QuTiP's original examples "focused on traditional models from cavity quantum electrodynamics, like Lindblad master equation simulations of the open Jaynes-Cummings, Rabi and Dicke models", with tools for the Wigner function and photonic g⁽²⁾(t). — [QuTiP 5 paper §1](https://arxiv.org/abs/2412.04705)
- The official v5 gallery has these sections: "Quantum Mechanics Lectures with QuTiP", "Time Evolution", "Optimal Control", "Permutational Invariant Lindblad Dynamics", "Hierarchical Equations of Motion", "Pulse-level Circuit Simulation", "Quantum Circuits and Algorithms", "Visualization" and "Miscellaneous". All notebooks target QuTiP 5. — [qutip-tutorials gallery](https://qutip.org/qutip-tutorials/)
  - **Lectures:**
    - L1 vacuum Rabi oscillations in the Jaynes–Cummings model
    - L2A two-qubit gate via a resonator coupler
    - L2B single-atom lasing
    - L3A Dicke model
    - L3B JC model in the ultrastrong-coupling regime
    - L4 correlation functions
    - L5 parametric amplifier
    - L6 quantum Monte Carlo trajectories
    - L7 iSWAP gate and process tomography
    - L8 adiabatic sweep
    - L9 squeezed states
    - L10 cavity QED in the dispersive regime
    - L11 charge qubits
    - L12 decay into squeezed vacuum
    - L13 resonance fluorescence
    - L14 Kerr nonlinearities
    - L15 cascaded systems
    - L16 gallery of Wigner functions

    — [gallery](https://qutip.org/qutip-tutorials/)
  - **Time evolution:**
    - 004 vacuum Rabi (mesolve)
    - 005 spin chain
    - 006 Monte Carlo "Birth and Death of Photons in a Cavity"
    - 007–010 Bloch–Redfield: two-level system; time-dependent operators; dissipative atom–cavity; phonon-assisted initialization
    - 011–012 Floquet
    - 013 non-Markovian Monte Carlo
    - 015–017 stochastic: heterodyne; inefficient detection; JC photocurrent
    - 018 "Stochastic vs. Monte-Carlo Solver: Cat states become coherent"
    - 019 steady state of an optomechanical system in the single-photon strong-coupling regime
    - 020 homodyned Jaynes–Cummings emission
    - 021 quasi-steady state of a periodically driven system
    - 023 Dysolve propagators
    - 022, 024, 025, 026: the "QuTiPv5 Paper Example" notebooks (smesolve homodyne; sesolve/mesolve; Floquet speed test; nm_mcsolve)

    — [gallery](https://qutip.org/qutip-tutorials/)
  - **HEOM:**
    - 1a–1e spin-boson model: introduction; very strong coupling; underdamped; Ohmic fitting; pure dephasing
    - 2 Fenna–Matthews–Olson (FMO) complex
    - 3 quantum heat transport
    - 4 dynamical decoupling of a non-Markovian environment
    - 5a fermionic single-impurity model
    - 5b discrete boson coupled to an impurity and fermionic leads

    — [gallery](https://qutip.org/qutip-tutorials/)
  - **Permutational-invariant (PIQS):**
    - overview
    - superradiant light emission
    - steady-state superradiance
    - open Dicke model
    - spin squeezing with noise
    - boundary time crystals
    - multiple spin ensembles
    - entropy and purity

    These links point to the legacy `qutip-notebooks` repo, not to v5-ported notebooks. — [gallery](https://qutip.org/qutip-tutorials/)
- QuTiP 5 paper §3.2.4 shows the value of microscopically derived solvers. A qubit's energies are adiabatically switched, H = (Δ/2) sin(ω_d t) σ_z. The Bloch–Redfield and HEOM solvers "capture this switching effect automatically", while "the naive approach, however, using a single constant collapse operator σ− fails and is insensitive to the drive" (Fig. 3b). — [QuTiP 5 paper](https://arxiv.org/abs/2412.04705)
- QuTiP 5 paper on HEOM: it is "a numerically exact method … under a minimal set of assumptions": a Gaussian bath, initially in a thermal equilibrium state, coupled linearly. It requires a multi-exponential decomposition of the bath correlation function (Matsubara or Padé with Nk terms), and the hierarchy is truncated at `max_depth`. The paper shows Ohmic spectral-density fitting (Fig. 12). It demonstrates the zero-temperature localization–delocalization transition of the spin-boson model, which "from a completely delocalized state" goes to "what appears to be a localized state" as coupling increases. — [QuTiP 5 paper §3.2.12](https://arxiv.org/abs/2412.04705)
- QuTiP-BoFiN paper (the basis of the HEOM solver; Phys. Rev. Research 5, 013181 (2023)) shows three things. It fits arbitrary spectral densities. It runs FMO energy transfer, "showing how a suitable non-Markovian environment can protect against pure dephasing". It benchmarks dynamical-decoupling strategies, finding "the Uhrig pulse-spacing scheme is less optimal than equally spaced pulses when the environment's spectral density is very broad". It also includes an integrable fermionic single-impurity benchmark. — [QuTiP-BoFiN, arXiv:2010.10806](https://arxiv.org/abs/2010.10806)
- **nm_mcsolve** (new in v5) handles time-local master equations whose rates γ_n(t) can go negative. It uses an "influence martingale" trajectory weighting and requires jump operators satisfying Σ_n A_n†A_n = α·1. The solver adds a zero-rate jump operator automatically if needed. The worked example is the damped Jaynes–Cummings model with a Lorentzian environment. — [QuTiP 5 paper §3.2.8](https://arxiv.org/abs/2412.04705)
- **smesolve** standard example: homodyne detection of a cavity, dρ = −i[H,ρ]dt + D[a]ρ dt + H[a]ρ dW. — [QuTiP 5 paper §3.2.11](https://arxiv.org/abs/2412.04705)
- **mcsolve** features:
  - 500 trajectories by default
  - `"map":"parallel"`, a `timeout`, and a `target_tol` that stops when the statistical error reaches a target
  - `improved_sampling`, which runs the no-jump trajectory only once and is useful when dissipation rates are very small
  - support for mixed initial states (new in v5)
  - planned work on waiting-time-distribution methods

  — [QuTiP 5 paper §3.2.7](https://arxiv.org/abs/2412.04705)
- **brmesolve** lets the secular approximation be "relaxed to a specified degree, leading to a so-called non-secular and non-Lindblad equation of motion". — [QuTiP 5 paper §3.2.9](https://arxiv.org/abs/2412.04705)
- **Correlation functions** (QuTiP 5 guide):
  - `correlation_2op_1t`, `correlation_2op_2t`, `correlation_3op_1t`, `correlation_3op_2t`, `correlation_3op` (new in 5.0, which also accepts e.g. `brmesolve`), and helpers `coherence_function_g1/g2`
  - `spectrum()`, which uses the exponential-series approach by default, and `spectrum_correlation_fft()`
  - worked examples: steady-state ⟨x(t)x(0)⟩ of a leaky cavity; the emission spectrum showing vacuum Rabi splitting in the JC model; g⁽¹⁾ of a coherent state decaying to thermal; g⁽²⁾ of coherent, thermal and Fock states

  — [QuTiP correlation guide](https://qutip.readthedocs.io/en/latest/guide/guide-correlation.html)
- **PIQS** (QuTiP's permutational-invariant solver): permutational invariance gives "an exponential reduction in the computational resources required to study the Lindblad dynamics of coupled spin-boson ensembles evolving under the effect of both local and collective noise". It is applied to spin squeezing, superradiance and quantum phase transitions under local dissipation. — [Shammah et al., PRA 98, 063815 (2018), arXiv:1805.05129](https://arxiv.org/abs/1805.05129)

**QuTiP family: optimal control and circuits**
- qutip-qtrl is the former `qutip.control` module, split out at QuTiP 5.0. It "offers support for both the CRAB and GRAPE methods". — [qutip-qtrl README](https://github.com/qutip/qutip-qtrl)
- qutip-qoc adds GOAT and JOPT, the latter using JAX autodiff through QuTiP 5's Diffrax support. It has a two-layer search: global with `dual_annealing`/`basinhopping` and local with any gradient-based `scipy.optimize.minimize`. It switches easily between GOAT, JOPT, GRAPE and CRAB, and its README mentions optional dependencies for an "RL (reinforcement learning) algorithm". — [qutip-qoc README](https://github.com/qutip/qutip-qoc)
- QuTiP 5 paper Table 10 gives a capability matrix:
  - GRAPE: no analytic controls; local search since v4; global since v5; no variable time; multi-objective since v5
  - CRAB: analytic and local since v4; global v5; no variable time
  - GOAT and JOPT: all features from v5, including optimizing the evolution time

  In the worked example, all four algorithms reach Hadamard-gate fidelity > 0.99 on a single qubit with optional dissipation (Fig. 18). — [QuTiP 5 paper §4.1](https://arxiv.org/abs/2412.04705)
- The gallery optimal-control notebooks are:
  - GRAPE (L-BFGS-B) for a Hadamard gate, a 2-qubit QFT, Lindbladian dynamics, symplectic dynamics and a CNOT
  - CRAB for 2-qubit state-to-state transfer and a QFT

  — [gallery](https://qutip.org/qutip-tutorials/)
- qutip-qip covers pulse-level noisy circuit simulation. Examples include Deutsch–Jozsa compiled on superconducting-qubit and spin-chain processors, cross-talk noise in an ion processor, and a Ramsey experiment with Lindblad dynamics. — [Li et al., Quantum 6, 630 (2022), arXiv:2105.09902](https://arxiv.org/abs/2105.09902)
- Gallery pulse-level notebooks include a 10-qubit QFT compiled and simulated, randomized benchmarking, and T-relaxation measurement with the idling gate. — [gallery](https://qutip.org/qutip-tutorials/)

**QuantumToolbox.jl**
- Solvers and features shown in the paper:
  - deterministic and stochastic dynamics (`sesolve`, `mesolve`, `mcsolve`, `ssesolve`, `smesolve`)
  - variational low-rank dynamics (based on Gravina & Savona, PRR 6, 023072 (2024))
  - HEOM through HierarchicalEOM.jl
  - `steadystate` ("direct, iterative, and eigenvalue-based methods")
  - `steadystate_fourier` for time-averaged steady states of periodically driven systems
  - correlation functions and spectra
  - GPU via CUDA arrays
  - distributed computing via Distributed.jl/SLURM
  - automatic differentiation

  — [QuantumToolbox.jl paper](https://arxiv.org/abs/2504.21440)
- The QuantumToolbox.jl paper presents two adaptive Fock-space algorithms. I found no equivalent in the other three libraries' API listings; that absence is my inference.
  - **Dynamical Fock dimension** (`dfd_mesolve`) "continuously monitors the population of the Fock states in time, increasing or decreasing the cutoff dimension of the Hilbert space accordingly. In this way, it is no longer required to check for the convergence of the results as a function of a fixed cutoff."
  - **Dynamical shifted Fock** (DSF) "efficiently simulates strongly driven systems by maintaining a low-dimensional Hilbert space" by monitoring the coherence α and applying displacement transformations at a threshold. It is demonstrated on a driven Jaynes–Cummings model with a Kerr nonlinearity and on coupled Kerr oscillators, compared against second-order cumulant expansion.

  — [QuantumToolbox.jl paper §3](https://arxiv.org/abs/2504.21440)
- The current API includes `mesolve_map`/`sesolve_map` for parameter sweeps; `spectrum` with `ExponentialSeries`, `PseudoInverse` and `Lanczos` solvers; `correlation_2op_1t/2op_2t/3op_1t/3op_2t`; `spectrum_correlation_fft`; `brmesolve` and `bloch_redfield_tensor`; `liouvillian_dressed_nonsecular`; `dfd_mesolve`, `dsf_mesolve` and `dsf_mcsolve`; `lr_mesolve`; and `steadystate_fourier`. — [QuantumToolbox.jl API docs source](https://github.com/qutip/QuantumToolbox.jl/blob/main/docs/src/resources/api.md)
- Steady-state methods: `SteadyStateDirectSolver`, `SteadyStateLinearSolver` (LinearSolve.jl, e.g. GMRES with ILU preconditioner, MKL Pardiso), `SteadyStateEigenSolver` and `SteadyStateODESolver`. — [QuantumToolbox.jl steady-state guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/steadystate)
- HierarchicalEOM.jl computes bosonic and fermionic spectra, stationary states and full auxiliary-density-operator (ADO) dynamics. It is exemplified on the single-impurity Anderson model and on an ultrastrongly coupled charge–cavity system with bosonic and fermionic reservoirs. — [Huang et al., Commun. Phys. 6, 313 (2023), arXiv:2306.07522](https://arxiv.org/abs/2306.07522)

**QuantumOptics.jl**
- Current API solver families:
  - `timeevolution`: `schroedinger`, `master`, `master_h`, `master_nh`, `master_dynamic`, `master_nh_dynamic`, `mcwf` (+`_h`, `_nh`, `_dynamic`), `master_bloch_redfield`, `bloch_redfield_tensor`
  - `stochastic`: `schroedinger`, `master`, `homodyne_carmichael`, plus `_semiclassical` and `_dynamic` variants
  - `semiclassical`: `schroedinger_dynamic`, `master_dynamic`, `mcwf_dynamic`
  - `steadystate`: `master`, `eigenvector`, `iterative`, `liouvillianspectrum`
  - `timecorrelations`: `correlation`, `spectrum`, `correlation2spectrum`, `correlation_dynamic`

  No HEOM appears in the API. — [QuantumOptics.jl API docs source](https://github.com/qojulia/QuantumOptics.jl/blob/master/docs/src/api.md)
- The 2018 paper's worked examples (v0.4.1) are a lossy Jaynes–Cummings model, a time-dependent JC model, the Gross–Pitaevskii / nonlinear Schrödinger equation (collision of two soliton-like wave packets), and a semiclassical model of cavity cooling. — [QuantumOptics.jl paper §4](https://arxiv.org/abs/1707.01060)
- Current example notebooks:
  - atom dephasing, cavity cooling, correlation spectrum, Doppler cooling, Jaynes–Cummings
  - lasing and cooling, many-body four-level system, N particles in a double well, optomechanical cooling
  - particle in a harmonic trap, particle into a barrier, pumped cavity, quantum kicked top, quantum Zeno effect
  - Raman, Ramsey, spin–orbit-coupled 1D BEC, superradiant laser, three-level maser
  - two-qubit entanglement, vortex, 2D wavepacket

  — [QuantumOptics.jl examples directory](https://github.com/qojulia/QuantumOptics.jl/tree/master/docs/examples/notebooks)
- The qojulia ecosystem handles large spin ensembles through QuantumCumulants.jl, "a Julia framework for generalized mean-field equations in open quantum systems". — [Plankensteiner, Hotter, Ritsch, Quantum 6, 617 (2022), arXiv:2105.01657](https://arxiv.org/abs/2105.01657)

**Dynamiqs**
- Solvers in the API reference:
  - general: `sesolve`, `mesolve`, `sepropagator`, `mepropagator`, `floquet`
  - stochastic: `jssesolve`, `dssesolve`, `jsmesolve`, `dsmesolve` (jump/diffusive SSE and SME)
  - methods: `Tsit5`, `Dopri5`, `Dopri8`, `Kvaerno3`, `Kvaerno5` (implicit), `Euler`, `EulerJump`, `EulerMaruyama`, `Rouchon1–3`, `Expm`, `Event`, `JumpMonteCarlo`, `DiffusiveMonteCarlo`, `LowRank`
  - gradient modes: `Direct`, `BackwardCheckpointed`, `Forward`, `HigherOrder`
  - optimal-control helpers: `snap_gate`, `cd_gate`

  The API index lists **no** steady-state, correlation-function/spectrum, Bloch–Redfield or HEOM routine. — [Dynamiqs Python API index source](https://github.com/dynamiqs/dynamiqs/blob/main/docs/python_api/index.md)
- A 2023 issue titled "Steadystate solvers" ("Useful but currently missing feature") was closed as completed on 2023-11-17, during the PyTorch era. The current JAX-era API index shows no steady-state function. — [dynamiqs #240](https://github.com/dynamiqs/dynamiqs/issues/240); [API index](https://github.com/dynamiqs/dynamiqs/blob/main/docs/python_api/index.md)
- Stated purpose: GPU-accelerated and differentiable simulation for "simulation of large quantum systems, gradient-based parameter estimation or quantum optimal control". It is sponsored by Alice & Bob (cat qubits) and was rewritten from PyTorch to JAX in early 2024. — [Dynamiqs "What is" page source](https://github.com/dynamiqs/dynamiqs/blob/main/docs/documentation/getting_started/whatis.md)
- Advanced examples: "Driven-dissipative Kerr oscillator" (constant-drive evolution and "Periodic revival of a coherent state"), "Continuous jump measurement" (`jssesolve`/`jsmesolve`, a qubit with σ± loss) and "Continuous diffusive measurement". Basics tutorials include batching, gradients, Floquet integration and time-dependent operators. — [Kerr example](https://www.dynamiqs.org/stable/documentation/advanced_examples/kerr-oscillator.html); [jump measurement](https://www.dynamiqs.org/stable/documentation/advanced_examples/continuous-jump-measurement.html); [diffusive measurement](https://www.dynamiqs.org/stable/documentation/advanced_examples/continuous-diffusive-measurement.html); [Floquet tutorial](https://www.dynamiqs.org/stable/documentation/basics/floquet-integration.html)
- Recent releases:
  - v0.3.1 (2025-02): diffusive SME `dsmesolve`
  - v0.3.2 (2025-03): `Event` jump-SSE method, "a major step towards a blazingly fast Monte Carlo solver"
  - v0.3.3 (2025-07): higher-order adaptive Rouchon methods; jump and diffusive Monte Carlo for Lindblad, "especially well on GPUs"
  - v0.3.4 (2025-11): Floquet tutorial; vectorized Lindblad
  - v0.3.6 (2026-06): low-rank `mesolve` and Runge–Kutta Kraus-map Rouchon schemes

  — [dynamiqs releases](https://github.com/dynamiqs/dynamiqs/releases)

### Inferences
- **Coverage ranking (inference from the API and docs evidence above):** QuTiP ≳ QuantumToolbox.jl > QuantumOptics.jl > Dynamiqs.
  - QuTiP is the only one with HEOM, nm_mcsolve, PIQS and the QOC/QIP family in one ecosystem.
  - QuantumToolbox.jl gets HEOM through HierarchicalEOM.jl and adds DFD/DSF, low-rank and parameter maps.
  - QuantumOptics.jl is the most physically diverse outside quantum optics: particles on grids/FFT, semiclassical hybrids, state-dependent dynamics.
  - Dynamiqs trades coverage for GPU speed and differentiability.
- The cavity-QED canon is the shared "lingua franca" of all four libraries: damped JC/vacuum Rabi, the driven-dissipative Kerr oscillator, cat-state decoherence, homodyne/photon-counting trajectories, and Ising chains for scaling. That makes it the natural ground for head-to-head comparison.
- For non-Markovian problems, QuTiP and QuantumToolbox.jl (via HierarchicalEOM.jl) are the only incumbents among the four with numerically exact machinery. QuantumOptics.jl and Dynamiqs would need Bloch–Redfield (QuantumOptics.jl only) or user-built pseudomode/chain models.

### Gaps
- I did not open each gallery notebook to extract its parameters. Only titles and links were verified.
- Some Dynamiqs tutorials may contain more sections (for example optimal control within the Kerr example). I only read the section headers.

---

## Q2. What scale and performance do they publish (Hilbert-space sizes, timings, GPU support, cross-library benchmarks)?

### Takeaway
The only recent benchmark that covers all four libraries is in the QuantumToolbox.jl paper (Quantum 2025), written by QuantumToolbox.jl's own authors.

On a driven-dissipative Kerr oscillator with Fock cutoff N = 50, QuantumToolbox.jl is fastest for mesolve, mcsolve and smesolve:
- **mesolve:** QuantumToolbox.jl ≈ 0.026 s, QuantumOptics.jl ≈ 0.045 s, Dynamiqs ≈ 0.046 s, QuTiP ≈ 0.21 s.

On a dissipative Ising chain:
- full master equations of 10 spins take 10³ s or more on one CPU in every library
- with a GPU (RTX 4090), 12 spins is the practical ceiling, again at ≈ 10³ s

QuTiP's own paper reports:
- a GPU crossover with up to 100× speed-ups
- a memory wall of **11 spins for mesolve and 22 for sesolve** on an 80 GB A100

QuantumOptics.jl's published benchmarks (2017–2020) are against QuTiP 4.2 and the Matlab QO Toolbox and are dated.

### Cited Findings

**QuantumToolbox.jl paper, §6 "Performance comparison with other packages" (Quantum 9, 1866, 2025)**
- Packages: QuTiP, QuantumOptics.jl, Dynamiqs and QuantumToolbox.jl. Solvers: `mesolve`, `mcsolve` and `smesolve`, plus AD of the master equation. Runs on CPU and GPU. Hardware: Intel i9-13900KF with 64 GB RAM; GPU an NVIDIA GeForce RTX 4090. "The GPU version of QuTiP was performed using its JAX backend, i.e., qutip-jax." — [QuantumToolbox.jl paper §6](https://arxiv.org/abs/2504.21440)
- Versions (paper Table 1, "latest version of each package as of Aug. 28, 2025"): Julia 1.11.6, QuantumToolbox.jl 0.34.1, QuantumOptics.jl 1.2.3; Python 3.13.5, QuTiP 5.2.1, QuTiP-JAX 0.1.1, dynamiqs 0.3.3. — [QuantumToolbox.jl paper](https://arxiv.org/abs/2504.21440)
- The paper's verdict: "QuantumToolbox.jl outperforms other packages in all scenarios. While the AD support in QuantumToolbox.jl should still be considered experimental, it is already competitive for backpropagation of the master equation." — [QuantumToolbox.jl paper Fig. 8 caption](https://arxiv.org/abs/2504.21440)
- **Benchmark problem 1, driven Kerr oscillator** (Fig. 8a–d): Ĥ = Δâ†â − Uâ†²â² + F(â + â†) with single-photon loss. — [QuantumToolbox.jl paper Eq. 31](https://arxiv.org/abs/2504.21440)
  - Script parameters: N = 50, Δ = 0.1, U = −0.05, F = 2, γ = 1, n_th = 0.2. Collapse operators are √(γ(1+n_th)) a and √(γ n_th) a†. t ∈ [0, 10] with 100 save points; `ntraj = 100` for mcsolve and smesolve. — [benchmark script qutip_benchmarks.py](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmarks.py)
  - The AD benchmark uses N = 100. — [benchmark script](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmarks.py)
- **Median timings** from the published result JSONs, converted from nanoseconds (my conversion; the raw values are BenchmarkTools.jl- and timeit-derived):

  | Solver (Kerr, N = 50) | QuantumToolbox.jl | QuantumOptics.jl | Dynamiqs | QuTiP |
  |---|---|---|---|---|
  | mesolve | 0.026 s | 0.045 s | 0.046 s | 0.206 s |
  | mcsolve, 100 trajectories | 0.018 s | 0.023 s | 6.8 s | 0.129 s |
  | smesolve, 100 trajectories | 1.0 s | 4.4 s | 11.4 s | 6.8 s |
  | reverse-mode gradient of mesolve (N = 100) | 12.8 s | n/a | ≈13.1–13.5 s | ≈16.7–17.3 s (qutip-jax) |

  Raw results: [QuTiP JSON](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmark_results.json), [dynamiqs JSON](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/dynamiqs_benchmark_results.json), [QuantumToolbox JSON](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/julia/quantumtoolbox_benchmark_results.json), [QuantumOptics JSON](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/julia/quantumoptics_benchmark_results.json), [autodiff JSONs](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/tree/main/src/benchmarks)
- Methodology caveats visible in the scripts:
  - The Dynamiqs script sets `set_device("cpu")` and **double** precision (`set_precision("double")`, "Set the same precision as the others") at module level, and uses `Tsit5(rtol=1e-6, atol=1e-8)` for mesolve. I did not check how the GPU runs override the device.
  - Its "mcsolve" entry is actually `jssesolve` with a **fixed-step `EulerJump(dt=1e-3)`**, not an adaptive Monte Carlo wave-function solver.
  - QuTiP's mcsolve used `"map":"parallel"` multiprocessing.

  — [dynamiqs_benchmarks.py](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/dynamiqs_benchmarks.py); [qutip_benchmarks.py](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmarks.py)
- **Benchmark problem 2, dissipative 1D Ising chain** (Fig. 8e–f): Ĥ = h_z Σσ_k^z + J_x Σσ_k^x σ_{k+1}^x with local σ_j⁻ dissipators. The script uses J_x = 25, h_z = 50, γ = 1. N runs over 2–10 spins on CPU and 2–12 on GPU; QuTiP-JAX on GPU only to N = 10. — [paper Eq. 32](https://arxiv.org/abs/2504.21440); [qutip_benchmarks.py](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmarks.py); [quantumtoolbox.jl script](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/julia/quantumtoolbox.jl)
- **Median mesolve times from the published N-scaling JSONs** (my conversion; full lists at [benchmarks dir](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/tree/main/src/benchmarks)):

  | Configuration | QuantumToolbox.jl | QuantumOptics.jl | Dynamiqs | QuTiP |
  |---|---|---|---|---|
  | CPU, N = 8 | 32.5 s | 53.6 s | 210 s | 70.6 s |
  | CPU, N = 10 | 1.13×10³ s | 2.2×10³ s | 2.09×10³ s | 1.9×10³ s |
  | GPU, N = 10 | 34.9 s | 110 s | 66.6 s | 144 s (qutip-jax) |
  | GPU, N = 12 | 1.16×10³ s | 2.15×10³ s | 2.14×10³ s | not run |

- **Distributed demo:** a 2D transverse-field Ising model on a 4×4 lattice (16 spins), solved with `mcsolve` with local decay. Trajectories are parallelized on a SLURM cluster with Distributed.jl/SlurmClusterManager.jl, using 10 nodes × 72 threads (720 threads) of Intel Xeon Platinum 8360Y. "Quantum trajectories are very useful here, as they avoid the need to compute the full Liouvillian of such a large system." — [QuantumToolbox.jl paper §4](https://arxiv.org/abs/2504.21440)
- Continuous benchmarking exists: the README links a "Benchmarks" CI page. — [QuantumToolbox.jl README](https://github.com/qutip/QuantumToolbox.jl); [benchmark page](https://qutip.org/QuantumToolbox.jl/benchmarks/)

**QuTiP 5 paper (own performance data; no cross-library comparison)**
- **GPU benchmark (Fig. 4):**
  - Problem: a 1D Ising spin chain. Hardware: NVIDIA A100 with 80 GB, against an AMD EPYC 7713 (64 cores).
  - "(a) shows a noiseless example for up to 22 spins (24 on CPU) using sesolve. Panel (b) shows the same problem in the presence of noise for up to 11 spins (12 on CPU) using mesolve."
  - "we observe a threshold in system size where the GPU outperforms the CPU calculation by up to two orders of magnitude"
  - "even with the jaxdia format and using a state-of-the-art graphics card with 80 gigabytes of RAM, the memory limit is reached already at 11 spins"

  — [QuTiP 5 paper §3.2.5](https://arxiv.org/abs/2412.04705)
- GPU rationale and caveat: "GPUs tend to shine when evaluating many small matrix-vector problems in parallel … However, when solving a single ODE of a large system … the potential advantage is less clear. The sequential nature of an ODE integrator makes it hard to parallelize." — [QuTiP 5 paper §3.2.5](https://arxiv.org/abs/2412.04705)
- **Floquet speed test:** `FloquetBasis`/`fsesolve` cost is "on average, independent of the time until which we want to evolve". A crossover time t_cross exists beyond which Floquet beats `sesolve`, studied versus Ising chain size N = 1–8 (Figs. 7–8). — [QuTiP 5 paper §3.2.10](https://arxiv.org/abs/2412.04705); [notebook 025](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/025_v5_paper-floquet-speed-test.ipynb)
- **Parallelism:**
  - v5 supports GPU through JAX data layers and an in-development qutip-cuquantum layer (multi-GPU, with NVIDIA)
  - multicore through `parallel_map` and trajectory solvers
  - MKL-backed linear algebra
  - MPI (new)
  - "Historically, we also supported OpenMP, but this support was removed with version 5"

  — [QuTiP 5 paper §6](https://arxiv.org/abs/2412.04705)
- Data layer: v5 adds `Dense` and `Dia` formats besides CSR, plus JAX `jax` (dense) and `jaxdia` formats through qutip-jax; "An experimental CSR format is available within JAX, but currently not yet supported by QuTiP." — [QuTiP 5 paper §1.1.1, §3.2.5](https://arxiv.org/abs/2412.04705)
- QuTiP maintains a separate benchmark repository; it was last pushed 2026-08-25 and I did not inspect it. — [qutip-benchmark](https://github.com/qutip/qutip-benchmark)

**QuantumOptics.jl (older data)**
- The 2018 paper benchmarks master-equation time evolution against QuTiP 4.2.0 and the Matlab QO Toolbox on three systems with different sparsity: a pumped cavity with decay (sparse H, dense ρ), a JC model with atom and cavity decay (sparse H, sparse ρ), and a particle in a harmonic trap (dense H, dense ρ). "Depending on the sparseness of the system, QuantumOptics.jl's flexible operator types can lead to considerable speed-ups compared to the purely sparse matrix approach in QuTiP and the QO Toolbox." Time-dependent versions were compared against QuTiP's Cython and pure-Python paths, excluding compile times. Runs were single-core on an Intel i7-5960X at 3.00 GHz. — [QuantumOptics.jl paper Figs. 6–7](https://arxiv.org/abs/1707.01060)
- The benchmark suite repository includes cavity, JC and particle variants of master-equation, time-dependent master, MCWF and Schrödinger evolution (including FFT particle), plus operator kernels: ptrace, expect, variance, Wigner, Q-function, coherent state, displace, dense/sparse arithmetic. It compares against QuTiP and the QO Toolbox. It was **last pushed 2020-11-24**, so the information is stale. — [QuantumOptics.jl-benchmarks](https://github.com/qojulia/QuantumOptics.jl-benchmarks)

**HierarchicalEOM.jl**
- It "achieves a significant speedup with respect to the corresponding method in … QuTiP". The benchmark (Fig. 7) covers building the HEOM Liouvillian, ADO time evolution and stationary states against QuTiP-BoFiN on the single-impurity-Anderson-model example. It shows single-thread gains plus multithreaded assembly; QuTiP-BoFiN did not support importance-threshold truncation at that time. — [HierarchicalEOM.jl paper, arXiv:2306.07522](https://arxiv.org/abs/2306.07522)

**Dynamiqs**
- Batching tutorial: the batched run takes 44.1 ms against 2.59 s for the unbatched loop, "about 60 times faster … even more significant for larger simulations, or when using a GPU". The device is not stated; the wording implies the timing is not from a GPU. — [Dynamiqs batching tutorial source](https://github.com/dynamiqs/dynamiqs/blob/main/docs/documentation/basics/batching-simulations.md)
- The README claims SME trajectories can be simulated "orders of magnitude faster" by batching over trajectories. It also offers a custom sparse diagonal (`dia`) format "offering substantial speedups for large systems", and parallelization across CPUs/GPUs. — [Dynamiqs README](https://github.com/dynamiqs/dynamiqs)
- The v0.3.0 release notes say the DIA sparse format "offers significant speedups for large systems" and link a November 2024 Alice & Bob blog post. — [dynamiqs v0.3.0 release](https://github.com/dynamiqs/dynamiqs/releases/tag/v0.3.0)
- **Unverified:** a search-engine snippet of that blog states v0.2.2 "is at least on par with QuTiP and features up to a 30x speed increase on the dissipative cat CNOT benchmark by leveraging GPU acceleration". The blog URL returned 404 on 2026-10-04. — [Alice & Bob blog (404)](https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/)
- Dynamiqs ships an internal benchmark suite (`task bench --tier physics`) for comparing runs; the README's "Coming soon" lists "Benchmark code to compare solvers and performance for different systems". — [Dynamiqs performance page source](https://github.com/dynamiqs/dynamiqs/blob/main/docs/documentation/basics/performance.md); [README](https://github.com/dynamiqs/dynamiqs)

### Inferences
- For exact Lindblad simulation, the practical ceiling in all four libraries on one CPU or GPU is **about 10–12 qubits**, i.e. Liouville dimension 4¹⁰–4¹² ≈ 10⁶–1.7×10⁷, with times of 10³ s. In single-mode bosonic problems the ceiling is Fock cutoffs of a few hundred. Beyond that, the incumbents switch to trajectories (mcsolve on clusters), permutation symmetry (PIQS), cumulants (QuantumCumulants.jl), low-rank methods (QuantumToolbox.jl and Dynamiqs `LowRank`) or HEOM pruning.
- The speed differences at small sizes (tens of milliseconds) are dominated by per-call overhead. At large sizes they are within roughly 2× across libraries on CPU. The large differences come from GPU use, batching and trajectory parallelism, not from the core Lindblad right-hand side. For pedagogical master's-level problems (N ≲ 100 Fock states, ≲ 8 qubits), every incumbent finishes in seconds, so raw speed is not decisive there.
- All cross-library numbers are produced by an interested party: QuantumToolbox.jl's own paper, QuantumOptics.jl's own paper, HierarchicalEOM.jl against BoFiN, and Dynamiqs's blog. The Dynamiqs-"mcsolve" entry in the QuantumToolbox.jl benchmark compares a fixed-step jump-SSE with adaptive MCWF. Its 6.8 s figure therefore says more about method choice than about the library.
- QuTiP issue #2733 ("smesolve more than 10x slower compared to v5.1.1", closed 2025-08-22) was filed with the QuantumToolbox.jl benchmark script as its reproducer. The smesolve number for QuTiP 5.2.1 may therefore reflect a performance regression that was being fixed at the time of the benchmark. This is uncertain: I did not determine whether 5.2.1 contains the fix. — [qutip #2733](https://github.com/qutip/qutip/issues/2733)

### Gaps
- I did not extract the QuTiP 5 paper's absolute GPU/CPU timings (Fig. 4 is only available as plotted data).
- I did not inspect the contents of the QuantumToolbox.jl continuous-benchmark page or the qutip-benchmark repository.
- No published Dynamiqs paper with benchmarks was found ("in preparation"). The Alice & Bob blog containing the 30× cat-CNOT claim could not be retrieved.
- I found no independent, third-party head-to-head benchmark of all four libraries.

---

## Q3. What limitations are documented in docs, issues and papers?

### Takeaway
The recurring documented limitations are:
1. Short pulses and discontinuous drives are silently skipped by adaptive integrators unless the user sets `max_step` (QuTiP) or declares discontinuity times (Dynamiqs).
2. Steady-state solvers are fragile: multiple steady states are unsupported in QuTiP, and iterative solvers give seed-dependent or spurious results in QuantumOptics.jl for nonlinear oscillators. Memory blows up for large sparse problems.
3. Precision: Dynamiqs defaults to single precision. QuTiP and JAX are locked to hardware floats. QuantumToolbox.jl now shows Float64 failing by about 10 orders of magnitude on exponentially small quantities (tunnelling splittings, Liouvillian gaps) and offers Double64/BigFloat.
4. Fock truncation is left to the user, except for QuantumToolbox.jl's DFD/DSF and QuTiP's ENR states.
5. Non-Markovian modelling needs exponential bath fits (HEOM) and fails for long-memory kernels (an open QuTiP ESPRIT issue).
6. Autodiff is first-class only in Dynamiqs and JAX-QuTiP; it is "experimental" in QuantumToolbox.jl and undocumented in QuantumOptics.jl.
7. GPU memory caps full Lindblad at about 11–12 qubits.

### Cited Findings

**Time-dependent Hamiltonians, pulses and discontinuities**
- QuTiP time dependence can be given "as a Python function, as a string (which will be compiled into machine code at the first usage), or as discrete time-dependent data which will be interpolated with cubic splines". — [QuTiP 5 paper §3.2.4](https://arxiv.org/abs/2412.04705)
- Array coefficients use a cubic spline by default (`order=3`); "When `order = 0`, the interpolation is step function". — [qutip/core/coefficient.py](https://github.com/qutip/qutip/blob/master/qutip/core/coefficient.py)
- QuTiP: "The max_step option is often important in time-dependent problems with periods of idling interspersed with short pulses; without setting a maximum time step for the solver to take, these short pulses might be ignored when the ODE solver takes too large time steps." — [QuTiP 5 paper §3.2.2](https://arxiv.org/abs/2412.04705)
- QuTiP integrator defaults: Adams (scipy zvode) is mesolve's default method, with `atol=1e-8`, `rtol=1e-6`, `nsteps=2500` and `max_step=0` (automatic). BDF, LSODA, DOP853 and Vern methods are also available. — [scipy_integrator.py](https://github.com/qutip/qutip/blob/master/qutip/solver/integrator/scipy_integrator.py); [mesolve.py](https://github.com/qutip/qutip/blob/master/qutip/solver/mesolve.py)
- Dynamiqs: "Adaptive step-size methods assume that the vector field is smooth. Whenever a Hamiltonian or a jump operator jumps discontinuously, the solver has to shrink the step size … rejecting many steps." Dynamiqs collects the discontinuity times of `dq.pwc()` and `clip()` operators automatically. For functions passed to `dq.modulated()`, the user must pass `discontinuity_ts=[...]`. — [Dynamiqs performance page](https://www.dynamiqs.org/stable/documentation/basics/performance.html)
- Dynamiqs `Event` jump-SSE method warns: "When using adaptive step size solvers for the no-click integration, you must specify either `dtmax` or `root_finder` to control the precision of the click times. Otherwise, the adaptive solver may choose to take very large step sizes, which results in imprecise click times." — [dynamiqs/method.py](https://github.com/dynamiqs/dynamiqs/blob/main/dynamiqs/method.py)
- QuantumOptics.jl time dependence:
  - user functions `H(t, psi)`, which may be state-dependent (the `_dynamic` solvers)
  - `TimeDependentSum`, e.g. `TimeDependentSum(1.0=>H_static, cos=>H_drive1, ...)`

  The docs page fetched does not mention discontinuities. — [QuantumOptics.jl time-dependent docs source](https://github.com/qojulia/QuantumOptics.jl/blob/master/docs/src/timeevolution/timedependent-problems.md)
- QuantumToolbox.jl time dependence: `QobjEvo(op, coef)` with `coef(p, t)`, built on SciMLOperators.jl. The page does not discuss discontinuities. — [QuantumToolbox.jl time-dependent guide source](https://github.com/qutip/QuantumToolbox.jl/blob/main/docs/src/users_guide/time_evolution/time_dependent.md)
- Historical QuTiP issues on time dependence (titles):
  - "Cannot use arrays for time-dependent control fields" (#932, 2018)
  - "QuTiP VS RK45: Which one gives the correct results for time-dependent systems?" (#1733, 2021)
  - "Mesolve function returns wrong result using default Ode method" (#1265, 2020)
  - "incorrect behaviour of correlation_2op_2t when using time-dependent Hamiltonian terms and collapse operators" (#1808, 2022)

  All are closed. — [#932](https://github.com/qutip/qutip/issues/932), [#1733](https://github.com/qutip/qutip/issues/1733), [#1265](https://github.com/qutip/qutip/issues/1265), [#1808](https://github.com/qutip/qutip/issues/1808)

**Stiffness and integrator choice**
- QuTiP offers BDF and LSODA alongside Adams. — [scipy_integrator.py](https://github.com/qutip/qutip/blob/master/qutip/solver/integrator/scipy_integrator.py)
- Dynamiqs offers implicit `Kvaerno3`/`Kvaerno5`. Its adaptive methods default to `rtol=1e-6`, `atol=1e-6` and `max_steps=100_000`. — [dynamiqs/method.py](https://github.com/dynamiqs/dynamiqs/blob/main/dynamiqs/method.py)
- QuantumToolbox.jl defaults to `DP5()` for mesolve, with ODE tolerances `abstol=1e-8`, `reltol=1e-6` and SDE tolerances `abstol=1e-3`, `reltol=2e-3`. Any DifferentialEquations.jl algorithm is accepted. — [time_evolution.jl](https://github.com/qutip/QuantumToolbox.jl/blob/main/src/time_evolution/time_evolution.jl); [mesolve.jl](https://github.com/qutip/QuantumToolbox.jl/blob/main/src/time_evolution/mesolve.jl)
- QuantumOptics.jl defaults to `DP5()` with `reltol=1e-6`, `abstol=1e-8`. — [timeevolution_base.jl](https://github.com/qojulia/QuantumOptics.jl/blob/master/src/timeevolution_base.jl)

**Numerical precision**
- Dynamiqs: "all objects are represented by default with single-precision floating-point numbers (float32 or complex64)". The docs warn about large numbers and about tight tolerances: simulations "may even get stuck. In such cases, it is recommended to switch to double-precision". They add that "Most GPUs do not have native support for double-precision … some recent NVIDIA GPUs (e.g. V100, A100, H100) do provide efficient support". — [Dynamiqs sharp bits](https://www.dynamiqs.org/stable/documentation/getting_started/sharp-bits.html)
- Dynamiqs can lower matrix-multiply precision to TF32/bfloat16 (`set_matmul_precision('high')`); "the speedup is not guaranteed. Always check the result against a 'highest' run". — [Dynamiqs performance page](https://www.dynamiqs.org/stable/documentation/basics/performance.html)
- Dynamiqs `EQX_ON_ERROR=off` removes runtime checks for "several tens of percent speedup", but "an invalid input silently produces a wrong result rather than an error". — [Dynamiqs performance page](https://www.dynamiqs.org/stable/documentation/basics/performance.html)
- **QuantumToolbox.jl arbitrary precision** (docs page added 2026-07-23, PR #745):
  - States, operators, superoperators, eigensolvers and time-evolution solvers are generic over the number type (`Double64`, `BigFloat`).
  - Worked example: the tunnelling splitting of a quartic double well. "`Float64` reports a splitting of ∼10⁻¹³ while the true value is ∼10⁻²³ — wrong by nearly ten orders of magnitude".
  - It names the Liouvillian gap of a driven-dissipative resonator near bistability as the open-system analogue.
  - It claims "Python QuTiP cannot do this at all. Its data layer (Dense and CSR) is compiled Cython with double complex hardcoded" and that JAX/XLA supports only hardware float types.

  — [QuantumToolbox.jl arbitrary-precision guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/arbitrary_precision); [commit history](https://github.com/qutip/QuantumToolbox.jl/commits/main/docs/src/users_guide/arbitrary_precision.md)
- Caveats stated on the same QuantumToolbox.jl page:
  - BigFloat costs "one to two orders of magnitude slowdown, and no multithreaded BLAS"; `Double64` is "only a few times slower than Float64" and about 10× faster than BigFloat in the double-well sweep
  - "No GPU support" for arbitrary precision
  - "At N = 150, even BigFloat returns a wrong splitting — the arithmetic is exact but the basis is not"
  - `eigsolve`'s default `tol = 1e-8` is hardcoded
  - only `sesolve`, `mesolve`, `eigenstates` and `eigsolve_al` are tested at high precision; `steadystate`, `mcsolve` and `spectrum` "are not currently tested at high precision"
  - "Extra precision does not generally make a time-evolution result better: the error of an ODE integration is dominated by the solver tolerances"

  — [QuantumToolbox.jl arbitrary-precision guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/arbitrary_precision)

**Steady states**
- QuTiP:
  - methods: `direct` (default), `power`, `eigenvalue` and `svd`
  - linear solvers: `solve`, `lstsq`, `spsolve`, `gmres`, `lgmres`, `bicgstab`, `mkl_spsolve`

  "In the situation where multiple steady states exist, the different solvers can often produce very different results; some may fail, others may produce linear combinations of possibilities … So far, QuTiP does not support an automated approach to dealing with this issue; therefore, the onus is on the user." For time-dependent drives the paper points to `steadystate_floquet`. — [QuTiP 5 paper §3.2.6](https://arxiv.org/abs/2412.04705)
- QuTiP issues:
  - `steadystate_floquet` renamed `steadystate_fourier` (#2632, 2025)
  - "None-CSR sparse formats tends to run out of memory when used with steadystate()+large problem sizes" (#2747, QuTiP 5.3.0.dev; closed 2026-04-15)
  - "Default steadystate method ('direct') fails with singular matrix in piqs.Dicke" (#2649, 2025)
  - still **open** since 2019: "qutip.piqs Dicke.liouvillian() gives wrong result for Liouvillian spectrum and nonlinear functions of steadystate" (#993)

  — [#2632](https://github.com/qutip/qutip/issues/2632), [#2747](https://github.com/qutip/qutip/issues/2747), [#2649](https://github.com/qutip/qutip/issues/2649), [#993](https://github.com/qutip/qutip/issues/993)
- QuantumOptics.jl:
  - open issue #428 (2024-12): for a quartic oscillator, `steadystate.master` and `steadystate.eigenvector` agree, "while the iterative solver fails to find the appropriate nullspace. The result of the iterative solver depends strongly on the random seed"
  - earlier #315 (2021): "Iterative solver for steadystate returning different results each time I run my code"

  — [QO.jl #428](https://github.com/qojulia/QuantumOptics.jl/issues/428); [QO.jl #315](https://github.com/qojulia/QuantumOptics.jl/issues/315)
- QuantumOptics.jl offers `steadystate.master`, `steadystate.eigenvector`, `steadystate.iterative` and `steadystate.liouvillianspectrum`. `steadystate.master` time-evolves until a steady state is reached and returns the time and state lists. — [QO.jl steady-state docs source](https://github.com/qojulia/QuantumOptics.jl/blob/master/docs/src/steadystate.md); [API](https://github.com/qojulia/QuantumOptics.jl/blob/master/docs/src/api.md)
- QuantumToolbox.jl on preconditioning: "The problem with precondioning is that it is only well defined for Hermitian matrices", while the Liouvillian is non-Hermitian. — [QuantumToolbox.jl steady-state guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/steadystate)
- QuantumToolbox.jl open issue "GPU support for `steadystate_fourier`" (#462); earlier "Request for Steady State Sparse GPU Solving" (#417, closed 2025-02). — [#462](https://github.com/qutip/QuantumToolbox.jl/issues/462); [#417](https://github.com/qutip/QuantumToolbox.jl/issues/417)
- Dynamiqs: no steady-state API (see Q1). — [API index](https://github.com/dynamiqs/dynamiqs/blob/main/docs/python_api/index.md)

**Correlation functions and spectra**
- QuTiP issue #1460 (2021): "correlation_2op_1t() and spectrum() giving the physically incorrect result unless solver='es'". The reporter asked to make `es` the default and to warn in the docs. — [qutip #1460](https://github.com/qutip/qutip/issues/1460)
- The QuTiP 5 guide says `spectrum()` uses the exponential-series approach by default. For FFT spectra: "we need to make sure to evaluate the correlation function for a sufficient long time and sufficiently high sampling rate so that the discrete Fourier transform (FFT) captures all the features". — [QuTiP correlation guide](https://qutip.readthedocs.io/en/latest/guide/guide-correlation.html)
- QuTiP issue "Bug in function spectrum_correlation_fft or wrong documentation" (#1537, 2020). — [#1537](https://github.com/qutip/qutip/issues/1537)

**Fock-space truncation and convergence**
- QuTiP "does not by default apply any approximation apart from truncating the individual systems' Hilbert spaces". ENR states can restrict the total excitation number. — [QuTiP 5 paper §3.3.1](https://arxiv.org/abs/2412.04705)
- QuantumToolbox.jl's DFD removes the need "to check for the convergence of the results as a function of a fixed cutoff", and DSF handles strongly driven systems. — [QuantumToolbox.jl paper](https://arxiv.org/abs/2504.21440)
- Rigorous, computable a-posteriori bounds for both Hilbert-space truncation and time-discretization errors of Lindblad simulations were published in 2026. They enable "fully adaptive simulations". The paper cites Dynamiqs, but there is no indication the method is implemented in any of the four libraries. — [Etienney, Robin, Rouchon, Quantum 10, 2031 (2026)](https://quantum-journal.org/papers/q-2026-03-16-2031/)

**Non-Markovian memory**
- HEOM requires bath correlation functions as sums of exponentials. At zero temperature these "are non-exponential … requiring us to apply a multi-exponential fit". Cost is controlled by `max_depth` and the number of exponents Nk. — [QuTiP 5 paper §3.2.4, §3.2.12](https://arxiv.org/abs/2412.04705)
- Open QuTiP issue #2950 (2026-07): "QuTiP's ESPRIT function does not provide the expected fitting results for bath correlation functions with strong non-Markovian effects, such as the giant-atom correlation function … does not capture the second peak at the delay time, no matter how many modes we consider. Using ESPIRA instead, the fit works, but it is considerably slower." — [qutip #2950](https://github.com/qutip/qutip/issues/2950)
- Active HEOM-scaling work in QuTiP:
  - "Low-rank factorization for correlated bath HEOM simulations" (#2765, closed 2025-10)
  - "PETSc Backend for HEOM Solver" (#2953, closed 2026-07)
  - "HEOM Solver Refactor with Direct Sum" (#2864, open)
  - "Backend support for heom solver" (#2976, open)

  — [#2765](https://github.com/qutip/qutip/issues/2765), [#2953](https://github.com/qutip/qutip/issues/2953), [#2864](https://github.com/qutip/qutip/issues/2864), [#2976](https://github.com/qutip/qutip/issues/2976)
- The HEOM method "usually results in time-consuming calculations and a large memory cost". — [HierarchicalEOM.jl paper abstract](https://arxiv.org/abs/2306.07522)
- QuTiP's planned solvers include a cumulant expansion, a "non-Markovian master equation expansion based on time-convolutionless truncation", a new Floquet master equation with flexible truncation, a Dyson-expansion solver for quickly driven systems, and Krylov solvers. — [QuTiP 5 paper §6.2](https://arxiv.org/abs/2412.04705)

**Measurement feedback**
- QuTiP 5 exposes feedback hooks to time-dependent coefficients: `MESolver.StateFeedback()`, `MCSolver.CollapseFeedback()`, and `SMESolver.WienerFeedback()`/`StateFeedback()`. For example, `QobjEvo([op, func], args={"W": SMESolver.WienerFeedback()})`. — [mesolve.py](https://github.com/qutip/qutip/blob/master/qutip/solver/mesolve.py); [mcsolve.py](https://github.com/qutip/qutip/blob/master/qutip/solver/mcsolve.py); [stochastic.py](https://github.com/qutip/qutip/blob/master/qutip/solver/stochastic.py)
- QuantumOptics.jl `_dynamic` solvers accept state-dependent `H(t, psi)`. — [QO.jl time-dependent docs source](https://github.com/qojulia/QuantumOptics.jl/blob/master/docs/src/timeevolution/timedependent-problems.md)

**Gradients and autodiff**
- QuTiP: autodiff only through the qutip-jax data layer with the Diffrax integrator. Paper examples are counting statistics and a qubit coupled to a mirror-terminated waveguide. — [QuTiP 5 paper §3.3.2](https://arxiv.org/abs/2412.04705)
- QuantumToolbox.jl: "preliminary support for automatic differentiation … this functionality is considered experimental and not all parts of the library are AD-compatible". It works with ForwardDiff, Zygote, Enzyme and Mooncake through SciMLSensitivity. — [QuantumToolbox.jl autodiff guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/autodiff)
- Dynamiqs: gradients are first-class:
  - several gradient methods, including Diffrax's "optimal online checkpointing"
  - "machine-precision accuracy"
  - derivatives with respect to evolution time
  - higher-order derivatives (Hessian)

  Sharp bit: for complex parameters, step along the conjugate of the JAX gradient. — [Dynamiqs README](https://github.com/dynamiqs/dynamiqs); [sharp bits](https://www.dynamiqs.org/stable/documentation/getting_started/sharp-bits.html)
- QuantumOptics.jl: no AD or GPU mention on the docs index, tutorial or time-evolution pages I checked. An open issue (2025-05) asks for "Correctness checks and unit tests for GPU arrays used as storage backend". — [QO.jl #445](https://github.com/qojulia/QuantumOptics.jl/issues/445)

**GPU and memory ceilings**
- QuTiP-JAX on an 80 GB A100 hits its memory limit at 11 spins for mesolve; distributing across GPU nodes "is a challenging task that we plan to explore". — [QuTiP 5 paper §3.2.5](https://arxiv.org/abs/2412.04705)
- Dynamiqs suggests unified memory (`TF_FORCE_UNIFIED_MEMORY=1`, `XLA_PYTHON_CLIENT_MEM_FRACTION=4`) to exceed GPU memory, "trad[ing] speed for capacity". — [Dynamiqs performance page](https://www.dynamiqs.org/stable/documentation/basics/performance.html)

**Parameter sweeps**
- Dynamiqs batching over Hamiltonians, initial states or jump operators gives about 60× on CPU in its tutorial. It warns against `for` loops in the sharp bits. — [batching tutorial](https://github.com/dynamiqs/dynamiqs/blob/main/docs/documentation/basics/batching-simulations.md); [sharp bits](https://www.dynamiqs.org/stable/documentation/getting_started/sharp-bits.html)
- QuantumToolbox.jl provides `mesolve_map`/`sesolve_map`. — [API](https://github.com/qutip/QuantumToolbox.jl/blob/main/docs/src/resources/api.md)
- QuTiP relies on `parallel_map` or MPI and solver reuse: "the speed-up can be significant if the solver is reused many times". — [QuTiP 5 paper §3.2.1, §6](https://arxiv.org/abs/2412.04705)

**API stability and maintenance**
- Dynamiqs: "under active development and some APIs and solvers are still finding their footing … new releases might introduce breaking changes". — [README](https://github.com/dynamiqs/dynamiqs)
- QuTiP: "there is an inevitability that some features that are not part of QuTiP-core may become abandoned or not well maintained". — [QuTiP 5 paper §6](https://arxiv.org/abs/2412.04705)
- QuTiP: the tensor-network data layer "exists currently in a very early alpha form". — [QuTiP 5 paper §6](https://arxiv.org/abs/2412.04705)
- QuantumOptics.jl's 2018 paper cites Julia compile time (LLVM against Cython/GCC) and says the package is "clearly not as feature-rich as other well-established frameworks such as QuTiP". This is a 2018 statement and partly outdated. — [QuantumOptics.jl paper §6](https://arxiv.org/abs/1707.01060)

### Inferences
- **Precision is now a contested axis.** QuantumToolbox.jl explicitly markets arbitrary precision, with a ten-orders-of-magnitude Float64 failure as its showcase, and claims QuTiP and JAX cannot follow. Its own caveats leave room for competition:
  - only the time-evolution and eigen-solvers are tested at high precision
  - steadystate, mcsolve and spectrum are untested
  - there is no GPU support
  - the default tolerances are Float64-hardcoded

  A competitor would have to beat it on robustness and verification (certified digits, convergence in Fock cutoff *and* precision), not just offer extended precision.
- **Discontinuous drives** (gates, pulse sequences, dynamical-decoupling trains) are a documented trap: QuTiP needs `max_step` and Dynamiqs needs `discontinuity_ts`. The Julia libraries leave it to DifferentialEquations.jl options, apparently undocumented in their user guides. Event-aware integration is a recognized pain point.
- **Steady-state reliability** for nonlinear, bistable or degenerate Liouvillians is weak across the board:
  - QuTiP: no multiple-steady-state handling; PIQS Liouvillian-spectrum bug open since 2019
  - QuantumOptics.jl: iterative solver seed-dependent
  - Dynamiqs: no steady-state solver

  That makes driven-dissipative bistability (Kerr, Dicke) a numerically decisive arena.
- Only Dynamiqs treats single precision as the default. Any head-to-head with Dynamiqs on default settings will mix precision effects with algorithm effects.

### Gaps
- I did not inspect QuantumToolbox.jl or QuantumOptics.jl issues on discontinuities or `tstops`, or docs pages on stiffness.
- QuantumOptics.jl's current GPU and AD status beyond issue #445 is undocumented in the pages checked.
- I found no explicit documentation of QuTiP's numerical-precision guarantees beyond the hard-wired double-complex data layer, and that claim comes from QuantumToolbox.jl's docs, not QuTiP's.
- I did not verify whether Dynamiqs ever shipped a steady-state solver in the PyTorch era.

---

## Q4. Which standard, named benchmark problems does each library use (for head-to-head comparison)?

### Takeaway
Six problems recur with published parameters and code:

| # | Problem | Main published source |
|---|---|---|
| 1 | Driven-dissipative Kerr oscillator | QuantumToolbox.jl four-library benchmark: N = 50 Fock states, explicit parameters and scripts |
| 2 | Dissipative transverse-field Ising chain | QuantumToolbox.jl and QuTiP 5 GPU scaling |
| 3 | Damped Jaynes–Cummings / vacuum Rabi | Present in all four; QuantumOptics.jl benchmark suite |
| 4 | Homodyne-monitored cavity | QuTiP, QuantumToolbox.jl and Dynamiqs SME examples |
| 5 | Spin-boson, FMO and Anderson impurity with HEOM | QuTiP-BoFiN and HierarchicalEOM.jl |
| 6 | Single-qubit Hadamard and two-qubit CNOT/QFT gate optimization | QuTiP QOC, QuantumToolbox.jl paper |

QuantumToolbox.jl's arbitrary-precision quartic double-well tunnelling splitting is a new precision benchmark with a published Float64 failure.

### Cited Findings

**Recognized benchmark problems with exact sources**

| Problem | Library / source | Exact location | Published parameters / outputs |
|---|---|---|---|
| Driven-dissipative Kerr oscillator: mesolve, mcsolve, smesolve, AD | QuantumToolbox.jl paper | §6, Fig. 8(a–d), Eq. (31) — [arXiv:2504.21440](https://arxiv.org/abs/2504.21440); scripts [qutip_benchmarks.py](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/python/qutip_benchmarks.py), [quantumtoolbox.jl](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures/blob/main/src/benchmarks/julia/quantumtoolbox.jl) | N = 50 (AD: N = 100), Δ = 0.1, U = −0.05, F = 2, γ = 1, n_th = 0.2, t ∈ [0, 10], 100 points, 100 trajectories; times in Q2 |
| Driven-dissipative Kerr oscillator (tutorial) | Dynamiqs | [Kerr example](https://www.dynamiqs.org/stable/documentation/advanced_examples/kerr-oscillator.html) | Time evolution, Wigner function, periodic revival of a coherent state |
| Kerr nonlinearity / cat states | QuTiP | [Lecture 14 Kerr](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/lectures/Lecture-14-Kerr-nonlinearities.ipynb); [018 cats become coherent (stochastic vs MC)](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/018_measures-trajectories-cats-kerr.ipynb) | Not extracted |
| Dissipative 1D Ising chain (mesolve CPU/GPU scaling) | QuantumToolbox.jl paper | §6, Fig. 8(e–f), Eq. (32) — [arXiv:2504.21440](https://arxiv.org/abs/2504.21440) | J_x = 25, h_z = 50, γ = 1, N = 2–10 (CPU), 2–12 (GPU) |
| 1D Ising chain, sesolve/mesolve on GPU | QuTiP 5 paper | §3.2.5, Fig. 4 — [arXiv:2412.04705](https://arxiv.org/abs/2412.04705) | A100 80 GB; to 22 spins (sesolve), 11 spins (mesolve) |
| Periodically driven Ising chain, Floquet vs sesolve | QuTiP 5 paper | §3.2.10, Figs. 7–8; [notebook 025](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/025_v5_paper-floquet-speed-test.ipynb) | Crossover time t_cross vs N (1–8); A = 2.5×2π, ω_d = 2π (qubit case) |
| 2D transverse-field Ising 4×4 with local decay (mcsolve, cluster) | QuantumToolbox.jl paper | §4 — [arXiv:2504.21440](https://arxiv.org/abs/2504.21440) | 720 threads, SLURM |
| Lossy Jaynes–Cummings model (and time-dependent JC) | QuantumOptics.jl paper | §4.1–4.2, Figs. 3, 6–7 — [arXiv:1707.01060](https://arxiv.org/abs/1707.01060); suite [`timeevolution_master_jaynescummings.jl` etc.](https://github.com/qojulia/QuantumOptics.jl-benchmarks) | vs QuTiP 4.2.0 / QO Toolbox |
| Pumped cavity; particle in harmonic trap | QuantumOptics.jl paper | Figs. 6–7 — [arXiv:1707.01060](https://arxiv.org/abs/1707.01060) | Sparse/dense contrasts |
| Vacuum Rabi oscillations, damped JC | QuTiP | [Lecture 1](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/lectures/Lecture-1-Jaynes-Cumming-model.ipynb); [time-evolution 004](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/004_rabi-oscillations.ipynb) | Not extracted |
| JC emission spectrum (vacuum Rabi splitting); g⁽¹⁾/g⁽²⁾ | QuTiP | [correlation guide](https://qutip.readthedocs.io/en/latest/guide/guide-correlation.html); [Lecture 4](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/lectures/Lecture-4-Correlation-Functions.ipynb) | Not extracted |
| Resonance fluorescence (Mollow) | QuTiP | [Lecture 13](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/lectures/Lecture-13-Resonance-flourescence.ipynb) | Not extracted |
| Correlation spectrum | QuantumOptics.jl | [correlation-spectrum notebook](https://github.com/qojulia/QuantumOptics.jl/tree/master/docs/examples/notebooks) | Not extracted |
| Homodyne-detected cavity (SME) | QuTiP 5 paper | §3.2.11, Eq. (34); [notebook 022](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/022_v5_paper-smesolve.ipynb) | Compared with mesolve |
| Homodyne detection (SSE/SME) | QuantumToolbox.jl | paper §3.5 ([arXiv](https://arxiv.org/abs/2504.21440)) | Not extracted |
| Jump / diffusive measurement (qubit, oscillator) | Dynamiqs | [jump](https://www.dynamiqs.org/stable/documentation/advanced_examples/continuous-jump-measurement.html), [diffusive](https://www.dynamiqs.org/stable/documentation/advanced_examples/continuous-diffusive-measurement.html) | Not extracted |
| Photon birth/death in a cavity (Monte Carlo, Haroche-type) | QuTiP | [time-evolution 006](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/006_photon_birth_death.ipynb); [Lecture 6](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/lectures/Lecture-6-Quantum-Monte-Carlo-Trajectories.ipynb) | Not extracted |
| Damped JC with non-Markovian (negative-rate) environment | QuTiP 5 paper | §3.2.8; [notebook 026](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/026_v5_paper-nm_mcsolve.ipynb) | Compared with exact |
| Driven qubit with Lindblad vs Bloch–Redfield vs HEOM; adiabatically switched qubit | QuTiP 5 paper | §3.2.4, Fig. 3 — [arXiv](https://arxiv.org/abs/2412.04705) | Naive Lindblad fails in 3b |
| Spin-boson HEOM (Drude–Lorentz, strong coupling, underdamped, Ohmic fit, pure dephasing); T = 0 localization | QuTiP | [HEOM 1a](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/heom/heom-1a-spin-bath-model-basic.ipynb)–[1e](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/heom/heom-1e-spin-bath-model-pure-dephasing.ipynb); paper §3.2.12 | Not extracted |
| FMO energy transfer; DD pulse-spacing comparison | QuTiP-BoFiN | [HEOM 2 FMO](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/heom/heom-2-fmo-example.ipynb); [HEOM 4 DD](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/heom/heom-4-dynamical-decoupling.ipynb); [BoFiN paper](https://arxiv.org/abs/2010.10806) | Uhrig vs equal spacing |
| Single-impurity Anderson model (fermionic HEOM) | QuTiP-BoFiN / HierarchicalEOM.jl | [HEOM 5a](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/heom/heom-5a-fermions-single-impurity-model.ipynb); [HierarchicalEOM.jl paper](https://arxiv.org/abs/2306.07522) | Speed vs BoFiN (Fig. 7) |
| Steady state, optomechanics in single-photon strong coupling | QuTiP | [notebook 019](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/time-evolution/019_optomechanical-steadystate.ipynb) | Large Fock space |
| Superradiance, steady-state superradiance, open Dicke, boundary time crystals | QuTiP PIQS | [superradiant emission](https://nbviewer.jupyter.org/github/qutip/qutip-notebooks/blob/master/examples/piqs-superradiant-light-emission.ipynb); [PIQS paper](https://arxiv.org/abs/1805.05129) | Legacy notebooks |
| Superradiant laser | QuantumOptics.jl | [notebook](https://github.com/qojulia/QuantumOptics.jl/tree/master/docs/examples/notebooks); QuantumCumulants.jl ([arXiv](https://arxiv.org/abs/2105.01657)) | Not extracted |
| Hadamard gate with GRAPE/CRAB/GOAT/JOPT | QuTiP-QOC | QuTiP 5 paper §4.1, Fig. 18; [notebook](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/miscellaneous/v5_paper-optimal-control.ipynb) | Fidelity > 0.99 |
| GRAPE Hadamard / QFT / CNOT / Lindbladian; CRAB state transfer / QFT | QuTiP (qtrl) | [OC notebooks 02–08](https://qutip.org/qutip-tutorials/) | Not extracted |
| Hadamard + CNOT pulses, 200 piecewise-constant parameters, Adam | QuantumToolbox.jl paper | §5 ([arXiv](https://arxiv.org/abs/2504.21440)) | Not extracted |
| 10-qubit QFT, randomized benchmarking, Deutsch–Jozsa at pulse level | qutip-qip | [10-qubit QFT](https://nbviewer.org/urls/qutip.org/qutip-tutorials/tutorials-v5/pulse-level-circuit-simulation/qip-10-qubit-QFT-algorithm.ipynb); [paper](https://arxiv.org/abs/2105.09902) | Not extracted |
| Quartic double-well tunnelling splitting (Float64 vs Double64 vs BigFloat) | QuantumToolbox.jl | [arbitrary-precision guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/arbitrary_precision) | Float64 ~10⁻¹³ vs true ~10⁻²³; N = 150 cutoff insufficient |
| Liouvillian gap near bistability | QuantumToolbox.jl | same page | Float64 noise floor |
| Dissipative cat CNOT | Dynamiqs (blog) | [blog, 404](https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/) | "up to 30x" claim, unverified |

- The QuTiP 5 paper's figures are reproducible through dedicated "QuTiPv5 Paper Example" notebooks (022, 024, 025, 026, the JAX and QOC notebooks) and a static GitHub code archive. — [QuTiP 5 paper abstract](https://arxiv.org/abs/2412.04705); [gallery](https://qutip.org/qutip-tutorials/)
- Steady-state methods for common quantum-optics problems were benchmarked comprehensively in P. D. Nation (2015), cited as ref. [34] in the QuTiP 5 paper and ref. [64] in the QuantumToolbox.jl paper. — [QuTiP 5 paper §3.2.6](https://arxiv.org/abs/2412.04705); [QuantumToolbox.jl paper references](https://arxiv.org/abs/2504.21440)

### Inferences
- The **driven-dissipative Kerr oscillator** is the best single head-to-head anchor. It has explicit parameters, open scripts and result JSONs for all four incumbents and appears in QuTiP's and Dynamiqs's own tutorials. Its physics (bistability, cat states, Liouvillian gap, switching) also extends naturally into regimes where truncation and precision decide the answer.
- The **dissipative Ising chain** is the standard *scaling* benchmark, but it rewards raw throughput and GPU memory more than numerical judgement.
- The **damped JC / vacuum Rabi spectrum** is the most "textbook-recognized" comparison. Its numerically delicate part is the emission spectrum, given the historical QuTiP `spectrum` bug (#1460) and the FFT-sampling caveats.
- Quoting a benchmark should include tolerances, precision (Dynamiqs is single precision by default) and the method actually used (e.g. fixed-step EulerJump against adaptive MCWF), because the published comparisons mix these.

### Gaps
- I did not extract the exact parameters of the QuTiP gallery notebooks (Kerr, Mollow, optomechanics, HEOM) or the Dynamiqs Kerr tutorial.
- I could not verify the Dynamiqs cat-CNOT benchmark definition.
- The QuantumToolbox.jl Liouvillian-gap example's parameters were not extracted (the section text was truncated in my read).

---

## Q5. Which of these problem classes is most likely to compel a master's-level quantum mechanics audience, and why?

### Takeaway
Everything in this section is inference grounded in the cited evidence. The strongest candidates are problems where the textbook model is familiar but numerics change the *qualitative* answer, and where incumbents show documented weak points:
1. Driven-dissipative Kerr bistability and the Liouvillian gap, including cat-state formation and decoherence.
2. Exponentially small tunnelling splittings, where Float64 is wrong by about 10 orders of magnitude.
3. Lindblad vs Bloch–Redfield vs HEOM on a driven or strongly coupled qubit, and spin-boson localization.
4. Damped Jaynes–Cummings / resonance-fluorescence spectra, computed from correlation functions.
5. Quantum trajectories under photon counting and homodyne detection.

Optimal control (GRAPE/CRAB) and GPU-scale spin chains matter to practitioners. For a master's audience, though, they are less about *why the numbers decide the physics*.

### Cited Findings (evidence underpinning the ranking)
- **Kerr/bistability:**
  - Kerr is the shared four-library benchmark. — [QuantumToolbox.jl paper §6](https://arxiv.org/abs/2504.21440)
  - QuTiP's gallery has a Kerr lecture and a "cats become coherent" trajectory notebook. — [gallery](https://qutip.org/qutip-tutorials/)
  - QuantumToolbox.jl docs state that the Liouvillian gap "becomes exponentially small near a bistability, and shrinks further with system size. Computed in Float64 it eventually stops being physics". — [arbitrary-precision guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/arbitrary_precision)
  - Steady-state solvers are fragile on nonlinear oscillators. — [QO.jl #428](https://github.com/qojulia/QuantumOptics.jl/issues/428)
  - QuTiP does not handle multiple steady states automatically. — [QuTiP 5 paper §3.2.6](https://arxiv.org/abs/2412.04705)
- **Tunnelling splitting:**
  - Float64 reports ~10⁻¹³ where the truth is ~10⁻²³; the Float64 curve "claims that beyond x₀ ≈ 3 the barrier stops mattering … That is not a small quantitative error; it is a qualitatively wrong statement".
  - The Fock-cutoff convergence trap: at N = 150 even BigFloat is wrong.

  — [QuantumToolbox.jl arbitrary-precision guide](https://qutip.org/QuantumToolbox.jl/stable/users_guide/arbitrary_precision)
- **Validity of master equations:**
  - The naive Lindblad model "fails and is insensitive to the drive" for the adiabatically switched qubit, while Bloch–Redfield and HEOM capture it.
  - HEOM shows the T = 0 spin-boson localization transition.

  — [QuTiP 5 paper §3.2.4, §3.2.12](https://arxiv.org/abs/2412.04705)
- **Spectra:**
  - The JC vacuum-Rabi-splitting spectrum and g⁽²⁾ statistics are QuTiP guide examples, with documented FFT sampling caveats. — [QuTiP correlation guide](https://qutip.readthedocs.io/en/latest/guide/guide-correlation.html)
  - Historical "physically incorrect" spectra unless `solver='es'`. — [qutip #1460](https://github.com/qutip/qutip/issues/1460)
- **Trajectories and measurement:**
  - MCWF unravelling "sheds light on the physical interpretation of the Lindblad master equation" and saves memory for large systems. — [QuTiP 5 paper §3.2.7](https://arxiv.org/abs/2412.04705)
  - Dynamiqs and QuantumToolbox.jl both feature jump and diffusive measurement examples. — [Dynamiqs jump measurement](https://www.dynamiqs.org/stable/documentation/advanced_examples/continuous-jump-measurement.html); [QuantumToolbox.jl paper §3.5](https://arxiv.org/abs/2504.21440)
- **Collective effects:** superradiance and the open Dicke model need permutational symmetry (PIQS) or cumulants to scale. — [PIQS paper](https://arxiv.org/abs/1805.05129); [QuantumCumulants.jl](https://arxiv.org/abs/2105.01657)

### Inferences
- **Ranking for a master's audience**, by how numerically decisive and pedagogically clear each class is:
  1. **Driven-dissipative Kerr oscillator: bistability, Liouvillian gap, switching time, cat states.** Students know the harmonic oscillator. The numerics decide whether a steady state is unique, how slow switching is (an exponentially small gap), and whether the Fock cutoff is converged. There is a published cross-library baseline to compare against, and a documented Float64 failure mode.
  2. **Exponentially small quantities in closed systems (double-well tunnelling splitting).** This is a pure demonstration that "the arithmetic decides the physics". A competitor now publishes the failure explicitly, so a head-to-head is on recognized ground. The convergence-in-basis vs convergence-in-precision lesson is valuable at master's level.
  3. **Which master equation? (Lindblad vs Bloch–Redfield vs HEOM).** The adiabatically switched qubit and the spin-boson model show approximations failing *qualitatively*. HEOM convergence (hierarchy depth, number of exponents, fit quality) is a rich numerical-judgement topic. The cost is high and requires HEOM machinery.
  4. **Spectra and photon statistics** (vacuum Rabi splitting, Mollow triplet, g⁽²⁾ antibunching). These connect to experiment and to the quantum regression theorem. Accuracy depends on the method (exponential series vs FFT sampling vs pseudo-inverse), with documented historical errors.
  5. **Quantum trajectories and continuous measurement** (photon birth/death, homodyne SME, cat decoherence seen through trajectories). These are conceptually powerful. The numerics centre on stochastic-integrator order, the number of trajectories and click-time precision (the Dynamiqs `Event` caveat).
  6. Lower priority for this audience: GRAPE/CRAB gate optimization (more engineering than physics); GPU-scale Ising chains (throughput rather than insight); PIQS superradiance (compelling physics, but it needs specialized symmetry machinery).
- Problems where incumbents' **documented weaknesses** coincide with recognized benchmarks are where a newcomer's comparison would be most informative:
  - steady states of bistable or degenerate Liouvillians
  - exponentially small gaps and splittings
  - pulse trains with discontinuities
  - spectra from correlation functions

### Gaps
- I found no survey of master's-level course usage (which QuTiP lectures are assigned in courses). The audience ranking above is reasoning, not data.
- I did not check whether QuTiP or QuantumOptics.jl publish worked Liouvillian-gap or metastability examples beyond QuantumToolbox.jl's arbitrary-precision page.
