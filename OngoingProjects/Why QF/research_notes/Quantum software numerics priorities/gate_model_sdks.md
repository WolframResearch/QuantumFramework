# Gate-model quantum SDKs: Qiskit (Aer, Dynamics, Experiments), Cirq + qsim, PennyLane. Numerical capabilities, featured problems, documented limits and named benchmarks (state as of 2026-10-04)

Method note: I verified facts on 2026-10-04 against primary sources: GitHub repositories (read through the `gh` API: releases, READMEs, release notes, PRs, doc sources and demo sources), official documentation pages and arXiv abstracts. Everything under "Cited Findings" comes from a linked source. Everything under "Inferences" is my reasoning. The Wolfram Quantum Framework is not evaluated here.

---

## 1. Which numerical problem classes do the official tutorials, demo galleries and papers feature, and how well are they served?

### Takeaway
The three ecosystems teach largely the same NISQ canon at toy scale (2–12 qubits): noisy-circuit simulation with Kraus/Pauli/thermal-relaxation noise, VQE on H2/LiH/H3+/HeH+ in STO-3G, QAOA MaxCut, gradients (parameter-shift, adjoint, backprop, SPSA, QNG), barren plateaus, ZNE, quantum volume/RB/XEB/tomography, and Trotterized spin dynamics. As of October 2026 the pulse-level and time-dependent-Hamiltonian layer is being abandoned. Qiskit removed Pulse in 2.0. Qiskit Dynamics is archived ("DEPRECATED"). qiskit-experiments dropped its calibration experiments. Aer is in "Reduced Maintenance Mode". PennyLane's `main` branch has removed `pulse`, `noise` (including ZNE), `qaoa`, all qutrit/qudit code and all continuous-variable code ahead of its next release.

### Cited Findings

#### Status snapshot (versions and maintenance, verified 2026-10-04)
- Qiskit SDK: latest release is 2.5.2 (2026-08-13), and 2.5.0 came out 2026-07-02 — [Qiskit releases](https://github.com/Qiskit/qiskit/releases)
- Qiskit 2.0 release note: "Qiskit Pulse has been completely removed in this release, following its deprecation in Qiskit v1.3." This covers pulse-based fake backends and calibrations in `QuantumCircuit`/`Target`. The note also says "Pulse migration to Qiskit Dynamics, as was the initial plan following the deprecation of Pulse, has been put on hold due to Qiskit Dynamics development priorities" — [Qiskit 2.0 note remove-pulse](https://github.com/Qiskit/qiskit/blob/stable/2.0/releasenotes/notes/2.0/remove-pulse-eb43f66499092489.yaml); [remove-pulse-calibrations](https://github.com/Qiskit/qiskit/blob/stable/2.0/releasenotes/notes/2.0/remove-pulse-calibrations-4486dc101b76ec51.yaml)
- Qiskit Aer: latest release is 0.17.2 (2025-09-17), with nothing since. The README carries a banner: "**Reduced Maintenance Mode** … Aer is currently operating under reduced maintenance. The maintainers are only able to address critical bug fixes. Feature requests and non-critical enhancements have been placed in the backlog and will be revisited as capacity allows." — [Aer releases](https://github.com/Qiskit/qiskit-aer/releases); [Aer README](https://github.com/Qiskit/qiskit-aer/blob/main/README.md)
- IBM's current Aer guide pins `qiskit[all]~=2.5.2` and `qiskit-aer~=0.17` — [IBM guide: Exact and noisy simulation with Qiskit Aer](https://quantum.cloud.ibm.com/docs/en/guides/simulate-with-qiskit-aer)
- Aer removed its own PulseSimulator in Aer 0.14 — [Aer 0.14 release note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.14/remove_pulse_simulator-f8de2f6d380f446a.yaml)
- Qiskit Dynamics: the GitHub repo is archived, with description "DEPRECATED - Tools for building and solving models of quantum systems in Qiskit". The README says "This repo is no longer being actively maintained." Release 0.6.0 (2025-10-07) "is the last release of the package for the foreseeable future" — [qiskit-dynamics repo](https://github.com/qiskit-community/qiskit-dynamics); [0.6.0 prelude](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/0.6.0-f7b07109b632757b.yaml)
- Dynamics 0.6 "requires `qiskit <= 1.3`, as a result of the removal of `qiskit.pulse` in Qiskit 2.0". Workflows that don't use `qiskit.pulse` "should still be compatible". `DynamicsBackend` "will now only work with Qiskit Experiments `0.8` releases" (`<0.9`), and `qiskit_ibm_runtime` is capped at 0.36.1 — [bound-qiskit-version note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/bound-qiskit-version-93a09ea5a3c4afbb.yaml); [runtime-bound note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/fix-qiskit-ibm-runtime-d3ff79ecceda0a65.yaml)
- qiskit-experiments: latest release is 0.14.2 (2026-08-25). PR #1476 "Deprecate pulse and restless related experiments and classes" merged 2024-10-25. PR #1511 "Remove pulse and calibration management related code" merged 2025-03-02, and 0.9.0 was released 2025-03-18 — [QE releases](https://github.com/qiskit-community/qiskit-experiments/releases); [PR 1476](https://github.com/qiskit-community/qiskit-experiments/pull/1476); [PR 1511](https://github.com/qiskit-community/qiskit-experiments/pull/1511)
- qiskit-experiments 0.14.0 highlights: new `PurityRB`, new `TiledExperiment` plus `partition_qubits` "for running experiments across a device", and removal of the `qiskit_ibm_experiment` dependency "for uploading to the discontinued IBM Experiment Service" (replaced by a local experiment service) — [0.14.0 prelude](https://github.com/qiskit-community/qiskit-experiments/blob/main/releasenotes/notes/prepare-0.14.0-faa0f9b0f22b4052.yaml)
- Qiskit algorithm libraries:
  - qiskit-algorithms README: "**Qiskit Algorithms is no longer officially supported by IBM**"; latest is 0.4.0 (2025-08-29).
  - qiskit-optimization is archived; last release 0.7.0 (2025-08-20).
  - qiskit-nature is still active; 0.8 released 2026-06-01.
  - Sources: [qiskit-algorithms](https://github.com/qiskit-community/qiskit-algorithms); [qiskit-optimization](https://github.com/qiskit-community/qiskit-optimization); [qiskit-nature releases](https://github.com/qiskit-community/qiskit-nature/releases)
- IBM-maintained numerical addons are active:
  - qiskit-addon-mpf 0.3.0, "Reducing the Trotter error of Hamiltonian dynamics with multi-product formulas".
  - qiskit-addon-aqc-tensor 0.3.1 (2026-07-16), "approximate quantum compilation with tensor networks".
  - qiskit-addon-sqd 0.13.1 (2026-08-12), "Sample-based Quantum Diagonalization".
  - Sources: [MPF](https://github.com/Qiskit/qiskit-addon-mpf); [AQC-Tensor](https://github.com/Qiskit/qiskit-addon-aqc-tensor); [SQD](https://github.com/Qiskit/qiskit-addon-sqd)
- Cirq: latest is v1.7.0 (2026-06-30); 1.5.0 was 2025-04-10 and 1.6.0 was 2025-07-23. qsim: v0.22.1 (2026-09-01) — [Cirq releases](https://github.com/quantumlib/Cirq/releases); [qsim releases](https://github.com/quantumlib/qsim/releases)
- Google and partner stack:
  - ReCirq's last tagged release is "v2020-10", though the repo was still pushed on 2026-10-02.
  - OpenFermion 1.8.1 (2026-07-17).
  - TensorFlow Quantum 0.7.6 (2026-02-25).
  - Mitiq 1.1.0 (2026-09-04).
  - Sources: [ReCirq](https://github.com/quantumlib/ReCirq/releases); [OpenFermion](https://github.com/quantumlib/OpenFermion/releases); [TFQ](https://github.com/tensorflow/quantum/releases); [Mitiq](https://github.com/unitaryfund/mitiq/releases)
- PennyLane: latest stable is 0.45.1 (2026-06-26), with 0.46.0b1 as a pre-release (2026-09-09). Lightning: 0.45.0 (2026-05-11) — [PennyLane releases](https://github.com/PennyLaneAI/pennylane/releases); [Lightning releases](https://github.com/PennyLaneAI/pennylane-lightning/releases)
- PennyLane 0.45: "The `qml` alias as in `import pennylane as qml` has been updated to `qp` in our source code and documentation." All 53 demo sources I sampled now use `import pennylane as qp` — [changelog-0.45.0](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-0.45.0.md); [demo sources](https://github.com/PennyLaneAI/qml/tree/master/demonstrations_v2)
- PennyLane development changelog, "Release 0.46.0 (development release)", breaking changes:
  - TensorFlow interfaces removed.
  - "`pennylane.pulse` module has been removed. This includes `ParametrizedHamiltonian` and `ParametrizedEvolution`, as well as the `stoch_pulse_grad` and `pulse_odegen` pulse-level gradient transforms."
  - "`pennylane.qaoa` module has been removed" (mixers, `maxcut` and other cost Hamiltonians, `cost_layer`/`mixer_layer`).
  - "`pennylane.noise` module has been removed, including `NoiseModel`, `add_noise`, `insert`, noise mitigation transforms (`mitigate_with_zne`, `fold_global`, `poly_extrapolate`, `richardson_extrapolate`, `exponential_extrapolate`), and `from_qiskit_noise`. Noise channels such as `AmplitudeDamping` are unaffected."
  - "All functionality related to qutrits/qudits has been removed. Qudit functionality in `pennylane.labs` still remains."
  - All continuous-variable (CV) code removed, including `DefaultGaussian`.
  - Python ≥3.12 required.
  - Source: [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md)
- PennyLane removal PRs:
  - PR #10238 "Remove the `qp.pulse` module" merged 2026-10-01. Its body says "The `qp.pulse` module is slated for removal in PL2."
  - PR #9867 "Remove everything related to qudits" merged 2026-07-24 ("All ops, devices, tests").
  - `default.qutrit`, `default.qutrit.mixed` and `default.gaussian` (source files `default_qutrit.py`, `default_qutrit_mixed.py`, `default_gaussian.py`) are present at tag v0.45.1 and absent on `main`.
  - Sources: [PR 10238](https://github.com/PennyLaneAI/pennylane/pull/10238); [PR 9867](https://github.com/PennyLaneAI/pennylane/pull/9867); [v0.45.1 devices](https://github.com/PennyLaneAI/pennylane/tree/v0.45.1/pennylane/devices); [main devices](https://github.com/PennyLaneAI/pennylane/tree/main/pennylane/devices)

#### Qiskit Aer (circuit-level simulation and noise)
- AerSimulator methods:
  - `automatic` (the default).
  - `statevector`: "can sample measurement outcomes from ideal circuits with all measurements at end of the circuit".
  - `density_matrix`: "may sample measurement outcomes from noisy circuits with all measurements at end of the circuit".
  - `stabilizer`: "can simulate noisy Clifford circuits if all errors in the noise model are also Clifford errors".
  - `extended_stabilizer`: approximate Clifford+T.
  - `matrix_product_state`.
  - `unitary`: "does not support measurement, reset, or noise".
  - `superop`: "does not support measurement".
  - `tensor_network`: "only available for GPU", accelerated by cuTensorNet.
  - Source: [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html)
- Official Aer tutorials: AerSimulator basics, device noise simulation, building noise models, custom gate noise, noise transformation, extended stabilizer, MPS method. There are how-tos for parallel and GPU execution — [Aer tutorials](https://qiskit.github.io/qiskit-aer/tutorials/index.html); [tutorial sources](https://github.com/Qiskit/qiskit-aer/tree/main/docs/tutorials)
- Noise-model toolkit: `QuantumError` describes "CPTP gate errors" applied after gate or reset instructions, and `ReadoutError` covers classical readout errors. Helpers include `kraus_error`, `mixed_unitary_error`, `coherent_unitary_error`, `pauli_error`, `depolarizing_error` and thermal relaxation — [Building noise models](https://qiskit.github.io/qiskit-aer/tutorials/3_building_noise_models.html)
- `NoiseModel.from_backend` composition:
  - Single-qubit gate errors: "a single qubit depolarizing error followed by a single qubit thermal relaxation error".
  - Two-qubit gate errors: "a two-qubit depolarizing error followed by single-qubit thermal relaxation errors on both qubits".
  - Readout errors per qubit.
  - Parameters are derived from T1, T2 and gate times.
  - The tutorial warns these "automatic models are only an *approximation* of the real errors that occur on actual devices, due to the fact that they must be build from a limited set of input parameters related to *average error rates* on gates".
  - The tutorial still uses the long-retired `ibmq_vigo` (outdated).
  - Source: [Device noise simulation](https://qiskit.github.io/qiskit-aer/tutorials/2_device_noise_simulation.html)
- Aer 0.10 added "simulation of circuits containing QASM 3.0 control-flow instructions": ForLoopOp, WhileLoopOp, ContinueLoopOp, BreakLoopOp, IfElseOp. It also added "relaxation noise on scheduled circuits in backend noise models", i.e. idle-time decoherence. Aer 0.13 added SwitchCaseOp — [Aer 0.10 note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.10/0-10-release-8c37dadcc1c82fcc.yaml); [Aer 0.13 note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.13/support_switch-41603d87cb8358fb.yaml)
- IBM's 2026 tutorial catalogue is built around hardware- and utility-scale workflows, with classical simulation in a supporting role. Contents by area:
  - Chemistry and diagonalization: SQD for chemistry and for a pooled nuclear Hamiltonian; sample-based Krylov; SqDRIFT; Krylov diagonalization of lattice Hamiltonians.
  - Optimization: QAOA, advanced QAOA, warm-start QAOA, multi-objective QAOA, Pauli-correlation encoding for MaxCut.
  - Many-body physics: neutron scattering with "AQC + Trotter dynamics"; "Ground-state energy estimation of the Heisenberg chain with VQE"; "Simulate time evolution of the transverse-field Ising model"; "Nishimori phase transition".
  - Dynamic circuits: long-range entanglement, "Simulation of kicked Ising Hamiltonian with dynamic circuits", cut-Bell-pair benchmarking.
  - Addons: MPF to reduce Trotter error, AQC-Tensor, operator backpropagation, circuit cutting, M3 readout mitigation.
  - Error mitigation: PEA, ZNE options in Estimator, PEC with shaded lightcones, propagated noise absorption, PEC with logical noise models.
  - Partner "Qiskit Functions", e.g. "Dissociation PES curves with Qunova HiVQE".
  - Source: [IBM Quantum tutorials index](https://quantum.cloud.ibm.com/docs/en/tutorials)

#### Qiskit Dynamics (ODE and pulse level; archived)
- JOSS paper by Puzzuoli, Wood, Egger, Rosand and Ueda (2023). Features: time-dependent Hamiltonians, Lindblad dynamics, rotating frames, RWA, multiple solvers, JAX autodiff/JIT/GPU, Dyson and Magnus perturbative solvers, and `DynamicsBackend` for pulse-level simulation — [JOSS 10.21105/joss.05853](https://joss.theoj.org/papers/10.21105/joss.05853)
- The README describes NumPy (default) or JAX array backends. JAX is "well-suited to tasks involving repeated evaluation of functions with different parameters; E.g. simulating a model of a quantum system over a range of parameter values, or optimizing the parameters of control sequence" — [README](https://github.com/qiskit-community/qiskit-dynamics)
- Tutorials and user guides (sources: [docs/tutorials](https://github.com/qiskit-community/qiskit-dynamics/tree/main/docs/tutorials), [docs/userguide](https://github.com/qiskit-community/qiskit-dynamics/tree/main/docs/userguide), rendered at [tutorials index](https://qiskit-community.github.io/qiskit-dynamics/tutorials/index.html)):
  - *Rabi oscillations*: a qubit, then a Lindblad equation with amplitude damping (Γ1) and pure dephasing (Γ2).
  - *Solving the Lindblad dynamics of a qubit chain*: steady state as a function of model parameters.
  - *Simulating backends at the pulse level with DynamicsBackend*: a 2-qubit fixed-frequency transmon model with `dim = 3` per transmon and explicit anharmonicities. It runs "rough calibrations for X and SX gates" through qiskit-experiments and characterizes the 2-qubit interaction with `CrossResonanceHamiltonian`.
  - *Gradient optimization of a pulse sequence*: an X gate with JAX gradients of the infidelity and a Gaussian-square parameterization.
  - *Simulating Qiskit Pulse schedules*.
  - *Systems modelling*: two transmons as Duffing oscillators (`dim=3`) with an exchange interaction.
  - User guides: rotating frame, RWA and sparse arrays for a transmon (`dim = 5`); Dyson/Magnus perturbative solvers on an oscillator (`dim = 10`); restricting operators to the low-energy subspace of a 3-transmon model.
- Several tutorials include code to "silence deprecation warnings from pulse" — [optimizing_pulse_sequence.rst](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/optimizing_pulse_sequence.rst); [dynamics_backend.rst](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/dynamics_backend.rst)
- `DynamicsBackend` returns multi-level (n-ary) outcomes, encoded as hex since 0.6, with `max_outcome_level` setting the number of levels. This supports leakage-level readout — [backend-hex-results note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/backend-hex-results-fd3762b9188cdd01.yaml)
- 0.6.0 added a `systems` module "for building abstract models of quantum systems … with the goal of minimizing the need for a user to explicitly build and work with arrays (an error-prone process)" — [systems-module note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/systems-module-227a233674a6100e.yaml)

#### qiskit-experiments (characterization and verification)
- Current experiment library, from the `main` docstring. No calibration section remains:
  - Verification: `StandardRB`, `InterleavedRB`, `PurityRB`, `LayerFidelity`, `LayerFidelityUnitary`, `StateTomography`, `ProcessTomography`, `MitigatedStateTomography`, `MitigatedProcessTomography`, `QuantumVolume`.
  - Single-qubit characterization: `T1`, `T2Hahn`, `T2Ramsey`, `Tphi`, `HalfAngle`, `FineAmplitude`/`FineXAmplitude`/`FineSXAmplitude`, `RamseyXY`, `FineFrequency`, `ReadoutAngle`, `FineDrag`/`FineXDrag`/`FineSXDrag`, `MultiStateDiscrimination`.
  - Two-qubit characterization: `ZZRamsey`.
  - Readout mitigation: `LocalReadoutError`, `CorrelatedReadoutError`.
  - Source: [library `__init__`](https://github.com/qiskit-community/qiskit-experiments/blob/main/qiskit_experiments/library/__init__.py)
- Manuals cover T1, T2 Hahn, T2* Ramsey, Tphi, readout mitigation, quantum volume, randomized benchmarking and state tomography — [manual sources](https://github.com/qiskit-community/qiskit-experiments/tree/main/docs/manuals); [rendered manuals](https://qiskit-community.github.io/qiskit-experiments/manuals/index.html)
- The T1 and Tphi manuals simulate with `AerSimulator.from_backend(FakePerth(), noise_model=NoiseModel.from_backend(FakePerth(), thermal_relaxation=True, gate_error=False, readout_error=False))`. The T1 manual also shows `MockIQBackend` producing kerneled (level-1) IQ data — [T1 manual](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/t1.html); [Tphi manual](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/tphi.html)
- The RB manual uses Aer with `depolarizing_error(5e-3, 1)` on `sx`/`x`, 0 on `rz` and `depolarizing_error(5e-2, 2)` on `cx`. It reports EPC, α (fit a·α^m+b) and EPG. Two-qubit EPC is corrected using EPGs from 1-qubit RB, and interleaved RB reports a systematic error — [RB manual](https://qiskit-community.github.io/qiskit-experiments/manuals/verification/randomized_benchmarking.html)

#### Cirq + qsim (+ ReCirq)
- Built-in Cirq simulators:
  - `cirq.Simulator`: pure state, handles "noisy evolutions that preserve state purity".
  - `cirq.DensityMatrixSimulator`: "the density matrix can represent all possible results of a noisy circuit".
  - External: qsim ("For most users we recommend qsim"), qsimh, quimb, and qFlex ("Not recommended - prefer qsim or quimb").
  - Source: [Cirq simulation](https://quantumai.google/cirq/simulate/simulation)
- cirq-core `sim` contains `clifford/`, `classical_simulator`, density-matrix, state-vector and product-state simulators. `contrib/quimb` adds `mps_simulator`, `density_matrix` and grid-circuit tensor-network contraction — [cirq-core/cirq/sim](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/sim); [contrib/quimb](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/contrib/quimb)
- Noise: Kraus channels (bit flip, amplitude damping, generalized amplitude damping, phase flip, phase damping, depolarizing, asymmetric depolarizing), the `cirq.kraus` protocol, and `cirq.NoiseModel` for whole circuits — [Representing noise](https://quantumai.google/cirq/noise/representing_noise); [Noisy simulation](https://quantumai.google/cirq/simulate/noisy_simulation)
- Quantum Virtual Machine: noisy qsim simulation driven by published device noise data for Willow, Weber and Rainbow. "In internal tests, the virtual and actual hardware are within experimental error of each other." For 25–32 qubits it recommends multinode parallelization — [Quantum Virtual Machine](https://quantumai.google/cirq/simulate/quantum_virtual_machine)
- Characterization (QCVV) docs: XEB theory, isolated XEB, parallel XEB, coherent vs incoherent XEB, error heatmaps; Google tutorials on echoes, spin echoes and identifying hardware changes — [Cirq docs tree](https://github.com/quantumlib/Cirq/tree/main/docs); [XEB theory](https://quantumai.google/cirq/noise/qcvv/xeb_theory)
- `cirq.experiments` helpers: `t1_decay`, `t2_decay`, `single_qubit_randomized_benchmarking`, `parallel_single_qubit_randomized_benchmarking`, `two_qubit_randomized_benchmarking`, `single_qubit_state_tomography`, `two_qubit_state_tomography`, n-qubit tomography, purity estimation, readout confusion matrix, XEB sampling/fitting/simulation, z-phase calibration — [cirq-core/cirq/experiments](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/experiments)
- Experiments gallery:
  - Base Cirq: textbook algorithms, Shor, VQE, quantum walks, Fourier checking, hidden linear function. The VQE example is a "2D +/- Ising model with transverse field", not a molecule.
  - ReCirq: QAOA (MaxCut, Ising, landscape and optimization analysis, hardware-grid circuits, tket routing, binary paintshop), Hartree–Fock VQE, Fermi–Hubbard spin-charge separation, OTOC scrambling, toric code, Rabi oscillations, KPZ with the Heisenberg chain, QC-QMC, lattice gauge theory, random circuit sampling, quantum chess.
  - Sources: [Cirq experiments](https://quantumai.google/cirq/experiments); [Cirq VQE tutorial](https://quantumai.google/cirq/experiments/variational_algorithm)
- qsim offers a "Schrödinger state-vector simulator designed to run on a single machine" and a new "truncated Matrix Product State simulator" (C++ only). It uses "AVX/FMA vector operations, OpenMP multithreading, and gate fusion" and cites the 2019 quantum-supremacy paper — [qsim overview](https://quantumai.google/qsim/overview); [qsim MPS docs](https://github.com/quantumlib/qsim/blob/main/docs/mps.md)
- qsim simulates noise only as quantum trajectories. Its noise tutorial says "a density matrix requires O(4^N) storage for N qubits. In qsimcirq, noisy circuits are instead simulated as 'trajectories'". A GPU trajectory binary, `qsim_qtrajectory_cuda`, "performs quantum trajactory simulations with amplitude damping and phase damping noise channels" — [Noisy qsimcirq](https://quantumai.google/qsim/tutorials/noisy_qsimcirq); [qsim usage](https://github.com/quantumlib/qsim/blob/main/docs/usage.md)
- Isakov et al. (2021) present a "delayed inner product algorithm for quantum trajectories which can result in an order of magnitude speedup for low noise". They also describe multinode trajectory simulation and compare noisy simulations with Google experiments — [arXiv:2111.02396](https://arxiv.org/abs/2111.02396)

#### PennyLane
- The demo source directory holds 216 entries — [qml/demonstrations_v2](https://github.com/PennyLaneAI/qml/tree/master/demonstrations_v2)
- Demos relevant to numerics, with first-publication date (they are actively maintained: most were last modified April–August 2026). Sources: [demo metadata/sources](https://github.com/PennyLaneAI/qml/tree/master/demonstrations_v2), each at `https://pennylane.ai/qml/demos/<name>`.
  - **VQE and chemistry:** *A brief overview of VQE* (2020-02); *Accelerating VQEs with quantum natural gradient* (2020-11); *VQE in different spin sectors* (2020-10); *How to implement VQD* (2024-08); *Adaptive circuits for quantum chemistry* (2021-09); *Givens rotations* (2021-06); *Building molecular Hamiltonians* (2020-08); *Differentiable Hartree-Fock* (2022-05); *Optimization of molecular geometries* (2021-06); *Modelling chemical reactions* (2021-07); *A quantum algorithm for vibronic dynamics* (2026-08).
  - **QAOA:** *Intro to QAOA* (2020-11); *QAOA for MaxCut* (2019-10); *Quantum natural SPSA optimizer* (2022-07).
  - **Gradients and optimization:** *Quantum natural gradient* (2019-10); *Adjoint Differentiation* plus supplementary benchmarking (2021-11); *Quantum gradients with backpropagation* (2020-08); *Generalized parameter-shift rules* (2021-08); *Optimization using SPSA* (2023-03); *Implicit differentiation of variational quantum algorithms* (2022-11).
  - **Barren plateaus:** *Barren plateaus in quantum neural networks* (2019-10); *Alleviating barren plateaus with local cost functions* (2020-09).
  - **Noise and mitigation:** *Noisy circuits* (2021-02); *How to use noise models* (2024-10); *How to import noise models from Qiskit* (2024-11); *Error mitigation with Mitiq and PennyLane* (2021-11); *Differentiating quantum error mitigation transforms* (2022-08); *Digital ZNE with Catalyst* (2024-11); *Is quantum computing useful before fault tolerance?* (2023-06); *Optimizing noisy circuits with Cirq* (2020-06).
  - **Pulse and hardware physics:** *Differentiable pulse programming with qubits* (2023-03); *Optimal control for gate compilation* (2023-08); *Gate calibration with reinforcement learning* (2024-04); *Pulse programming on OQC Lucy* (2023-10); *Quantum computing with superconducting qubits* (2022-03).
  - **Dynamics and many-body:** *Exploring Trotterization* (2026-07-28); *Seeing quantum phase transitions with quantum computers* (2025-10); *The Quantum Graph Recurrent Neural Network* (2020-07, Ising dynamics); *g-sim: Lie-algebraic classical simulations* (2024-06).
  - **Benchmarks and simulators:** *Quantum volume* (2020-12); *Beyond classical computing with qsim* (2020-11); *Introducing matrix product states* (2024-09); *How to simulate quantum circuits with tensor networks using DefaultTensor* (2024-07); *Efficient Simulation of Clifford Circuits* (2024-04); *Classical shadows* (2021-06).
  - **Mid-circuit measurement and qudits:** *Introduction to mid-circuit measurements* (2024-05); *How to create dynamic circuits with mid-circuit measurements* (2024-05); *Qutrits and quantum algorithms* (2023-05).
- Simulator devices used in the sampled demos:
  - Native: `lightning.qubit`, `default.qubit`, `default.mixed`, `default.tensor`, `default.clifford`, `default.qutrit`, `default.gaussian`.
  - Plugin devices: `cirq.qsim`, `cirq.mixedsimulator`, `qiskit.aer`, `qulacs.simulator`, `qrack.simulator`.
  - Source: [demo sources](https://github.com/PennyLaneAI/qml/tree/master/demonstrations_v2)
- Demos that still call modules removed on `main`, found by grep of the sources:
  - Pulse/`evolve`: `oqc_pulse`, `tutorial_pulse_programming101`, `tutorial_optimal_control`, `tutorial_rl_pulse`.
  - Noise/ZNE transforms: `tutorial_diffable-mitigation`, `tutorial_error_mitigation`, `tutorial_how_to_use_noise_models`, `tutorial_how_to_import_qiskit_noise_models`, `tutorial_mitigation_advantage`.
  - `qaoa` module: `tutorial_qaoa_intro`, `qnspsa`.
  - `default.qutrit`: `tutorial_qutrits_bernstein_vazirani`.
  - `default.gaussian`: `tutorial_sc_qubits`.
  - Source: [demo sources](https://github.com/PennyLaneAI/qml/tree/master/demonstrations_v2); compare with [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md)
- Reference paper: Bergholm et al., "PennyLane: Automatic differentiation of hybrid quantum-classical computations" — [arXiv:1811.04968](https://arxiv.org/abs/1811.04968) (not re-read in this session).

### Inferences
- **Coverage matrix.** My synthesis of the findings above; ✔ = featured in official material, ~ = partial or indirect, ✘ = absent, ⚠ = present but deprecated/removed/unmaintained.

| Problem class | Qiskit (Aer / Dynamics / Experiments / IBM tutorials) | Cirq + qsim (+ReCirq) | PennyLane |
|---|---|---|---|
| Noisy circuits, noise models | ✔ Aer Kraus/Pauli/thermal/readout, `from_backend` (⚠ Aer reduced maintenance) | ✔ channels, NoiseModel, QVM with device data; qsim trajectories | ✔ `default.mixed` channels (≤23 wires); ⚠ `NoiseModel`/`add_noise` removed on `main` |
| Pulse-level transmon, leakage, DRAG | ⚠ Dynamics only (archived; pinned to Qiskit ≤1.3); Pulse removed in 2.0; FineDrag survives as gate-level characterization | ✘ | ⚠ `qp.pulse` (transmon drive/interaction, JAX ODE) removed on `main` |
| T1/T2/Ramsey/Hahn, RB, tomography, QV | ✔ qiskit-experiments (strongest; active) | ✔ `cirq.experiments` (T1/T2, RB, tomography, XEB) | ~ QV demo; channels only |
| VQE molecules (H2, LiH) | ~ shifted to SQD/Krylov; qiskit-nature active; qiskit-algorithms unsupported | ~ OpenFermion/ReCirq HF-VQE; core VQE is an Ising toy | ✔ broadest (datasets, differentiable HF, geometry opt.) |
| QAOA MaxCut | ✔ hardware-scale tutorials; qiskit-optimization archived | ✔ ReCirq QAOA | ✔ demos, ⚠ `qaoa` module removed on `main` |
| Gradients, QNG | ~ via unsupported qiskit-algorithms / qiskit-machine-learning | ✘ in core (via TFQ or PennyLane plugin) | ✔ core strength |
| Barren plateaus | ✘ (not featured) | ~ (historically TFQ) | ✔ |
| ZNE / error mitigation | ✔ PEA/ZNE/PEC/PNA/M3 on hardware via runtime | ~ via Mitiq (external) | ✔ demos, ⚠ ZNE transforms removed on `main` |
| Trotterized many-body dynamics | ✔ TFIM, kicked Ising, MPF, AQC-Tensor | ✔ Fermi–Hubbard, KPZ/Heisenberg, OTOC | ✔ *Exploring Trotterization* (2026) |

- **Judgement on fit.**
  - Qiskit serves characterization (qiskit-experiments) and Markovian noisy-circuit simulation (Aer) best. Aer's frozen state and IBM's pivot to hardware-scale workflows mean simulator numerics are no longer IBM's teaching focus.
  - Cirq is strongest on characterization (XEB), hardware-calibrated noise (QVM), qudits and classical control. It lacks autodiff and pulse physics in core.
  - PennyLane has the deepest pedagogy for variational and differentiable numerics, but its own main branch is now stripping the physics-flavoured numerics: pulses, noise models/ZNE, qutrits and CV.
- **Direction of travel.** All three are moving toward utility/FTQC workflows: SQD/Krylov/addons for IBM; Catalyst, resource estimation, QSVT and QRAM demos for PennyLane. Open-system, time-dependent and multi-level simulation is being orphaned. This is the clearest gap in the incumbents' coverage.

### Gaps
- I did not open each IBM tutorial to check whether it runs on Aer or only on hardware. The IBM index lists titles, not simulator usage.
- Qiskit Nature tutorial content (whether H2/LiH dissociation tutorials still exist in 0.8) was not checked.
- I did not find an official statement of PennyLane's "PL2" scope or timeline beyond the PR text, nor whether 0.46.0 final will ship all `main`-branch removals.
- TensorFlow Quantum's barren-plateau tutorial and gradient support for Cirq were not re-verified, so their current status is unknown.

---

## 2. What scale and performance are published (statevector limits, MPS reach, GPU acceleration, gradient cost)?

### Takeaway
Published statevector scale is memory-bound and essentially the same everywhere: about 29–32 qubits per node or GPU, 36 qubits on 8×A100-80GB, and 41 qubits multi-node for PennyLane Lightning. MPS and stabilizer methods go to 50–1000+ qubits only for low-entanglement or Clifford(+few-T) circuits. Density-matrix and noisy simulation is much smaller. Documented performance claims:

| Claim | Source |
|---|---|
| Noisy simulation "fewer than 18 qubits" on an 8 GB laptop | qsim |
| Hard cap of 23 wires on `default.mixed` | PennyLane |
| Parameter-shift cost grows linearly with the number of parameters | PennyLane docs |
| Adjoint method: "significantly lower memory usage and a similar runtime" vs backprop | PennyLane docs |
| Dyson/Magnus solvers: about 2–4× (solution) and 10–60× (gradient) faster than standard ODE solvers on GPU for transmon gates | Qiskit Dynamics paper |

No official cross-SDK simulator benchmark exists.

### Cited Findings
- **qsim hardware guidance**:
  - On a laptop with ≥8 GB: "Noiseless simulations that use fewer than 29 qubits" and "Noisy simulations that use fewer than 18 qubits"; "The qubit upper bounds in this chart are not technical limits".
  - Memory rule of thumb: "8 · 2^N bytes for an N-qubit circuit".
  - "GPU hardware starts to outperform CPU hardware significantly (up to 15x faster) for circuits with more than 20 qubits."
  - "On an NVIDIA A100 GPU with 40GB of memory, the maximum number of qubits is 32." With "eight 80-GB NVIDIA A100 GPUs (640GB of total GPU memory), you can simulate up to 36 qubits", which needs cuQuantum Appliance for multi-GPU. The cuStateVec backend "does not enable multi-GPU support".
  - Benchmarks compared an A100 with a c2-standard-4 CPU, including a noisy case with a "phase damping channel (p=0.01)".
  - Gate-fusion advice: keep fusion size low "up to 22 qubits", higher above that.
  - Source: [qsim: Choosing hardware](https://quantumai.google/qsim/choose_hw)
- **Conflicting qsim memory figures:**
  - qsim overview: "Circuits of up to 30 qubits can be simulated in qsim with ~16GB of RAM".
  - Cirq simulation page: qsim is "Recommended for deep circuits with up to 30 qubits (consumes 8GB RAM)".
  - Sources: [qsim overview](https://quantumai.google/qsim/overview); [Cirq simulation](https://quantumai.google/cirq/simulate/simulation)
- **qsim other platforms:** AMD Instinct GPUs are supported via HIP/ROCm (build from source). There is a 32-qubit, depth-14 circuit tutorial and an HTCondor multinode tutorial — [AMD GPU](https://github.com/quantumlib/qsim/blob/main/docs/tutorials/amd_gpu.md); [q32d14 tutorial](https://quantumai.google/qsim/tutorials/q32d14); [multinode](https://quantumai.google/qsim/tutorials/multinode)
- **Aer GPU:**
  - Requires CUDA ≥ 11.2. The GPU package "is only available on x86_64 Linux". GPU methods are statevector, density matrix and unitary, plus `tensor_network` per the API page.
  - Multi-GPU and multi-node runs use cache-blocking "chunks": e.g. "30-qubits circuit is distributed into 2^10 chunks with 20-qubits", subject to `sizeof(complex)*2^(blocking_qubits+4) < size of the smallest memory space`.
  - Sources: [Aer README](https://github.com/Qiskit/qiskit-aer/blob/main/README.md); [Running with multiple GPUs/nodes](https://qiskit.github.io/qiskit-aer/howtos/running_gpu.html)
- **Aer options:** `precision` is 'single' or 'double' (default double). `max_memory_mb` defaults to system memory. `batched_shots_gpu` batching applies only up to `batched_shots_gpu_max_qubits` (default 16). `runtime_parameter_bind_enable` binds parameters at runtime. MPS truncation threshold defaults to 1e-16, with an optional max bond dimension — [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html)
- **Aer MPS:** a 50-qubit EPR state runs in 0.31 s, and the tutorial mentions 500- and 1000-qubit runs as feasible but slower. "In the worst case, the tensors may grow exponentially". Converting to a state vector loses the benefit — [Aer MPS tutorial](https://qiskit.github.io/qiskit-aer/tutorials/7_matrix_product_state_method.html)
- **Aer extended stabilizer:** "able to handle circuits of up to 63 qubits". The example is a 40-qubit, 60-gate random circuit. Default approximation error is 0.05, and Markov-chain mixing is repeated per shot — [Extended stabilizer tutorial](https://qiskit.github.io/qiskit-aer/tutorials/6_extended_stabilizer_tutorial.html)
- **PennyLane Lightning HPC paper** (Asadi et al., 2024): "up to 30 qubits on a single device or node, and up to 41 qubits using multiple nodes". It covers CPUs (SIMD, multithreading) and NVIDIA and AMD GPUs, batched multi-GPU execution, and "distributed forward and gradient-based quantum circuit executions across multiple nodes". Workloads were QAOA, VQE and synthetic circuits — [arXiv:2403.02512](https://arxiv.org/abs/2403.02512)
- **Lightning changelog:** MPI-distributed simulation added to `lightning.kokkos`. Mid-circuit measurements via the tree-traversal algorithm are supported. `lightning.tensor` uses cuTensorNet (worksize tunable). There is a `lightning.amdgpu` build, and CUDA 13 builds for `lightning.gpu`/`lightning.tensor` — [Lightning CHANGELOG](https://github.com/PennyLaneAI/pennylane-lightning/blob/main/.github/CHANGELOG.md)
- **PennyLane `default.mixed` cap:** "This device does not currently support computations on more than 23 wires" — [default_mixed.py @ v0.45.1](https://github.com/PennyLaneAI/pennylane/blob/v0.45.1/pennylane/devices/default_mixed.py)
- **PennyLane tensor-network demos:**
  - The DefaultTensor how-to sweeps `num_qubits` from 50 to 200 (step 50) with `max_bond_dim` 50, cutoff = complex128 machine epsilon and `contract="auto-mps"`.
  - The MPS demo simulates the H6 (STO-3G, bond length 1.3) VQE circuit with `default.tensor` at bond dimension 30. It then extrapolates in bond dimension ("finite-size scaling … run the same simulation with an increasing bond dimension and check that it saturates").
  - Sources: [DefaultTensor how-to](https://pennylane.ai/qml/demos/tutorial_How_to_simulate_quantum_circuits_with_tensor_networks); [MPS demo](https://pennylane.ai/qml/demos/tutorial_mps)
- **PennyLane `default.clifford`** uses stim and targets "circuits scaling up to more than thousands of qubits" — [Clifford demo](https://pennylane.ai/qml/demos/tutorial_clifford_circuit_simulations)
- **PennyLane gradient cost:**
  - parameter-shift: "Number of circuit executions required scales linearly with the number of trainable circuit parameters".
  - adjoint: "significantly lower memory usage and a similar runtime" compared with backprop.
  - Source: [PennyLane gradients and training](https://docs.pennylane.ai/en/stable/introduction/interfaces.html)
- **PennyLane adjoint demos:** the adjoint demo cites Jones & Gacon (arXiv:2009.02823) and concludes adjoint "scales nicely in time without excessive memory requirements", whereas backprop "requires more and more memory as the circuit gets larger". The supplementary benchmark compares adjoint (Lightning) with parameter-shift and backprop (`default.qubit`), sweeping wires at 6 layers and layers at 12 wires — [Adjoint Differentiation](https://pennylane.ai/qml/demos/tutorial_adjoint_diff); [benchmarking supplement](https://pennylane.ai/qml/demos/adjoint_diff_benchmarking)
- **PennyLane SPSA demo:** compares executions per optimization step: 2 for SPSA vs 2 × (number of parameters) for parameter-shift (each multiplied by 15 Hamiltonian terms for the H2 cost) — [SPSA demo](https://pennylane.ai/qml/demos/tutorial_spsa)
- **Qiskit Dynamics perturbative solvers** (Puzzuoli et al., J. Comput. Phys. 2023): "Dyson and Magnus-based solvers provide a speedup over traditional ODE solvers, ranging from roughly 2× to 4× for a solution and 10× to 60× for a gradient", benchmarked on "a two-transmon entangling gate" on GPU, with accuracy shown on "a single transmon". The user guide says that "once compiled, the Dyson-based solver has a significant speed advantage" — [arXiv:2210.11595](https://arxiv.org/abs/2210.11595); [Dyson/Magnus user guide source](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/userguide/perturbative_solvers.rst)
- **PennyLane pulse ODE (now removed on `main`):** "the step sizes between t0 and t1 are chosen adaptively to stay within a given error tolerance. We can backpropagate through this ODE solver and obtain the gradient via `jax.grad`". The optimal-control demo optimizes 15P = 45 parameters for P = 3 and reports wall-clock CPU time — [Pulse programming demo](https://pennylane.ai/qml/demos/tutorial_pulse_programming101); [Optimal control demo](https://pennylane.ai/qml/demos/tutorial_optimal_control)
- **qsim trajectories:** "an order of magnitude speedup for low noise" with the delayed-inner-product algorithm — [arXiv:2111.02396](https://arxiv.org/abs/2111.02396)

### Inferences
- **Precision arithmetic:**
  - At double precision (Aer default, complex128 = 16 B per amplitude), 30 qubits need 16 GiB and 34 qubits need 256 GiB.
  - qsim's 8 B per amplitude means single precision (complex64). This fits the A100-40GB → 32 qubits figure (32 GiB) and is consistent with the `MPSStateSpace<For, float>` template in its MPS docs. qsim's docs don't state the precision explicitly.
- **Density-matrix arithmetic:** an n-qubit density matrix costs as much as a 2n-qubit statevector, so about 15 qubits fit in 16 GiB at double precision. PennyLane's 23-wire cap (≈1 PB at complex128) is therefore never the binding limit in practice. qsim's laptop guidance of "noisy < 18 qubits" refers to trajectory sampling, not density matrices.
- **What the bar is and isn't:** raw qubit count is commoditized at about 30 statevector qubits per node across all incumbents and is not where the bar is interesting. The meaningful bars are accuracy controls (MPS truncation and bond-dimension extrapolation, trajectory convergence) and differentiable time evolution. That second capability is now orphaned in Qiskit (Dynamics archived) and PennyLane (pulse removed).

### Gaps
- No official, current cross-SDK simulator benchmark (Aer vs qsim vs Lightning on identical circuits) was found.
- Aer's docs publish no qubit-count or throughput table.
- Lightning's per-device documentation pages were not fetched; only the 2024 paper abstract and the changelog were.
- The Isakov et al. paper body (qubit counts, hardware) was not read beyond the abstract.
- Gradient-cost numbers for Qiskit's (unsupported) qiskit-algorithms gradient classes were not found.

---

## 3. What limitations are documented (qudits, time-dependent Hamiltonians, pulse fidelity, precision, gradient restrictions, non-Markovian noise, mid-circuit measurement)?

### Takeaway
The documented limits cluster as follows:
1. **Qubit-only assumptions.** Aer, qsim and PennyLane after its qudit removal are qubit-only. Cirq is the exception with first-class qudits. Dynamics handles arbitrary levels but is archived.
2. **Time-dependent Hamiltonians and pulse simulation** have been removed or orphaned in both Qiskit and PennyLane. Cirq never offered them.
3. **Precision defaults are a trap.** Cirq defaults to complex64, and the JAX-based tools need explicit x64.
4. **Gradient restrictions.** Adjoint and backprop are analytic-only and statevector-only; parameter-shift cost is linear in the number of parameters.
5. **Noise is always per-instruction Markovian CPTP,** built from average error rates. No SDK documents non-Markovian noise.
6. **Mid-circuit measurement and feed-forward are supported everywhere,** but they forfeit the "sample once at the end" shortcut.

### Cited Findings

#### Qubits vs qudits
- Cirq: "In Cirq, qudits work exactly like qubits except they have a `dimension` attribute different than 2, and they can only be used with gates specific to that dimension." The page shows a qutrit gate in a circuit — [Cirq qudits](https://quantumai.google/cirq/build/qudits)
- PennyLane shipped `default.qutrit` and `default.qutrit.mixed` through v0.45.1, and its qutrit demo uses `default.qutrit` and cites Galda et al. on "Ternary Decomposition of the Toffoli Gate on Fixed-Frequency Transmon Qutrits". All qutrit/qudit functionality has since been removed on `main` (PR #9867, merged 2026-07-24) — [Qutrits demo](https://pennylane.ai/qml/demos/tutorial_qutrits_bernstein_vazirani); [PR 9867](https://github.com/PennyLaneAI/pennylane/pull/9867); [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md)
- qsim: the memory model is stated "for an N-qubit circuit". qsimcirq raises `NotImplementedError` for matrix gates on more than 6 qubits ("only up to 6-qubit gates are supported") and for controlled gates with more than 4 targets — [qsim choose_hw](https://quantumai.google/qsim/choose_hw); [qsim_circuit.py](https://github.com/quantumlib/qsim/blob/main/qsimcirq/qsim_circuit.py)
- Aer: every documented method is qubit-based (statevector, density matrix, stabilizer, MPS, unitary, superop). An issue search for "qudit" in qiskit-aer returned no issues — [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html); [issue search](https://github.com/Qiskit/qiskit-aer/issues?q=qudit)
- Qiskit Dynamics supports arbitrary subsystem dimensions (`dim = 3`, 4, 5 and 10 appear in its docs) and multi-level `DynamicsBackend` outcomes — [docs sources](https://github.com/qiskit-community/qiskit-dynamics/tree/main/docs); [hex-results note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/backend-hex-results-fd3762b9188cdd01.yaml)

#### Time-dependent Hamiltonians, pulse level, leakage, DRAG
- Qiskit 2.0 removed Pulse. Dynamics is archived and pinned to `qiskit <= 1.3` for pulse workflows, qiskit-experiments 0.8 and qiskit-ibm-runtime ≤0.36.1. Aer removed PulseSimulator in 0.14. qiskit-experiments removed "pulse and calibration management related code" (rough calibrations such as Rabi/DRAG calibrations are gone). Only gate-level "Fine" amplitude and DRAG characterizations remain — [Qiskit 2.0 note](https://github.com/Qiskit/qiskit/blob/stable/2.0/releasenotes/notes/2.0/remove-pulse-eb43f66499092489.yaml); [Dynamics bound note](https://github.com/qiskit-community/qiskit-dynamics/blob/main/releasenotes/notes/bound-qiskit-version-93a09ea5a3c4afbb.yaml); [Aer 0.14 note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.14/remove_pulse_simulator-f8de2f6d380f446a.yaml); [QE PR 1511](https://github.com/qiskit-community/qiskit-experiments/pull/1511); [QE library](https://github.com/qiskit-community/qiskit-experiments/blob/main/qiskit_experiments/library/__init__.py)
- PennyLane's `pennylane.pulse` (ParametrizedHamiltonian/Evolution, transmon drive and interaction builders, `stoch_pulse_grad`, `pulse_odegen`) is removed on `main` ("slated for removal in PL2") — [PR 10238](https://github.com/PennyLaneAI/pennylane/pull/10238); [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md)
- PennyLane's optimal-control demo states explicitly: "we do not consider sources of noise in the system, such as leakage, dephasing" — [Optimal control demo](https://pennylane.ai/qml/demos/tutorial_optimal_control)
- Cirq's documentation tree has no pulse-level or ODE-simulation guide; its sections are build, simulate, noise/QCVV, transform, hardware and Google — [Cirq docs tree](https://github.com/quantumlib/Cirq/tree/main/docs)

#### Numerical precision
- Cirq: "the simulator uses numpy's `float32` precision (which is `complex64` for complex numbers) by default"; `np.complex128` can be requested — [Cirq simulation](https://quantumai.google/cirq/simulate/simulation)
- Aer: precision option 'single' or 'double', default 'double' — [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html)
- The Qiskit Dynamics tutorials have to call `jax.config.update("jax_enable_x64", True)` explicitly — [dynamics_backend.rst](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/dynamics_backend.rst); [systems_modelling.rst](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/systems_modelling.rst)
- The qsim MPS API is shown as `MPSStateSpace<For, float>`; the MPS is C++-only — [qsim MPS docs](https://github.com/quantumlib/qsim/blob/main/docs/mps.md); [qsim overview](https://quantumai.google/qsim/overview)
- Historical Aer MPS numerics bugs (both closed): "MPS: invariant that sum of probabilities = 1 is not preserved" (2021) and "Kraus noise crashes MPS simulator" (2020) — [Aer #1205](https://github.com/Qiskit/qiskit-aer/issues/1205); [Aer #991](https://github.com/Qiskit/qiskit-aer/issues/991)

#### Gradient-method restrictions
- PennyLane:
  - Simulation-based methods are "Reverse accumulation; not hardware compatible; statevector simulators only".
  - backprop is analytic only (shots=None); adjoint is analytic only and "raises error when shots>0".
  - Hardware-compatible methods are parameter-shift (finite-difference fallback), finite-diff, Hadamard-test variants and SPSA.
  - `"best"` prefers device gradients, then backprop, then parameter-shift, then finite differences.
  - Source: [PennyLane gradients](https://docs.pennylane.ai/en/stable/introduction/interfaces.html)
- PennyLane's TensorFlow interface is removed on `main`, as are the pulse gradient transforms — [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md)
- Open PennyLane issues touching numerics:
  - "Slow performance when compiling (tf and jax-jit) gradients through `mitigate_with_zne`" (open since 2022).
  - "[BUG] Trotter error estimate should add up expectation values of nested commutators, not states" (open, 2025).
  - Sources: [#2801](https://github.com/PennyLaneAI/pennylane/issues/2801); [#8152](https://github.com/PennyLaneAI/pennylane/issues/8152)
- Qiskit's algorithm library, which carries VQE/QAOA and the gradient classes after Qiskit 1.0, is "no longer officially supported by IBM" — [qiskit-algorithms README](https://github.com/qiskit-community/qiskit-algorithms)

#### Noise-model scope (Markovian, approximate)
- Aer:
  - Noise is a `QuantumError` (CPTP) attached after gates or resets, plus a `ReadoutError`. Idle relaxation is applied on scheduled circuits.
  - Backend-derived models are "only an approximation … built from … average error rates".
  - The stabilizer method requires Clifford errors.
  - Sources: [Building noise models](https://qiskit.github.io/qiskit-aer/tutorials/3_building_noise_models.html); [Device noise](https://qiskit.github.io/qiskit-aer/tutorials/2_device_noise_simulation.html); [Aer 0.10 note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.10/0-10-release-8c37dadcc1c82fcc.yaml); [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html)
- Cirq: the pure-state simulator handles only noise "that preserve[s] state purity", and the density-matrix simulator handles general channels — [Cirq simulation](https://quantumai.google/cirq/simulate/simulation)
- qsim: trajectory results are "stochastic in nature"; "Simulating many repetitions of a noisy circuit requires executing the entire circuit once for each repetition" — [Noisy qsimcirq](https://quantumai.google/qsim/tutorials/noisy_qsimcirq)
- PennyLane: channels include AmplitudeDamping, GeneralizedAmplitudeDamping, PhaseDamping, DepolarizingChannel and ThermalRelaxationError(t1, t2, tg) on `default.mixed` (≤23 wires) — [Noisy circuits](https://pennylane.ai/qml/demos/tutorial_noisy_circuits); [How to use noise models](https://pennylane.ai/qml/demos/tutorial_how_to_use_noise_models)
- Qiskit Dynamics: dissipation enters in Lindblad form (`static_dissipators` with jump operators) — [Rabi tutorial source](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/Rabi_oscillations.rst); [Lindblad tutorial source](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/Lindblad_dynamics_simulation.rst)

#### Mid-circuit measurement and feed-forward
- Aer supports control flow (if/else, for, while, break, continue since 0.10; switch since 0.13). Its fast-sampling paths assume "all measurements at end of the circuit" — [Aer 0.10 note](https://github.com/Qiskit/qiskit-aer/blob/main/releasenotes/notes/0.10/0-10-release-8c37dadcc1c82fcc.yaml); [AerSimulator API](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html)
- qsim: "Multiple samples can be generated with minimal additional cost for circuits with no intermediate measurements" — [qsim overview](https://quantumai.google/qsim/overview)
- Cirq: `ClassicallyControlledOperation` with key, bitmask and sympy conditions. "The Cirq built-in simulators provide support for classical control, but caution should be exercised when exporting these circuits to other environments" — [Classical control](https://quantumai.google/cirq/build/classical_control)
- PennyLane has mid-circuit-measurement demos (introduction, dynamic circuits, collecting statistics), and Lightning supports MCMs via tree traversal — [MCM intro](https://pennylane.ai/qml/demos/tutorial_mcm_introduction); [dynamic MCM how-to](https://pennylane.ai/qml/demos/tutorial_how_to_create_dynamic_mcm_circuits); [Lightning CHANGELOG](https://github.com/PennyLaneAI/pennylane-lightning/blob/main/.github/CHANGELOG.md)

### Inferences
- **Non-Markovian noise.** None of the consulted documentation offers non-Markovian noise; every noise API is per-operation CPTP or GKSL Lindblad. Modelling memory effects would require explicitly simulating environment modes, e.g. a bosonic bath in a Lindblad or Hamiltonian model. Among the incumbents only Qiskit Dynamics (archived) could do that.
- **Leakage-aware simulation** (transmon as a qutrit plus DRAG) is now effectively unsupported in maintained incumbent code paths. The only maintained route would be Cirq qudits with hand-built gates and channels, which has no time-dependent Hamiltonian solver.
- **Precision pitfalls are documented but left to the user** (complex64 default in Cirq; x64 opt-in in JAX). A tool with arbitrary-precision or verified arithmetic would address a real, documented pain point. That is an inference about demand, not a cited complaint.

### Gaps
- I found no documented statement anywhere on non-Markovian noise. That absence is my inference, not a cited limitation.
- qsim's handling of non-qubit `cirq.Qid`s was not tested; only the absence of qudit handling in the Python bridge was observed.
- Lightning's support for mixed states or noise channels was not verified.
- PennyLane's `default.qubit` precision defaults were not verified.
- The status of Cirq's `CliffordSimulator` was not separately confirmed (only the `clifford/` module exists).

---

## 4. Which standard named benchmark problems does each use (exact links for head-to-head comparison)?

### Takeaway
Recognized, citable benchmark problems recur across the SDKs. Each has an exact tutorial that a later study can reproduce:
- H2 (and LiH, H3+, HeH+, H6) ground-state energies in STO-3G against FCI and chemical accuracy (1.6 mHa).
- MaxCut QAOA on small graphs, with Farhi et al. as the canonical reference.
- Quantum natural gradient (Stokes et al.) and barren-plateau gradient-variance scaling (McClean et al.).
- ZNE (Temme et al.) and the 127-qubit kicked-Ising "utility" experiment (Kim et al. 2023).
- Quantum volume heavy-output test (Cross et al.).
- RB EPC/EPG (Magesan et al.), T1/T2*/T2-echo fits and XEB / random circuit sampling (Arute et al. 2019).
- TFIM, Heisenberg and Fermi–Hubbard Trotter dynamics with commutator-scaling Trotter error (Childs et al. 2021).
- Transmon Rabi and Lindblad T1/T2 dynamics, and X/CNOT/Toffoli pulse optimization.

### Cited Findings

#### Chemistry (STO-3G molecules)
- **H2 (PennyLane):**
  - *A brief overview of VQE* loads the H2 dataset Hamiltonian and runs a `DoubleExcitation` ansatz on Hartree–Fock in `lightning.qubit`. Its commented geometry is ±0.70108983 (bohr) — [tutorial_vqe](https://pennylane.ai/qml/demos/tutorial_vqe)
  - *Differentiable Hartree-Fock* optimizes H2 geometry jointly with circuit parameters to "-1.1373060483 Ha" — [tutorial_differentiable_HF](https://pennylane.ai/qml/demos/tutorial_differentiable_HF)
  - *How to implement VQD* uses H2 at `bondlength=0.742` in STO-3G (excited states) — [tutorial_vqe_vqd](https://pennylane.ai/qml/demos/tutorial_vqe_vqd)
  - *Accelerating VQEs with QNG* uses H2 at `bondlength=0.7` and compares `GradientDescentOptimizer` with `QNGOptimizer(approx="block-diag")` (max 500 iterations), reporting accuracy vs FCI in Ha and kcal/mol — [tutorial_vqe_qng](https://pennylane.ai/qml/demos/tutorial_vqe_qng)
  - *VQE in different spin sectors* uses H2 at ±0.6614 bohr, singlet and triplet — [tutorial_vqe_spin_sectors](https://pennylane.ai/qml/demos/tutorial_vqe_spin_sectors)
  - The SPSA demo uses H2 at ±0.6614 bohr and runs on the PennyLane-Qiskit `qiskit.aer` / `qiskit.ibmq` devices — [tutorial_spsa](https://pennylane.ai/qml/demos/tutorial_spsa)
- **H2 (IBM):** the Qunova HiVQE partner tutorial produces "Dissociation PES curves" — [IBM tutorial: qunova-hivqe](https://quantum.cloud.ibm.com/docs/en/tutorials/qunova-hivqe)
- **LiH (PennyLane):** *Adaptive circuits for quantum chemistry*, LiH in STO-3G: "the exact energy of the ground electronic state of LiH which is -7.8825378193 Ha", reached "within chemical accuracy" — [tutorial_adaptive_circuits](https://pennylane.ai/qml/demos/tutorial_adaptive_circuits)
- **H3+ and reactions (PennyLane):**
  - *Optimization of molecular geometries* for H3+: "within chemical accuracy (0.0016 Ha)" — [tutorial_mol_geo_opt](https://pennylane.ai/qml/demos/tutorial_mol_geo_opt)
  - *Modelling chemical reactions*: an H2 bond-stretch scan along r, then a linear three-hydrogen system (H at 0, r and 4.0 bohr; "three electrons in six spin orbitals" in STO-3G) as a reaction coordinate — [tutorial_chemical_reactions](https://pennylane.ai/qml/demos/tutorial_chemical_reactions)
- **H6 (PennyLane):** STO-3G, bond length 1.3, MPS bond dimension 30 with bond-dimension extrapolation — [tutorial_mps](https://pennylane.ai/qml/demos/tutorial_mps)
- **HeH+ (PennyLane):** STO-3G, bond length 1.5, ctrl-VQE on a coupled-transmon pulse model, converging to chemical accuracy. This relies on `qp.pulse`, now removed on `main` — [tutorial_pulse_programming101](https://pennylane.ai/qml/demos/tutorial_pulse_programming101)
- **H2O (PennyLane):** Hamiltonian construction in STO-3G — [tutorial_quantum_chemistry](https://pennylane.ai/qml/demos/tutorial_quantum_chemistry)
- **Cirq / ReCirq:** Hartree–Fock VQE (molecular data pages) and QC-QMC — [ReCirq HF-VQE](https://quantumai.google/cirq/experiments/hfvqe); [QC-QMC](https://quantumai.google/cirq/experiments/qcqmc)
- **IBM SQD / Krylov:** chemistry, nuclear and lattice Hamiltonians — [SQD tutorial](https://quantum.cloud.ibm.com/docs/en/tutorials/sample-based-quantum-diagonalization); [Krylov tutorial](https://quantum.cloud.ibm.com/docs/en/tutorials/krylov-quantum-diagonalization)
- **Canonical references** (cited in the demo sources above): Peruzzo et al. VQE [arXiv:1304.3061](https://arxiv.org/abs/1304.3061), cited by [Cirq VQE](https://quantumai.google/cirq/experiments/variational_algorithm); Kandala et al. 2017 hardware-efficient VQE [arXiv:1704.05018](https://arxiv.org/abs/1704.05018), cited by [phase-transitions demo](https://pennylane.ai/qml/demos/tutorial_quantum_phase_transitions)

#### Optimization (QAOA MaxCut)
- *QAOA for MaxCut*: a 4-vertex graph on 4 wires. "For `n_layers=1`, we find an objective function value of around C=3", while `n_layers=2` "recover[s] the optimal" cut. It cites Farhi et al. [arXiv:1411.4028](https://arxiv.org/abs/1411.4028) — [tutorial_qaoa_maxcut](https://pennylane.ai/qml/demos/tutorial_qaoa_maxcut)
- *Intro to QAOA* uses the graph with edges [(0,1),(1,2),(2,0),(2,3)] — [tutorial_qaoa_intro](https://pennylane.ai/qml/demos/tutorial_qaoa_intro)
- *Quantum natural SPSA*: QAOA on a 4-node, 4-edge random graph — [qnspsa](https://pennylane.ai/qml/demos/qnspsa)
- ReCirq QAOA reproduces Google's hardware QAOA, including MaxCut, Ising and landscape analysis — [ReCirq QAOA](https://quantumai.google/cirq/experiments/qaoa); [QAOA MaxCut](https://quantumai.google/cirq/experiments/qaoa/qaoa_maxcut)
- IBM QAOA tutorials — [IBM QAOA](https://quantum.cloud.ibm.com/docs/en/tutorials/quantum-approximate-optimization-algorithm); [Advanced QAOA](https://quantum.cloud.ibm.com/docs/en/tutorials/advanced-techniques-for-qaoa)

#### Gradients, natural gradient and barren plateaus
- QNG: [tutorial_quantum_natural_gradient](https://pennylane.ai/qml/demos/tutorial_quantum_natural_gradient) uses the block-diagonal Fubini–Study approximation and cites Stokes et al. [arXiv:1909.02108](https://arxiv.org/abs/1909.02108)
- Barren plateaus: [tutorial_barren_plateaus](https://pennylane.ai/qml/demos/tutorial_barren_plateaus) "partly reproduce[s] some of the findings in McClean et. al., 2018", computing gradient variance over random circuits as qubit number grows. Local-cost mitigation: [tutorial_local_cost_functions](https://pennylane.ai/qml/demos/tutorial_local_cost_functions), citing Cerezo et al. [arXiv:2001.00550](https://arxiv.org/abs/2001.00550)
- Adjoint differentiation: [tutorial_adjoint_diff](https://pennylane.ai/qml/demos/tutorial_adjoint_diff), citing Jones & Gacon [arXiv:2009.02823](https://arxiv.org/abs/2009.02823)

#### Error mitigation and the utility experiment
- *Error mitigation with Mitiq and PennyLane*: ZNE with `fold_global` and Richardson extrapolation on `default.mixed`, 4 wires — [tutorial_error_mitigation](https://pennylane.ai/qml/demos/tutorial_error_mitigation)
- *Differentiating quantum error mitigation transforms*: `DepolarizingChannel` with p = 0.05 on 4 wires — [tutorial_diffable-mitigation](https://pennylane.ai/qml/demos/tutorial_diffable-mitigation)
- *Digital ZNE with Catalyst* cites Temme et al. [arXiv:1612.02058](https://arxiv.org/abs/1612.02058) — [zne_catalyst](https://pennylane.ai/qml/demos/zne_catalyst)
- *Is quantum computing useful before fault tolerance?* rescales Kim et al. 2023 (127-qubit 2D transverse-field-Ising dynamics). It uses "only 9 [qubits], placed on a 3 × 3 grid", `DepolarizingChannel` p = 0.005, 10 Trotter layers and observable Z₄ on `default.mixed`. It discusses the MPS and isometric-TNS classical simulations — [tutorial_mitigation_advantage](https://pennylane.ai/qml/demos/tutorial_mitigation_advantage)
- IBM counterparts:
  - [Combine error mitigation options (ZNE etc.)](https://quantum.cloud.ibm.com/docs/en/tutorials/combine-error-mitigation-techniques)
  - [Utility-scale PEA](https://quantum.cloud.ibm.com/docs/en/tutorials/probabilistic-error-amplification)
  - [Simulate a kicked Ising model with the TEM function](https://quantum.cloud.ibm.com/docs/en/tutorials/simulate-kicked-ising-tem)
  - [Kicked Ising with dynamic circuits](https://quantum.cloud.ibm.com/docs/en/tutorials/dc-hex-ising)

#### Verification and characterization benchmarks
- **Quantum volume (Cross et al., PRA 100, 032328):**
  - qiskit-experiments `QuantumVolume` — [QV manual](https://qiskit-community.github.io/qiskit-experiments/manuals/verification/quantum_volume.html)
  - PennyLane *Quantum volume* (num_qubits = 5; heavy-output definition; it now uses a fake Lima device because the original "has since been retired") — [quantum_volume](https://pennylane.ai/qml/demos/quantum_volume); [Cross et al. DOI](https://doi.org/10.1103/PhysRevA.100.032328)
- **Randomized benchmarking:**
  - qiskit-experiments StandardRB/InterleavedRB/PurityRB/LayerFidelity, with the documented noise parameters given in §1 — [RB manual](https://qiskit-community.github.io/qiskit-experiments/manuals/verification/randomized_benchmarking.html)
  - Cirq single- and two-qubit RB helpers — [cirq/experiments](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/experiments)
- **T1, T2*, T2-echo, Tphi:**
  - qiskit-experiments manuals on FakePerth-derived Aer noise — [T1](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/t1.html); [T2 Ramsey](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/t2ramsey.html); [T2 Hahn](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/t2hahn.html); [Tphi](https://qiskit-community.github.io/qiskit-experiments/manuals/characterization/tphi.html)
  - Cirq `t1_decay` / `t2_decay` — [cirq/experiments](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/experiments)
- **State and process tomography:** [qiskit-experiments state tomography manual](https://qiskit-community.github.io/qiskit-experiments/manuals/verification/state_tomography.html); Cirq single-, two- and n-qubit tomography helpers ([cirq/experiments](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/experiments))
- **XEB and random circuit sampling (Arute et al. 2019):**
  - Cirq XEB docs — [XEB theory](https://quantumai.google/cirq/noise/qcvv/xeb_theory)
  - ReCirq RCS — [RCS demonstration](https://quantumai.google/cirq/experiments/random_circuit_sampling/rcs_experiment_demonstration)
  - PennyLane *Beyond classical computing with qsim* (`cirq.qsim` device; cites [doi:10.1038/s41586-019-1666-5](https://doi.org/10.1038/s41586-019-1666-5)) — [qsim_beyond_classical](https://pennylane.ai/qml/demos/qsim_beyond_classical)

#### Many-body dynamics and Trotterization
- *Exploring Trotterization* (2026-07-28): H = αX + βZ with coefficients [0.2, 1.3], first- and second-order (Strang) and Suzuki formulas, Trotter error vs exact evolution, and a T-count vs Trotter-error trade-off with phase-gradient rotations. It cites Childs et al. [PRX 11, 011020 (2021)](https://doi.org/10.1103/PhysRevX.11.011020) — [exploring_trotterization](https://pennylane.ai/qml/demos/exploring_trotterization)
- *The Quantum Graph Recurrent Neural Network*: learning Ising dynamics with `trotter_step = 0.01` — [tutorial_qgrnn](https://pennylane.ai/qml/demos/tutorial_qgrnn)
- IBM:
  - [Simulate time evolution of the TFIM](https://quantum.cloud.ibm.com/docs/en/tutorials/time-evolution)
  - [Multi-product formulas to reduce Trotter error](https://quantum.cloud.ibm.com/docs/en/tutorials/multi-product-formula)
  - [AQC for time-evolution circuits](https://quantum.cloud.ibm.com/docs/en/tutorials/approximate-quantum-compilation-for-time-evolution)
  - [Heisenberg chain VQE](https://quantum.cloud.ibm.com/docs/en/tutorials/spin-chain-vqe)
  - [Compilation methods for Hamiltonian simulation](https://quantum.cloud.ibm.com/docs/en/tutorials/compilation-methods-for-hamiltonian-simulation-circuits)
- ReCirq:
  - [Fermi–Hubbard spin-charge separation](https://quantumai.google/cirq/experiments/fermi_hubbard)
  - [KPZ and the Heisenberg spin chain](https://quantumai.google/cirq/experiments/kpz/kpz)
  - [OTOC scrambling](https://quantumai.google/cirq/experiments/otoc)
  - [Lattice gauge theory](https://quantumai.google/cirq/experiments/lattice_gauge/lattice_gauge)

#### Transmon / pulse-level problems (all now archived or removed paths)
- Qiskit Dynamics tutorials:
  - [Rabi oscillations with T1/T2 Lindblad](https://qiskit-community.github.io/qiskit-dynamics/tutorials/Rabi_oscillations.html)
  - [Lindblad qubit-chain steady state](https://qiskit-community.github.io/qiskit-dynamics/tutorials/Lindblad_dynamics_simulation.html)
  - [Two-transmon DynamicsBackend with X/SX calibration and CR Hamiltonian tomography](https://qiskit-community.github.io/qiskit-dynamics/tutorials/dynamics_backend.html)
  - [Gradient X-gate pulse optimization](https://qiskit-community.github.io/qiskit-dynamics/tutorials/optimizing_pulse_sequence.html)
  - [Dyson/Magnus solvers](https://qiskit-community.github.io/qiskit-dynamics/userguide/perturbative_solvers.html)
- PennyLane:
  - [Optimal control for CNOT and Toffoli](https://pennylane.ai/qml/demos/tutorial_optimal_control)
  - [RL gate calibration on coupled transmons](https://pennylane.ai/qml/demos/tutorial_rl_pulse)
  - [OQC Lucy transmon pulses](https://pennylane.ai/qml/demos/oqc_pulse)
  - [Superconducting-qubit primer (transmon–cavity readout)](https://pennylane.ai/qml/demos/tutorial_sc_qubits)

### Inferences
- **Most reproducible head-to-head targets.** These come with exact numbers in incumbent tutorials:
  - H2/STO-3G at 0.742 Å: the coordinate ±0.70108983 bohr is 1.402 bohr ≈ 0.742 Å by my arithmetic, against the published energy −1.1373060483 Ha.
  - LiH/STO-3G exact energy −7.8825378193 Ha.
  - H3+ to chemical accuracy (0.0016 Ha).
  - H6 MPS bond-dimension convergence.
  - QAOA MaxCut on the 4-node graph (C ≈ 3 at p = 1; optimum at p = 2).
  - The 9-qubit (3×3) kicked-Ising ZNE reproduction of Kim et al.
  - The qiskit-experiments RB noise model (5e-3 and 5e-2 depolarizing) with EPC/EPG outputs.
  - FakePerth T1/T2 fits.
  - The Trotter-error study of H = 0.2X + 1.3Z.
- **The 127-qubit Kim et al. experiment is a community classical-simulation benchmark.** It is known from background literature, not re-verified here: Tindall et al., PRX Quantum 2024, belief-propagation tensor networks. It is far beyond small-scale study, but its 3×3 rescaling (the PennyLane demo) is a recognized small analogue.
- **Additional canonical references, background only, not re-fetched:** Magesan et al. RB (PRL 106, 180504, 2011), Motzoi et al. DRAG (PRL 103, 110501, 2009) and McClean et al. barren plateaus (Nat. Commun. 9, 4812, 2018; arXiv:1803.11173).

### Gaps
- No incumbent tutorial was found that publishes a full H2 or LiH dissociation curve against FCI with tabulated numbers. The chemical-reactions demo plots curves, but I didn't extract the tabulated values. Qiskit Nature tutorials were not checked.
- The exact graph and depth used in IBM's QAOA tutorial, and the ReCirq QAOA problem sizes, were not extracted.
- No incumbent tutorial was found for DRAG leakage suppression with explicit leakage numbers. qiskit-experiments' FineDrag manual page was not opened; only the class listing was.

---

## 5. Which of these problem classes would be most compelling to a master's-level QM audience, and why?

### Takeaway
The most compelling studies are those where the numerics decide the physics answer and where the incumbents' support is now weakest:
1. Driven-dissipative multilevel transmon dynamics: Rabi/Ramsey/T1/T2, leakage, DRAG, via time-dependent Lindblad equations.
2. Characterization as an inference problem: T1/T2/RB/tomography fits with physical constraints.
3. Molecular dissociation curves (H2, LiH) against FCI and chemical accuracy.
4. Trotter error vs exact propagation in TFIM/Heisenberg chains.
5. Open-system noise and ZNE bias: density-matrix vs trajectory convergence.

Barren plateaus and QAOA are standard but more about optimization than quantum mechanics.

### Cited Findings
- The pulse/ODE layer is orphaned:
  - Qiskit Pulse was removed and migration to Dynamics was "put on hold".
  - Dynamics is archived and pinned to qiskit ≤1.3.
  - PennyLane's `qp.pulse` is removed on `main`.
  - PennyLane's optimal-control demo ignores "leakage, dephasing".
  - Sources: [Qiskit 2.0 note](https://github.com/Qiskit/qiskit/blob/stable/2.0/releasenotes/notes/2.0/remove-pulse-eb43f66499092489.yaml); [Dynamics README](https://github.com/qiskit-community/qiskit-dynamics); [PR 10238](https://github.com/PennyLaneAI/pennylane/pull/10238); [optimal control demo](https://pennylane.ai/qml/demos/tutorial_optimal_control)
- The incumbents' own transmon material treats transmons as multilevel Duffing oscillators (`dim = 3`, 5 or 10), with Lindblad Γ1/Γ2 dissipators and gate optimization by autodiff through ODE solvers — [Dynamics tutorials](https://qiskit-community.github.io/qiskit-dynamics/tutorials/index.html); [systems_modelling.rst](https://github.com/qiskit-community/qiskit-dynamics/blob/main/docs/tutorials/systems_modelling.rst)
- Characterization is the one Qiskit numerics area still under active development (qiskit-experiments 0.14.x, 2026). It covers T1/T2/Tphi/RB/tomography/QV with documented fit models (a·α^m + b, EPC/EPG). Cirq offers parallel helpers and XEB — [QE releases](https://github.com/qiskit-community/qiskit-experiments/releases); [RB manual](https://qiskit-community.github.io/qiskit-experiments/manuals/verification/randomized_benchmarking.html); [cirq/experiments](https://github.com/quantumlib/Cirq/tree/main/cirq-core/cirq/experiments)
- Chemistry benchmarks with exact reference numbers and a chemical-accuracy criterion are the most widely shared problem class across galleries — [tutorial_adaptive_circuits](https://pennylane.ai/qml/demos/tutorial_adaptive_circuits); [tutorial_mol_geo_opt](https://pennylane.ai/qml/demos/tutorial_mol_geo_opt); [tutorial_differentiable_HF](https://pennylane.ai/qml/demos/tutorial_differentiable_HF)
- Trotter-error numerics are current:
  - PennyLane published *Exploring Trotterization* on 2026-07-28.
  - PennyLane has an open bug in its Trotter-error estimator.
  - IBM ships an MPF addon specifically to reduce Trotter error.
  - Sources: [exploring_trotterization](https://pennylane.ai/qml/demos/exploring_trotterization); [PennyLane #8152](https://github.com/PennyLaneAI/pennylane/issues/8152); [qiskit-addon-mpf](https://github.com/Qiskit/qiskit-addon-mpf)
- Noise and mitigation numerics:
  - PennyLane's `NoiseModel` and ZNE transforms are removed on `main`.
  - Aer is in reduced maintenance, and its device models are "approximations" from average error rates.
  - qsim offers only trajectory noise.
  - Sources: [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/main/doc/releases/changelog-dev.md); [Aer README](https://github.com/Qiskit/qiskit-aer/blob/main/README.md); [Aer device noise](https://qiskit.github.io/qiskit-aer/tutorials/2_device_noise_simulation.html); [noisy qsimcirq](https://quantumai.google/qsim/tutorials/noisy_qsimcirq)

### Inferences
Ranked for a master's-level QM audience. The ranking is my judgement; the supporting facts are cited above.

1. **Transmon as a driven, damped anharmonic oscillator (qutrit or more).** Covers Rabi and Ramsey fringes, T1/T2 via Lindblad, leakage out of the qubit subspace, DRAG derivation and verification, and the validity of the RWA and rotating frame. Numerics decide the answer: leakage population and gate infidelity depend on pulse width, anharmonicity and RWA validity, which have no closed form. Master's QM students know the harmonic and anharmonic oscillator, perturbation theory and master equations, so this connects textbook QM to hardware. The incumbent bar is the archived Qiskit Dynamics tutorials and the (now removed) PennyLane pulse demos. Directly comparable targets: DynamicsBackend two-transmon `dim = 3`, the X-gate pulse optimization, and the Dyson/Magnus 2–4× and 10–60× speed claims.
2. **Characterization experiments as numerical inverse problems.** T1 and T2*/T2-echo fits, Ramsey detuning, RB decay to EPC/EPG, and tomography with positivity constraints, simulated from a Lindblad or Kraus model and then fitted. Statistics meets open-system QM. qiskit-experiments is the active bar with concrete noise parameters (FakePerth thermal relaxation; RB depolarizing 5e-3 and 5e-2), so results can be compared directly.
3. **Molecular dissociation (H2 → LiH) to chemical accuracy.** Covers static correlation at stretched bonds, Hartree–Fock failure, FCI vs VQE ansätze, and Jordan–Wigner qubit counts. It is the most common cross-SDK benchmark with exact reference energies (−1.1373060483 Ha for H2; −7.8825378193 Ha for LiH; 1.6 mHa accuracy). It is less novel because PennyLane's coverage is deep.
4. **Trotterized many-body dynamics vs exact propagation** (TFIM, Heisenberg, kicked Ising 3×3). Commutator-scaling error bounds, first- vs second-order formulas, MPF extrapolation and light-cone or entanglement growth. Numerics decide step-size and accuracy trade-offs. The incumbents' 2026 material (PennyLane demo, IBM MPF/AQC tutorials) shows the topic is current.
5. **Open-system noise and its mitigation.** Covers exact density-matrix evolution vs stochastic trajectories (convergence ∝ 1/√N), ZNE extrapolation bias, and Markovian vs non-Markovian effects. The last is an explicit blind spot: no incumbent documents non-Markovian noise.
6. **Lower priority for this audience:**
   - Barren plateaus: concentration of measure is compelling, but it is QML-centric.
   - QAOA MaxCut: combinatorics more than QM.
   - Quantum volume and XEB: protocol definitions, though QV's heavy-output threshold of 2/3 is a crisp, checkable numeric criterion.

### Gaps
- No survey data or curriculum evidence on master's-course interest in these topics was gathered. The ranking is reasoned from physics content and incumbent coverage, not measured demand.
- I didn't check whether incumbent courses (e.g. IBM Quantum Learning, PennyLane Codebook) use these problems pedagogically; only docs, galleries and repos were examined.
