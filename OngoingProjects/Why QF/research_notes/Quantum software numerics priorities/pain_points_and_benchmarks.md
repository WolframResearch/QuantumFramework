# Quantum simulation software: numerical pain points, independent benchmarks, and teaching use

Compiled 2026-10-04. The evidence comes from:

- arXiv and journal papers.
- GitHub issue trackers, read through the GitHub API (issue bodies, dates, states, archive flags).
- The Stack Exchange API.
- Project forums (PennyLane Discourse, the QuTiP Google Group, Julia Discourse).
- Official documentation.

Conventions used throughout:

- "Older" means before 2023.
- "Title only" means the issue body was not read.
- "Snippet" means the claim comes from a search-engine excerpt, not the fetched page.
- When one reporter filed several issues, they are counted once.
- Developer- or vendor-authored comparisons are labelled as such. Only a few studies are independent.

## 1. What do comparison papers, benchmark studies and reviews (2022–2026) report about the frameworks' relative strengths and weaknesses?

### Takeaway
Independent (non-developer) benchmarks exist almost only for gate-model state-vector simulators. They find run-time spreads of more than 100× that depend on task, size, precision and hardware, plus installation, precision-labelling and GPU-correctness defects.

For open-system, bosonic, photonic, HEOM and stabilizer tools, almost every comparison is written by the developers of the winning package. These comparisons agree on one picture:

- QuTiP is the universal baseline. It pays Python per-call overhead and relies on the qutip-jax plug-in for GPU and autodiff.
- Julia tools (QuantumToolbox.jl, HierarchicalEOM.jl) and JAX tools (Dynamiqs) win on speed, batching and gradients.
- None of these comparisons measures accuracy or robustness: tolerance, truncation or precision sensitivity.

### Cited Findings

**A. Independent gate-model simulator studies** (authors with no stake in the compared packages)

- **Jamadagni, Läuchli & Hempel, SciPost Phys. Core 7, 075 (2024)** — [arXiv 2401.09076](https://arxiv.org/abs/2401.09076); [SciPost](https://scipost.org/SciPostPhysCore.7.4.075); [HTML v1](https://arxiv.org/html/2401.09076v1)
  - Containerized benchmark on an HPC cluster, about 24 package configurations: Braket, Cirq, cuQuantum, HiQ, HybridQ, Intel-QS, myQLM, PennyLane, ProjectQ, Qcgpu, Qibo/Qibojit, Qiskit, QPanda, Qrack, qsimcirq, QuEST, Qulacs, Quantum++, SV-Sim, Yao.
  - Three tasks:
    - Trotterized XYZ-Heisenberg chain dynamics.
    - Sycamore-style random circuits with fSim gates.
    - Quantum Fourier transform (QFT).
  - "Depending on the simulation task, the computational performance of packages can differ by more than two orders of magnitude at both small and large problem sizes."
  - Hardware and parallelism:
    - GPUs gave up to about 500× at large N but carry overhead at small N.
    - Multithreading gave 5–30× for N > 20.
    - Most packages have constant overhead until about 15–20 qubits.
  - Fastest packages:
    - qsimcirq is best in the large-N, single-precision limit.
    - qiskit and qpanda are among the fastest in double precision.
    - yao and qrack have low small-N overhead.
  - Defects found:
    - Qiskit and PennyLane ran in double precision despite documenting single-precision support.
    - PennyLane on GPU gave inconsistent random-circuit results.
    - Qcgpu misbehaved at N=32 with no error.
    - QPanda timed out on GPU.
    - HybridQ is capped at 32 qubits.
    - Qiskit and myQLM had installation or compile failures.

- **Cerrudo-Herrera et al., *SIMULATION* (online 22 Apr 2026)**, abstract-level only — [DOI](https://doi.org/10.1177/00375497261438515); [code](https://github.com/Cerrudoxx/ComputaexQuantumBenchmark)
  - Seven state-vector simulators: Qiskit, Qulacs, Qibo, Qsimov, Cirq, PennyLane, Intel-QS.
  - Circuits: Grover, QFT and quantum volume, at 3–30 qubits on one HPC node.
  - Qulacs is fastest below 22 qubits; Qiskit is fastest beyond that.
  - Qiskit and Qulacs scale best across cores.
  - Intel-QS uses the least memory below 24 qubits but is slow on Grover.

**B. Developer- or vendor-authored comparisons (self-reported)**

- **VQCSim, Firc et al., ICCAD 2026** — [arXiv 2607.11985](https://arxiv.org/abs/2607.11985)
  - The authors benchmark their own simulator against generic frameworks; affiliations were not checked.
  - In hybrid ML workflows, "framework dispatch and orchestration overhead" dominates run time.
  - Their compile-once PyTorch state-vector simulator gave pooled median speedups of 4.49× (inference) and 26.78× (training) on MQT Bench circuits.

- **Malarchick, arXiv 2603.18052 (Mar 2026, rev. Sep 2026)** — [arXiv](https://arxiv.org/abs/2603.18052)
  - A single author compares their own C, GPU and FPGA kernels with QuTiP. They have no stake in QuTiP, but the comparison is not neutral either.
  - Dense Lindblad propagation at transmon-relevant dimensions d = 3, 9, 27.
  - Custom C code beat QuTiP 5.2.3 by up to two orders of magnitude at d=3; the advantage inverted by d=27.
  - The dense LU solve inside Padé propagator construction took more than 90% of run time at the largest d.
  - This quantifies Python per-call overhead for small problems.

- **QuantumToolbox.jl, Mercurio et al., *Quantum* 9, 1866 (2025)** — [Quantum](https://quantum-journal.org/papers/q-2025-09-29-1866/); [arXiv 2504.21440](https://arxiv.org/abs/2504.21440); [benchmark code](https://github.com/albertomercurio/QuantumToolbox.jl-Paper-Figures)
  - Compared against QuTiP 5.2.1 (GPU via qutip-jax 0.1.1), QuantumOptics.jl 1.2.3 and dynamiqs 0.3.3. Versions are as of 28 Aug 2025. Hardware: i9-13900KF CPU and RTX 4090 GPU.
  - Benchmark models:
    - A driven Kerr oscillator, H = Δa†a − U a†²a² + F(a + a†), with single-photon loss. Used for mesolve, mcsolve, smesolve and autodiff.
    - A 1D transverse-field Ising chain with local σ⁻ dissipators, up to N=12 spins on GPU.
  - Results:
    - "QuantumToolbox.jl consistently outperforms the other libraries across all tested scenarios."
    - Its own autodiff is "experimental" but competitive.
    - The autodiff panel marks at least one package "N/A".
  - Critique of Python packages:
    - They are tied to a single backend (SciPy or JAX).
    - Adding GPU and autodiff "required major restructuring."
    - JAX "restricts computations to JAX-compatible operations, limiting flexibility in data structures and control flows."
  - Self-acknowledged Julia weaknesses: a young ecosystem, compile latency and a smaller community.
  - Conflict note: the authors overlap with the QuTiP team (Nori), and the package lives in the qutip GitHub organisation.

- **Dynamiqs (Alice & Bob, ENS, Sherbrooke)** — [IEEE QCE'24 poster / MQSF 2025 pitch](https://www.cda.cit.tum.de/images/mqsf2025/pitches/02_MQSF2025_Dynamiqs_Gouzien.pdf)
  - Built on JAX and Diffrax, with a QuTiP-like API, GPU support, batching and end-to-end differentiability.
  - The poster states the gradient trade-offs: autodiff is "fast and reliable, but large memory", the adjoint method is "low memory, but slower", and recursive checkpointing is recommended.
  - Vendor blog claims (snippet only; the page returned 404 on 2026-10-04) — [blog](https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/):
    - On CPU, Dynamiqs is "on par with QuTiP", with speedups at large batch sizes.
    - v0.2.2 is up to 30× faster.
    - 100,000 qubit SME trajectories took 10 min in QuTiP on CPU, 20 s in Dynamiqs on CPU and 0.7 s in Dynamiqs on GPU.

- **HierarchicalEOM.jl, *Commun. Phys.* (2023)**, compared with QuTiP-BoFiN (*PRR* 5, 013181, 2023) — [arXiv 2306.07522](https://arxiv.org/abs/2306.07522); [Commun. Phys.](https://www.nature.com/articles/s42005-023-01427-2); [BoFiN](https://arxiv.org/abs/2010.10806)
  - Faster HEOM matrix construction: multithreaded, while BoFiN is single-threaded.
  - Faster time evolution and steady states.
  - The two teams overlap.

- **PQLS, Simanovskis, Adve & Rostampoor, arXiv 2609.12309 (Sep 2026)** — [arXiv](https://arxiv.org/abs/2609.12309)
  - Steady-state Lindblad solver, benchmarked against QuTiP, QuTiP-JAX and RydIQule on CPU and GPU.
  - Claims "substantial computational speedups" from batched tensor processing that cuts Python overhead in parameter sweeps.

- **Piquasso, Kolarovszki et al., *Quantum* 9, 1708 (2025)**, compared with Strawberry Fields 0.23.0 and The Walrus 0.21.0 — [arXiv 2403.04006](https://arxiv.org/abs/2403.04006); [PDF](https://quantum-journal.org/wp-content/uploads/2025/04/q-2025-04-15-1708.pdf)
  - Strawberry Fields uses a per-mode ("local") Fock cutoff, so storage scales as c^d. Piquasso uses a total-photon ("global") cutoff.
  - With a global cutoff, passive elements cause no further probability loss after truncation; with a local cutoff they generally do.
  - On a 10-mode, 4-layer CV neural network, Piquasso is faster than Strawberry Fields at every level of probability loss.
  - "We were forced to stop the calculation with Strawberry Fields due to excessive memory usage" (cutoff 10, 12 × 16 GB RAM).

- **Leonteva et al., ACM (2025)** — [arXiv 2504.14027](https://arxiv.org/abs/2504.14027); [DOI](https://doi.org/10.1145/3776567)
  - Authors are at QPerfect SAS, which makes MIMIQ, one of the compared emulators. This is not independent.
  - Seven emulators: Qiskit-MPS, Quimb-MPS, QMatchaTea-MPS, MIMIQ-MPS, Quimb-TN, Pyqrack, MQT-DDS. Circuits: 13 MQT Bench circuits from 4 to 1,024 qubits, with a 300 s limit.
  - Problems solved at 100 qubits:

    | Emulator | Solved at 100 qubits |
    |---|---|
    | MIMIQ-MPS | 12 |
    | QMatchaTea-MPS | 10 |
    | Quimb-MPS | 9 |
    | Qiskit-MPS | 8 |
    | MQT-DDS | 3 |
    | Quimb-TN | 3 |
    | Pyqrack | 1 |

  - Quimb-TN needed manual hyperparameter tuning.
  - MIMIQ was improved during the study.

- **Mazumder et al., arXiv 2607.09882 (Jul 2026, rev. Sep 2026)** — [arXiv](https://arxiv.org/abs/2607.09882)
  - Approximate simulators: BlueQubit, AWS Braket, Quantum Rings, Qiskit pauli-prop and PauliPropagation.jl. The authors' affiliations were not checked for conflicts.
  - Pauli-path simulation of IBM's 127-qubit kicked-Ising benchmark ran up to about 1,700× faster on GPU at fine truncation.
  - Only GPUs reached truncation thresholds δ < 10⁻⁵.

- **QuaSARQ, Osama, Thanos & Laarman, arXiv 2603.14641 (2026)**, GPU stabilizer simulator — [arXiv](https://arxiv.org/abs/2603.14641)
  - Compared with Stim, Qiskit-Aer, Qibo, Cirq and PennyLane.
  - Up to 180,000 qubits at depth 1,000 (about 130 million gates).
  - Claims up to 105× speedup and more than 80% energy reduction, mainly in many-shot sampling.
- **Stim (Gidney, *Quantum* 2021, older)** set the stabilizer-simulation baseline over Qiskit's and Cirq's Clifford simulators — [arXiv 2103.02202](https://arxiv.org/abs/2103.02202)
- **QuantumClifford.jl** claims low-level parity with Stim: about 15 µs for an in-place multiplication of 1M-qubit Paulis — [repo](https://github.com/QuantumSavory/QuantumClifford.jl)

- **Benchpress, Nation et al.**, all authors IBM — [arXiv 2409.08844](https://arxiv.org/abs/2409.08844); [*Nature Computational Science* 2025](https://www.nature.com/articles/s43588-025-00792-y)
  - Covers circuit construction, manipulation and transpilation only, not simulation. Seven SDKs: Braket 1.86.1, BQSKit 1.1.2, Cirq 1.4.1, Tket 1.31.0, Qiskit 1.2.0, QTS 0.4.8, Staq 3.5.
  - Construction tests: Qiskit 2.0 s, Tket 14.2 s with one failure, BQSKit 50.9 s with two failures.
  - Cirq was 55× faster than Qiskit on the Hamiltonian-simulation construction test.
  - Transpilation failures: BQSKit 19%, Tket 8%, QTS 2%.
  - Tket beat Qiskit on two-qubit-gate depth overall.

- **Other developer claims in machine-learning simulators**
  - TensorCircuit (older, 2022) claims speedups over Qiskit of "nearly a million times" for gradients of moderate circuits — [arXiv 2205.10091](https://arxiv.org/abs/2205.10091)
  - TQml (2025) claims up to about 10× over PennyLane-Torch on single-threaded CPU, with PennyLane overtaking it on GPU at the largest sizes (snippet) — [arXiv 2506.04891](https://arxiv.org/html/2506.04891v3)

**C. Reviews and qualitative comparisons**

- **QuTiP 5, *Physics Reports* 1153 (Jan 2026), 62 pp.** — [arXiv 2412.04705](https://arxiv.org/abs/2412.04705)
  - QuTiP has been "at the forefront of open-source quantum software for the past 13 years."
  - QuTiP 5 adds a pluggable data layer (JAX, CuPy), plus qutip-qip and qutip-qoc.
- **Wayo et al., photonic-simulator review, arXiv 2502.05245 (Feb 2025)** — [arXiv](https://arxiv.org/abs/2502.05245)
  - Covers Strawberry Fields, Piquasso, QuTiP, SimulaQron, Perceval and QuantumOptics.jl.
  - Non-Gaussian states need tensor networks or GPUs.
  - A 2026 *EPJ Quantum Technology* review with a similar scope exists. Its relation to the arXiv paper was not verified. — [EPJ QT](https://link.springer.com/article/10.1140/epjqt/s40507-026-00496-w)
- **Larssen, Peng & Markidis, arXiv 2605.23034 (May 2026)** — [arXiv](https://arxiv.org/abs/2605.23034)
  - Compares effective two-level, three-mode Duffing and circuit-based transmon models for flux-tunable qubits.
  - Benchmark suite: flux-dependent spectra, two-qubit interaction parameters, driven dynamics, CZ gate, leakage, runtime.
  - Multilevel models reveal driven-dynamics effects that two-level models miss.
  - The authors recommend circuit-based models "for high-fidelity reference simulation or detailed leakage analysis."
- **Díaz-Camacho et al., arXiv 2512.01807 (Dec 2025)** — [arXiv](https://arxiv.org/abs/2512.01807)
  - Distributed-QC emulators (Qiskit Aer, SquidASM, Interlin-q, SQUANCH) on a distributed inverse QFT.
  - "Many platforms either lacked support for teleportation protocols or required complex workarounds."
- **QuantumOptics.jl, Krämer et al., *CPC* 2018 (older)** — [arXiv 1707.01060](https://arxiv.org/abs/1707.01060); [benchmark repo](https://github.com/qojulia/QuantumOptics.jl-benchmarks)
  - Compared usability and performance with QuTiP and the MATLAB QOToolbox.
  - Its benchmark repository's last commit was 2020-11-24, so it is stale.

**D. Ecosystem changes that alter any 2026 comparison** (verified through the GitHub API on 2026-10-04 unless noted)

- **Archived:**
  - Strawberry Fields is archived; its last release was v0.23.0-post1 on 2022-06-02 — [repo](https://github.com/XanaduAI/strawberryfields)
  - MrMustard is archived — [repo](https://github.com/XanaduAI/MrMustard)
  - The Walrus is archived — [repo](https://github.com/XanaduAI/thewalrus)
  - On the PennyLane forum (Feb 2026), Strawberry Fields is called deprecated and users are pointed to PennyLane — [forum](https://discuss.pennylane.ai/t/fock-backend-truncation/3545)
- **IBM stack:**
  - qiskit-dynamics is archived. Its README says "This repo is no longer being actively maintained." — [repo](https://github.com/qiskit-community/qiskit-dynamics)
  - Qiskit Pulse was deprecated for removal in Qiskit 2.0 as "not in scope for the direction of the project anymore." — [Qiskit #13063](https://github.com/Qiskit/qiskit/issues/13063)
  - qiskit-algorithms, home of Qiskit's gradient code, says it "is no longer officially supported by IBM" (not archived) — [repo](https://github.com/qiskit-community/qiskit-algorithms)
- **Qudits:**
  - PennyLane PR #9867 "Remove everything related to qudits" was merged 2026-07-24 and will ship in 0.46 — [PR](https://github.com/PennyLaneAI/pennylane/pull/9867)
  - Cirq issue #7461 "Deprecate qudit support?" was opened 2025-07-02 and closed as not-planned on 2025-07-03 — [issue](https://github.com/quantumlib/Cirq/issues/7461)
- **TensorFlow Quantum** is not archived and was revived in 2025–26 (v0.7.5 and v0.7.6). It still lags Keras 3 and current Cirq — [repo](https://github.com/tensorflow/quantum); [#1090](https://github.com/tensorflow/quantum/issues/1090)

### Inferences
- The independent literature answers "which state-vector simulator is fastest for N qubits." It does not answer which open-system or bosonic tool is accurate, robust or convenient. All open-system speed claims are self-reported, on models chosen by the authors, and none reports error against a converged reference.
- Across the open-system comparisons, the competitive axes are GPU use, batching over parameters, autodiff, and low per-call overhead. That pattern suggests incumbents compete on throughput, not on correctness guarantees such as convergence in truncation, tolerance or precision.
- 2025–26 consolidation has removed maintained tools from three niches:
  - Fock-space photonics: the Xanadu stack is archived.
  - Pulse-level simulation in the IBM stack: Qiskit Pulse removed, qiskit-dynamics archived.
  - Qudits in PennyLane: removed.
- Benchmarks also mislabel or silently change numerical precision. Single precision was claimed but not used (Jamadagni), and MIMIQ changed during the study (Leonteva). Reproducible, version-pinned comparisons are rare.

### Gaps
- No independent benchmark of QuTiP against QuantumOptics.jl, QuantumToolbox.jl and Dynamiqs from 2022–2026 was found. All are developer-authored.
- No cross-package accuracy study for Lindblad solvers was found, such as error versus tolerance or truncation on a common reference problem.
- Three papers were read only at abstract level:
  - Cerrudo-Herrera 2026 (SAGE paywall).
  - The ACM version of Leonteva (HTTP 403); the arXiv PDF was read instead.
  - The *Nature Computational Science* version of Benchpress, which may differ from the arXiv HTML.
- Not found: independent tensor-network library comparisons (quimb, ITensor, TeNPy, cuTensorNet) from 2022–2026, or independent stabilizer-simulator comparisons beyond developer papers.
- The Dynamiqs vendor speedups could not be re-verified because the blog page now returns 404.

## 2. Which numerical pain points recur in GitHub issues, Stack Exchange and forums?

### Takeaway
The most persistent and best-documented problem is adaptive ODE steppers missing or under-resolving time structure: thin or delayed pulses and fast lab-frame oscillations. In QuTiP this has been reported by at least six distinct people from 2020 to 2026. Each time the fix is manual `max_step`/`rtol`/`atol` tuning, and in 2026 a user asked for automatic pulse detection.

Close behind:

- **Truncation convergence is left to the user.** Examples: scqubits ("8 months" to notice), Strawberry Fields and MrMustard cutoffs, QuTiP forum advice.
- **Silent wrong answers in near-singular linear algebra.** Fidelity above 1, steady states that are wrong or depend on the initial guess.
- **Normalization failures from single precision combined with fixed tolerances.** At least 8 reports across 4 circuit frameworks.
- **Qudit and mixed-dimension support is shrinking** in circuit SDKs.
- **Mid-circuit measurement and feedforward** cost exponential memory or exact gradients.
- **Autodiff paths are fragile.**
- **Reproducibility breaks** mainly through dependency drift and non-deterministic transpilation, not RNG seeds.

### Cited Findings

**2.1 Fock-space truncation and convergence checks**

- QuTiP Google Group (Aug 2022, older), coupled driven oscillators gave different populations for cutoffs N = 2..8 — [thread](https://groups.google.com/g/qutip/c/mZdhAHvSmTk)
  - N. Lambert: results converge around N = 5–6, the basis must cover the drive-induced displacement (about 2K/ω₀), and convergence should be checked by trying several cutoffs.
- **scqubits #274 (2025-07-30)**, verified in the issue body — [issue](https://github.com/scqubits/scqubits/issues/274)
  - The default charge cutoff of 5 in custom circuits is "too small for devices with EJ > EC … it took me about 8 months to realize this was a problem."
  - Maintainer: the default "does not have any significance"; the cutoff needed for convergence depends on the circuit and the parameter regime.
- The scqubits paper (*Quantum* 2021, older) leaves truncation convergence to the user — [arXiv 2107.08552](https://arxiv.org/abs/2107.08552). There is a docs page on "Bases, truncation, and convergence" — [docs](https://scqubits.readthedocs.io/en/v4.1/guide/circuit/ipynb/custom_circuit_bases.html)
  - scqubits #306 (2026-06) reviews a new "Convergence Diagnostics" module (`scq.check_convergence`). Its release status is unknown. — [issue](https://github.com/scqubits/scqubits/issues/306)
- **Strawberry Fields**
  - #753 (2025-02-25), "Cubic phase inaccuracies due to high Fock cutoff requirement" — [issue](https://github.com/XanaduAI/strawberryfields/issues/753)
    - The reporter said it was blocking their work.
    - Xanadu staff replied (2025-02-27): "we won't be able to make improvements to the cubic phase gate."
  - PennyLane forum "Fock Backend truncation" (2023-10) — [forum](https://discuss.pennylane.ai/t/fock-backend-truncation/3545)
    - Cutoff 2 gave unnormalised probabilities.
    - Staff: gates are built recursively, so elements are exact inside the cutoff but operators are not unitary, and high cutoffs can become numerically unstable.
  - The docs advise checking that `state.trace()` stays near 1 (snippet) — [docs](https://strawberryfields.readthedocs.io/en/latest/introduction/circuits.html)
  - Older reports — [#354](https://github.com/XanaduAI/strawberryfields/issues/354), [#488](https://github.com/XanaduAI/strawberryfields/issues/488), [#670](https://github.com/XanaduAI/strawberryfields/issues/670):
    - [#354](https://github.com/XanaduAI/strawberryfields/issues/354) (2020): homodyne measurement probabilities below 0 or above 1.
    - [#488](https://github.com/XanaduAI/strawberryfields/issues/488) (2020): "pure states are not pure."
    - [#670](https://github.com/XanaduAI/strawberryfields/issues/670) (2022): random failures in `measure_fock`.
- **MrMustard automatic cutoff selection (2023–24)** — [issues](https://github.com/XanaduAI/MrMustard/issues); [#250](https://github.com/XanaduAI/MrMustard/issues/250)
  - [#250](https://github.com/XanaduAI/MrMustard/issues/250) "Autocutoffs not being respected."
  - [#238](https://github.com/XanaduAI/MrMustard/issues/238): probabilities do not update when the cutoff is raised.
  - [#396](https://github.com/XanaduAI/MrMustard/issues/396) and [#462](https://github.com/XanaduAI/MrMustard/issues/462)/#490 (AUTOSHAPE).
  - [#481](https://github.com/XanaduAI/MrMustard/issues/481): MemoryError.
  - The repository is now archived.
- **Stack Exchange**
  - "Weird behavior when simulating vacuum Rabi oscillation on QuTip" (2023-09-18), unanswered — [QC.SE 34202](https://quantumcomputing.stackexchange.com/questions/34202)
    - A Hamiltonian with a (n+n†)⁴ nonlinearity gives anomalous dynamics only at 15×15 and 30×30 truncations; other sizes behave as expected.
  - "Truncated Qumode States and Support" (2024) — [QC.SE 35477](https://quantumcomputing.stackexchange.com/questions/35477)
  - "truncation before matrix exponential: how to do it right?" (2017, older) — [Phys.SE 317182](https://physics.stackexchange.com/questions/317182)
- Piquasso documents that local-cutoff truncation keeps leaking probability under passive gates and blows up memory (see 1B) — [PDF](https://quantum-journal.org/wp-content/uploads/2025/04/q-2025-04-15-1708.pdf)
- **QuantumToolbox.jl (2025)** added `dfd_mesolve`, a callback that grows or shrinks the cutoff during evolution. The paper says that with it, "it is no longer required to check for the convergence of the results as a function of a fixed cutoff." It also added `dsf_mesolve` (dynamical shifted Fock) for strongly driven systems whose states "occupy a significant portion of the Hilbert space." — [arXiv 2504.21440](https://arxiv.org/abs/2504.21440)
- **Dynamiqs benchmark notes (2026)** — [benchmarks README](https://github.com/dynamiqs/dynamiqs/blob/main/benchmarks/README.md) and `systems.py` in the same folder
  - The driven Kerr case is "kept at moderate n because the Kerr spectrum grows as n², which bounds the explicit-RK step size." Truncation therefore causes stiffness.
  - The vectorized n²×n² Liouvillian "sets a memory ceiling."
  - The Liouvillian propagator scales as O(n⁶).
- SQcircuit [#40](https://github.com/stanfordLINQS/SQcircuit/issues/40) (2026-08, open): `remove_dependent_columns()` silently drops independent columns. It may be the root cause of [#32](https://github.com/stanfordLINQS/SQcircuit/issues/32) and [#39](https://github.com/stanfordLINQS/SQcircuit/issues/39). — [issue](https://github.com/stanfordLINQS/SQcircuit/issues/40)
- Weak signal: QuTiP issue titles containing "truncation" or "cutoff" return 0 hits. Truncation trouble shows up as usage questions in forums, not as bug reports.

**2.2 Mixed-dimension and qudit support**

- **Qiskit has no native qudits.**
  - The only qudit-specific issue is #12068 (2024, open) — [issue](https://github.com/Qiskit/qiskit/issues/12068)
  - A QAMP mentorship issue says Qiskit has no proper means for qudit circuits (2022, older) — [QAMP #39](https://github.com/qiskit-advocate/qamp-fall-22/issues/39)
  - Four Stack Exchange questions, 2020–2024:
    - [QC.SE 17666](https://quantumcomputing.stackexchange.com/questions/17666) (1,101 views).
    - [QC.SE 37561](https://quantumcomputing.stackexchange.com/questions/37561) (2024, "How to create a qudit of any dimension d in qiskit?").
    - [QC.SE 9905](https://quantumcomputing.stackexchange.com/questions/9905).
    - [QC.SE 34112](https://quantumcomputing.stackexchange.com/questions/34112).
- **Cirq #7461 (2025-07)** proposed deprecating qudits because they get in the way of performance work. It was closed after Infleqtion objected — [issue](https://github.com/quantumlib/Cirq/issues/7461)
  - Infleqtion: "We get a lot of mileage out of qudit support… to model 2+ energy levels of atomic systems in order to simulate leakage errors."
  - Infleqtion wrote its own qudit Kraus channel.
  - The maintainer noted that Cirq's Y gate lacks qudit support.
  - Qudit bugs were still being fixed in 2026 — [PR #8138](https://github.com/quantumlib/Cirq/pull/8138), [PR #8318](https://github.com/quantumlib/Cirq/pull/8318)
- **PennyLane**
  - Open requests and bugs:
    - #2190 "Provide a native way to simulate qudit" (2022, open) — [issue](https://github.com/PennyLaneAI/pennylane/issues/2190)
    - #5223, mixed qubit/qutrit/qudit interaction (2024, open) — [issue](https://github.com/PennyLaneAI/pennylane/issues/5223)
    - #4367, `qml.matrix` gives the wrong size on qutrit circuits (2023, open) — [issue](https://github.com/PennyLaneAI/pennylane/issues/4367)
  - Forum (2026-03): staff said there are "no current plans to support GPU acceleration of qutrit simulations." — [forum](https://discuss.pennylane.ai/t/qutrits-simulator-doubts/9286)
  - All qudit functionality was then removed (PR #9867, 2026-07-24). The changelog says "Qudit functionality in pennylane.labs still remains." — [changelog-dev](https://github.com/PennyLaneAI/pennylane/blob/master/doc/releases/changelog-dev.md)
- **Open-system tools: dims and tensor-structure errors**
  - QuTiP #2386 (2024) — [issue](https://github.com/qutip/qutip/issues/2386)
    - The excitation-number-restricted basis broke steady-state solvers in 5.0.0 through a dims mismatch.
    - After the fix, results differed from QuTiP 4.7.
  - More QuTiP dims bugs: [#2595](https://github.com/qutip/qutip/issues/2595) and [#2612](https://github.com/qutip/qutip/issues/2612) (2025); older segfaults on incompatible dims, [#1456](https://github.com/qutip/qutip/issues/1456) and [#1782](https://github.com/qutip/qutip/issues/1782).
  - Dynamiqs [#1037](https://github.com/dynamiqs/dynamiqs/issues/1037) (2025-12): dims with `sprepost`/`vectorize`.
  - QuantumToolbox.jl #447 (2025-04): cannot build an operator with generalized dims of different length — [issue](https://github.com/qutip/QuantumToolbox.jl/issues/447)
  - QuantumOptics.jl [#377](https://github.com/qojulia/QuantumOptics.jl/issues/377) (2023-12): `schroedinger()` and `master()` behave differently on sum bases.
  - scqubits [#256](https://github.com/scqubits/scqubits/issues/256) (2024-12): `identity_wrap` is wrong for callables.
- Yao.jl added qudits to YaoToEinsum (PR #530, 2024-11) — [PR](https://github.com/QuantumBFS/Yao.jl/pull/530)

**2.3 Time-dependent Hamiltonians, discontinuous drives and event handling**

This is the most recurrent failure in the open-system tools.

- **QuTiP v5 docs, "Working with pulses"** — [docs](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-time.html)
  - Variable-step methods can miss thin pulses. `max_step` "should be set to under half the pulse width to be certain they are not missed."
  - The docs' own example misses a pulse on t ∈ [0.7, 0.75].
- **Independent QuTiP reports of skipped or mis-resolved pulses, 2020–2026 (six or more reporters)**
  - #1265 (2020, older): a π/2–wait–π/2 sequence gives Pe = 0.5 instead of 1 once the separation exceeds about 23 — [issue](https://github.com/qutip/qutip/issues/1265)
  - #1945 (2022, older), a drive that is zero before t = 100 is skipped — [issue](https://github.com/qutip/qutip/issues/1945)
    - Maintainer: "adaptive solver skips the real part of the pulse. You need to set max_step."
    - The user proposed an automatic max_step = span/20.
  - #2552 (2024-11), a Ramsey simulation shows discontinuities — [issue](https://github.com/qutip/qutip/issues/2552)
    - Maintainer: "a pulse is too thin and missed… always set max_step."
  - #2355 (2024-03): an ultrafast pulse leaves a ladder system unexcited; unanswered — [issue](https://github.com/qutip/qutip/issues/2355)
  - #2814 (2026-01-30, open; verified), "Detect pulses in mesolve to adjust time steps of the Integrator" — [issue](https://github.com/qutip/qutip/issues/2814)
    - Pulses after long delays are skipped.
    - The user wants automatic "time barriers" so the solver need not take nanosecond steps over the whole wait.
  - Title only:
    - [#2813](https://github.com/qutip/qutip/issues/2813) (2026-01): "mevolve sometimes ignores hamiltonian time coefficients if last coefficient is 0."
    - [Discussion #2653](https://github.com/qutip/qutip/discussions/2653) (2025-03): Hahn-echo T2 comes out too short.
- **Interpolated (array) coefficients behave inconsistently** — [#2063](https://github.com/qutip/qutip/issues/2063), [#2253](https://github.com/qutip/qutip/issues/2253)
  - QuTiP [#2063](https://github.com/qutip/qutip/issues/2063) (2023): a QobjEvo interpolates differently when called directly than inside a solver.
  - [#2253](https://github.com/qutip/qutip/issues/2253) (2023-10): `[ops, func]` and `[ops, ndarray]` give different results.
- **Cython compilation of string coefficients: six reports**
  - [#954](https://github.com/qutip/qutip/issues/954) (2019, still open): runtime compilation on Windows.
  - #2162 (2023): without `filelock`, QuTiP silently falls back to uncompiled `eval` — [issue](https://github.com/qutip/qutip/issues/2162)
  - #2639 (2025-02), the string `'np.heaviside(t,1)'` fails because Cython types it as complex — [issue](https://github.com/qutip/qutip/issues/2639)
    - A second user (2025-11) saw it on Windows/AMD but not macOS/Intel with the same notebook.
  - [#2770](https://github.com/qutip/qutip/issues/2770) (2025-10): compilation fails in an editable install.
- **Dynamiqs**
  - #1085 (2026-05): the adaptive Rouchon3 method crashes whenever H has discontinuities, caused by a Diffrax API change — [issue](https://github.com/dynamiqs/dynamiqs/issues/1085)
  - Its benchmark suite includes `sesolve_pwc`, an "adaptive stepper crossing piecewise-constant discontinuities," and reports `nrej` (rejected steps) as the cost — [README](https://github.com/dynamiqs/dynamiqs/blob/main/benchmarks/README.md)
- **QuantumToolbox.jl**
  - #597 (2025-11): `mesolve` fails when a collapse operator has the form A + f(t)B — [issue](https://github.com/qutip/QuantumToolbox.jl/issues/597)
  - #504 (2025-07): times and states mismatch in `sesolve` with callbacks — [issue](https://github.com/qutip/QuantumToolbox.jl/issues/504)
  - On the other hand, its 2025 paper uses solver callbacks, triggered "by a condition that is verified along the dynamics," for adaptive truncation — [arXiv 2504.21440](https://arxiv.org/abs/2504.21440)
- **The IBM stack:** with Qiskit Pulse removed and qiskit-dynamics archived (1D), IBM appears to have no maintained pulse-level or time-dependent Hamiltonian simulator. This is an inference; Qiskit Aer's pulse support was not separately checked.
  - Qiskit Experiments 0.13 removed its Dynamics-based PulseBackend because it "could no longer be maintained" — [release notes](https://qiskit-community.github.io/qiskit-experiments/release_notes.html)
  - PennyLane #6100 (2024-08): `qml.evolve` with a Rydberg pulse under JAX gives probabilities summing to more than 1 — [issue](https://github.com/PennyLaneAI/pennylane/issues/6100)

**2.4 Stiff systems, tolerances and solver failures**

- QuTiP "ODE integration error: Try to increase the allowed number of substeps by increasing the nsteps parameter." Four reports, all older (2021–22): [#1605](https://github.com/qutip/qutip/issues/1605), [#1896](https://github.com/qutip/qutip/issues/1896), and two Google Group threads ([1](https://groups.google.com/g/qutip/c/lfYrgaxbOsQ), [2](https://groups.google.com/g/qutip/c/RPzVHvJ-R38)).
  - [#1896](https://github.com/qutip/qutip/issues/1896) occurred in tutorial code.
  - No titled 2023+ GitHub issue was found. That may mean the error became rarer after QuTiP 5 (inference).
- **Default tolerances fail silently on large frequencies (five reporters)**
  - #2560 (2024-11, verified): with H = diag(0, 5000), P(0) drifts. Maintainer: "You need to change the options rtol." — [issue](https://github.com/qutip/qutip/issues/2560)
  - #2258 / Discussion #2255 (2023-11): a 5 GHz oscillator drifts; fixed with `max_step=1/(wr*100), atol=rtol=1e-9` — [issue](https://github.com/qutip/qutip/issues/2258)
  - #2051 (2022, older): "decay" in purely unitary evolution — [issue](https://github.com/qutip/qutip/issues/2051)
  - #2549 (2024-10): transmon results depend on the number of tlist points; the cause was output aliasing — [issue](https://github.com/qutip/qutip/issues/2549)
- QuantumOptics.jl #516 (2023-12, open): the default DP5 solver and tolerances are undocumented, and the reporter calls DiffEq's defaults (reltol 1e-3, abstol 1e-6) too loose. Maintainer: "This is indeed poorly documented." — [issue](https://github.com/qojulia/QuantumOptics.jl/issues/516)
- Dynamiqs #1168 (2026-09-30) — [issue](https://github.com/dynamiqs/dynamiqs/issues/1168)
  - The stiff Kvaerno3/5 methods hit the maximum step count on a constant two-level `mesolve`, even with 10⁶ steps.
  - It works with `assume_hermitian=False` or `vectorized=True`.
  - The report says it was drafted with AI help; the reporter reproduced the failure.
- QuantumToolbox.jl — [#444](https://github.com/qutip/QuantumToolbox.jl/issues/444), [#617](https://github.com/qutip/QuantumToolbox.jl/issues/617)
  - [#444](https://github.com/qutip/QuantumToolbox.jl/issues/444) (2025-04): `mesolve` is very slow or fails where the equivalent QuTiP script works.
  - [#617](https://github.com/qutip/QuantumToolbox.jl/issues/617) (2025-12): the vectorized sparse-Liouvillian right-hand side is slow; the maintainers plan a non-vectorized master equation.
- Julia Discourse: a QuantumOptics.jl Liouvillian is slow at 4096×4096 — [thread](https://discourse.julialang.org/t/liouvillian-superoperator-way-cannot-run-by-timeevolution-master-dynamic-of-quantumoptics-jl/113650)
- Stack Exchange (older): "Fastest numerical method to solve Lindblad Master Equation?" (2020, score 10, 3,647 views) — [Phys.SE 593284](https://physics.stackexchange.com/questions/593284)
  - "Are there simple ways to numerically solve the time-dependent Schrödinger equation?" (2014, score 43, 11.5k views) — [SciComp.SE 10876](https://scicomp.stackexchange.com/questions/10876)

**2.5 Floating-point precision, positivity, underflow and steady states**

*Open-system tools*

- QuTiP #3024 (2026-09-30, open) — [issue](https://github.com/qutip/qutip/issues/3024)
  - `fidelity(rho, rho)` = 1.0000000102 for a rank-2 density matrix.
  - Diagnosis: the eigenvalue filter compares to 0 instead of a tolerance.
- QuTiP #3026 (2026-09-30, open; verified title; same reporter as #3024) — [issue](https://github.com/qutip/qutip/issues/3026)
  - `steadystate(method="eigen")` "silently returns a completely wrong, non-physical state" for a Jaynes–Cummings model with N = 4.
  - It works with `sparse=False`.
- Older negative-eigenvalue reports — [#1919](https://github.com/qutip/qutip/issues/1919), [#646](https://github.com/qutip/qutip/issues/646):
  - [#1919](https://github.com/qutip/qutip/issues/1919) (2022): PIQS (permutation-invariant module) Dicke entropy is −Inf.
  - [#646](https://github.com/qutip/qutip/issues/646) (2017): negative eigenvalues give negative entropy.
- QuTiP #2175 (2023-06, open), non-unique steady states — [issue](https://github.com/qutip/qutip/issues/2175)
  - N. Lambert: depending on method it "can either fail, return one of the possibilities, or some linear combination… The default one (direct) tends to fail."
- More QuTiP steady-state failures:
  - [#2649](https://github.com/qutip/qutip/issues/2649) (2025-03): the default direct solver fails with a singular matrix on `piqs.Dicke`.
  - [#2747](https://github.com/qutip/qutip/issues/2747) (2025-09): out of memory with non-CSR sparse formats.
  - #2999 (2026-09, open): gmres "Tolerance was not reached" recurs in QuTiP's own CI — [issue](https://github.com/qutip/qutip/issues/2999)
- QuantumOptics.jl #428 (2024-12, open), quartic oscillator — [issue](https://github.com/qojulia/QuantumOptics.jl/issues/428)
  - `steadystate.iterative` gives spurious results that depend on the random initial guess.
  - The `master` and `eigenvector` methods agree with each other.
- Bloch–Redfield does not guarantee positivity — [arXiv 2402.06354](https://arxiv.org/abs/2402.06354); [docs](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-bloch-redfield.html)
  - The QuTiP docs warn of negativity in the non-secular form (snippet).
- HierarchicalEOM.jl #298 (2026-06), confirmed and fixed in v2.15.1 — [issue](https://github.com/qutip/HierarchicalEOM.jl/issues/298)
  - The sparsity pattern changed abruptly at large dimension.
  - A basis-independent quantity jumped at the same point.
- SQcircuit [#20](https://github.com/stanfordLINQS/SQcircuit/issues/20) (2023): limited precision of physical constants.

*Circuit tools*

- **PennyLane "probabilities do not sum to 1": four independent GitHub bugs (2023–2026)**
  - [#3894](https://github.com/PennyLaneAI/pennylane/issues/3894) and [#5444](https://github.com/PennyLaneAI/pennylane/issues/5444).
  - [#6100](https://github.com/PennyLaneAI/pennylane/issues/6100): pulses under JAX.
  - [#9000](https://github.com/PennyLaneAI/pennylane/issues/9000) (2026-01): 32-bit JAX exceeds a "hardcoded cutoff value of 1e-7."
  - Fix PR [#7076](https://github.com/PennyLaneAI/pennylane/pull/7076) forces normalization.
  - Older forum threads from 2020–22, including Strawberry Fields Gaussian boson sampling: [1](https://discuss.pennylane.ai/t/numpy-says-probabilities-do-not-sum-to-1/487), [2](https://discuss.pennylane.ai/t/gbs-for-molecular-vibronic-spectra-probabilities-do-not-sum-to-1/1307), [3](https://discuss.pennylane.ai/t/probabilities-of-a-circuit-arent-added-to-one/2153)
- **Cirq simulates in complex64 by default; normalization failures in three or more reports**
  - [#6402](https://github.com/quantumlib/Cirq/issues/6402) (2024). The maintainers renormalized the output: "there will always be numbers that can't be represented exactly leading to an accumulation of errors."
  - [#6228](https://github.com/quantumlib/Cirq/issues/6228) (2023) and [#5268](https://github.com/quantumlib/Cirq/issues/5268) (2022).
  - PR [#6672](https://github.com/quantumlib/Cirq/pull/6672) uses complex128 for cross-entropy benchmarking (XEB) circuits.
  - PR [#8298](https://github.com/quantumlib/Cirq/pull/8298) (2026) removes the renormalization again.
- **qsim** advertises single precision — [README](https://github.com/quantumlib/qsim)
  - [#754](https://github.com/quantumlib/qsim/issues/754): NumPy 2 breaks expectation tests (2025).
  - [PR #1014](https://github.com/quantumlib/qsim/pull/1014): dtype changes to fix numerical comparison failures (2026).
  - Older: [PR #422](https://github.com/quantumlib/qsim/pull/422) sets flush-to-zero and denormals-are-zero, trading subnormal underflow for speed (2021).
- **Qiskit Aer**
  - [#2462](https://github.com/Qiskit/qiskit-aer/issues/2462) (2026-09, open): single precision with `save_statevector` crashes from 14 qubits up.
  - [#1870](https://github.com/Qiskit/qiskit-aer/issues/1870): the `precision=` option does not behave as expected.
  - [#2340](https://github.com/Qiskit/qiskit-aer/issues/2340): multi-GPU sampling returns unexpected results.
  - [#2410](https://github.com/Qiskit/qiskit-aer/issues/2410): matrix-product-state `save_amplitudes` is wrong above 64 qubits.
  - [#2276](https://github.com/Qiskit/qiskit-aer/issues/2276): for rx(π⁵⁰), consecutive floats are 2³⁰ apart.
- Older circuit reports:
  - Qibo [#517](https://github.com/qiboteam/qibo/issues/517): probabilities stop summing to 1 at 26 or more qubits in single precision.
  - Azure Quantum [#287](https://github.com/microsoft/azure-quantum-python/issues/287).
- PennyLane #10185 (2026-09, open): gate fusion produces a NaN Euler angle on `H; H` — [issue](https://github.com/PennyLaneAI/pennylane/issues/10185)

**2.6 Large parameter sweeps, batching and memory**

*Open-system tools*

- QuTiP #2819 (2026-02-10, open; verified title), repeated `mesolve` in a loop leaks memory — [issue](https://github.com/qutip/qutip/issues/2819)
  - Only `gc.collect()` keeps usage flat.
  - The maintainer's workaround example uses Fock cutoff 120 and a coherent state with α = 8.
- QuTiP #3010 (2026-09, open): `parallel=True` in stochastic solvers fails with lambda coefficients because multiprocessing must pickle them — [issue](https://github.com/qutip/qutip/issues/3010)
  - Older: [#1202](https://github.com/qutip/qutip/issues/1202) (2020), spawn-based multiprocessing freezes `mcsolve`.
- scqubits — [#318](https://github.com/scqubits/scqubits/issues/318)
  - [#318](https://github.com/scqubits/scqubits/issues/318) (2026-06, open): `ParameterSweep.energy_by_dressed_index(subtract_ground=True)` mutates cached eigenvalues, silently corrupting later lookups.
  - [#241](https://github.com/scqubits/scqubits/issues/241) (2024): GPU diagonalisation (CuPy) is not optimised for sweeps.
- GPU memory — [Dynamiqs issues](https://github.com/dynamiqs/dynamiqs/issues)
  - QuantumToolbox.jl [#533](https://github.com/qutip/QuantumToolbox.jl/issues/533) (2025-08): CUDA allocation grows with integration time.
  - Dynamiqs [#1132](https://github.com/dynamiqs/dynamiqs/issues/1132)/#1133 (2026-07): the LowRank method densifies.
  - Dynamiqs [#961](https://github.com/dynamiqs/dynamiqs/issues/961) (2025): `asqarray` densifies sparse arrays.

*Circuit tools*

- Qiskit parameter binding:
  - [#7676](https://github.com/Qiskit/qiskit/issues/7676) (2022, older): "takes too long time."
  - Still being optimised in [PR #14987](https://github.com/Qiskit/qiskit/pull/14987) (2025).
- Aer [#2438](https://github.com/Qiskit/qiskit-aer/issues/2438) (2026-05): EstimatorV2 has an O(N²) loop and a fork-unsafe executor.
- PennyLane Lightning [#1195](https://github.com/PennyLaneAI/pennylane-lightning/issues/1195) (2025-07): batched execution on `lightning.gpu` is slower than a sequential loop.
- PennyLane forum threads on parallelizing over inputs or parameters (2023–25): [1](https://discuss.pennylane.ai/t/is-there-a-way-to-parallelize-the-same-circuit-for-multiple-input-data/2516), [2](https://discuss.pennylane.ai/t/parallel-circuit-execution-during-optimisation/2998)

*Overhead studies*

- VQCSim (2026) and Malarchick (2026) both measure per-call and dispatch overhead, each by comparing their own code with the frameworks. The same overhead is the stated motivation for PQLS and for Dynamiqs batching (see section 1).

**2.7 Non-Markovian dynamics**

- QuTiP #2950 (2026-07-25, open), ESPRIT fails to fit long-memory bath correlation functions (a giant atom with delay) at any number of modes — [issue](https://github.com/qutip/qutip/issues/2950)
  - QuTiP developer: "the bug is on me"; the implementation does not match the algorithm.
- OQuPy #150 (2025-06, open), the default `scipy.integrate.quad` fails on the oscillatory integrand `CustomSD.eta_function` at long memory times — [issue](https://github.com/tempoCollaboration/OQuPy/issues/150)
  - Usually it warns; in niche cases it fails silently.
- HierarchicalEOM.jl [#298](https://github.com/qutip/HierarchicalEOM.jl/issues/298): see 2.5.
- QuTiP [#2580](https://github.com/qutip/qutip/issues/2580) (2024-12, open), title only: "Environments should work with more solvers."
- No direct "HEOM out of memory" user report was found. Memory-scaling concerns come from design papers only.

**2.8 Measurement feedback and feedforward**

- QuTiP #1571 (2021-06, still open): "Iterable access to solver results and possibility of feedback to solvers" — [issue](https://github.com/qutip/qutip/issues/1571)
  - Google Group: a user wants to update H on the fly from ⟨a⟩ during `smesolve` — [thread](https://groups.google.com/d/topic/qutip/WrwLh_iLKo8)
- Stochastic-solver correctness and performance problems:
  - QuTiP #2733 (2025-08): `smesolve` is more than 10× slower in 5.2.0 than in 5.1.1 — [issue](https://github.com/qutip/qutip/issues/2733)
  - Title only, all from one reporter, so one source:
    - Dynamiqs [#1092](https://github.com/dynamiqs/dynamiqs/issues/1092): `dssesolve` is wrong for some `nsteps_per_save`.
    - Dynamiqs [#1113](https://github.com/dynamiqs/dynamiqs/issues/1113): biased click statistics.
    - Dynamiqs [#1155](https://github.com/dynamiqs/dynamiqs/issues/1155): possible Rouchon1 bug.
    - Dynamiqs [#1162](https://github.com/dynamiqs/dynamiqs/issues/1162): `dsmesolve` fails under JIT.
- **PennyLane mid-circuit-measurement methods (docs v0.45)** — [docs](https://docs.pennylane.ai/en/stable/introduction/dynamic_quantum_circuits.html)

  | Method | Memory | Trade-off |
  |---|---|---|
  | Deferred measurement | O(2^n_MCM) | "requires an additional qubit for each mid-circuit measurement" |
  | One-shot | O(1) | no analytic mode; finite-difference gradients only |
  | Tree-traversal | — | only on some devices |

  - User reports:
    - "Running circuits with many mid-circuit measurements" (2024): a developer called deferred measurement "a known scaling issue" — [forum](https://discuss.pennylane.ai/t/running-circuits-with-many-mid-circuit-measurements/4355)
    - "Mid-circuit measurement: slow and memory explosion" (2024-11) — [forum](https://discuss.pennylane.ai/t/mid-circuit-measurement-slow-and-memory-explosion/7452)
    - [#5443](https://github.com/PennyLaneAI/pennylane/issues/5443): mid-circuit measurements break broadcasting.
    - [#6541](https://github.com/PennyLaneAI/pennylane/issues/6541): tree-traversal fails under autograd.
- **Stim** — [README](https://github.com/quantumlib/Stim); [#1017](https://github.com/quantumlib/Stim/issues/1017); [gates.md](https://github.com/quantumlib/Stim/blob/main/doc/gates.md)
  - Feedback is limited to classically controlled Paulis on measurement records, such as `CX rec[-3] 0`.
  - "There is no support for non-Clifford operations, such as T gates and Toffoli gates." The maintainer repeated this in 2025-12.
  - Related Stack Exchange question: "Simulating flag qubits and conditional branches using Stim" (2021, older, 1,005 views) — [QC.SE 22281](https://quantumcomputing.stackexchange.com/questions/22281)
  - The stimcirq bridge gained feedback gates only in Apr 2026 — [#1053](https://github.com/quantumlib/Stim/issues/1053)
- Qiskit #14968 (2025-08, open): "Poor routing in transpiler for dynamic circuits" — [issue](https://github.com/Qiskit/qiskit/issues/14968)

**2.9 Automatic-differentiation gaps**

- qutip-jax — [#95](https://github.com/qutip/qutip-jax/issues/95)
  - [#95](https://github.com/qutip/qutip-jax/issues/95) (2025-08, open): the second derivative in the QuTiP 5 paper's own JAX tutorial fails since JAX 0.7.0.
  - [#62](https://github.com/qutip/qutip-jax/issues/62): not jit-compatible.
  - [#93](https://github.com/qutip/qutip-jax/issues/93): vectorisation breaks linear combinations of states.
- QuantumToolbox.jl #547 (2025-09, open): forward-mode autodiff fails for QobjEvo with time-dependent coefficients, because a ComplexF64 ScalarOperator cannot hold a Dual — [issue](https://github.com/qutip/QuantumToolbox.jl/issues/547)
  - The paper itself calls its autodiff "experimental."
- QuantumOptics.jl — [#376](https://github.com/qojulia/QuantumOptics.jl/issues/376)
  - [#357](https://github.com/qojulia/QuantumOptics.jl/issues/357) (2023, open): add ForwardDiff support to all solvers.
  - [#376](https://github.com/qojulia/QuantumOptics.jl/issues/376) (2023-12, open): ForwardDiff fails with TimeDependentOperator.
  - [#432](https://github.com/qojulia/QuantumOptics.jl/issues/432) (2025): breakage after the upstream `promote_dual` removal.
- Dynamiqs Hessians are a requested feature: [#1021](https://github.com/dynamiqs/dynamiqs/issues/1021) (2025) and a unitaryHACK bounty, [#1082](https://github.com/dynamiqs/dynamiqs/issues/1082) (2026).
- The Xanadu differentiable continuous-variable stack (Strawberry Fields TensorFlow backend, MrMustard, The Walrus) is archived (1D).
- Qiskit — [qiskit-algorithms #216](https://github.com/qiskit-community/qiskit-algorithms/issues/216)
  - Gradients live in the unsupported qiskit-algorithms package.
  - [#216](https://github.com/qiskit-community/qiskit-algorithms/issues/216) (2025): controlled-rotation gradients fail.
  - [#137](https://github.com/qiskit-community/qiskit-algorithms/issues/137): no V2 gradient primitives.
- PennyLane — [interfaces.rst](https://github.com/PennyLaneAI/pennylane/blob/master/doc/introduction/interfaces.rst)
  - The adjoint method "reverses through the circuit… applying the inverse (adjoint) gate," so it needs reversible unitary circuits.
  - Mid-circuit-measurement gradients are finite-difference only (2.8).
  - Lightning [#808](https://github.com/PennyLaneAI/pennylane-lightning/issues/808): `qml.state` is unsupported with adjoint.
  - Catalyst [#384](https://github.com/PennyLaneAI/catalyst/issues/384): differentiating `QubitUnitary` core-dumps.
- quimb — [#355](https://github.com/jcmgray/quimb/issues/355)
  - [#355](https://github.com/jcmgray/quimb/issues/355) (2026-04, open): NaN in gradient.
  - [#147](https://github.com/jcmgray/quimb/issues/147) (2022, older): JAX autodiff incompatible with SVD.

**2.10 Reproducibility and version drift**

*Circuit tools*

- Qiskit #16490 (2026-06-26, fixed 2026-08-13), "seed-independent" non-determinism in transpilation — [issue](https://github.com/Qiskit/qiskit/issues/16490); [fix PR](https://github.com/Qiskit/qiskit/pull/16732)
  - A fixed circuit, backend and transpiler seed gave run-to-run differences in gate count and depth.
  - The fix sums VF2Layout errors "in a deterministic order."
- Qiskit [#14514](https://github.com/Qiskit/qiskit/issues/14514) (2025): the circuit unitary changes across optimization levels.
  - Older: [#7770](https://github.com/Qiskit/qiskit/issues/7770), equivalent circuits give different output distributions.
- Seed handling:
  - Older: Aer [#585](https://github.com/Qiskit/qiskit-aer/issues/585), `seed_simulator` does not apply to a list of circuits.
  - Stack Exchange: QAOA optimization fails once a seed is set in AerSampler (2023) — [QC.SE 33616](https://quantumcomputing.stackexchange.com/questions/33616)

*Open-system tools: dependency drift*

- scqubits [#320](https://github.com/scqubits/scqubits/issues/320) and [#327](https://github.com/scqubits/scqubits/issues/327) (2026, two reporters): QuTiP 5.3 breaks `HilbertSpace`.
- SQcircuit [#23](https://github.com/stanfordLINQS/SQcircuit/issues/23) (2024): incompatible with QuTiP 5.
- Dynamiqs:
  - [#1011](https://github.com/dynamiqs/dynamiqs/issues/1011) (2025): incompatible with QuTiP 5.2.
  - [#996](https://github.com/dynamiqs/dynamiqs/issues/996) (2025) and [#1140](https://github.com/dynamiqs/dynamiqs/issues/1140) (2026): new JAX versions break it.
- QuantumOptics.jl [#432](https://github.com/qojulia/QuantumOptics.jl/issues/432) and [#441](https://github.com/qojulia/QuantumOptics.jl/issues/441): upstream changes broke autodiff and precompilation.
- QuTiP version and platform effects:
  - [#2386](https://github.com/qutip/qutip/issues/2386): results differ from QuTiP 4.7.
  - [#2639](https://github.com/qutip/qutip/issues/2639): platform-dependent coefficient failure.
  - [#2300](https://github.com/qutip/qutip/issues/2300) and [#2316](https://github.com/qutip/qutip/issues/2316) (2024): QuTiP 4.7 breaks with SciPy 1.12.
  - [#2733](https://github.com/qutip/qutip/issues/2733): performance regression (2.8).
- No QuTiP seed-reproducibility bug reports were found for `mcsolve`. This is a weak negative signal.

**2.11 Where the signal lives**

- A Stack Exchange API sample (2026-10-04) shows that QC.SE and Phys.SE carry little of this traffic — [QC.SE 5000](https://quantumcomputing.stackexchange.com/questions/5000); [QC.SE 6418](https://quantumcomputing.stackexchange.com/questions/6418)
  - The top-voted QuTiP questions are about usage, such as plotting Bloch spheres and measuring qubits, and mostly date from 2017–2021.
  - Numerics threads are few and often unanswered, e.g. QC.SE 34202.
  - Most of the evidence above comes from GitHub issues and project forums.

### Inferences
- Discontinuity- and event-aware time integration is the clearest unmet need in open-system tools. Evidence:
  - The same manual `max_step` workaround has been given for six years.
  - QuTiP's documentation example itself misses a pulse.
  - An open 2026 request asks for automatic pulse detection.
  - Dynamiqs benchmarks count rejected steps at discontinuities.
  - Every incumbent leaves pulse alignment to the user.
- Truncation convergence is a user-owned, poorly tooled task. Automated aids are either new (QuantumToolbox.jl's dynamical Fock dimension in 2025; scqubits diagnostics under review in 2026) or gone (the archived MrMustard autocutoff). Worked examples that show a convergence check change the answer would address a real, documented failure: "8 months" lost in scqubits [#274](https://github.com/scqubits/scqubits/issues/274).
- Silent wrong answers cluster where conditioning matters:
  - eigen- and square-root-based functionals of nearly singular density matrices;
  - steady-state solvers;
  - fitting bath correlation functions;
  - single-precision state vectors checked against fixed absolute tolerances.

  Tolerance control relative to the problem's own scale, or extended precision, would change the result in these regimes. This is inferred from the issue diagnoses, not stated by maintainers.
- In circuit SDKs, multilevel physics (leakage, qutrits) is being dropped. Qiskit never had it, PennyLane removed it in 2026, and Cirq keeps it for one industrial user's leakage modelling. Mixed-dimension work is pushed back to open-system tools.
- Reproducibility failures come mostly from fast-moving dependency stacks (JAX, Diffrax, SciPy, QuTiP minor versions) and order-dependent floating-point sums, not from seeding.

### Gaps
- Stack Exchange coverage is thin. Fork searches of Phys.SE, QC.SE and SciComp.SE found no 2023+ threads on solver stiffness, negative eigenvalues or HEOM memory.
- The QuTiP 5 seed API for `mcsolve`/`smesolve` and any concrete seed-reproducibility complaint were not verified.
- No direct user reports of HEOM memory exhaustion were found.
- Several 2026 Dynamiqs stochastic-solver bugs and QuTiP [#2813](https://github.com/qutip/qutip/issues/2813) are title-only.
- It is not established when the Strawberry Fields TensorFlow backend was deprecated, or whether the scqubits convergence module has shipped.
- The Qiskit Slack and other non-public channels could not be sampled.
- The Bloch–Redfield positivity warning and the Strawberry Fields `trace()` advice come from search snippets of the docs, though forum answers corroborate the latter.

## 3. Which numerical problems are standard in master's-level QM, quantum optics and quantum information teaching, and which do instructors use to show numerics are necessary?

### Takeaway
Five problem families recur across three or more independent courses or texts:

- Driven two-level dynamics: Rabi, Ramsey, detuning and the breakdown of the rotating-wave approximation (RWA).
- The Lindblad master equation for a damped or driven qubit or cavity, compared with quantum-jump trajectories.
- Wigner and Husimi phase-space pictures of Fock, cat and squeezed states.
- Grid solutions of the Schrödinger equation: Numerov or shooting for eigenproblems, split-operator or Crank–Nicolson for time evolution.
- Exact diagonalization of spin chains.

Textbook circuit algorithms dominate quantum-information courses, where simulation mostly illustrates cases that can be done by hand.

Problems where numerics clearly decide the answer:

- RWA breakdown and ultrastrong coupling.
- Driven-dissipative steady states and spectra (Mollow triplet).
- Finite-time Landau–Zener sweeps.
- Anharmonic potentials.
- Surface-code thresholds by Monte Carlo.

### Cited Findings

**QuTiP-based course material**

- **QuTiP "in Education" page** lists courses that use QuTiP — [page](https://qutip.org/education)
  - Courses listed:
    - MIT HSSP 2024. To our knowledge, not stated on the page, HSSP is a high-school programme, so this entry is not master's level.
    - Yale CHEM 584.
    - Caltech "Applied Physics: Quantum Electronics" (2019 syllabus).
    - U. Tokyo 量子技術序論.
    - ETH Zurich "Classical and Quantum Parametric Phenomena" (Winter 2024).
  - It describes QuTiP as "a user-friendly Python library for learning and teaching quantum mechanics through simulations."
- **J.R. Johansson, qutip-lectures** (about 2013, maintained for QuTiP 5; reportedly written for an invited Chalmers course) — [repo](https://github.com/jrjohansson/qutip-lectures)
  - Jaynes–Cummings model, including ultrastrong coupling.
  - Cavity–qubit gates and single-atom lasing.
  - Dicke model.
  - Correlation functions.
  - Parametric amplifier.
  - Monte Carlo trajectories.
  - iSWAP gate.
  - Adiabatic quantum computing.
  - Squeezed states.
  - Dispersive circuit QED.
  - Charge qubits.
  - Kerr nonlinearities.
  - Wigner-function gallery.
- **QuTiP 5 tutorials (current)** — [tutorials](https://qutip.org/qutip-tutorials/)
  - 26 time-evolution notebooks: Schrödinger, master equation, Monte Carlo, Bloch–Redfield, Floquet, stochastic, steady state.
  - 8 optimal-control notebooks.
  - Iterative maximum-likelihood tomography.
  - 8 permutation-invariant notebooks: superradiance, spin squeezing, time crystals.
  - 10 HEOM notebooks: spin–bath, FMO complex, heat transport, fermions.
  - Pulse-level circuit simulation: 10-qubit QFT, randomized benchmarking, Deutsch–Jozsa.
  - Lectures 0–16, including resonance fluorescence.
- **Trinity College Dublin, MSc Quantum Science & Technology, "Open Quantum Systems"** (master's; year not stated) — [student repo 1](https://github.com/Ilumirnau/Open-Quantum-Systems); [repo 2](https://github.com/a24l/Open-Quantum-Systems)
  - Worksheet 1: time-dependent driving, the limits of RWA validity, XXY/Ising open chains.
  - Worksheet 2: thermal states, decoherence from Gaussian frequency noise, Kraus maps, purity.
  - Worksheet 3: driven and coupled oscillators, Husimi-Q of cat and Fock states.
  - Worksheet 4: Lindblad dephasing, quantum-jump trajectories, a cavity with a nonlinear medium.
- **"Numerical Methods for Quantum Optics and Open Quantum Systems" lectures** (instructor, institution and year not identified; recent) — [site](https://daniele-pybectn.github.io/Computational_quantum_optics_lectures/)
  - Python ODE numerics.
  - Liouville equation in phase space.
  - QuTiP closed and open dynamics.
  - Jaynes–Cummings ("proof of field quantization") and its emission spectrum.
  - Resonance fluorescence and the Mollow triplet.
  - Wigner functions.
- **UMass ECE 550/650 (Niffenegger, 2025; mixed undergraduate/graduate)** — [repo](https://github.com/UMassIonTrappers/Introduction-to-Quantum-Computing)
  - QuTiP labs: time-dependent Schrödinger equation and Rabi oscillations, RWA with detuning, Ramsey, noise and composite pulses, two-qubit gates.
  - Qiskit labs:
    - T1 and T2.
    - GHZ parity.
    - Quantum volume.
    - Deutsch–Jozsa and Bernstein–Vazirani.
    - Shor code.
    - Grover, QFT, QPE, Shor.
- **Campaioli, Cole & Hapuarachchi, tutorial on quantum master equations, *PRX Quantum* 5, 020202 (2024)** (graduate) — [arXiv 2303.16449](https://arxiv.org/abs/2303.16449); [code](https://github.com/frnq/qme)
  - Lindblad, Redfield, Floquet, Suzuki–Trotter and sparse solvers, with 34 code examples.
- **Dawes, QuTiP in undergraduate QM (2019, older)**: "from two-state systems to optical interaction with multilevel atoms" — [arXiv 1909.13651](https://arxiv.org/abs/1909.13651)

**Textbooks and lecture notes**

- **Giannozzi, "Numerical Methods in Quantum Mechanics"** (Univ. Udine Laurea Magistrale, i.e. master's) — [notes](https://www.fisica.uniud.it/~giannozz/Corsi/MQ/LectureNotes/mq.pdf)
  - Numerov on a grid for 1D and spherical potentials.
  - The full table of contents could not be retrieved.
- **Izaac & Wang, *Computational Quantum Mechanics* (Springer 2018, advanced undergraduate, older)** — [listing](https://www.abebooks.com/9783319999296/Computational-Quantum-Mechanics-Undergraduate-Lecture-331999929X/plp)
  - Covers precision, eigenvalue problems, Fourier methods, Schrödinger equations in 1D and higher dimensions, time propagation, central potentials, multi-electron systems.
  - Blurb: "Quantum mechanics undergraduate courses mostly focus on systems with known analytical solutions; the finite well, simple Harmonic, and spherical potentials. This textbook introduces the numerical techniques required…"
- **Steck, *Quantum and Atom Optics*** (graduate): a "Split-Operator Methods" chapter with "numerical tests for fourth and sixth order methods" (seen via the search index) — [PDF mirror](http://atomoptics.uoregon.edu/~tbrown/files/relevant_papers/quantum-optics-notes.pdf)
- **Newman, *Computational Physics*, Exercise 9.8**, as used in CUNY City Tech PHYS 4150 (undergraduate) — [course page](https://openlab.citytech.cuny.edu/phys4150/lab-exercises/partial-differential-equations/)
  - Crank–Nicolson for an electron in a 10⁻⁸ m box with N = 1000 slices.
  - Other Newman quantum exercises (shooting method, spectral method) were not verified.
- **ETH Zurich "Computational Quantum Physics" (Troyer, FS 2011; graduate, older)** — [cqp.pdf](https://edu.itp.phys.ethz.ch/fs11/cqp/cqp.pdf)
  - One-body and many-body Schrödinger equations.
  - Path-integral and quantum Monte Carlo.
  - Exact many-body methods.
- **"Quantumandu" school notes (ICTP / Tribhuvan, Jul 2024)**: exact-diagonalization code for spin chains plus tensor networks — [arXiv 2503.03564](https://arxiv.org/abs/2503.03564)

**Landau–Zener as a numerics-required teaching problem** (two course-level sources plus student theses)

- **Guttieres, Petrovic & Freericks, *Am. J. Phys.* (arXiv Jun 2023)**, "Computational projects with the Landau-Zener problem in the quantum mechanics classroom" — [arXiv 2306.11633](https://arxiv.org/abs/2306.11633)
  - Calls it "an excellent model quantum system for a computational project."
  - Teaches accuracy, discretization and extrapolation.
- **Washington State University Physics 555 "Quantum Technologies" (graduate, Fall 2024)** — [LZ notes](https://physics-555-quantum-technologies.readthedocs.io/en/fall2024/Notes/LandauZener.html)
  - Numerically integrates the time-dependent Schrödinger equation and compares with analytic predictions.
  - Because H(t) does not commute with itself at different times, "Solving the differential equation numerically is usually the easiest approach."
- **KTH student theses**, "Simulating the Landau-Zener Problem" and a second Landau–Zener thesis. These are bachelor-level; the year was not verified. — [DiVA](https://kth.diva-portal.org/smash/get/diva2:1879424/FULLTEXT01.pdf)

**Quantum-information and quantum-error-correction teaching**

- **NordIQuEst Quantum Autumn School 2024, surface-code threshold tutorial** — [notebook](https://nordiquest.net/application-library/training-material/qas2024/notebooks/surface_code_threshold.html)
  - Tools: Stim and PyMatching.
  - Repetition code at d = 3, 5, 7.
  - Rotated surface-code memory with p swept over 0.002–0.009 and 10,000 shots per point.
  - Computes the suppression factor Λ.
- **Kastoryano QEC course (Cologne WS 2018/19, older)**: thresholds are done analytically (counting bound about 0.037; toric-code optimum 0.113 via a stat-mech mapping), with no coding — [notes](https://www.thp.uni-koeln.de/kastoryano/ExSheets/Notes_v9.pdf)
- **Watrous / IBM Quantum Learning (2024–25)**, theory-focused, 16 lessons — [arXiv 2507.11536](https://arxiv.org/abs/2507.11536)
  - Density matrices, channels, distance measures, the stabilizer formalism, toric and surface codes.
- **PennyLane Codebook (self-study)**: modules through Hamiltonian simulation — product formulas, linear combinations of unitaries (LCU), qubitization, quantum signal processing (QSP), quantum singular value transformation (QSVT) — [Codebook](https://pennylane.ai/codebook)
- **Johns Hopkins graduate CS class (Kashani & Zaret, 2023)**: Yao.jl for superposition, Bell and GHZ states — [arXiv 2302.12889](https://arxiv.org/abs/2302.12889)

**Explicit "why numerics" statements**

- QuTiP 5 paper: the open-system problem "is often impossible to solve in practice when the number of constituent systems and the dimensionality of the Hilbert space become too large." It calls QuTiP "a research, teaching, and industrial tool." — [arXiv 2412.04705](https://arxiv.org/html/2412.04705v1)
- The Izaac & Wang blurb (above) frames numerics as what lets students go beyond the few analytically solvable potentials.
- The WSU Phys 555 quote above makes the same point for Landau–Zener and any non-commuting H(t).
- QuantumToolbox.jl (2025) on strongly driven nonlinear systems: these need states that "occupy a significant portion of the Hilbert space," and semiclassical cumulant approximations "may fail to accurately characterize quantum fluctuations." — [arXiv 2504.21440](https://arxiv.org/abs/2504.21440)

### Inferences
- Problems where numerics materially decide the answer and that also appear in teaching:
  - RWA breakdown and Bloch–Siegert-type shifts (TCD worksheet 1; ultrastrong Jaynes–Cummings in qutip-lectures).
  - Driven-dissipative steady states and emission spectra (Mollow triplet; Jaynes–Cummings spectrum).
  - Trajectories versus master equation (TCD, QuTiP, Campaioli).
  - Finite-time Landau–Zener (AJP 2023 project paper, WSU Phys 555, student theses).
  - Anharmonic and double-well potentials (Giannozzi, Izaac & Wang, Steck).
  - Exact diagonalization beyond hand size (ETH, Quantumandu).
  - Monte Carlo QEC thresholds (NordIQuEst).
- Recent quantum-optics and open-systems teaching runs almost entirely on QuTiP. Julia tools appear occasionally, and no 2022–2026 master's quantum-optics course built on Mathematica was found.
- None of the syllabi or tables of contents examined lists Fock-space truncation or solver-tolerance convergence checks as a topic. Combined with section 2.1, the convergence step seems to be learned the hard way, which is a teaching gap. This is based on what was absent, not confirmed.
- In circuit-based quantum-information courses (UMass Qiskit labs, IBM/Watrous, PennyLane Codebook), simulation mostly illustrates small, exactly solvable cases. Numerics rarely decides the answer there, except in noise, threshold and leakage studies.

### Gaps
- Not found in any verified syllabus:
  - Jaynes–Cummings collapse and revival as an assignment.
  - g²(τ) antibunching.
  - The Bloch–Siegert shift by name.
  - Kerr-cat generation.
  - SDP-based quantities.
  - scqubits, QuantumToolbox.jl, Dynamiqs, Strawberry Fields or Cirq used in course assignments.
- Missing course details:
  - The instructor, institution and year of the "Numerical Methods for Quantum Optics" course.
  - The year of the TCD module.
  - Companion code for García-Ripoll's 2022 superconducting-circuits book.
- Could not be retrieved: Giannozzi's full table of contents, Newman's full exercise list (HTTP 403), Steck's PDF, and the TU Delft "Applications of Quantum Mechanics" Landau–Zener page (page too large to fetch).
- Instructor quotes tying a specific problem to the need for numerics are rare; most "why numerics" statements come from book and paper prefaces.

## 4. Which named problems recur across frameworks as de facto standard benchmarks?

### Takeaway
- **Gate model:** QFT, GHZ/W states, Grover, quantum volume, Sycamore-style random circuits and Trotterized Heisenberg/transverse-field Ising dynamics recur across independent studies and the MQT Bench and SupermarQ suites. IBM's 127-qubit kicked-Ising experiment has become the reference problem for approximate classical simulation.
- **Open systems:** the Jaynes–Cummings model, the driven-damped cavity, the driven Kerr oscillator, the dissipative transverse-field Ising chain and qubit stochastic-master-equation trajectories recur. HEOM tools use spin-boson and FMO.
- **QEC:** Stim's generated repetition, surface and color-code memory circuits.
- **Superconducting circuits:** transmon, fluxonium and 0-π.
- **Benchmark design gap:** the open-system benchmarks time a fixed configuration rather than measuring accuracy against a converged reference. The exceptions are photonic probability loss versus cutoff, approximate-simulator truncation δ, and Dynamiqs' rejected-step counts.

### Cited Findings

**Gate-model circuits**

- QFT:
  - Jamadagni 2024 — [arXiv 2401.09076](https://arxiv.org/abs/2401.09076)
  - Cerrudo-Herrera 2026 — [DOI](https://doi.org/10.1177/00375497261438515)
  - Díaz-Camacho 2025, inverse QFT — [arXiv 2512.01807](https://arxiv.org/abs/2512.01807)
  - Leonteva 2025, MQT Bench `qft`/`qftentangled`; easiest circuits for MPS — [arXiv 2504.14027](https://arxiv.org/abs/2504.14027)
  - QuTiP tutorials, 10-qubit pulse-level QFT — [tutorials](https://qutip.org/qutip-tutorials/)
- Grover and quantum volume: Cerrudo-Herrera 2026; quantum volume also in Benchpress circuit construction ([arXiv 2409.08844](https://arxiv.org/abs/2409.08844)) and the UMass Qiskit labs.
- Random circuits: Sycamore-style fSim random circuits (Jamadagni). MQT Bench `random` is the "very hard" class, and the only one without evidence of polynomial-time simulability (Leonteva).
- Hamiltonian simulation of spin chains:
  - Trotterized XYZ-Heisenberg (Jamadagni).
  - SupermarQ "Hamiltonian simulation" of the 1D transverse-field Ising model — [arXiv 2202.11045](https://arxiv.org/abs/2202.11045); [Quantum Benchmark Zoo](https://quantumbenchmarkzoo.org/content/benchmarking-tools/benchmark-suites/Supermarq)
  - Benchpress Hamiltonian-simulation construction.
  - AppQSim, a 2025 application-oriented Hamiltonian-simulation suite (title only) — [arXiv 2503.04298](https://arxiv.org/abs/2503.04298)
- **SupermarQ (HPCA 2022)**: GHZ, Mermin–Bell, bit-flip and phase-flip repetition codes, QAOA (Sherrington–Kirkpatrick and MaxCut), VQE, transverse-field Ising Hamiltonian simulation — [arXiv 2202.11045](https://arxiv.org/abs/2202.11045)
- **MQT Bench**, Quetschlich, Burgholzer & Wille, *Quantum* 7, 1062 (2023) — [arXiv 2204.13719](https://arxiv.org/abs/2204.13719); [repo](https://github.com/munich-quantum-toolkit/bench)
  - More than 70,000 circuits from 2 to 130 qubits at four abstraction levels.
  - Families in the 2026 repo:
    - Basic algorithms: ae, bv, dj, ghz, ghz_dynamic, graphstate, grover, hhl, qaoa, qft, qftentangled, qnn, qpeexact, qpeinexact, qwalk, randomcircuit, shor, wstate.
    - VQE ansätze: real_amp, su2, two_local.
    - Arithmetic: several adders and multipliers.
    - QEC: Steane and Shor codes.
    - Dynamic circuits: dynamic_qft.

**Approximate and utility-scale classical simulation**

- IBM's 127-qubit kicked-Ising "utility" experiment (*Nature* 2023) became a reference problem for classical methods:
  - "Efficient tensor network simulation of IBM's Eagle kicked Ising experiment" (belief-propagation tensor networks; reports about 10⁻¹⁴ accuracy for average magnetization at 5 Trotter steps, per a search summary) — [arXiv 2306.14887](https://arxiv.org/abs/2306.14887)
  - "Fast and converged classical simulations of evidence for the utility of quantum computing before fault tolerance" — [Sci. Adv.](https://www.science.org/doi/10.1126/sciadv.adk4321)
  - Pauli-path GPU benchmarks (Mazumder 2026) — [arXiv 2607.09882](https://arxiv.org/abs/2607.09882)

**Stabilizer and QEC**

- Stim's built-in `stim.Circuit.generated` families (verified in the API reference):
  - `repetition_code:memory`
  - `surface_code:rotated_memory_x/z`
  - `surface_code:unrotated_memory_x/z`
  - `color_code:memory_xyz`

  Source: [Stim API](https://github.com/quantumlib/Stim/blob/main/doc/python_api_reference_vDev.md). The NordIQuEst threshold tutorial uses these circuits (section 3).
- Random Clifford circuits at scale: QuaSARQ, up to 180,000 qubits — [arXiv 2603.14641](https://arxiv.org/abs/2603.14641). Random Clifford construction is also a Benchpress test.

**Open-system and bosonic solvers**

- **QuTiP's own benchmark suite** (qutip-benchmark, `qutip_benchmark/benchmarks/bench_solvers.py`; repo active, last commit 2026-08-25) — [repo](https://github.com/qutip/qutip-benchmark)
  - mesolve and mcsolve on "Jaynes-Cummings", "Cavity" and "Qubit Spin Chain", for Hilbert dimensions 4–128.
  - steadystate on Jaynes–Cummings and Cavity.
- **QuantumOptics.jl benchmarks** (older; last commit 2020-11-24) — [repo](https://github.com/qojulia/QuantumOptics.jl-benchmarks)
  - Elementary operations: coherent state, displacement, expectation values and variances, partial trace, Q-function, Wigner function.
  - Schrödinger, master-equation and Monte Carlo wave-function (MCWF) evolution for:
    - a cavity;
    - the Jaynes–Cummings model;
    - a particle in a potential (including an FFT variant);
    - time-dependent cavity, Jaynes–Cummings and particle variants.
- **QuantumToolbox.jl paper (2025)**: driven Kerr oscillator with photon loss for mesolve, mcsolve, smesolve and autodiff; dissipative transverse-field Ising chain for CPU/GPU scaling — [arXiv 2504.21440](https://arxiv.org/abs/2504.21440)
- **Dynamiqs benchmark suite (2026)** — [README](https://github.com/dynamiqs/dynamiqs/blob/main/benchmarks/README.md)
  - Driven transmon with a DRAG-like pulse.
  - Transverse-field Ising chain.
  - Piecewise-constant drive crossing discontinuities.
  - Driven-damped cavity, "the canonical large open bosonic system."
  - Cat qubit stabilized by two-photon dissipation, batched over drive amplitude.
  - Cross-resonance gate between two transmons, batched as a calibration sweep.
  - Reverse-mode gradient through mesolve, "the pulse-optimization workload."
  - Floquet driven Kerr resonator.
  - Jump and diffusive SSE/SME unravelings of a cavity.
- **Small transmon Lindblad propagation** at d = 3, 9, 27 with GRAPE-style optimisation — [arXiv 2603.18052](https://arxiv.org/abs/2603.18052)
- **Qubit SME trajectories**: Dynamiqs vendor claim, 100,000 trajectories — [blog snippet](https://alice-bob.com/blog/dynamiqs-gpu-opensource-quantum-simulation-library/)

**Non-Markovian (HEOM)**

- The Fenna–Matthews–Olson light-harvesting complex and spin-boson/spin–bath models are the showcase problems of QuTiP-BoFiN, HierarchicalEOM.jl and the QuTiP HEOM tutorials — [BoFiN](https://arxiv.org/abs/2010.10806); [HierarchicalEOM.jl](https://arxiv.org/abs/2306.07522); [tutorials](https://qutip.org/qutip-tutorials/)

**Photonic**

- Timings of the hafnian, loop hafnian and torontonian (Piquasso vs The Walrus).
- Gaussian boson sampling and boson sampling (Piquasso vs Perceval).
- Multi-layer CV neural networks under a Fock cutoff (Piquasso vs Strawberry Fields), with probability loss against the cutoff as an accuracy axis.

Source: [arXiv 2403.04006](https://arxiv.org/abs/2403.04006)

**Superconducting circuits**

- scqubits ships transmon, fluxonium, flux qubit, 0-π (reduced and full), cos2φ qubit, oscillator and generic or custom circuits — [repo](https://github.com/scqubits/scqubits)
- The 2026 model-adequacy benchmark suite (Larssen et al.) covers flux-dependent spectra, two-qubit couplings, driven dynamics, CZ gate, leakage and runtime — [arXiv 2605.23034](https://arxiv.org/abs/2605.23034)

### Inferences
- A small set of problems appears in teaching (section 3), in incumbent benchmarks and in GitHub bug reports alike:
  - the Jaynes–Cummings model;
  - the driven-damped cavity;
  - the driven Kerr oscillator and cat states;
  - the dissipative or Trotterized transverse-field Ising/Heisenberg chain;
  - the transmon and a two-transmon gate;
  - QFT/GHZ;
  - surface-code memory.

  Examples of the bug-report overlap: QuTiP [#3026](https://github.com/qutip/qutip/issues/3026) involves a Jaynes–Cummings steady state, [#2549](https://github.com/qutip/qutip/issues/2549) a transmon/Duffing model, and [#2819](https://github.com/qutip/qutip/issues/2819) a Fock cutoff of 120 with α = 8. Worked studies on these problems can be compared directly with published incumbent numbers.
- Open-system benchmarks reward throughput at fixed settings. None scores error against a converged reference, robustness to discontinuous drives, or sensitivity to truncation and tolerance. Such accuracy-at-cost benchmarks are a gap. Dynamiqs' `nrej` metric and Piquasso's probability-loss-versus-cutoff plot are the closest precedents.
- Pulse-level transmon and cross-resonance workloads, Floquet drives and piecewise-constant controls are appearing in 2026 benchmark suites. This matches the IBM stack's withdrawal from pulse simulation (section 1D) and the QuTiP pulse-skipping reports (2.3).

### Gaps
- QASMBench and other circuit suites were not verified in this pass.
- No consensus open-system benchmark suite shared across packages was found. Each package benchmarks its own chosen models.
- HEOM benchmark parameters (bath cutoffs, hierarchy depth) used across packages were not compared.
- Tensor-network circuit-contraction benchmarks (e.g. Sycamore contraction with quimb/cotengra) were not separately verified for 2022–2026.
