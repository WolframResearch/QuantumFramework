# Bosonic and continuous-variable photonic simulation: the bar set by Strawberry Fields, The Walrus, MrMustard and Piquasso (QuTiP and QuantumOptics.jl for comparison), as of October 2026

Conventions: "checked 2026-10-04" marks facts read directly from the GitHub REST API, PyPI JSON or HTTP headers on that date. Quotes come from primary sources (repository files, issue threads, papers, staff forum posts). Inferences and the absence-based claims are kept in the Inferences sections.

## 1. Which numerical problem classes do these tools feature in tutorials, demos and papers?

### Takeaway
The four photonic incumbents mostly feature four kinds of problem: Gaussian boson sampling (GBS) and hafnian, torontonian and permanent combinatorics; compiling and decomposing Gaussian circuits; differentiable state and gate learning (CV quantum neural networks); and heralded non-Gaussian state preparation with photon-number-resolving (PNR) detection (cat, GKP, cubic-phase and Fock states). Kerr dynamics and Wigner-function galleries mainly live in QuTiP. None of the four photonic tools showcases a worked lossy-interferometry or quantum-Fisher-information study (NOON versus coherent versus squeezed light under loss). Several mention metrology as a motivation.

### Cited Findings

#### Strawberry Fields (SF, Xanadu)
- SF describes itself as "a full-stack Python library for designing, simulating, and optimizing continuous-variable quantum optical circuits". Its README advertises "graph and network optimization, machine learning, and chemistry" and an "end-to-end differentiable TensorFlow backend". — [SF README](https://github.com/XanaduAI/strawberryfields)
- SF has four backends. `fock` simulates "in a truncated Fock basis using NumPy". `gaussian` uses "the Gaussian formalism". `tf` is a truncated Fock basis "using TensorFlow". `bosonic` represents "states as linear combinations of Gaussian functions in phase space". — [SF backends docs](https://strawberryfields.readthedocs.io/en/stable/code/sf_backends.html)
- The example scripts still in the repository are IQP, boson_sampling, bosonic_tutorial_sampling, custom_operation, gate_teleportation, gaussian_boson_sampling, gaussian_cloning, gbs_data_visualization, gkp_phase_gate, hamiltonian_simulation, measurement_based_squeezing, optimization, quantum_neural_network and teleportation, plus an X8 hardware job file (checked 2026-10-04). — [SF examples folder](https://github.com/XanaduAI/strawberryfields/tree/master/examples)
- The documentation sidebar lists `sf.apps` (the GBS applications layer), `sf.tdm` (time-domain multiplexing) and "GBS datasets". — [SF docs landing page](https://strawberryfields.readthedocs.io/en/stable/)
- SF asks users to cite two papers: the platform paper and the applications paper on GBS-based algorithms. — [Killoran et al., Quantum 3, 129 (2019)](https://quantum-journal.org/papers/q-2019-03-11-129/); [Bromley et al., QST 5, 034010 (2020)](https://iopscience.iop.org/article/10.1088/2058-9565/ab8504/meta)
- The canonical GBS applications are dense-subgraph search and molecular vibronic spectra. — [Arrazola & Bromley, PRL 121, 030503 (2018)](https://doi.org/10.1103/PhysRevLett.121.030503); [Huh et al., Nat. Photon. 9, 615 (2015)](https://doi.org/10.1038/nphoton.2015.153)
- The bosonic backend paper simulates GKP, cat and Fock states as linear combinations of Gaussians. Its applications are Gaussian channels, threshold and PNR measurement, and gate teleportation, "with levels of accuracy that are not feasible with existing methods". — [Bourassa et al., PRX Quantum 2, 040315 (2021), arXiv:2103.05530](https://arxiv.org/abs/2103.05530)
- Optimization demos follow Arrazola et al. (state preparation and gate synthesis by machine learning) and Killoran et al. (CV quantum neural networks). — [Arrazola et al., QST 4, 024004 (2019)](https://doi.org/10.1088/2058-9565/aaf59e); [Killoran et al., PRR 1, 033063 (2019)](https://doi.org/10.1103/PhysRevResearch.1.033063)
- The website demos no longer exist. A Jan 2026 forum user reported that `strawberryfields.ai/photonics/demos/run_quantum_neural_network.html` was missing. Staff replied that the website "and all associated demos" had been pulled. — [PennyLane forum, 2026-02-03](https://discuss.pennylane.ai/t/strawberryfields-tutorials-missing/9223)
- In Aug 2025 staff pointed bosonic-backend users to a "realistic bosonic qubits" demo that shows "the effect that having finite energy has in the result". — [PennyLane forum, 2025-08-04](https://discuss.pennylane.ai/t/numerical-accuracy-of-the-simulator-and-the-fidelity-of-the-model-relative-to-a-physical-implementation/8855)

#### The Walrus (Xanadu; Gupt, Izaac, Quesada)
- The Walrus offers fast hafnians, loop hafnians and torontonians, including for structured matrices. It also provides a loop-hafnian interface for Gaussian-state calculations, sampling algorithms, classical approximation of hafnians of non-negative matrices, and multidimensional Hermite polynomials. — [The Walrus README](https://github.com/XanaduAI/thewalrus); [JOSS 4, 1705 (2019)](https://doi.org/10.21105/joss.01705)
- The gallery has four benchmark tutorials (basics, hafnian, permanent, torontonian), each with a timing plot. — [gallery.rst](https://github.com/XanaduAI/thewalrus/blob/master/docs/gallery/gallery.rst)
- The hafnian tutorial times random symmetric matrices for n = 2…30. — [hafnian_tutorial.ipynb](https://github.com/XanaduAI/thewalrus/blob/master/docs/gallery/hafnian_tutorial.ipynb)
- A "Non-Gaussian states gallery" builds Fock, kitten, cubic-phase, photon-added, cat, four-cat and GKP states from "Gaussian circuits and photon-number-resolved measurements". It credits Sabapathy et al. and Su et al. — [gallery.rst](https://github.com/XanaduAI/thewalrus/blob/master/docs/gallery/gallery.rst); [Sabapathy et al., PRA 100, 012326 (2019)](https://doi.org/10.1103/PhysRevA.100.012326); [Su, Myers & Sabapathy, PRA 100, 052301 (2019)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.100.052301)
- The companion paper on realistic heralded preparation is "Simulating realistic non-Gaussian state preparation". — [Quesada et al., PRA 100, 022341 (2019)](https://doi.org/10.1103/PhysRevA.100.022341)

#### MrMustard (Xanadu; Miatto, Yao, Quesada and others)
- MrMustard calls itself a "differentiable simulator with a sophisticated built-in optimizer". It "operates seamlessly across phase space and Fock space" with NumPy and JAX backends. — [MrMustard README](https://github.com/XanaduAI/MrMustard)
- Its representations are Bargmann, phase space, characteristic function, quadrature and Fock. It offers "Riemannian optimization on symplectic/unitary/orthogonal groups". — [MrMustard README](https://github.com/XanaduAI/MrMustard)
- Its components include a Mach-Zehnder, a "Sauron" state, `PNRSampler`, `ThresholdSampler` and `HomodyneSampler`. Projectors and POVMs are implemented as dual states. — [MrMustard README](https://github.com/XanaduAI/MrMustard)
- The README quick start builds a four-lobed cat state, plots Wigner functions and projects onto a Fock state. It then claims the Fock amplitudes at `fock_array(shape=(100, 4))` are "exact down to machine precision". — [MrMustard README](https://github.com/XanaduAI/MrMustard)
- The flagship paper showcases a 216-mode interferometer optimization (Borealis) and cat- and cubic-phase-state preparation with 2- and 3-mode circuits and Fock measurements. — [Yao, Miatto & Quesada, SciPost Phys. 17, 082 (2024), arXiv:2209.06069](https://arxiv.org/abs/2209.06069)
- A second paper covers noisy linear-optical circuits with PNR detectors. It lists GBS, GKP, cat and NOON state preparation, quantum computing and quantum metrology as applications. — [De Prins, Yao, Apte & Miatto, Quantum 7, 1097 (2023), arXiv:2303.08879](https://arxiv.org/abs/2303.08879)
- Features added in 2025–2026:
  - `stellar_roots()` and `plot_stellar_roots()`;
  - a preliminary "Wormhole algorithm" for Fock amplitudes "with complexity lineary in the total photon number";
  - a numerically stable `fock_diagonals`;
  - a `Kgate` (Kerr) "diagonal in the Fock basis" (merged 2026-04-21).

  Sources: [MrMustard CHANGELOG](https://github.com/XanaduAI/MrMustard/blob/develop/.github/CHANGELOG.md); [commit list](https://github.com/XanaduAI/MrMustard/commits/develop)
- The repository has no tutorial notebooks, and the docs are API reference only (checked 2026-10-04). — [MrMustard docs](https://mrmustard.readthedocs.io/en/latest/); [repo tree](https://github.com/XanaduAI/MrMustard)

#### Piquasso (Budapest Quantum Computing Group)
- Piquasso is a "full-stack open-source software platform for the simulation and programming of photonic quantum computers" with "optional high-performance C++ backends". — [Kolarovszki et al., Quantum 9, 1708 (2025)](https://quantum-journal.org/papers/q-2025-04-15-1708/)
- It has PureFock, Fock, Gaussian and Sampling simulators, with TensorFlow and JAX for differentiation. — [arXiv:2403.04006v3 (HTML)](https://arxiv.org/html/2403.04006v3)
- The paper's showcased benchmarks are:
  - CVQNN state learning against SF;
  - GBS against SF;
  - hafnian, loop-hafnian and torontonian kernels against The Walrus;
  - boson sampling against Perceval.

  Source: [arXiv:2403.04006v3](https://arxiv.org/html/2403.04006v3)
- v7.0.0 added a tutorial on dense k-subgraphs via GBS and mid-circuit measurement. Its example heralds a single-photon-like state from two-mode squeezing with `r = log(1+sqrt(2))`. — [Piquasso v7.0.0 release](https://github.com/Budapest-Quantum-Computing-Group/piquasso/releases/tag/v7.0.0)
- v7.2.0 added a SNAP gate, a "Boson Sampling as a hardware accelerator for Monte Carlo integration" tutorial, and parity expectation values for Gaussian states. — [v7.2.0 release](https://github.com/Budapest-Quantum-Computing-Group/piquasso/releases/tag/v7.2.0)
- v8.0.0 (2026-06-30) added:
  - a linear cross-entropy benchmarking (LXEB) module;
  - a CUDA permanent kernel via JAX FFI;
  - `Kerr`/`CrossKerr` in `PassiveSimulator`;
  - partial distinguishability;
  - imperfect PNR detectors;
  - a `UniformLoss` channel;
  - automatic simulator selection.

  Source: [v8.0.0 release](https://github.com/Budapest-Quantum-Computing-Group/piquasso/releases/tag/v8.0.0)

#### QuTiP and QuantumOptics.jl (comparison only)
- QuTiP's tutorial set covers parametric amplifiers (Lecture 5), squeezed oscillator states (Lecture 9), decay into a squeezed vacuum (Lecture 12), Kerr nonlinearities (Lecture 14) and a gallery of Wigner functions (Lecture 16). It also includes single-photon interference, "measures-trajectories-cats-kerr" and an optomechanical steady state. — [qutip-tutorials repository](https://github.com/qutip/qutip-tutorials); [Lecture 14, Kerr](https://github.com/qutip/qutip-tutorials/blob/main/tutorials-v5/lectures/Lecture-14-Kerr-nonlinearities.md)
- QuantumOptics.jl's example notebooks are cavity-QED and atom–cavity problems: Jaynes–Cummings, pumped cavity, optomechanical and cavity cooling, lasing, superradiant laser, Ramsey, and others. None is a photonic CV circuit, GBS or metrology example. — [QuantumOptics.jl-examples](https://github.com/qojulia/QuantumOptics.jl-examples)

### Inferences
- How well each class is served (inference from the findings above):
  - **GBS and hafnians:** very well served (The Walrus, Piquasso C++/CUDA kernels, SF apps).
  - **Gaussian circuits and decompositions:** well served (all four tools).
  - **Heralded cat, GKP and cubic-phase preparation and optimization:** well served by The Walrus's gallery and MrMustard's papers, but MrMustard ships no tutorials.
  - **GKP and cat states without a cutoff:** SF's bosonic backend only, now deprecated.
  - **Kerr and Wigner negativity:** mostly QuTiP. SF and Piquasso support Kerr gates; MrMustard added one only in 2026.
  - **Lossy interferometry and QFI versus N:** not showcased by any incumbent. MrMustard's papers name metrology and NOON states as applications without a worked lossy-interferometry or QFI example. This absence is inferred from the example lists, the gallery, the README, the changelog and the release notes reviewed. It is not proof of absence.
- The incumbents' showcase problems are aimed at photonic quantum-computing research audiences: advantage experiments, compilers and ML-style circuit training. Fewer of them are framed as textbook quantum-optics questions answered by numerics.

### Gaps
- I could not see the full list of SF's former website demos. strawberryfields.ai now redirects, and the Wayback Machine could not be fetched from this environment. Only the repository example scripts and the demos named in forum posts are verified.
- I did not verify whether PennyLane's current CV demos (`gbs`, `qonn`, `quantum_neural_net`, `tutorial_gaussian_transformation`, `tutorial_mbqc`, `tutorial_photonics`; see [PennyLaneAI/qml](https://github.com/PennyLaneAI/qml)) still execute on current PennyLane releases.
- Other bosonic tools were outside scope and not surveyed: Dynamiqs, Bosonic Qiskit, Perceval (discrete-variable linear optics), Gabs.jl.

## 2. How are Fock-space truncation and cutoff convergence handled, what pain points do users report, and how far do cutoffs and mode counts go?

### Takeaway
There are two truncation semantics:
- **Recurse, then truncate** (SF, MrMustard, The Walrus): exact matrix elements of the infinite-dimensional operator are computed by recursion and then truncated. Norm leaks, and a trace below 1 is the user's warning sign.
- **Truncate, then exponentiate** (QuTiP and QuantumOptics.jl by default): the truncated generator is exponentiated. The result stays unitary and normalized, but amplitudes near the cutoff are wrong.

Documented failures cluster at cutoffs of a few tens to a few hundred photons:
- cubic-phase gate errors, with overflow above cutoff 70;
- GKP Fock "probabilities" below 0 and above 1 near cutoff 40;
- Wigner instability around cutoff 200;
- memory saturation at about 50 photons per mode on a desktop.

The mitigations are MrMustard's automatic Fock-shape selection (99.999% probability, max 50), Piquasso's global photon-number cutoff and SF's cutoff-free Gaussian-sum backend.

### Cited Findings

#### Two truncation semantics
- SF staff (Oct 2023) explained that "taking the matrix exponential of truncated infinite dimensional operators is not the same as taking the matrix exponential of the infinite operators then truncating". SF computes operators "exactly using recursive methods" (citing arXiv:2004.11002), so "density matrices may not be fully normalized if the cutoff is not high enough". — [PennyLane forum, "Fock Backend truncation"](https://discuss.pennylane.ai/t/fock-backend-truncation/3545)
- In the same thread staff wrote: "Having low cutoff dimension means your state representation is incomplete and, as such, probabilities will not sum up to one." — [forum thread 3545](https://discuss.pennylane.ai/t/fock-backend-truncation/3545)
- SF's own guidance is only qualitative: "increasing the cutoff results in higher accuracy at a cost of increased memory consumption". — [SF doc/introduction/circuits.rst](https://github.com/XanaduAI/strawberryfields/blob/master/doc/introduction/circuits.rst)
- The recursive Fock-amplitude technique for Gaussian gates comes from Miatto and Quesada. — [Miatto & Quesada, Quantum 4, 366 (2020)](https://quantum-journal.org/papers/q-2020-11-30-366/)
- QuTiP's `coherent(N, alpha, method='operator')` (the default) displaces vacuum "using the displacement operator defined in the truncated Hilbert space", which "guarantees that the resulting state is normalized". The `'analytic'` method "does not guarantee that the state is normalized if truncated to a small number of Fock states, but would in that case give more accurate coefficients". — [qutip/core/states.py](https://github.com/qutip/qutip/blob/master/qutip/core/states.py)
- According to a 2020 issue that is still open, QuTiP's `displace` uses matrix exponentiation ("as qutip does now"). The issue asks for an analytic `method='analytical'` based on Laguerre polynomials. — [qutip#1293](https://github.com/qutip/qutip/issues/1293)
- In QuantumOptics.jl, `displace` and `squeeze` are "computed as the matrix exponential of finite-dimensional (truncated) creation and annihilation operators". A separate `displace_analytical` gives exact matrix elements (borrowed from IonSim.jl). — [QuantumOpticsBase.jl src/fock.jl](https://github.com/qojulia/QuantumOpticsBase.jl/blob/master/src/fock.jl)
- A QuantumOptics.jl user found `exp(dense(5im*create(b)))` differing from QuTiP by about 1.37e12 at dimension 200, "and higher it is, worse is the error". — [QuantumOptics.jl#281](https://github.com/qojulia/QuantumOptics.jl/issues/281)

#### Strawberry Fields pain points
- The cubic-phase gate docstring warns: "The cubic phase gate has lower accuracy than the Kerr gate at the same cutoff dimension." Non-Gaussian gates "can only be used in the Fock backends". — [strawberryfields/ops.py](https://github.com/XanaduAI/strawberryfields/blob/master/strawberryfields/ops.py)
- Open issue #753 (2025-02-25) reports for the cubic-phase gate on vacuum:
  - the trace is "less than one for all tested cutoff dimensions";
  - the largest eigenvalue below 1 suggests "spurious correlations with the second qumode";
  - errors accumulate so that "the errors prevent accurate simulation of even a single quartic phase gate";
  - "When increasing the cutoff dimension past 70 … overflow errors prevent the simulation from even running."

  Source: [SF#753](https://github.com/XanaduAI/strawberryfields/issues/753)
- A Sept 2026 preprint gives the exact Fock amplitudes of the cubic-phase gate on any pure Gaussian input (Airy functions). It notes the usual route "build[s] the cubic quadrature in a truncated Fock space and exponentiate[s] it, which alters the generator". Benchmarked against the exact result, the truncated construction's error "at cutoffs in common use is found to be of the same size as the quantity being computed". — [Cordshooli, arXiv:2609.13077](https://arxiv.org/abs/2609.13077)
- A user simulating two GKP logical-0 states (ε = 0.125) with cutoffs above 40 saw Fock "probabilities" both negative and above 1 near the cutoff. Staff called it "a known issue stemming from instabilities on Gaussian operations" in The Walrus. Follow-ups in Feb 2026 found it unresolved as SF was deprecated. — [PennyLane forum thread 2044](https://discuss.pennylane.ai/t/fock-probabilities-greater-than-1-in-circuits-with-gkp-states/2044)
- Wigner-function instability:
  - In 2022 a user saw instability "for the 200 and 500 cutoff dimension" with highly squeezed states, while "the data in the state is fine".
  - In 2024 another user saw it with 6 dB squeezing (r = −0.7) "at around a cutoff of 200".
  - Staff: "It's still unclear why you get this instability", and suggested switching to MrMustard.

  Source: [PennyLane forum thread 2293](https://discuss.pennylane.ai/t/wigner-calculation-unstable-for-the-fock-backend/2293)
- In the Fock backend, `measure_homodyne` on two cat states built probability vectors with p < 0 or p > 1 in "30-40% of the time". The issue is now closed. — [SF#354](https://github.com/XanaduAI/strawberryfields/issues/354)
- Older correctness bugs found by users, now closed: two-mode squeezing with zero parameters was not the identity at cutoff 4 (2019); "reduced_dm not correct for subsystems of more than 1 mode"; "pure states are not pure". — [SF#186](https://github.com/XanaduAI/strawberryfields/issues/186); [SF#469](https://github.com/XanaduAI/strawberryfields/issues/469); [SF#488](https://github.com/XanaduAI/strawberryfields/issues/488)
- Open memory issue: "Additional memory creation from partial trace during measurements" (2020). — [SF#464](https://github.com/XanaduAI/strawberryfields/issues/464)
- The bosonic backend paper sets out the Fock-cutoff problem explicitly:
  - "high-quality and therefore high-energy cat and GKP states require a high photon-number cutoff, incurring large memory loads and processing times";
  - "displacements and squeezing can quickly push Fock distributions beyond the energy cutoff";
  - loss forces a density matrix, "squaring the number of elements";
  - Fock density-matrix elements scale "like |α|^{4N}".

  Source: [Bourassa et al., arXiv:2103.05530](https://arxiv.org/abs/2103.05530)
- In that paper's comparison, SF's Fock backend was run with "a photon number cutoff of 50 photons per mode—going beyond this value saturates the memory of the standard desktop terminal". The authors add that "although n = 32 is comfortably within the cutoff of 50 photons, the Fock representation cannot accurately construct the Wigner function for α = 4". The naive cutoff estimate is "an underestimation of the required photon number cutoff". — [Bourassa et al., arXiv:2103.05530](https://arxiv.org/abs/2103.05530)
- SF's answer is the bosonic backend, where "a `cutoff_dim` is not needed … since it doesn't use Fock states" and "cat states can be approximated to arbitrary precision". — [PennyLane forum, 2025-08-04](https://discuss.pennylane.ai/t/numerical-accuracy-of-the-simulator-and-the-fidelity-of-the-model-relative-to-a-physical-implementation/8855)

#### MrMustard
- Default settings, which choose the Fock array shape automatically:
  - `AUTOSHAPE_PROBABILITY = 0.99999`, `AUTOSHAPE_MAX = 50`, `AUTOSHAPE_MIN = 1`;
  - `STABLE_FOCK_CONVERSION = True` ("more stable, but slower");
  - an internal 128-bit precision setting for Hermite polynomials.

  Source: [mrmustard/settings/settings.py](https://github.com/XanaduAI/MrMustard/blob/develop/mrmustard/settings/settings.py)
- Changelog fixes and changes:
  - "Fixed conflated use of cutoff (max photon number) and shape (array size)";
  - "Made `fock_diagonals` numerically stable";
  - "Wigner computation now uses only Clenshaw";
  - a quadrature bug that forced "the number of quadrature points to equal the Fock cutoff and silently giving wrong results".

  Source: [MrMustard CHANGELOG](https://github.com/XanaduAI/MrMustard/blob/develop/.github/CHANGELOG.md)
- Cutoff-related issues: "AUTOSHAPE for fock states" (open since 2024-09); "MemoryError" (2024-09); "Autocutoffs not being respected" (2023); "Wigner plotting unstable" (2023); "probability vectors do not update upon cutoff increase" (2023). — [MrMustard#490](https://github.com/XanaduAI/MrMustard/issues/490); [#481](https://github.com/XanaduAI/MrMustard/issues/481); [#250](https://github.com/XanaduAI/MrMustard/issues/250); [#252](https://github.com/XanaduAI/MrMustard/issues/252); [#238](https://github.com/XanaduAI/MrMustard/issues/238)
- For noisy circuits with PNR detection, the cost falls from O(M² ∏_{D∪U} C_i²) to O(M² ∏_U C_i² ∏_D C_i), where D and U are the detected and undetected modes. This allows "twice the number of modes" with the same resources. — [De Prins et al., arXiv:2303.08879](https://arxiv.org/abs/2303.08879)

#### Piquasso
- Piquasso uses a "global cutoff" on total photon number, so the state dimension is C(d+c−1, c−1). SF uses a "local cutoff" per mode, giving c^d. — [arXiv:2403.04006v3](https://arxiv.org/html/2403.04006v3)
- The paper reports, against SF v0.23.0:
  - a 10-mode, 4-layer CVQNN (Piquasso global cutoff 3–16 versus SF local cutoff 3–6): "For any probability loss, the Piquasso Fock simulator … has lower execution time";
  - gradient timings that "match until 8 modes", after which Piquasso pulls ahead;
  - SF runs "halted due to excessive memory usage".

  Source: [arXiv:2403.04006v3](https://arxiv.org/html/2403.04006v3)
- In 2025 Piquasso added a validation error when "specified cutoff is not high enough". In v7.0.0 the Fock simulator's "cutoff auto-defaults, can be set via Config". — [piquasso#472](https://github.com/Budapest-Quantum-Computing-Group/piquasso/issues/472); [v7.0.0 release](https://github.com/Budapest-Quantum-Computing-Group/piquasso/releases/tag/v7.0.0)

#### QuTiP (phase-space functions)
- `wigner()` has the methods 'clenshaw' (the signature default), 'iterative', 'laguerre' and 'fft'. "The 'clenshaw' method is the preferred method for dealing with density matrices that have a large number of excitations (>~50)." The same docstring also says "The 'iterative' method is default", which contradicts the signature. — [qutip/wigner.py](https://github.com/qutip/qutip/blob/master/qutip/wigner.py)
- A July 2026 bug: "`qfunc()` cannot accurately calculate the Husimi-Q function of coherent states with large average photon numbers". It is now closed. — [qutip#2932](https://github.com/qutip/qutip/issues/2932)

#### Published scale and performance figures (GBS and hafnians)
- The fastest exact hafnian algorithms run in O(n³ 2^{n/2}). They were benchmarked up to 56×56 on the Titan supercomputer, and one would need "the 288000 CPUs of this machine for about a month and a half to compute the hafnian of a 100 × 100 matrix". — [Björklund, Gupt & Quesada, ACM JEA (2019), arXiv:1805.12498](https://arxiv.org/abs/1805.12498)
- Bulmer et al. simulated GBS "with up to 100 modes and 92 photons" on about 100,000 cores. They cut estimated runtimes from "600 million years" to "several months". — [Bulmer et al., Sci. Adv. 8, eabl9236 (2022), arXiv:2108.01622](https://arxiv.org/abs/2108.01622)
- Borealis has 216 squeezed modes. Exact sampling would take ">9,000 years" against "36 μs" on the device. — [Madsen et al., Nature 606, 75 (2022)](https://www.nature.com/articles/s41586-022-04725-x)
- A tensor-network classical sampler "can simulate the ideal distribution better than the experiment can". — [Oh et al., Nat. Phys. 20, 1461 (2024)](https://www.nature.com/articles/s41567-024-02535-8)
- Jiuzhang 4.0 used 1,024 squeezed states in 8,176 modes, with up to 3,050 photons detected. — [Nature (2026)](https://www.nature.com/articles/s41586-026-10523-6); [arXiv:2508.09092](https://arxiv.org/html/2508.09092v1)
- An April 2026 classical algorithm's output-count statistics are "closer to the exact solution than current experiments up to 1152 modes". — [Goodman et al., arXiv:2604.12330](https://arxiv.org/abs/2604.12330)

### Inferences
- The incumbents differ in where truncation goes wrong:
  - SF and MrMustard keep exact amplitudes but lose normalization. Trace loss is a free diagnostic.
  - QuTiP and QuantumOptics.jl keep normalization but distort the amplitudes, and nothing flags it unless the user knows to use `'analytic'` or `displace_analytical`.
- Neither group offers built-in, user-facing cutoff-convergence studies (error versus cutoff curves) as a documented workflow. MrMustard's autoshape is the closest, a probability-mass heuristic capped at 50 by default. A tool that made convergence and error certification explicit would beat the current bar. This is inferred from the docs and issues reviewed.
- Practical ceilings implied by the sources:
  - single-mode non-Gaussian work in SF degrades between cutoff 40 and 200 (GKP, cubic phase, Wigner);
  - multimode mixed-state Fock simulation hits desktop memory at about 50 photons per mode or about 8–10 modes at small local cutoffs.

  Arithmetic check: 6^10 ≈ 6.0×10^7 complex amplitudes is about 1 GB in complex128 for a pure state, and squared for a density matrix.
- Exact GBS beyond about 100 modes or about 90 photons is supercomputer territory and is dominated by specialist kernels (The Walrus, Piquasso C++/CUDA, Bulmer et al.). It is not a credible numerical arena for a general-purpose tool.

### Gaps
- I found no published SF or QuTiP table of maximum cutoff × mode count against runtime. The figures above are scattered single data points.
- I could not verify the exact Piquasso benchmark timings (absolute seconds). The numbers come from a summary of the arXiv HTML and should be re-read from the paper's figures before quoting.
- I did not check whether The Walrus's 2025 fixes (#403 Takagi, #404 `hermite_renorm`) resolve the GKP near-cutoff instability.
- I did not fetch the body of the MrMustard MemoryError issue (#481) beyond its expected-behaviour statement.

## 3. What is the maintenance status of each tool in 2025-2026, and what does Xanadu now recommend?

### Takeaway
Xanadu's open-source photonic simulation stack is effectively frozen:
- **Strawberry Fields** is officially "fully deprecated": website and demos were pulled in early 2026, the repository is archived, and the last release dates from 2022.
- **The Walrus and MrMustard** public repositories are archived. MrMustard's last public commits describe syncing with a private repository, and PyPI still serves the 2024 version.

Xanadu staff now point users to PennyLane, whose CV support is limited to the Gaussian `default.gaussian` device. Piquasso, QuTiP and QuantumOptics.jl are actively maintained.

### Cited Findings
- **Strawberry Fields:**
  - GitHub `archived: true`; last release v0.23.0 on 2022-06-01, plus a post1 the next day; PyPI latest 0.23.0 (checked 2026-10-04). — [SF releases](https://github.com/XanaduAI/strawberryfields/releases); [PyPI](https://pypi.org/project/StrawberryFields/)
  - The README states support for "Python version 3.8, 3.9, or 3.10". — [SF README](https://github.com/XanaduAI/strawberryfields)
  - A 2025-05-20 commit "Updates for scipy, numpy 2, and Python 3.12 (#757)" was never released. — [SF commits](https://github.com/XanaduAI/strawberryfields/commits/master)
  - PR #761 (2025-12-19) made the README say "We are not accepting Pull Requests for new features at this time" and deleted CHALLENGES.md. — [SF#761](https://github.com/XanaduAI/strawberryfields/pull/761)
  - PR #763 (2026-01-16), "Decommissioning the quantum cloud", added warnings that "cloud access is no longer available". — [SF#763](https://github.com/XanaduAI/strawberryfields/pull/763)
  - Staff, 2026-02-03: "StrawberryFields has been under maintenance mode for several years and it has now been fully deprecated so we've decided to pull down the StrawberryFields website and all associated demos. The last version of Strawberry Fields (v0.23) is still available on PyPI but it's no longer maintained. I encourage you to use PennyLane instead." — [PennyLane forum thread 9223](https://discuss.pennylane.ai/t/strawberryfields-tutorials-missing/9223)
  - `https://strawberryfields.ai/` and its demo URLs return HTTP 301 to xanadu.ai (checked 2026-10-04). The readthedocs API docs remain online, labelled "0.24.0-dev", with no demos. — [strawberryfields.ai](https://strawberryfields.ai/); [SF docs](https://strawberryfields.readthedocs.io/en/stable/)
  - Users kept posting SF questions through Feb 2026, for example "Cannot import strawberryfields module" (Oct 2025) and "Hybrid backends: qubits and bosons" (Feb 2026). — [SF forum category](https://discuss.pennylane.ai/c/photonic-software/strawberry-fields/11)
- **PennyLane–SF bridge:**
  - The plugin "will not be supported in newer versions of PennyLane. It is compatible with versions of PennyLane up to and including 0.29". Its last release is v0.29.1 (2023-05-05). — [pennylane-sf README](https://github.com/PennyLaneAI/pennylane-sf)
  - In July 2023 staff said "If you need to do CV quantum computing … use Strawberryfields". — [forum thread 3142](https://discuss.pennylane.ai/t/whats-the-main-difference-between-pennylane-and-strawberry-fields/3142)
  - The PennyLane 0.45.1 deprecations page lists no deprecation of `default.gaussian`; the only CV item is that `qml.X`/`qml.P` were renamed `QuadX`/`QuadP` (removed in v0.33). — [PennyLane deprecations](https://docs.pennylane.ai/en/stable/development/deprecations.html); [default.gaussian](https://pennylane.ai/devices/default-gaussian)
- **Recommendation drift:** in May 2024 Xanadu staff recommended MrMustard for an SF Fock-backend Wigner problem. By Feb 2026 the recommendation for SF users was PennyLane. — [forum thread 2293](https://discuss.pennylane.ai/t/wigner-calculation-unstable-for-the-fock-backend/2293); [forum thread 9223](https://discuss.pennylane.ai/t/strawberryfields-tutorials-missing/9223)
- **The Walrus:**
  - GitHub `archived: true`; last push 2026-07-24; last release v0.22.0 (2025-05-15; PyPI 2025-05-16). — [The Walrus releases](https://github.com/XanaduAI/thewalrus/releases); [PyPI](https://pypi.org/project/thewalrus/)
  - Last commits (Oct–Nov 2025): "Faster hermite_renorm using strides (#404)" and "Fix computation of Takagi unitary (#403)". The docs read 0.23.0-dev. — [commits](https://github.com/XanaduAI/thewalrus/commits/master); [docs](https://the-walrus.readthedocs.io/en/latest/)
- **MrMustard:**
  - GitHub `archived: true`; last push 2026-07-24; latest GitHub release v1.0.0a1 (2025-07-22). — [MrMustard releases](https://github.com/XanaduAI/MrMustard/releases)
  - PyPI latest is still 0.7.3 (2024-03-28), with `requires_python` "<3.12,>=3.9" (checked 2026-10-04). — [PyPI](https://pypi.org/project/mrmustard/)
  - The last merged PR (2026-06-01) is described as "Syncing the public MrMustard repository to the private". — [MrMustard#654](https://github.com/XanaduAI/MrMustard/pull/654)
  - The unreleased 1.0.0a2 changelog removes TensorFlow entirely (backends are now NumPy and JAX) and drops Python 3.10. — [CHANGELOG](https://github.com/XanaduAI/MrMustard/blob/develop/.github/CHANGELOG.md)
- **Piquasso:** not archived; v8.0.1 (2026-07-01); commits through 2026-09-24; PyPI supports Python ≥3.10, <3.15 (checked 2026-10-04). Peer-reviewed in Quantum (April 2025). — [Piquasso repo](https://github.com/Budapest-Quantum-Computing-Group/piquasso); [PyPI](https://pypi.org/project/piquasso/); [Quantum 9, 1708](https://quantum-journal.org/papers/q-2025-04-15-1708/)
- **QuTiP:** v5.3.1 (2026-08-04), active (checked 2026-10-04). — [QuTiP repo](https://github.com/qutip/qutip)
- **QuantumOptics.jl:** v1.2.10 (2026-09-11), with a "Release v1.3 with QuantumOpticsBase visualizations" commit on 2026-09-12 (checked 2026-10-04). — [QuantumOptics.jl repo](https://github.com/qojulia/QuantumOptics.jl)
- **Popularity proxies** (GitHub stars, checked 2026-10-04): SF 854, The Walrus 109, MrMustard 102, Piquasso 61, QuTiP 2,081, QuantumOptics.jl 624. — GitHub API for each repo linked above
- **Xanadu's current public focus** is hardware. It published an integrated photonic GKP source in Nature in June 2025. — [Nature 642, 587 (2025)](https://www.nature.com/articles/s41586-025-09044-5)

### Inferences
- In 2026 no actively maintained, publicly released Xanadu tool covers non-Gaussian CV (Fock-space) photonic simulation:
  - SF is deprecated and stuck at Python ≤3.10 on PyPI;
  - MrMustard's PyPI release is from 2024 (Python <3.12) and the public repo is frozen;
  - PennyLane only offers Gaussian CV.

  Piquasso is the main live open-source photonic simulator. QuTiP and QuantumOptics.jl are live general-purpose Fock-space tools.
- Users of SF tutorials (still cited in many courses and papers) have lost their reference material. This is an opening for well-documented worked studies, though it is only inferred from the demo removal and the continuing forum demand.
- The timing (SF cloud decommission in Jan 2026, The Walrus and MrMustard archived around July 2026, MrMustard synced to private) suggests Xanadu moved photonic-simulation development in-house. This is speculation and not confirmed by any announcement found.

### Gaps
- GitHub does not expose archive dates. The archive dates for The Walrus and MrMustard are inferred from `pushed_at` (2026-07-24) and are not verified.
- I found no official Xanadu blog or announcement explaining The Walrus or MrMustard archival, and no statement naming a successor library.
- I did not test whether `pip install strawberryfields` or `mrmustard` works on current Python and NumPy versions.

## 4. Which standard named benchmark problems exist for head-to-head comparison?

### Takeaway
Recognised, citable problems exist for each class:
- GBS: hafnian timing curves and experimental instances (Jiuzhang, Borealis), with The Walrus and Piquasso kernels as reference implementations;
- CVQNN state learning, the Piquasso-versus-SF benchmark;
- heralded non-Gaussian state preparation (Sabapathy, Su, Quesada; The Walrus gallery);
- GKP and cat states in the Gaussian-sum formalism (Bourassa 2021);
- the lossy Mach–Zehnder phase-estimation family, with a canonical figure at η = 0.9 and the bound Δφ ≥ √((1−η)/(ηN));
- Kerr collapse–revival and cat formation (Yurke–Stoler, Kirchmair);
- a 2026 exact cubic-phase-gate benchmark that exposes cutoff error.

### Cited Findings

#### GBS, hafnians and boson sampling
- The original GBS proposal is Hamilton et al. — [PRL 119, 170501 (2017)](https://doi.org/10.1103/PhysRevLett.119.170501)
- Boson sampling's complexity basis is Aaronson & Arkhipov. — [Theory of Computing 9 (2013)](https://doi.org/10.4086/toc.2013.v009a004)
- Experimental instances:
  - Jiuzhang 1.0 — [Zhong et al., Science 370, 1460 (2020)](https://doi.org/10.1126/science.abe8770)
  - Borealis, 216 modes — [Madsen et al., Nature 606, 75 (2022)](https://www.nature.com/articles/s41586-022-04725-x)
  - Jiuzhang 4.0, 8,176 modes and up to 3,050 photons — [Nature (2026)](https://www.nature.com/articles/s41586-026-10523-6)
- Classical-simulation benchmarks:
  - hafnian timing up to 56×56 — [Björklund et al., arXiv:1805.12498](https://arxiv.org/abs/1805.12498)
  - "Quadratic speed-up for simulating GBS" — [Quesada et al., PRX Quantum 3, 010306 (2022)](https://doi.org/10.1103/PRXQuantum.3.010306)
  - Bulmer et al., 100 modes and 92 photons — [arXiv:2108.01622](https://arxiv.org/abs/2108.01622)
  - Oh et al. — [Nat. Phys. (2024)](https://www.nature.com/articles/s41567-024-02535-8)
  - Goodman et al. — [arXiv:2604.12330](https://arxiv.org/abs/2604.12330)
- Software reference points:
  - The Walrus hafnian, permanent and torontonian benchmark tutorials — [gallery](https://github.com/XanaduAI/thewalrus/blob/master/docs/gallery/gallery.rst)
  - Piquasso against The Walrus v0.21.0 (hafnian, loop hafnian, torontonian), SF v0.23.0 (full GBS circuits) and Perceval v0.11.2 (boson sampling) — [arXiv:2403.04006v3](https://arxiv.org/html/2403.04006v3)
  - Piquasso's LXEB module — [v8.0.0](https://github.com/Budapest-Quantum-Computing-Group/piquasso/releases/tag/v8.0.0)
- Canonical GBS application benchmarks are dense subgraphs and molecular vibronic spectra. SF's applications layer and paper target GBS applications in general; I did not check its module-by-module coverage. — [Arrazola & Bromley 2018](https://doi.org/10.1103/PhysRevLett.121.030503); [Huh et al. 2015](https://doi.org/10.1038/nphoton.2015.153); [Bromley et al. 2020](https://iopscience.iop.org/article/10.1088/2058-9565/ab8504/meta)

#### Optimization and state preparation
- CVQNN state learning — [Killoran et al., PRR 1, 033063 (2019)](https://doi.org/10.1103/PhysRevResearch.1.033063); [Arrazola et al., QST 4, 024004 (2019)](https://doi.org/10.1088/2058-9565/aaf59e)
- The quantitative head-to-head is Piquasso against SF on 4-layer CVQNNs, at matched probability loss of about 10⁻². — [arXiv:2403.04006v3](https://arxiv.org/html/2403.04006v3)
- Heralded non-Gaussian states via PNR detection:
  - [Sabapathy et al., PRA 100, 012326 (2019)](https://doi.org/10.1103/PhysRevA.100.012326)
  - [Su et al., PRA 100, 052301 (2019)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.100.052301)
  - [Quesada et al., PRA 100, 022341 (2019)](https://doi.org/10.1103/PhysRevA.100.022341)
  - The Walrus gallery: Fock, kitten, cubic, photon-added, cat, four-cat, GKP — [gallery.rst](https://github.com/XanaduAI/thewalrus/blob/master/docs/gallery/gallery.rst)
  - MrMustard: cat and cubic-phase preparation, and the 216-mode Borealis interferometer — [arXiv:2209.06069](https://arxiv.org/abs/2209.06069)
- GKP and cat states: original GKP encoding; cutoff-free Gaussian-sum simulation, including an explicit comparison against SF Fock at cutoff 50. — [Gottesman, Kitaev & Preskill, PRA 64, 012310 (2001)](https://doi.org/10.1103/PhysRevA.64.012310); [Bourassa et al., arXiv:2103.05530](https://arxiv.org/abs/2103.05530)
- Cubic-phase gate exact Fock amplitudes as a benchmark against truncated exponentiation. — [arXiv:2609.13077](https://arxiv.org/abs/2609.13077)

#### Lossy interferometry and phase sensitivity (the metrology family)
- Shot-noise versus squeezed-vacuum injection — [Caves, PRD 23, 1693 (1981)](https://doi.org/10.1103/PhysRevD.23.1693)
- Heisenberg limit with coherent plus squeezed vacuum — [Pezzé & Smerzi, PRL 100, 073601 (2008)](https://doi.org/10.1103/PhysRevLett.100.073601)
- Optimal states under loss — [Dorner et al., PRL 102, 040403 (2009)](https://doi.org/10.1103/PhysRevLett.102.040403); [Demkowicz-Dobrzański et al., PRA 80, 013825 (2009)](https://doi.org/10.1103/PhysRevA.80.013825); [Kołodyński & Demkowicz-Dobrzański, PRA 82, 053804 (2010)](https://doi.org/10.1103/PhysRevA.82.053804)
- General noisy bounds — [Escher et al., Nat. Phys. 7, 406 (2011)](https://doi.org/10.1038/nphys1958); [Demkowicz-Dobrzański, Kołodyński & Guţă, Nat. Commun. 3, 1063 (2012)](https://doi.org/10.1038/ncomms2067)
- Laser-powered optimum — [Lang & Caves, PRL 111, 173601 (2013)](https://doi.org/10.1103/PhysRevLett.111.173601)
- NOON review — [Dowling, Contemp. Phys. 49, 125 (2008)](https://doi.org/10.1080/00107510802091298)
- Canonical benchmark figure: the review's Fig. 12 shows phase-estimation precision with "equal losses in both arms (η = 0.9)". In it:
  - optimal N-photon states saturate the asymptotic limit √((1−η)/(ηN));
  - "The NOON states … achieve nearly optimal precision only for low N (≤10) and rapidly diverge becoming out-performed by classical strategies";
  - coherent plus squeezed vacuum (Caves 1981) "in the presence of loss also saturates the asymptotic quantum limit".

  Source: [Demkowicz-Dobrzański, Jarzyna & Kołodyński, Prog. Opt. 60 (2015), arXiv:1405.7703](https://arxiv.org/abs/1405.7703)
- The same review gives the lossy bound "∆ϕ ≥ √((1 −η)/(ηN)), where η is the overall power transmission". For coherent plus squeezed input it gives Δφ ≈ √(e^{−2r} + f(η))/√⟨N⟩ (Eqs. 164–165). The quantum gain under loss is "bound to a constant factor". — [arXiv:1405.7703](https://arxiv.org/abs/1405.7703)
- Experimental anchors:
  - unconditional sub-shot-noise photonic phase estimation — [Slussarenko et al., Nat. Photon. 11, 700 (2017)](https://doi.org/10.1038/s41566-017-0011-5)
  - squeezed light in LIGO — [LIGO Scientific Collaboration, Nat. Photon. 7, 613 (2013)](https://doi.org/10.1038/nphoton.2013.177)

#### Kerr nonlinearity and Wigner negativity
- Kerr-generated cat states — [Yurke & Stoler, PRL 57, 13 (1986)](https://doi.org/10.1103/PhysRevLett.57.13)
- Collapse and revival from the single-photon Kerr effect — [Kirchmair et al., Nature 495, 205 (2013)](https://doi.org/10.1038/nature11902)
- Tutorial coverage — [QuTiP Lecture 14](https://github.com/qutip/qutip-tutorials/blob/main/tutorials-v5/lectures/Lecture-14-Kerr-nonlinearities.md)
- Negativity volume as the standard non-classicality measure — [Kenfack & Życzkowski, J. Opt. B 6, 396 (2004)](https://doi.org/10.1088/1464-4266/6/10/003)

### Inferences
- The metrology family is the only one with an extensive literature of closed-form limits and canonical comparison plots and no incumbent software implementation to beat. A "reproduce Fig. 12 of arXiv:1405.7703 and extend it to finite N and realistic η with certified cutoff convergence" study would be a recognisable head-to-head target against the literature rather than against a library.
- For GBS, a credible head-to-head is small-scale correctness: exact photon-number probabilities and threshold-click statistics for small instances, checked against `thewalrus` or Piquasso values. Competing on speed against The Walrus, Piquasso or Bulmer-style kernels looks unrealistic.
- For CVQNN learning, the Piquasso paper's protocol (10 modes, 4 layers, matched probability loss of about 10⁻², local versus global cutoff) is a ready-made, published head-to-head.

### Gaps
- I found no community-maintained benchmark suite (like qubit-world QASMBench or MQT Bench) for CV or bosonic simulators. Benchmarks are paper-specific.
- I did not extract numeric data points from Fig. 12 of the review or from the Dorner et al. optimal-state curves. They would need digitizing or recomputing for a quantitative comparison.

## 5. Which problem class is most compelling for a master's-level audience, and why?

### Takeaway
Inference, supported by the cited material: lossy interferometric phase estimation is the most compelling. It compares quantum Fisher information or phase sensitivity versus photon number for coherent, squeezed-plus-coherent and NOON inputs under loss. It is textbook content with exact analytic anchors. Finite-N, realistic-η answers (crossovers, optimal states) need numerics that are sensitive to Fock truncation. No incumbent photonic tool showcases it. Kerr cat dynamics with Wigner negativity is a strong second, and GBS is famous but a poor arena for numerical competition.

### Cited Findings
- The metrology problem has textbook-level analytic limits to validate numerics:
  - shot noise and the squeezed-light improvement — [Caves 1981](https://doi.org/10.1103/PhysRevD.23.1693)
  - the Heisenberg limit — [Pezzé & Smerzi 2008](https://doi.org/10.1103/PhysRevLett.100.073601)
  - the lossy bound √((1−η)/(ηN)) — [arXiv:1405.7703](https://arxiv.org/abs/1405.7703)
- The comparison has a qualitative twist that numerics decides: at η = 0.9, NOON states are near-optimal only for N ≤ 10 and then lose to classical strategies, while coherent plus squeezed light saturates the lossy bound. — [arXiv:1405.7703, Fig. 12](https://arxiv.org/abs/1405.7703)
- Optimal N-photon input states for lossy interferometers are the subject of dedicated papers, "Optimal Quantum Phase Estimation" and "Quantum phase estimation with lossy interferometers". The review plots them as a separate curve from NOON and Caves-type inputs. — [Dorner et al. 2009](https://doi.org/10.1103/PhysRevLett.102.040403); [Demkowicz-Dobrzański et al. 2009](https://doi.org/10.1103/PhysRevA.80.013825); [arXiv:1405.7703, Fig. 12](https://arxiv.org/abs/1405.7703)
- There are real-world anchors: squeezed light in LIGO and photonic sub-shot-noise experiments. — [Nat. Photon. 7, 613 (2013)](https://doi.org/10.1038/nphoton.2013.177); [Slussarenko et al. 2017](https://doi.org/10.1038/s41566-017-0011-5)
- Kerr dynamics:
  - cat formation at fractional revival times and collapse–revival are classic results — [Yurke & Stoler 1986](https://doi.org/10.1103/PhysRevLett.57.13); [Kirchmair et al. 2013](https://doi.org/10.1038/nature11902)
  - QuTiP already teaches it — [QuTiP Lecture 14](https://github.com/qutip/qutip-tutorials/blob/main/tutorials-v5/lectures/Lecture-14-Kerr-nonlinearities.md)
  - high-excitation Wigner functions are numerically delicate: QuTiP prefers Clenshaw above about 50 excitations, and SF's Wigner function was unstable at cutoff about 200 — [qutip/wigner.py](https://github.com/qutip/qutip/blob/master/qutip/wigner.py); [forum thread 2293](https://discuss.pennylane.ai/t/wigner-calculation-unstable-for-the-fock-backend/2293)
- GKP and cat states are physically topical (Xanadu's 2025 on-chip GKP source) but numerically demanding. With a 50-photon cutoff, SF's Fock representation "cannot accurately construct the Wigner function for α = 4", which motivated a separate Gaussian-sum formalism. — [Nature 642, 587 (2025)](https://www.nature.com/articles/s41586-025-09044-5); [arXiv:2103.05530](https://arxiv.org/abs/2103.05530)
- GBS's numerical frontier is supercomputer-scale hafnian combinatorics. — [arXiv:1805.12498](https://arxiv.org/abs/1805.12498); [arXiv:2108.01622](https://arxiv.org/abs/2108.01622)

### Inferences
- **Lossy interferometry, ranked first.** It uses exactly the master's-level toolkit: beam splitters, phase shifts, coherent, squeezed and Fock states, loss as a channel or Kraus map, and Fisher information. The answer ("which input state wins at this N and η?") genuinely depends on computation, and closed-form limits give built-in correctness checks. It also has a known pedagogical trap, truncation, that a careful study can make explicit with convergence curves. That is the weakness documented for SF, QuTiP and QuantumOptics.jl in Q2.
- **Kerr evolution and Wigner negativity, ranked second.** It is visually compelling (negativity islands, revival times). Its numerics are sensitive to cutoff and Wigner-evaluation method, and it ties to a landmark experiment. QuTiP already sets a tutorial-level bar, so differentiation would have to come from accuracy and convergence certification.
- **Heralded cat and GKP preparation by PNR detection, ranked third.** It is physically current, but conceptually heavier (Gaussian-state Fock amplitudes, hafnians, optimization). The Walrus gallery and MrMustard papers set a sophisticated research bar, although both tools are now archived.
- **GBS, ranked least suitable as a numerics showcase.** It is the most famous, but credible small-instance demonstrations reduce to hafnian evaluation, and the frontier is owned by specialist kernels. It works for teaching complexity (HOM interference to hafnians), not for "numerics materially determines the answer" at a scale that beats incumbents.

### Gaps
- I found no survey of master's-level quantum-optics curricula or instructor preferences. The ranking is a reasoned inference, not an evidence-based measure of audience appeal.
- I did not verify which lossy-interferometry quantities (QFI versus classical Fisher information for photon counting versus homodyne) have closed forms at finite N. A study would need to settle that before claiming where numerics are indispensable.
