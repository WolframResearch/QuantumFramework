# QEC, stabilizer and large-scale circuit simulation incumbents: Stim + sinter + PyMatching, quimb + cotengra, TensorCircuit(-NG); Qiskit Aer MPS as a comparison point

Scope note: these notes cover only what the incumbents do and how well, as of 2026-10-04. They do not evaluate QF. "Verified" means read in the primary source during this session (paper text, docs page, source file, or PyPI metadata). Lines marked *(summarized)* come from a page that a fetch tool summarized, not from raw text I read myself. Lines marked *(search snippet)* were not opened beyond the search result.

## Which numerical problem classes do these tools feature?

### Takeaway
The Stim, sinter and PyMatching stack features one dominant problem class: Monte Carlo estimates of logical error rates for stabilizer-code memory experiments under circuit-level Pauli noise. That covers thresholds, logical error rate against distance, qubit footprint, and decoders configured automatically from detector error models (DEMs). quimb and cotengra feature tensor-network contraction of large circuits and many-body networks: amplitudes, local expectation values reduced by the reverse light cone, contraction paths for random circuit sampling, MPS/TEBD dynamics (light cone, entanglement growth), and, recently, PEPS/PEPO and belief propagation for 127-qubit "utility" circuits. TensorCircuit and its successor TensorCircuit-NG feature differentiable, GPU-accelerated and distributed variational workflows. Their newer additions are qudits, a Stim-backed stabilizer engine (used for measurement-induced phase transitions), MPS, analog dynamics and fermionic Gaussian states. Qiskit Aer's MPS method is a general-purpose simulator for low-entanglement circuits.

### Cited Findings

#### Versions in force (PyPI upload dates, checked 2026-10-04)
- **stim 1.16.0** (2026-05-22). Earlier: 1.15.0 (2025-05-07), 1.14.0 (2024-09-24), 1.13.0 (2024-03-18), 1.12.0 (2023-08-22). Dev builds 1.17.dev* were uploaded through 2026-09-22. **sinter** and **stimcirq** ship in lockstep, also at 1.16.0 — [PyPI stim](https://pypi.org/project/stim/), [PyPI sinter](https://pypi.org/project/sinter/). (A summarized fetch of the GitHub releases page gave "May 7, 2026" for v1.15.0. PyPI shows 2025-05-07, so the summary date is wrong.)
- **PyMatching 2.4.0** (2026-05-22). 2.3.0 (2025-08-10) introduced correlated matching; 2.3.1 (2025-09-25) — [PyPI PyMatching](https://pypi.org/project/PyMatching/).
- Other decoders sinter can call: **fusion-blossom 0.2.13** (2025-02-03), **ldpc 2.4.1** (2025-12-08, BP+OSD family), **beliefmatching 0.2.0** (2025-06-28) — [PyPI fusion-blossom](https://pypi.org/project/fusion-blossom/), [PyPI ldpc](https://pypi.org/project/ldpc/), [PyPI beliefmatching](https://pypi.org/project/beliefmatching/).
- **quimb 1.15.0** (2026-08-10); **cotengra 0.8.2** (2026-06-22) — [PyPI quimb](https://pypi.org/project/quimb/), [PyPI cotengra](https://pypi.org/project/cotengra/).
- **tensorcircuit 0.12.0** (2024-03-15). This is the last release of the original package; its repository is archived and points to TensorCircuit-NG — [PyPI tensorcircuit](https://pypi.org/project/tensorcircuit/), [archived repo](https://github.com/refraction-ray/tensorcircuit). **tensorcircuit-ng 1.10.0** was released 2026-10-02. Releases are frequent: 1.6.0 (2026-03-30), 1.7.0 (06-14), 1.8.0 (07-18), 1.9.0 (08-11), 1.10.0 (10-02) — [PyPI tensorcircuit-ng](https://pypi.org/project/tensorcircuit-ng/).
- **qiskit-aer 0.17.2** (2025-09-17); the online docs show 0.17.1 — [PyPI qiskit-aer](https://pypi.org/project/qiskit-aer/), [AerSimulator docs](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html).

#### Stim, sinter and PyMatching (QEC stack)
- Stim describes itself as a tool "for high-performance simulation and analysis of quantum stabilizer circuits", aimed at QEC. Its headline features are fast sampling through `compile_sampler()`, automatic decoder configuration through `stim.Circuit.detector_error_model()`, and stabilizer utilities (`PauliString`, `Tableau`, `TableauSimulator`) — [Stim README](https://github.com/quantumlib/Stim).
- The official getting-started notebook defines the canonical workflow, in this order: build a circuit; add detector annotations; generate example QEC circuits; decode with PyMatching; estimate a **repetition-code threshold** by Monte Carlo; streamline with sinter; estimate a **surface-code threshold and footprint** — [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
- `stim.Circuit.generated` ships these built-in code tasks: `repetition_code:memory`, `surface_code:rotated_memory_x`, `surface_code:rotated_memory_z`, `surface_code:unrotated_memory_x`, `surface_code:unrotated_memory_z`, `color_code:memory_xyz`. The noise knobs are `after_clifford_depolarization`, `before_round_data_depolarization`, `before_measure_flip_probability` and `after_reset_flip_probability` — [Stim Python API reference](https://github.com/quantumlib/Stim/blob/main/doc/python_api_reference_vDev.md), [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
- The notebook shows a DEM for a d=9 repetition code with `before_round_data_depolarization=0.04` and `before_measure_flip_probability=0.01`. It contains lines such as `error(0.0266667) D0 D1` and `error(0.01) D0 D8` — [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
- Stim's noise instruction set, from the gate reference on main: `DEPOLARIZE1`, `DEPOLARIZE2`, `E`/`CORRELATED_ERROR`, `ELSE_CORRELATED_ERROR`, `HERALDED_ERASE`, `HERALDED_PAULI_CHANNEL_1`, `I_ERROR`, `II_ERROR`, `PAULI_CHANNEL_1`, `PAULI_CHANNEL_2`, `X_ERROR`, `Y_ERROR`, `Z_ERROR`. Collapsing operations are M/MX/MY, MR*, R*, pair measurements `MXX`/`MYY`/`MZZ`, and `MPP`. Feedback takes the form `CX rec[-k] q` — [Stim gates.md](https://github.com/quantumlib/Stim/blob/main/doc/gates.md).
- Features by release *(summarized from the GitHub releases page)* — [Stim releases](https://github.com/quantumlib/Stim/releases):
  - v1.12 (2023): heralded erasure (`HERALDED_ERASE`, `HERALDED_PAULI_CHANNEL_1`), `MXX`/`MYY`/`MZZ`, `stim.FlipSimulator`.
  - v1.13 (2024): stabilizer flows (`stim.Flow`, `has_flow`, `detecting_regions`), `SPP`, and `sinter.FusionBlossomCompiledDecoder`.
  - v1.14 (2024): `flow_generators`, `time_reversed_for_flows`, and `sinter.Sampler` for generalized Monte Carlo.
  - v1.15 (2025): free-form instruction **tags** (for example `CX[adiabatic] 0 1`) and Pauli targets in `OBSERVABLE_INCLUDE`.
  - v1.16 (2026): the "stimflow" glue library for building QEC circuits, `stim.CliffordString`, `missing_detectors`, and loop-folding acceleration in `stim sample`.
- sinter "takes Stim circuits annotated with noise, detectors, and logical observables", samples them with Stim, decodes them, and records logical-failure statistics. It runs in parallel over cores with multiprocessing, offers `max_shots`/`max_errors` stopping, resumable collection (`--save_resume_filepath`), and `sinter.plot_error_rate` with uncertainty bands — [sinter README](https://github.com/quantumlib/Stim/blob/main/glue/sample/README.md), [sinter API](https://github.com/quantumlib/Stim/blob/main/doc/sinter_api.md).
- sinter's built-in decoder table on main: `vacuous`, `pymatching`, `pymatching-correlated`, `fusion_blossom`, `hypergraph_union_find`, `mw_parity_factor`, `perfectionist`. Custom decoders plug in through `custom_decoders` (`sinter.Decoder` / `sinter.Sampler`) — [sinter built-in decoders source](https://github.com/quantumlib/Stim/blob/main/glue/sample/src/sinter/_decoding/_decoding_all_built_in_decoders.py), [sinter API](https://github.com/quantumlib/Stim/blob/main/doc/sinter_api.md).
- PyMatching implements minimum-weight perfect matching (MWPM) with the exact "sparse blossom" algorithm. It loads Stim DEMs, parity-check matrices, or NetworkX/rustworkx graphs, and supports weighted graphs with or without boundaries and batch decoding. v2.3+ adds **two-pass correlated matching** (`enable_correlations=True`, sparse blossom run twice with reweighting) — [PyMatching README](https://github.com/oscarhiggott/PyMatching), [PyMatching releases](https://github.com/oscarhiggott/PyMatching/releases).
- Featured real-world application: Google's Willow below-threshold experiment. The real-time decoder used "a specialized version of the Sparse Blossom algorithm". The neural-network decoder was pretrained on synthetic SI1000 (p=0.4%) data "generated by Stim". The released QEC datasets include processing examples "with open-source tools such as Stim and PyMatching" — [Google Quantum AI, arXiv:2408.13687 / Nature 638, 920 (2025)](https://arxiv.org/abs/2408.13687).

#### quimb and cotengra
- quimb's circuit guide is built around `Circuit`, an exact tensor network of arbitrary geometry. It offers:
  - `amplitude`;
  - `local_expectation` with "reverse lightcone" cancellation, so that only gates in the causal cone are kept;
  - `partial_trace` for reduced density matrices;
  - three sampling schemes: marginal-based, `sample_gate_by_gate`, and `sample_chaotic`;
  - `*_rehearse` methods that return the contraction width W and the cost C (log10 FLOPs) before contracting;
  - a simplification pipeline with default `'ADCRS'`.
  The classes `CircuitMPS`, `CircuitPermMPS` and `CircuitDense` are also documented — [quimb circuit guide](https://quimb.readthedocs.io/en/latest/circuit/circuit.html) *(summarized)*.
- quimb's example gallery includes "MPS Evolution with TEBD", "Real time simple update", "Random Unitary Evolution", "Quenching Random Product State" (exact), "Periodic DMRG", "MERA", "TRG", "Circuit Sampling Methods", "Training Quantum Circuits", "Bayesian Optimizing QAOA Circuit Energy" and "Converting Circuit to MPO" — [quimb docs index](https://quimb.readthedocs.io/en/latest/) *(summarized)*.
- The TEBD example evolves the isotropic Heisenberg chain from a product state with two flipped spins, L=44 and t∈[0,80]. It plots local ⟨Sz⟩ (a ballistic **light cone**), **block entanglement entropy growth**, and the Schmidt gap, and compares with an MBL variant — [quimb TEBD example](https://quimb.readthedocs.io/en/latest/examples/ex_TEBD_evo.html) *(summarized)*.
- quimb 1.15.0 (2026-08-10) added several circuit engines: `CircuitPEPSSimpleUpdate` (PEPS of arbitrary geometry), `CircuitPEPOSimpleUpdate` (a Heisenberg-picture, backward PEPO evolution), `CircuitMPSLazy`, a shared `CircuitBase`, OpenQASM 3 parsing, an infinite-PEPS module, and gauging by belief propagation (`gauge_d2bp`) — [quimb changelog](https://quimb.readthedocs.io/en/latest/changelog.html) *(summarized)*.
- cotengra is a drop-in replacement for einsum/ncon on large networks. It provides:
  - a hyper-optimizer that samples contraction trees and tunes meta-parameters (greedy, random-greedy, simulated annealing, KaHyPar partitioning);
  - dynamic **slicing**, subtree reconfiguration, and explicit `ContractionTree` objects;
  - hyper-edge and arbitrary einsum support;
  - backend-agnostic contraction through autoray.
  — [cotengra docs](https://cotengra.readthedocs.io/en/latest/) *(summarized)*.
- cotengra 0.8.0 (May 2026) made subtree reconfiguration the default in `HyperOptimizer`, turned on parallelism in the `"auto"`/`"auto-hq"` presets, and added peak-memory tree reordering, a benchmarking example, and a `pytblis` contraction backend — [cotengra changelog](https://github.com/jcmgray/cotengra/blob/main/docs/changelog.md).
- Featured application: Begušić, Gray and Chan simulated the 127-qubit kicked-Ising "utility" observables on IBM Eagle with sparse Pauli dynamics and tensor networks. They report that "All tensor network simulations were performed using quimb and cotengra" — [Begušić, Gray, Chan, arXiv:2308.05077 / Sci. Adv. 10, eadk4321 (2024)](https://arxiv.org/abs/2308.05077).

#### TensorCircuit and TensorCircuit-NG
- The original TensorCircuit (2023 paper) is "based on tensor network contraction". It supports automatic differentiation (AD), JIT compilation, vmap and GPU, is "especially suited to variational algorithms", and offers `tc.DMCircuit` (density matrix), Monte Carlo trajectory noise with `tc.Circuit`, `cond_measure`/`conditional_gate`/`post_select` for mid-circuit logic, MPS inputs, Qiskit import/export, and integration with the cotengra path finder — [TensorCircuit, arXiv:2205.10091 / Quantum 7, 912 (2023)](https://arxiv.org/abs/2205.10091).
- TensorCircuit-NG (paper Feb 2026) bills itself as "universal, composable, scalable". Its engines are `Circuit`, `DMCircuit`, `MPSCircuit`, `StabilizerCircuit` (backed by Stim), `QuditCircuit` (d≥3), `AnalogCircuit` (time-dependent Hamiltonians), `FGSCircuit` (fermionic Gaussian states) and `U1Circuit`. Backends are JAX, TensorFlow, PyTorch, NumPy and cupy. It adds distributed data parallelism and model-parallel slicing of tensor-network contractions — [TensorCircuit-NG, arXiv:2602.14167](https://arxiv.org/abs/2602.14167), [tensorcircuit-ng README](https://github.com/tensorcircuit/tensorcircuit-ng) *(README summarized)*.
- TC-NG's featured stabilizer application is the **measurement-induced phase transition (MIPT) in random Clifford circuits**. The workflow is random 2-qubit Cliffords, projective measurement with probability p, and half-chain entanglement from the tableau, swept over p and L to find the crossing p_c — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167).
- TC-NG also features quench dynamics with fermionic Gaussian states ("melting of a Néel state", entanglement asymmetry, subsystem information capacity), time evolution by ODE/Krylov/Trotter/Chebyshev, and noise models with readout error mitigation — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167), [tensorcircuit-ng README](https://github.com/tensorcircuit/tensorcircuit-ng).

#### Qiskit Aer MPS (comparison point only)
- `method="matrix_product_state"` is "a tensor-network statevector simulator that uses a Matrix Product State (MPS) representation… with or without truncation… The default behaviour is no truncation." Its options are `matrix_product_state_max_bond_dimension` (default None), `matrix_product_state_truncation_threshold` (default 1e-16), `mps_sample_measure_algorithm`, `mps_log_data`, `mps_swap_direction`, `mps_omp_threads`, `mps_lapack`, and others — [AerSimulator docs (0.17.1)](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html).
- The MPS state implementation accepts these operation types: `kraus` (stochastic per-shot application), `measure`, `reset`, `initialize`, classical control (`bfunc`/`jump`/`mark`), `roerror`, `save_densmat`/`save_expval`/`save_mps`/`save_probs`. It has **no** `superop` — [Aer MPS source](https://github.com/Qiskit/qiskit-aer/blob/main/src/simulators/matrix_product_state/matrix_product_state.hpp).

### Inferences
- Stim releases from 2024 to 2026 have gone mostly into constructing and verifying circuits (flows, tags, stimflow, `missing_detectors`), not into new physics (non-Pauli noise, qudits). The physics envelope of the Stim stack is therefore stable: Clifford circuits plus Pauli or heralded-Pauli noise. Users who need non-Pauli noise or leakage move to forks or new tools (see the limitations question).
- quimb's 2026 direction (PEPO in the Heisenberg picture, BP gauging, PEPS circuits) follows the "utility-experiment" problem class: local observables of large, shallow-to-moderate-depth 2D circuits. Its stated limit is that exact contraction is exponential in treewidth.
- TensorCircuit-NG competes on breadth and on GPU/AD throughput rather than on QEC depth. Its stabilizer engine is a front-end to Stim, so it adds no QEC capability beyond Stim's.
- Problem classes the user listed that are directly featured by an incumbent: threshold and logical error rate against distance (Stim notebook); DEMs and matching (Stim/PyMatching); RCS path optimization (cotengra paper); TEBD light cone and entanglement growth (quimb example); large-circuit expectation values (quimb `local_expectation`, Begušić et al.); quenches (quimb exact example, TC-NG FGS); MIPT (TC-NG). I found no Lieb-Robinson *bound* demonstration as such in any of the tools; the light cone appears only as an emergent feature of the TEBD example.

### Gaps
- The quimb circuit guide and TEBD example were read through a summarizing fetch, not raw text. Exact cell outputs (for example the TEBD runtime) should be re-checked before head-to-head use.
- I did not find an official list of TensorCircuit-NG tutorials by topic (for example whether a Lieb-Robinson or QEC tutorial exists). The README claims "100+ example scripts, 40+ tutorial notebooks" without topic detail in what I read.

## What scale and performance are published?

### Takeaway
The published bars are as follows:
- **Stim:** analyses a distance-100 surface code circuit (≈20k qubits, 8M gates, 1M measurements) in 15 s, then samples full shots at about 1 kHz (2021 laptop numbers).
- **Sparse blossom:** decodes d=17 surface code circuits at 0.1% noise in under 1 µs per round per core, and d=29 in 3.5 µs per round.
- **cotengra:** found paths giving a >10,000× cheaper Sycamore simulation than Google's original estimate. Measured single-amplitude contractions run from sub-second (Bristlecone/7×7 at depth 32) to an extrapolated ~10^11 s for Sycamore m=20 on a 5 GB consumer GPU (2021). It reached under 10 minutes per depth-20 Sycamore sample on NVIDIA's Selene supercomputer.
- **TensorCircuit:** 600-qubit, 7-layer 1D VQE at 18 s per step on an A100 (2023).
- **TensorCircuit-NG:** exact differentiable 40-qubit, 20-layer circuits at ~18 min per step on 8×H200 (2026).
- **Qiskit Aer MPS:** only a 50-qubit GHZ example in 0.31 s (old tutorial).

### Cited Findings

#### Stim
- From the paper abstract: "With no foreknowledge, Stim can analyze a distance 100 surface code circuit (20 thousand qubits, 8 million gates, 1 million measurements) in 15 seconds and then begin sampling full circuit shots at a rate of 1 kHz." The three techniques credited are linear-time deterministic measurement (tracking the inverse tableau, instead of quadratic CHP), cache-friendly 256-bit SIMD, and batched Pauli-frame sampling against a single reference sample — [Gidney, Stim, arXiv:2103.02202 / Quantum 5, 497 (2021)](https://arxiv.org/abs/2103.02202). **Older (2021) information.**
- Benchmark hardware in the paper was a ThinkPad laptop (Intel i7-8650U @ 1.90 GHz, 16 GB). The comparison set was chp (Aaronson-Gottesman), Qiskit `method='stabilizer'`, `cirq.CliffordSimulator` and graphsim. The conclusion: "at larger sizes Stim outperforms everything else by one or more orders of magnitude, and when collecting thousands of samples Stim outperforms everything else by many orders of magnitude". The exception is multi-level distillation, "where graphsim shines" — [arXiv:2103.02202](https://arxiv.org/abs/2103.02202).
- The README claims "a few seconds of analysis and then… sample shots at kilohertz rates" for "thousands of qubits and millions of operations", and that "Stim can multiply Pauli strings with 100 billion terms in one second" (256-bit AVX) — [Stim README](https://github.com/quantumlib/Stim).
- The paper's stated goal was not reached: "my attempt at simulating surface code circuits with tens of thousands of qubits and millions of operations in a tenth of a second… that goal wasn't reached" — [arXiv:2103.02202](https://arxiv.org/abs/2103.02202).

#### sinter and PyMatching
- sinter has been "tested… on 2 core machines, 4 core machines, and 96 core machines… generally achieves good resource utilization". It publishes no throughput figures — [sinter README](https://github.com/quantumlib/Stim/blob/main/glue/sample/README.md). v1.12 "doubled" sinter performance on high-core machines — [Stim releases](https://github.com/quantumlib/Stim/releases) *(summarized)*.
- Shot budgets used in the official notebook:
  - repetition-code sweep: 10,000 shots per point, done by hand;
  - sinter repetition-code sweep: `max_shots=100_000`, `max_errors=500`;
  - surface-code threshold sweep: `max_shots=1_000_000`, `max_errors=5_000`;
  - footprint sweep at p=0.1%: `max_shots=5_000_000`, `max_errors=100`.
  — [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb). The notebook pins `sinter~=1.14`, so its examples date from late 2024.
- Sparse blossom processes "both X and Z bases of distance-17 surface code circuits in less than one microsecond per round of syndrome extraction on a single core" at 0.1% circuit-level depolarizing noise. At "distance 29 with the same noise model (more than sufficient to achieve 10^-12 logical error rates), PyMatching takes 3.5 microseconds per round". Fitted runtime scaling against node count N for comparators: NetworkX ∝ N^2.15 and PyMatching v0.7 ∝ N^1.19 at p=0.1%; N^2.77 and N^1.21 in a second, noisier benchmark — [Higgott & Gidney, Sparse Blossom, arXiv:2303.15933 / Quantum 9, 1600 (2025)](https://arxiv.org/abs/2303.15933).
- The PyMatching README claims v2 is "100-1000x faster than previous versions" and "over 100,000x faster than NetworkX", with runtime roughly linear in the number of graph nodes — [PyMatching README](https://github.com/oscarhiggott/PyMatching) *(summarized)*.
- Hardware context: Willow's real-time decoder (a specialized Sparse Blossom) achieved "an average decoder latency of 63 µs at distance-5 up to a million cycles, with a cycle time of 1.1 µs". Λ = 2.14 ± 0.02 for the d=7 memory (0.143% error per cycle) — [arXiv:2408.13687](https://arxiv.org/abs/2408.13687).

#### Non-Clifford extensions of the Stim ecosystem (for scale context)
- Tsim (Apr 2026) samples "in linear time in the number of Clifford gates and exponentially only in the number of non-Clifford gates", implements the Stim API and circuit format with added T and arbitrary rotations, uses GPU acceleration, and "for low-magic circuits… can match the sampling performance of Stim" — [Haenel, Luo, Zhao, arXiv:2604.01059](https://arxiv.org/abs/2604.01059).
- Simulating d=3 magic-state-cultivation circuits reached "nearly 4×10^6 shots per second… on a laptop" at p=0.0005, about 1.1× slower than a Stim simulation of the all-Clifford proxy (T replaced by S). d=5 needs only about 8 Clifford ZX terms on average — [Wan, Zhong, Zapirain, Quantum 10, 2134 (2026)](https://quantum-journal.org/papers/q-2026-06-12-2134/) *(summarized)*.

#### cotengra and quimb
- Gray & Kourtis, abstract: "we estimate a speed-up of over 10,000× compared to the original expectation for the classical simulation of the Sycamore 'supremacy' circuits" — [arXiv:2002.01935 / Quantum 5, 410 (2021)](https://arxiv.org/abs/2002.01935). **Older (2021).**
- Measured single-amplitude contractions (Table 2) used single precision on an NVIDIA Quadro P2000 (5 GB, 3.031 TFLOPs), paths sliced to width W_s=27, quimb for the tensor network and JAX for the contraction — [arXiv:2002.01935](https://arxiv.org/abs/2002.01935):

  | Circuit | Time | Sliced cost C_s (FLOPs) | Slicing overhead |
  |---|---|---|---|
  | Bristlecone-70 (1+32+1) | 0.418 s | 4.91×10^10 | — |
  | Bristlecone-70 (1+40+1) | 277 s | 3.14×10^13 | — |
  | Rectangular-7×7 (1+48+1) | ≈9.4×10^4 s (est.) | 1.20×10^16 | — |
  | Sycamore-53, m=12 | 574 s | 1.80×10^14 | — |
  | Sycamore-53, m=20 | ≈9.74×10^10 s (est.) | 3.10×10^22 | 6410× |
  | Sycamore-53*, m=20 (χ=2 swapped fSim decomposition) | ≈7.17×10^9 s (est.) | 1.50×10^21 | 431× |

  FLOP efficiency reached up to ~85% of peak.
- Slicing trade-off: for Sycamore m=20, "the required memory can be brought down by a factor of ∼16,000 whilst keeping the FLOPs increase < 10"; "W_S ∼ 27 is required to fit a contraction on a standard consumer GPU" — [arXiv:2002.01935](https://arxiv.org/abs/2002.01935).
- With NVIDIA, cotengra achieved "state-of-the-art simulation of the Sycamore chip with cutensor on the Selene supercomputer, producing a sample from a circuit of depth 20 in less than 10 minutes" — [cotengra docs](https://cotengra.readthedocs.io/en/latest/).
- quimb guide example: an **80-qubit GHZ** circuit (randomly ordered CNOTs) took 449 ms for the first sample and ~105 ms for later batches of 8 (cached paths). Simplification reduced a 668-tensor network to 87 tensors — [quimb circuit guide](https://quimb.readthedocs.io/en/latest/circuit/circuit.html) *(summarized)*.
- quimb TEBD example: L=44, t up to 80, 101 output times, cutoff 1e-12, Trotter error target 1e-3 (estimated 9.94×10^-4). Maximum bond dimension reached 15 for Heisenberg and 19 for the MBL variant; energy 8.75 was conserved; runtime ≈2 min — [quimb TEBD example](https://quimb.readthedocs.io/en/latest/examples/ex_TEBD_evo.html) *(summarized; runtime is hardware-dependent)*.
- Utility benchmark (127 qubits, kicked Ising): the most accurate method combined a mixed Schrödinger/Heisenberg tensor network with the BP Bethe free entropy, at "effective wavefunction-operator sandwich bond dimension >16,000,000" (16,777,216). It achieved "absolute accuracy, without extrapolation, in the observables of <0.01", "orders of magnitude faster than the quantum experiment" — [arXiv:2308.05077](https://arxiv.org/abs/2308.05077).

#### TensorCircuit and TensorCircuit-NG
- The original paper claims "up to 600 qubits with moderate circuit depth and low-dimensional connectivity". Table 9 runs a 1D critical TFIM VQE with 7 ladder layers (split ZZ gates at bond dimension 2, MPO Hamiltonian, cotengra with subtree reconfiguration) on an A100 40G with TensorFlow — [arXiv:2205.10091](https://arxiv.org/abs/2205.10091). **Older (2023).**

  | Qubits n | Time per step (energy + gradient) | Accuracy |
  |---|---|---|
  | 200 | 5.7 s | 99.6% |
  | 400 | 11.8 s | 99.5% |
  | 600 | 18.2 s | 99.4% |

- Table 6 times the value plus gradient of a TFIM hardware-efficient ansatz (JAX, V100 GPU) for three circuit sizes — [arXiv:2205.10091](https://arxiv.org/abs/2205.10091):

  | Framework | n=10, d=3 | n=16, d=16 | n=22, d=11 |
  |---|---|---|---|
  | TensorCircuit (GPU) | 0.0026 s | 0.023 s | 0.19 s |
  | TensorFlow Quantum | 0.005 s | 0.026 s | 0.68 s |
  | PennyLane (CPU) | 0.012 s | 0.31 s | 24.84 s |
  | Qibo (GPU) | 0.033 s | 0.198 s | OOM |

- Table 5 compares contraction paths for n=40, d=6: the default greedy contractor gives log10 FLOPs 7.373, while cotengra with subtree reconfiguration gives about 7.000 — [arXiv:2205.10091](https://arxiv.org/abs/2205.10091).
- TC-NG distributed results on 8× NVIDIA H200 (141 GB each, NVLink) — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167):
  - 32-qubit, 16-layer TFIM VQE: 17.86 s → 2.38 s per step (7.5× on 8 GPUs);
  - scaling "up to N = 40 qubits and L = 20 layers", 11,700 parameters, "approximately 18 minutes per step";
  - runtime grows as T ∝ 2^(1.1N) (rendered "21.1N" in the extracted PDF text);
  - slicing `target_size` was 2^29 elements, with cotengra `max_repeats=640`.
- TC-NG MIPT on a 20-qubit, 40-layer Haar-random circuit (complex64, batch 1000): 84 ms per trajectory on an H200 against 1.097 s on an Apple M4 Pro CPU. Sparse TFIM Hamiltonian construction at L=24: 0.059 s on an H200 — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167).
- Vendor claims in the TC-NG README: "10 to 10^6+ times acceleration compared to TensorFlow Quantum, Pennylane, CuQuantum, TorchQuantum or Qiskit" and a "1000+ qubits 1D VQE" converging within 0.5% error — [tensorcircuit-ng README](https://github.com/tensorcircuit/tensorcircuit-ng) *(summarized; not independently verified)*.

#### Qiskit Aer MPS
- The tutorial example is a 50-qubit GHZ state in 0.31 s — [Aer MPS tutorial](https://qiskit.github.io/qiskit-aer/tutorials/7_matrix_product_state_method.html) *(summarized)*. The tutorial outputs appear to come from an old Aer version (0.8.x), so treat the timing as **older information**.
- MPS has no GPU support in Aer's GPU table: GPU is "No" for stabilizer, `matrix_product_state` and `extended_stabilizer`, and the separate `tensor_network` method is "GPU only" through cuTensorNet — [AerSimulator docs](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html).

### Inferences
- The incumbent QEC bar has two parts: sampling throughput (kHz full shots at 10^4 qubits; ≥10^6 shots per data point is routine in the notebook) and decoding throughput (µs per round). A head-to-head on *large-distance* surface-code thresholds would be judged on these, so it plays to the incumbents' strengths.
- The tensor-network bar is set less by "number of qubits" than by contraction width and cost (log2 of the largest tensor, log10 FLOPs). Any head-to-head should report W and C (as quimb's `rehearse` does), not wall-clock alone, because the published wall-clocks are tied to specific GPUs (Quadro P2000, A100, H200).
- Several headline numbers (Stim 2021, cotengra 2021, TensorCircuit 2023) are 3–5 years old on old hardware. Current versions are probably faster, but I found no refreshed official benchmarks.

### Gaps
- No 2025–2026 official Stim throughput benchmark was found (only the qualitative "kilohertz" in the README), and no sinter throughput figure exists.
- The cotengra 0.8.0 "benchmarking example" (`ex-benchmarking`) was not opened. Its numbers would be the most current cotengra performance data.
- Gray & Chan's "Hyper-optimized compressed contraction of tensor networks with arbitrary geometry" (arXiv:2206.07044, PRX 2024) was not retrieved, so its compressed-contraction performance claims are unverified here.
- I found no Aer MPS benchmark newer than the old tutorial.

## What limitations are documented (non-Clifford gates, non-Pauli noise, leakage, qudits, mid-circuit feedback, correlated noise)?

### Takeaway
Stim documents three hard limits: no non-Clifford gates, only Pauli noise channels ("no amplitude decay"), and only single-control Pauli feedback. Its Pauli-frame method cannot represent relaxation or leakage. In 2025–2026, leakage, non-Pauli noise and T gates were handled outside Stim: Riverlane's Deltakit-Stim fork, Tsim, approximate thermal-relaxation schemes, and Clifford-proxy tricks for magic-state cultivation. PyMatching requires graphlike (weight ≤2) error mechanisms. quimb has no built-in noise channels, per a 2023 maintainer reply on an issue that is still open. TensorCircuit-NG supports noise (density matrix and Monte Carlo trajectories) and qudits, but a noisy tensor-network engine is only proposed. Aer MPS handles Kraus noise only through stochastic trajectories and has no GPU.

### Cited Findings

#### Stim (verbatim limits and supporting sources)
- The README states the main limitations:
  1. "There is no support for non-Clifford operations, such as T gates and Toffoli gates. Only stabilizer operations are supported."
  2. "`stim.Circuit` only supports Pauli noise channels (eg. no amplitude decay). For more complex noise you must manually drive a `stim.TableauSimulator`."
  3. "`stim.Circuit` only supports single-control Pauli feedback. For multi-control feedback, or non-Pauli feedback, you must manually drive a `stim.TableauSimulator`."
  — [Stim README](https://github.com/quantumlib/Stim).
- From the paper: "In a Pauli frame simulation, noise processes must be Pauli channels… dephasing and depolarization are Pauli channels but relaxation to the ground state and leaking outside the computational basis aren't" — [arXiv:2103.02202](https://arxiv.org/abs/2103.02202).
- The maintainer (Strilanc) wrote in Dec 2025: "Stim doesn't support simulating non-Clifford gates." — [Stim issue #1017](https://github.com/quantumlib/Stim/issues/1017).
- **Leakage and erasure inside Stim:** only heralded erasure (`HERALDED_ERASE`: "a 1 is recorded and the target qubit is erased to the maximally mixed state by applying X_ERROR(0.5) and Z_ERROR(0.5)") and `HERALDED_PAULI_CHANNEL_1`. When converted to a DEM, the channel "is split into multiple potential effects". `I_ERROR`/`II_ERROR` "have no effect… [they] may be interpreted in arbitrary ways by external tools", acting as hooks for custom noise in downstream tools — [Stim gates.md](https://github.com/quantumlib/Stim/blob/main/doc/gates.md).
- **Correlated noise:** Stim supports correlated *Pauli* errors (`E`/`CORRELATED_ERROR`, `ELSE_CORRELATED_ERROR`, `PAULI_CHANNEL_2`, `DEPOLARIZE2`). It has nothing for correlated non-Pauli or coherent noise, consistent with the Pauli-only limit — [Stim gates.md](https://github.com/quantumlib/Stim/blob/main/doc/gates.md), [Stim README](https://github.com/quantumlib/Stim).
- **Qudits:** the gate reference is qubit-only, and a keyword search of the Stim issue tracker for "qudit" returned no issues — [Stim gates.md](https://github.com/quantumlib/Stim/blob/main/doc/gates.md), [Stim issues](https://github.com/quantumlib/Stim/issues).
- **Feedback:** open issues include "Add `stim.Circuit.with_feedback_ensuring_flows`" (#860, 2024) and "Fix time_reversed_for_flows handling feedback incorrectly" (#764, 2024) — [Stim issue #860](https://github.com/quantumlib/Stim/issues/860), [Stim issue #764](https://github.com/quantumlib/Stim/issues/764).

#### Workarounds outside Stim (2025–2026)
- **Deltakit-Stim** (Riverlane) is a Stim fork that adds `LEAKAGE(p)`, `RL` (leakage reset) and `HERALD_LEAKAGE_EVENT`, plus an "adaptive" DEM carrying metadata for "non-computational errors (leakage, erasure, atom loss)". Riverlane claims its Local Clustering Decoder "can reduce physical qubit requirements by a factor of 4 compared to standard non-adaptive decoding" — [Deltakit-Stim GitHub](https://github.com/Deltakit/deltakit-stim) *(summarized)*, [Riverlane announcement](https://www.riverlane.com/news/introducing-deltakit-stim) *(search snippet)*.
- An exact thermal-relaxation channel "is inherently non-Pauli and cannot be represented as a stochastic Pauli channel, so approximations are required to incorporate it into the current stabilizer-based QEC simulation toolchain" — [arXiv:2512.09189](https://arxiv.org/html/2512.09189) *(search snippet)*.
- **Non-Clifford QEC circuits:** the original magic-state-cultivation work simulated most d=5 circuits "by replacing the T(T†) gates entirely with S(S†) gates" (a Clifford proxy) — reported by [Wan, Zhong, Zapirain, Quantum 10, 2134 (2026)](https://quantum-journal.org/papers/q-2026-06-12-2134/) *(search snippet)*. Tsim and stabilizer-rank/ZX methods now target this gap — [arXiv:2604.01059](https://arxiv.org/abs/2604.01059).

#### PyMatching
- It requires graphlike error mechanisms: "each error must produce exactly one or two detection events". Stim's `decompose_errors=True` is the standard workaround. Correlated (hyperedge) information is used only through the v2.3 two-pass reweighting — [PyMatching README](https://github.com/oscarhiggott/PyMatching) *(summarized)*, [PyMatching releases](https://github.com/oscarhiggott/PyMatching/releases).
- Robustness example: a segfault with correlations enabled for pure-Y data noise (`PAULI_CHANNEL_1(0, p, 0)`) was fixed in v2.3.1 — [PyMatching issue #175](https://github.com/oscarhiggott/PyMatching/issues/175).

#### quimb and cotengra
- Noise: asked about noisy simulation, the maintainer replied "there are not built in ways to do this… You could for instance add some unitary noise simply by making the parameters noisy" (May 2023; the issue is still open). **Older information.** It is not verified whether later versions added channels — [quimb issue #182](https://github.com/jcmgray/quimb/issues/182).
- Exponential cost: the docs warn of an "exponential slow-down expected for generic quantum circuits" (only the prefactor can be reduced) and say "you are unfortunately quite unlikely to achieve the best performance without some tweaking" — [quimb circuit guide](https://quimb.readthedocs.io/en/latest/circuit/circuit.html) *(summarized)*.
- Slicing overhead can be severe at the frontier: 6410× for Sycamore-53 m=20 when squeezed into 5 GB. FLOP efficiency drops on hyper-edge networks ("the performance extracted from the GPU via JAX is not great… hyper-edges resulting in pairwise contractions that do not dispatch to matrix-matrix multiplication") — [arXiv:2002.01935](https://arxiv.org/abs/2002.01935).
- A quimb issue about mid-circuit measurement, "measurement didn't collapse the state correct" (#344, Jan 2026, closed), indicates measurement/collapse support exists but is a newer, less-exercised path *(title only; content not read)* — [quimb issue #344](https://github.com/jcmgray/quimb/issues/344).

#### TensorCircuit and TensorCircuit-NG
- Density-matrix simulation via `tc.DMCircuit` costs as much as "twice as many qubits" of pure-state simulation. The alternative is Monte Carlo trajectories with `tc.Circuit` (general Kraus, plus vmap-batched trajectories). Mid-circuit measurement and feedback use `cond_measure`/`conditional_gate`; `post_select` returns unnormalized states — [arXiv:2205.10091](https://arxiv.org/abs/2205.10091).
- TC-NG states that its stabilizer formalism allows "thousands of qubits", but generic phase transitions (for example MIPT in Haar-random circuits) require full simulation. The MPS engine "is inherently approximate" — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167).
- An open issue (#15, June 2025) proposes a "Locally Purified Density Operator Engine for Noisy Quantum Circuit Simulation", so noisy *tensor-network* simulation is not yet a built-in engine. Qudit support arrived in 2025 (#10, closed Feb 2025). A sampling integer overflow "for large qudit/qubit systems" was fixed in Aug 2025 (#41) — [TC-NG issue #15](https://github.com/tensorcircuit/tensorcircuit-ng/issues/15), [#10](https://github.com/tensorcircuit/tensorcircuit-ng/issues/10), [#41](https://github.com/tensorcircuit/tensorcircuit-ng/issues/41).

#### Qiskit Aer (comparison)
- The `stabilizer` method handles noisy Clifford circuits only "if all errors in the noise model are also Clifford errors". `extended_stabilizer` is approximate for Clifford+T, with cost growing with T count. MPS has no GPU and applies Kraus noise stochastically per shot; it has no superoperator op — [AerSimulator docs](https://qiskit.github.io/qiskit-aer/stubs/qiskit_aer.AerSimulator.html), [Aer MPS source](https://github.com/Qiskit/qiskit-aer/blob/main/src/simulators/matrix_product_state/matrix_product_state.hpp).

### Inferences
- Taken together, these limits leave a structural gap that no incumbent fills natively for QEC. It contains:
  - exact (density-matrix or channel-level) treatment of amplitude damping/T1, coherent over-rotations, and leakage modeled as a genuine third level (a qutrit), on small codes;
  - multi-control or non-Pauli feedback.
  Users currently use Pauli-twirled approximations, forks (Deltakit-Stim), or new tools (Tsim). The size of the error from Pauli-twirling a non-Pauli channel at small distance is a question where numerics decide the answer, and where exact small-scale simulation is the ground truth.
- For tensor networks, the documented failure mode is not "too many qubits" but treewidth and entanglement growth, plus slicing overhead. Problems with bounded light cones (shallow circuits, short-time quenches) are where tensor-network tools shine, and deep generic circuits are where they fail.
- PyMatching's graphlike requirement means Y errors and hyperedges need decomposition or correlated reweighting. Decoder choice (MWPM, correlated MWPM, BP+OSD, union-find, NN) therefore changes published thresholds, which is a source of non-comparability across studies.

### Gaps
- I found no Stim, sinter or PyMatching documentation describing coherent-error support beyond the general "Pauli noise only" statement. The Willow paper injects coherent errors on hardware (line "We inject a range of coherent errors…"), but I did not extract how they were simulated.
- I could not open Deltakit-Stim release and version data, or the base Stim version it forks.
- quimb's current (1.15) support for noise channels and qudits in the `Circuit` classes was not verified; the general tensor-network machinery supports arbitrary index dimensions, but circuit-level qudit gates were not confirmed.
- I found no documented limitation for cotengra itself regarding qudits or noise; it is agnostic to the physics.

## Which standard named benchmark problems does each use (with exact links)?

### Takeaway
The recognised, reproducible benchmarks are:
- **Stim notebook:** the repetition-code memory threshold (phenomenological noise) and the rotated surface-code memory threshold (≈1% under uniform circuit noise), plus the footprint extrapolation at p=0.1%.
- **Noise-model standards:** SD6 and SI1000 (Gidney et al. 2021).
- **Decoder benchmark:** d×d×d surface-code memory at p=0.1% circuit-level depolarizing noise, timed per round.
- **Stim paper suite:** five sampling tasks.
- **Tensor networks:** Sycamore-53 / Bristlecone-70 / Rectangular-7×7 random circuits (single amplitude); the IBM Eagle 127-qubit kicked-Ising utility circuit; the quimb TEBD Heisenberg light-cone example.
- **TensorCircuit:** the 1D critical TFIM VQE at 200–600 qubits.
- **TensorCircuit-NG:** 32–40-qubit distributed VQE and the 20-qubit, 40-layer MIPT.
- **Aer MPS:** the 50-qubit GHZ example.

### Cited Findings

#### QEC benchmarks
1. **Repetition-code memory threshold (phenomenological)**. `stim.Circuit.generated("repetition_code:memory", rounds=3d, distance=d, before_round_data_depolarization=p)`, d∈{3,5,7,(9)}, p∈{0.05…0.5}, decoded with PyMatching. The notebook remarks that the performance looks "amazingly good" partly because of the phenomenological model and because depolarizing errors are Z errors one third of the time, to which the code is immune — [getting_started.ipynb §7–8](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
2. **Rotated surface-code memory threshold under uniform circuit noise**. `"surface_code:rotated_memory_z"`, rounds=3d, d∈{3,5,7}, with `after_clifford_depolarization = after_reset_flip_probability = before_measure_flip_probability = before_round_data_depolarization = p` and p∈{0.008…0.012}. Result: "the threshold of the surface code is roughly 1%" (per-round plot via `failure_units_per_shot_func`) — [getting_started.ipynb §9](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
3. **Surface-code footprint ("trillion rounds")**. At p=0.1%, d∈{3,5,7,9}, fit log error rate against d (`scipy.stats.linregress`: slope −1.19, r = −0.9989). Conclusion: "a distance 20 patch would be sufficient to survive a trillion rounds. That's a surface code with around 800 physical qubits". Stated caveats: one circuit realization, one noise model, one decoder, Z basis only, and extrapolation far beyond the data — [getting_started.ipynb §9](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
4. **SD6 and SI1000 circuit noise models**, defined in [Gidney, Newman, Fowler, Broughton, "A Fault-Tolerant Honeycomb Memory", arXiv:2108.10457 / Quantum 5, 605 (2021)](https://arxiv.org/abs/2108.10457):

   | Element | SD6 | SI1000 |
   |---|---|---|
   | Two-qubit gate | CX(p) | CZ(p) |
   | Single-qubit Clifford | p | p/10 |
   | Initialization (InitZ) | p | 2p |
   | Measurement (MZ) | p | 5p |
   | Idle | p | p/10 |
   | Resonator idle | — | 2p |
   | Cycle | — | ≈1000 ns |

   Reported thresholds, per d-round block:

   | Code | SD6 | SI1000 |
   |---|---|---|
   | Surface code | 0.5–0.7% | 0.3–0.5% |
   | Honeycomb code | 0.2–0.3% | 0.1–0.15% |

   SI1000 at p=0.4% was later used by Google to generate Stim pretraining data for the Willow neural decoder — [arXiv:2408.13687](https://arxiv.org/abs/2408.13687).
5. **Decoder-speed benchmark**: d×d×d rotated surface-code memory, circuit-level depolarizing noise p=0.1% (and a second, noisier setting). Metric: decoding time per round, single core; comparators are PyMatching v0.7 and NetworkX — [arXiv:2303.15933, Figs. 10–12](https://arxiv.org/abs/2303.15933).
6. **The Stim paper's five sampling tasks**, compared against chp, Qiskit stabilizer, Cirq Clifford and graphsim — [arXiv:2103.02202 §6, Figs. 1, 5–8](https://arxiv.org/abs/2103.02202):
   - time to the 1000th sample on d×d×d surface-code memory;
   - time to the first sample on d×d×d unrotated surface-code memory;
   - d×d×d Bacon-Shor memory;
   - random circuits with n qubits and n layers (random H/S/I on each qubit, CNOTs on a random pairing, random-basis measurement of 5% of qubits per layer);
   - nested 7-to-1 S-state distillation (5 levels give 2401 independent 15-qubit pieces).
7. **Real-hardware reference data**: Willow surface-code memories at d=3, 5, 7 (Λ = 2.14 ± 0.02; ε7 = 0.143% per cycle) and repetition codes up to d=29. The datasets are released with Stim/PyMatching processing examples — [arXiv:2408.13687](https://arxiv.org/abs/2408.13687).
8. **Non-Clifford QEC benchmark (emerging)**: magic-state cultivation circuits at d=3 and d=5 (Gidney, Shutty, Jones 2024), now used to benchmark simulators beyond Stim — [Quantum 10, 2134 (2026)](https://quantum-journal.org/papers/q-2026-06-12-2134/), [Tsim arXiv:2604.01059](https://arxiv.org/abs/2604.01059).

#### Tensor-network and large-circuit benchmarks
9. **Random circuit sampling (single amplitude, perfect fidelity)** on Sycamore-53 (m=12…20 cycles), Bristlecone-70 (1+32+1 … 1+40+1) and Rectangular-7×7 (1+32+1 … 1+48+1). Metrics: contraction width W, cost C, sliced cost C_s, wall time. Also contraction width and cost of random regular graphs (k=3,4,5) and lattices, compared with tree-decomposition baselines (QuickBB, FlowCutter) and with optimal widths where known — [Gray & Kourtis, arXiv:2002.01935, Table 2, Figs. 4–10](https://arxiv.org/abs/2002.01935).
10. **IBM Eagle 127-qubit kicked-Ising "utility" circuits** (the heavy-hex lattice experiment of Kim et al., Nature 2023), up to 20 Trotter steps. Converged classical reference values come from quimb and cotengra — [Begušić, Gray, Chan, arXiv:2308.05077](https://arxiv.org/abs/2308.05077).
11. **quimb TEBD light cone**: Heisenberg chain, L=44, two flipped spins, t≤80; observables are ⟨Sz⟩(x,t), block entropy and Schmidt gap; MBL variant included — [quimb ex_TEBD_evo](https://quimb.readthedocs.io/en/latest/examples/ex_TEBD_evo.html).
12. **quimb circuit examples**: 80-qubit GHZ sampling ([circuit guide](https://quimb.readthedocs.io/en/latest/circuit/circuit.html)), QAOA energy with Bayesian optimization ([ex_tn_qaoa_energy_bayesopt](https://quimb.readthedocs.io/en/latest/examples/ex_tn_qaoa_energy_bayesopt.html)), circuit sampling methods ([ex_tn_circuit_sample_explore](https://quimb.readthedocs.io/en/latest/examples/ex_tn_circuit_sample_explore.html)), quench of a random product state (exact; [ex_quench](https://quimb.readthedocs.io/en/latest/examples/ex_quench.html)) — [quimb docs index](https://quimb.readthedocs.io/en/latest/).
13. **TensorCircuit 1D TFIM VQE at the critical point, open boundaries**, n=200/400/600, 7 ladder layers; reference energy from quimb two-site DMRG; script `examples/vqe_extra_mpo.py`. **PQC value-and-gradient benchmark** for (n,d)=(10,3), (16,16), (22,11) against PennyLane, TFQ and Qibo; scripts `examples/vqeh2o_benchmark.py` and `examples/vqetfim_benchmark.py` — [arXiv:2205.10091 §7.4, Tables 6 & 9](https://arxiv.org/abs/2205.10091).
14. **TensorCircuit-NG**: distributed 1D TFIM VQE (N=32…40, L=N/2, SU(4) ladder); MIPT in Haar-random circuits (20 qubits, 40 layers); Clifford MIPT via `StabilizerCircuit` (script `VB_stabilizer_mipt.py`); MPS VQE (`VC_mps_vqe.py`); qudit VQE (`IVB_qudit_simulation.py`) — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167). The TC-NG issue tracker shows a 2026 program of "Reproduce arXiv:…" issues (for example 1709.01662, 2207.05612, 2503.01966, 2510.08344) — [TC-NG issues](https://github.com/tensorcircuit/tensorcircuit-ng/issues).
15. **Qiskit Aer MPS**: 50-qubit GHZ tutorial — [Aer MPS tutorial](https://qiskit.github.io/qiskit-aer/tutorials/7_matrix_product_state_method.html).

### Inferences
- The most "recognised" QEC head-to-head is benchmark 2 or 3 (the Stim notebook's rotated surface code with sinter and PyMatching). It is fully specified in code, so an independent implementation can reproduce exactly the same circuits (via exported `.stim` text) and compare logical error rates point by point. A tougher, more hardware-relevant variant uses SI1000 (benchmark 4).
- The "≈1%" (Stim notebook, uniform noise, per round) and "0.5–0.7%" (SD6, per d-round block) surface-code thresholds differ because of the noise-model definition and the normalization (per round against per block). Any head-to-head must pin both, or the comparison is meaningless.
- For tensor networks, the most recognised named problems are Sycamore RCS (benchmark 9) and the IBM 127-qubit kicked Ising (benchmark 10). Both sit at the research frontier (10^14–10^22 FLOPs; bond dimension 1.6×10^7). The quimb TEBD light cone (benchmark 11) and the TensorCircuit 1D TFIM VQE (benchmark 13) are recognised, reproducible and laptop-to-GPU-scale.

### Gaps
- I did not retrieve the exact Zenodo or dataset URL for Google's released Willow data (the paper cites it as its Ref. [12]).
- I did not retrieve the original Kim et al. (Nature 2023) utility paper or Arute et al. (Nature 2019) Sycamore paper; the benchmark definitions above come from the papers that use them.
- There is no single "official" benchmark suite for quimb, cotengra or TC-NG comparable to Stim's notebook. The examples serve that role informally.

## Which of these problem classes is most compelling to a master's-level audience, and why?

### Takeaway
This section is inference. The strongest classes are those where a single, physically meaningful number comes out only of computation, and where that number moves with modelling choices a QM student understands:
1. **Surface-code threshold and logical error rate against distance**, including how the answer shifts between noise models (uniform, SD6, SI1000) and between Pauli and non-Pauli noise.
2. **Light cone and entanglement growth after a quench** (TEBD/MPS), showing why some dynamics are classically tractable and others are not.
3. **The measurement-induced entanglement transition**, where a critical measurement rate p_c appears only from finite-size numerics.

Random circuit sampling is compelling as a complexity story but is less of a quantum-mechanics lesson.

### Cited Findings
- The Stim tutorial frames thresholds as a purely numerical object: "trying a bunch of physical error rates and code distances… see where the curves cross… That's the threshold physical error rate". It then moves on to footprint: "the threshold isn't the only metric you care about… What you really want to estimate is the quality and corresponding quantity of qubits", ending with "distance 20… around 800 physical qubits" at p=0.1% — [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
- The answer depends on the noise model: surface-code thresholds are 0.5–0.7% (SD6) against 0.3–0.5% (SI1000), and the ranking of honeycomb against surface code changes with the model — [arXiv:2108.10457](https://arxiv.org/abs/2108.10457). The Stim notebook gets "roughly 1%" with its own uniform model — [getting_started.ipynb](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).
- There is a direct experimental anchor: Willow's Λ = 2.14 ± 0.02 per +2 distance, "beyond break-even… by a factor of 2.4 ± 0.3", limited by "rare correlated error events occurring approximately once every hour" — [arXiv:2408.13687](https://arxiv.org/abs/2408.13687).
- Pauli-frame simulators cannot represent relaxation or leakage — [arXiv:2103.02202](https://arxiv.org/abs/2103.02202) — and the field treats exact thermal relaxation as needing approximation inside stabilizer toolchains — [arXiv:2512.09189](https://arxiv.org/html/2512.09189) *(search snippet)*.
- The TEBD example ties a textbook concept (ballistic light cone of local magnetization, linear entanglement growth) to a numerical method whose cost is set by the bond dimension the entanglement requires: 15 for Heisenberg, 19 for MBL at L=44 — [quimb TEBD example](https://quimb.readthedocs.io/en/latest/examples/ex_TEBD_evo.html) *(summarized)*. quimb's `local_expectation` exploits the same causal structure ("reverse lightcone") to make large-circuit observables cheap — [quimb circuit guide](https://quimb.readthedocs.io/en/latest/circuit/circuit.html).
- MIPT is pitched by TC-NG as requiring simulation "far beyond the reach of state-vector methods" (Clifford) or batched trajectories (Haar), with p_c found from crossings of entanglement curves — [arXiv:2602.14167](https://arxiv.org/abs/2602.14167).
- The utility/RCS frontier shows the stakes of classical simulation claims: 127-qubit observables converged to <0.01 "orders of magnitude faster than the quantum experiment" — [arXiv:2308.05077](https://arxiv.org/abs/2308.05077); a >10,000× reduction in the estimated Sycamore simulation cost — [arXiv:2002.01935](https://arxiv.org/abs/2002.01935).

### Inferences
- **Most compelling overall: threshold, logical error rate against distance, and footprint for small surface and repetition codes, extended to non-Pauli noise.**
  - Students already know stabilizers, Pauli errors and Kraus operators. The threshold and Λ cannot be obtained analytically for circuit-level noise. The incumbents' own sources show the answer moving by a factor of about 2–3 with the noise model, and Pauli-frame tools cannot even pose the amplitude-damping, coherent-error or leakage versions.
  - A study comparing the Pauli-twirled (Stim-representable) model against the exact channel at d=3–5 would be pedagogically rich: Kraus operators, twirling, leakage as a qutrit level.
  - It is also numerically decisive (exact density matrices at small d) and squarely in the documented gap of Stim and PyMatching.
  - Head-to-head anchor: reproduce the Stim notebook's d=3 and d=5 `rotated_memory_z` points (same circuits, same p) as the Pauli-noise baseline, then depart from it.
- **Second: quench dynamics, light cone and entanglement growth** (Heisenberg or TFIM chain, product or Néel initial state).
  - Concepts at master's level: Lieb-Robinson velocity, area-law against volume-law growth, why MPS cost grows with time.
  - The quimb TEBD example (L=44) is a recognised reference. Exact state-vector results at L≈20–26 give ground truth for checking TEBD truncation error, which is a concrete numerical question with a definite answer.
- **Third: the Clifford MIPT.** It is conceptually striking (a phase transition driven by measurement, with entanglement as the order parameter), needs only a stabilizer engine plus a tableau entanglement formula, and p_c comes only from finite-size numerics. It is a natural "stabilizer engine" showcase that does not require decoders.
- **Less suited as a master's QM study**: Sycamore RCS contraction-path optimization and 127-qubit utility circuits. They are very instructive about computational complexity (treewidth, slicing, FLOPs), but the physics content per unit of effort is lower, and the incumbent bar (GPU clusters, 10^14–10^22 FLOPs, bond dimension 1.6×10^7) is far out of reach for a teaching-scale comparison. A *small* RCS instance is still a good way to show cost scaling with depth through W and C.
- **Avoid as a head-to-head**: large-distance threshold sweeps under Pauli noise (d≥9, ≥10^6 shots per point). This is exactly where Stim, sinter and PyMatching are strongest (kHz sampling, µs-per-round decoding), and where the gap is likely orders of magnitude (inference from the published throughput, not measured here).

### Gaps
- I found no published survey of teaching usage or student preference for these problem classes. The ranking above is reasoned inference from what the sources show about numerical decisiveness and conceptual accessibility.
- I found no incumbent-published, ready-made benchmark for "exact non-Pauli noise on small surface codes" against which a QF study could be compared directly. The comparison would have to be constructed from the Stim baseline plus the definitions of SD6 and SI1000.
