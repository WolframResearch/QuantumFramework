# Related papers: what a quantum-trajectory simulator must support (2024-26 audit)

Collected 2026-08-08 for the 2D/3D measurement engine. Layout mirrors the parent corpus in
`../../related-papers/` (`eprint/<id>.tar.gz` originals, `tex/<id>/` extracted sources). This
folder is the increment that the parent corpus was thin on: the **simulator software** baseline,
the **trajectory-method cost/representation** literature, the **estimation/filter** methods, and
the explicit **multi-axis (3D) + feedback** experiments. The parent corpus (levitated expansion,
squeezing, two-particle entanglement, rotors, collapse tests) is not duplicated here; read the two
together. All main texts were read on collection day (the software/methods cluster in full, the
experimental papers' modeling sections).

One intended paper, 2506.01837 (microparticle-geometry design for 6-DOF detection efficiency), is
**PDF-only on arXiv (no TeX source)** and only tangential, so it is not saved here; it is noted as a
pointer for the per-axis efficiency point.

## Simulator software: the baseline to meet or beat
- **2412.04705** Lambert et al. *QuTiP 5: The Quantum Toolbox in Python.* The reference open-source
  stack. Ships stochastic-master-equation and stochastic-Schrödinger solvers (diffusive and jump
  unravelings), Monte Carlo wavefunction, a swappable data layer (dense / sparse / JAX / CuPy) so the
  same model runs on CPU or GPU, automatic differentiation via QuTiP-JAX + Diffrax, and control
  tooling (QuTiP-QOC). This is the feature floor any new simulator is measured against.
- **dynamiqs** (Alice&Bob / Gouzien et al.; software, arXiv paper in preparation, repo
  github.com/dynamiqs/dynamiqs). JAX, GPU-accelerated, differentiable; batches SME trajectories
  "orders of magnitude" faster by vectorizing over the noise; QuTiP-compatible API. Not saved (no
  arXiv TeX); recorded here as the differentiable-GPU baseline.

## Trajectory methods: cost and representation
- **2606.13779** Sander, Cichy, Eigel, Eisert, Fröhlich, Peham, Wille. *Computational regimes in
  matrix-product-state-based quantum trajectory simulations.* The key accounting result: total cost
  of a fixed-accuracy trajectory simulation splits into three channels, per-trajectory memory,
  per-trajectory runtime, and sampling effort (number of noise realizations), and the choice of
  "unraveling" (how the averaged dynamics is turned into random trajectories) redistributes cost
  between them rather than reducing it. Two dimensionless inflation factors (bond-dimension alpha,
  sampling kappa) and a hardware-aware decision map. For many-body / many-mode problems; the
  single-particle grid is the trivial (exact-state) end of this same spectrum.
- **2510.07211** *Efficient tensor-network simulations of weakly-measured quantum circuits.* Samples
  measurement outcomes with a Markov chain and contracts the network along space; area-law states
  reach hundreds of qubits, volume-law tens. Relevant if the engine ever grows to many coupled modes
  or spins rather than one continuous particle.

## Estimation and filtering: the tracker side
- **2510.16754** *Post-processed estimation of quantum state trajectories* (experiment on a
  continuously measured nanomechanical resonator). Quantum-state **smoothing**: using the future part
  of the record to sharpen the estimate of the past state, beyond real-time filtering. Compensates
  for gaps in the record and for information the detector never caught. Establishes smoothing /
  retrodiction as a first-class simulator output, not just forward filtering.
- **2501.13885** *Quantum model reduction for continuous-time quantum filters.* Exact reduced-order
  filters in Belavkin form: only the part of the filter state needed for the target observables is
  kept, with a principled construction (minimal realization + non-commutative conditional
  expectation). The rigorous version of "carry a handful of numbers, not the whole state," and it
  stays a physical filter.
- **2603.05468** *Kraus-Constrained Sequence Learning for Quantum Trajectories from Continuous
  Measurement.* Names the pain point directly: standard SME solvers "require exact model
  specification, known system parameters, and are sensitive to parameter mismatch." Learns the
  conditional-state map from records with a Kraus-structured output layer that stays completely
  positive and trace preserving by construction, and handles parameter drift. Motivates
  model-mismatch and per-shot-drift support.

## Multi-axis (3D) and feedback experiments: the reality the simulator must match
- **2205.10193** Pontin, Fu, et al. *Simultaneous cooling of all six degrees of freedom of an
  optically levitated nanoparticle by elliptic coherent scattering.* Three translational plus three
  librational modes cooled at once; detection channels see mixtures of several DOF, so the modes are
  coupled in both dynamics and readout. The canonical "the object is 6 coupled numbers" reference.
- **2506.21341** Melo, Veldhuizen, Tomassi, Meyer, Quidant (ETH). *Cooling via measurement-free
  coherent feedback.* All-optical feedback with a real delay (1.3 km fiber, ~6.3 us) and phase
  (colored) noise; detects x, y, z simultaneously by balanced detection; outlook explicitly names
  "unique three-dimensional cooling ... through multidirectional detection schemes" as the frontier.
  The feedback loop, the delay, and the non-white noise are all requirements this paper makes concrete.
- **2505.10157** *Quantum feedback cooling of a trapped nanoparticle by using a low-pass filter.*
  Compares low-pass-filter feedback against cold damping, delayed feedback, and linear-quadratic-
  Gaussian (LQG) control at detection efficiency eta; the feedback filter is a colored, causal object
  and eta is the deciding knob. The control-theory side of the requirements.
- **2307.06765** Pontin, Ulbricht, et al. *Sensing directional noise baths in levitated
  optomechanics.* A directed stochastic force shows up only in the x-y **cross-correlation** spectrum,
  not in the per-axis spectra: a per-axis 1D model is structurally blind to it. Direct evidence that
  coupled multi-axis modeling with cross-correlation readouts is required, not optional. "All
  calculations and numerics are fully 3D."
- **2602.07531** *Enhanced ground-state cooling of a levitated micromagnet via quantum interference in
  a cavity-magnomechanical system.* Simultaneous internal (magnon) and external (motion) cooling; a
  reminder that the "modes" are not always three translations, and coupling to auxiliary systems is
  part of the picture.

## Multi-axis readout, cross-talk, and feedback-without-filtering (added from the API sweep)
- **2603.11976** *Interference-Based 3D Optical Cold Damping of a Levitated Nanoparticle* (2026).
  A single interferometric scheme that cools all three translational axes at once; makes 3D cold
  damping a concrete, current experiment rather than an outlook.
- **2506.17172** *Feedback cooling scheme for an optically levitated oscillator with controlled
  cross-talk* (2025). The readout of one axis leaks into another (cross-talk), and the scheme turns
  that coupling into a control resource. Direct evidence that per-axis-independent modeling is not
  faithful to the hardware.
- **2409.08827** *Three-Dimensional and Selective Displacement Sensing of a Levitated Nanoparticle
  via Spatial Mode Decomposition* (2024). How x, y, z are actually separated at the detector, with
  per-axis efficiency; the measurement side of the coupled-multi-axis picture.
- **2603.06370** *Quantum Feedback Cooling without State Filtering* (2026). A contrast point: cooling
  by acting directly on the record without running a state estimator. Useful for scoping which
  problems actually need the tracker and which do not.

## What the audit concludes (see ../SIMULATOR-REQUIREMENTS.md)
The frontier runs closed-loop, imperfectly-monitored, coupled multi-mode motion under realistic
(thermal, recoil, phase, shot-to-shot) noise, and verifies with smoothing. The linear/Gaussian part
of that is dimension-agnostic and cheap; the expensive grid is needed only where a coordinate goes
non-Gaussian (double well, nonlinearity, cat/Fock preparation), and there it is usually one or two
coordinates, not all three. The missing capabilities, in priority: feedback actuation, detector
inefficiency (mixed conditional states), realistic/colored/quenched noise, coupled multi-mode +
cross-correlation readouts, smoothing, time-dependent protocols, non-Gaussian where it matters,
model mismatch, and batched/GPU/differentiable execution.
