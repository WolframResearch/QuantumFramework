# Lanyon AI: product and technical method for numerical solver verification

Researched 2026-10-04. Labels used below:
- **[Fact]**: read directly in a primary source, or checked by this researcher against public artifacts. The method is stated each time.
- **[Company claim]**: a statement by Lanyon or its founders that I could not check independently.
- **[Third party]**: a source independent of the company.

Identity check: Lanyon AI (Lanyon AI Inc., Princeton NJ, lanyon.ai, GitHub `lanyonai`) is a real company doing exactly this work. It was founded in 2026 by Jonathan Gorard, Ammar Hakim and James (Jimmy) Juno. No disambiguation from the unrelated "Lanyon" entities was needed.

## 1. What does the product do, concretely? Inputs, outputs, solver domains, languages and ecosystems

### Takeaway
Lanyon is a proprietary AI agent, also called "Lanyon", that *generates* new solvers; it does not audit existing code. The user types natural-language prompts through slash commands (`/specify`, `/deepthought`, `/verify`, `/compile`, `/simulate`). An LLM writes a short spec in a Racket (Lisp) DSL. Deterministic symbolic code then expands that spec into Lean 4 proofs (over Mathlib `Real`) and double-precision C kernels. Finally the agent writes a simulation driver and runs it in a container.

Every public demo is a finite-volume or DG solver for a first-order hyperbolic or advection-diffusion PDE system on structured grids, in plasma physics, fluids and electromagnetism. There is no public evidence that it ingests user code, or supports Fortran, Python, Julia, MATLAB or the Wolfram Language.

### Cited Findings
- **[Fact]** Homepage tagline (2026): "Formal Verification for a Computable Universe" and "Specification → Implementation + Proof". The pitch has three steps:
  - "Our agent reasons in a highly information-dense, scientifically-aware, formal domain-specific language".
  - "Our neurosymbolic compiler generates optimized implementations and proofs, deterministically, from the same DSL source".
  - "Machine-checkable proofs certify that the implementation is correct, with no possibility of misformalization."

  [lanyon.ai](https://lanyon.ai/)
- **[Fact]** Workflow as documented by CTO Ammar Hakim (July 2026) — [Advection-Diffusion note](https://lanyon.ai/research/advection-diffusion/):
  - Step 1, `/specify` (or `/deepthought` for complex systems): Lanyon emits "a Lisp DSL fragment that describes the system of equations".
  - Step 2, `/verify`: "Lanyon will refuse to build the actual simulation code if /verify does not pass."
  - Step 3, `/compile`: produces the C code.
  - `/simulate` then takes a prose description of the run. "Lanyon will now create the simulation driver and run the simulation inside a container."
- **[Fact]** Example prompts are deliberately terse, for instance "/specify Create a 1D advection-diffusion solver". The note adds: "you do not need to specify the scheme, the limiter or any other parameters. Lanyon will choose appropriate defaults for you." `/simulate` prompts can point at a web page or paper for initial conditions, and "the RAG pipeline ... will attempt to use the initial conditions from the specified webpage." — [Advection-Diffusion note, Jul 2026](https://lanyon.ai/research/advection-diffusion/)
- **[Fact]** `/specify` uses Lanyon's "default fast model"; `/deepthought` engages its "deep reasoning model". The underlying LLM(s) are not named. Even for general-relativistic Maxwell, "Lanyon still takes less than 2 minutes and less than a few thousand tokens to complete each prompt". — [GR Maxwell note, Jul 2026](https://lanyon.ai/research/gr-maxwell/)
- **[Fact]** Outputs published per equation family, as a DSL spec (`specifications/*.rkt`), Lean proofs (`proofs/*.lean`), C kernels (`implementations/*.c`) and screen captures (`screencaps/*.mov`, stored as LFS pointers):
  - [AdvectionDiffusion](https://github.com/lanyonai/AdvectionDiffusion): linear advection; isotropic and anisotropic advection-diffusion; 1D to 3D; FV, DG p=1 and DG p=2.
  - [BurgersEquation](https://github.com/lanyonai/BurgersEquation): inviscid and viscous Burgers.
  - [CompressibleEuler](https://github.com/lanyonai/CompressibleEuler): isothermal and full Euler.
  - [MaxwellEquations](https://github.com/lanyonai/MaxwellEquations): Maxwell and perfectly hyperbolic Maxwell.
  - [GeneralRelativisticMaxwell](https://github.com/lanyonai/GeneralRelativisticMaxwell): Maxwell in curved spacetime (Kerr demos).
  - [IdealMHD](https://github.com/lanyonai/IdealMHD): with and without hyperbolic divergence cleaning.
  - [ElectrostaticVlasov](https://github.com/lanyonai/ElectrostaticVlasov): phase spaces 1x1v through 3x3v (2D to 6D).
- **[Fact]** The DSL is plain Racket. For example, `maxwell_1d.rkt` begins `#lang racket` and defines a hash with `'state`, `'parameters-assumptions`, `'fluxes` and `'wavespeeds` as quasi-quoted S-expressions, plus a `lax-friedrichs-1d` flux block. — [maxwell_1d.rkt](https://github.com/lanyonai/MaxwellEquations/blob/main/specifications/maxwell_1d.rkt)
- **[Fact]** The C output is a set of header-style kernels: flux, wave-speed, Lax-Friedrichs fluctuation, minmod reconstruction, and validity and consistency checks, all on `double`. I inspected all 72 published C files on 2026-10-04:
  - Every file ends with an empty `main` whose body is `// Insert simulation drivers here.`
  - None contains time-stepping, boundary-condition, I/O or `malloc` code.

  So the time integrator, boundary and initial conditions, and driver come from `/simulate` and are not in the published, "verified" artifacts. — [linear_advection_1d.c](https://github.com/lanyonai/AdvectionDiffusion/blob/main/implementations/linear_advection_1d.c)
- **[Fact]** Discretizations described:
  - Finite volume with second-order reconstruction, a minmod limiter, and the wave-propagation (fluctuation) form; Lax/Rusanov fluxes by default, with Roe averaging discussed.
  - Modal DG with orthonormal bases on [-1,1]^d and recovery-based diffusion.
  - "In this deep-dive Lanyon will assume simple, rectangular domains and so will use a uniform, rectangular mesh."

  Sources: [Advection-Diffusion note](https://lanyon.ai/research/advection-diffusion/); [Maxwell note](https://lanyon.ai/research/maxwell-equations/); [Ideal MHD note, Aug 2026](https://lanyon.ai/research/ideal-mhd/)
- **[Fact]** Published demonstrations:
  - Brio-Wu shock (1024 cells), Orszag-Tang vortex (1024×1024), Mach 2 flow over a conducting cylinder — [Ideal MHD](https://lanyon.ai/research/ideal-mhd/)
  - 1D/2D Riemann problems and Mach 2/Mach 4 flows over cylinders — [Euler](https://lanyon.ai/research/euler-equations/)
  - Wald and Blandford-Znajek setups around a Kerr black hole (spin 0.9, 256×256 in r-θ) — [GR Maxwell](https://lanyon.ai/research/gr-maxwell/)
  - Plane wave and a pulse in a metal box — [Maxwell](https://lanyon.ai/research/maxwell-equations/)
- **[Fact]** A "Godunov-Peshkov-Romenskii Unified Model of Continuum Mechanics" note is listed as "In preparation". — [Research index](https://lanyon.ai/research/)
- **[Fact]** Proof assistant: "Though Lanyon presently uses Lean as its default proof assistant, we intend to support other theorem-provers such as Rocq and Agda in the future." — [Advection-Diffusion note, Jul 2026](https://lanyon.ai/research/advection-diffusion/)
- **[Fact]** The research notes' computer algebra is done in **Maxima**, not Mathematica. `LanyonScripts` holds `.mac` files for Burgers, Euler, ideal MHD, isothermal Euler, Maxwell and shallow water, using the Maxima `clifford` package. — [LanyonScripts](https://github.com/lanyonai/LanyonScripts); [burgers.mac](https://github.com/lanyonai/LanyonScripts/blob/main/eqnsys/burgers.mac)
- **[Fact]** Lanyon forked the Gkeyll plasma code on 2026-08-24 (GitHub API). The founders' 2025 paper says the generated C "may either be used to run standalone simulations, or integrated into a larger computational multiphysics framework such as Gkeyll." — [lanyonai GitHub](https://github.com/lanyonai); [arXiv:2503.13877](https://arxiv.org/abs/2503.13877)

### Inferences
- Lanyon is a **correct-by-construction solver generator**: prompt → DSL → (C + Lean). It is not a verifier of third-party or legacy solver code. Nothing public suggests it can take an existing Fortran, C++ or Python CFD/FEA code and check it.
- Domain coverage today is narrow: explicit FV/DG for first-order hyperbolic and advection-diffusion systems on structured meshes, with a plasma, astrophysics and gas-dynamics flavour that matches the founders' Gkeyll background. Nothing published covers:
  - FEA or structural mechanics, elliptic solvers, implicit time integration, or linear algebra;
  - unstructured meshes, climate or finance.
- Because the driver and time integrator are generated outside the published "verified" artifacts, the guarantee (whatever its strength) covers only the spatial kernel building blocks, not the full simulation.

### Gaps
- The LLM(s) behind `/specify` and `/deepthought`, the RAG corpus, and the DSL's full grammar are undisclosed.
- How users access the product (CLI, web app or API) is unknown. The slash-command interface and "container" wording suggest a hosted agent harness, but this is unconfirmed.
- Whether an ElectrostaticVlasov research note exists: it is not on the research index as of 2026-10-04. The repo exists.

## 2. What is the technical method, and which sense of "verification"?

### Takeaway
The method is **neurosymbolic code-and-proof generation**:
- An LLM proposes a spec in a proprietary Racket/Lisp DSL.
- An in-house Lisp symbolic rewriting prover (descended from the founders' 2025 "Shock with Confidence" pipeline) emits Lean 4 theorems and proofs about Real-valued definitions of the scheme's building blocks.
- The same spec is emitted as C.

"Verification" here means **formal verification in the computer-science sense**, applied to algebraic properties of a discretization:
- Hyperbolicity, wave-speed bounds, interface consistency, flux-jump/conservation identities, reconstruction exactness and symmetry.

This partly overlaps the "code verification" category of ASME V&V and Roache / Oberkampf & Roy. It does not use the standard tools of that category: no method of manufactured solutions and no order-of-accuracy or grid-convergence testing. Solution verification (discretization-error estimation) and validation (comparison with experiment) are absent.

The technique is not Herbie/FPTaylor-style floating-point error analysis, interval or validated numerics, a posteriori error estimation, or differential testing.

### Cited Findings
- **[Company claim]** "the LLM proposes a formal specification in our own proprietary domain-specific language (DSL), and then purely symbolic algorithms expand that specification into code and proofs simultaneously. If the specification cannot be rigorously proven to be correct, the code simply doesn't generate." — [Announcing Lanyon AI, J. Gorard, Jul 2026](https://lanyon.ai/blog/welcome/)
- **[Company claim]** "Lanyon's internal (Lisp-based) symbolic theorem-prover builds upon the earlier work of Gorard and Hakim (2025) on formal verification of PDE solvers within finite precision arithmetic. As a consequence, the symbolic expressions that get translated into Lean are always parenthesized in such a way that any algebraic manipulations are guaranteed to be consistent with the IEEE-754 axioms. In cases where there exists any possibility for ambiguity, we use simp only rather than simp, restrict uses of ring_nf and field_simp only to cases where a commutative (semi)ring or field structure can be safely assumed". — [Linear PDE benchmarking, J. Juno, Jul 2026](https://lanyon.ai/research/linear-benchmarking/)
- **[Fact]** Precursor paper, "Shock with Confidence: Formal Proofs of Correctness for Hyperbolic Partial Differential Equation Solvers" (Gorard & Hakim, arXiv 2503.13877, submitted 2025-03-18, "prepared for submission to ACM"):
  - It introduces "a new formal verification pipeline for such algorithms in Racket" that generates C and "formal proofs of various mathematical and physical correctness properties ... including L^2 stability, flux conservation, and physical validity".
  - It uses "a custom-built theorem-proving and automatic differentiation framework that fully respects the algebraic structure of floating-point arithmetic".

  [arXiv:2503.13877](https://arxiv.org/abs/2503.13877)
- **[Fact]** The paper body is more hedged than its abstract:
  - "we have tried wherever possible to admit only those algebraic transformations that are permitted under the IEEE 754 standard"; the rewrite rules are "somewhat complex and ad hoc".
  - Each generated proof is "a symbolic piece of Racket code", that is, executed Racket rather than a Lean or Coq kernel.
  - The conclusions report "full proofs" for several solvers, "conditional/partial proofs" for others (incomplete superbee and van Leer limiter proofs; conditional isothermal Euler proofs), and "notable limitations".

  [arXiv:2503.13877 (PDF)](https://arxiv.org/pdf/2503.13877)
- **[Fact]** That pre-Lanyon Racket pipeline is open source:
  - `gkylcas` contains 45 `.rkt` files under `provable-algorithms/finite_volume/`, with Lax/Roe code generators and limiter tests.
  - Gkeyll contains executable Racket proofs such as `proof_inviscid_burgers_lax_cfl_stability.rkt` and `proof_isothermal_euler_mom_x_lax_strict_hyperbolicity.rkt`.
  - Gorard's last public commits to `provable-algorithms` are dated 2026-02-12 (GitHub API, checked 2026-10-04).

  [gkeyllorg/gkylcas](https://github.com/gkeyllorg/gkylcas); [gkeyllorg/gkeyll](https://github.com/gkeyllorg/gkeyll)
- **[Fact]** A second founders' paper, "BEACONS: Bounded-Error, Algebraically-Composable Neural Solvers for Partial Differential Equations" (Gorard, Hakim, Juno; arXiv 2602.14853, 2026-02-16), builds "formally-verified neural network solvers for PDEs". It derives "rigorous extrapolatory bounds on the worst-case L^inf errors of shallow neural network approximations". It has a Racket DSL, a code generator, and "a bespoke automated theorem-proving system for producing machine-checkable certificates of correctness". Demos are linear advection, inviscid Burgers and compressible Euler, in 1D and 2D. — [arXiv:2602.14853](https://arxiv.org/abs/2602.14853)
- **[Fact]** The Lean proofs are stated over Mathlib reals. Every published proof file starts `import Mathlib`, and state, parameter and flux fields are `Real` (for example `structure State where f : Real`). Across all 72 published Lean files there is no `Float` or IEEE type. — inspection on 2026-10-04 of [linear_advection_1d.lean](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean) and the [lanyonai repos](https://github.com/lanyonai)
- **[Fact]** The theorems come from a fixed template; the same named properties recur in every repository.
  - Finite-volume files carry 13 kinds per direction: `Hyperbolicity`, `WaveStability`, `DiffusiveFluxConsistency`, `WaveConsistency`, `WaveJumpCondition`, `Left/RightFluctuationsConsistent`, `FluxConservative`, `Left/RightReconstructionConsistent`, `Left/RightReconstructionLinearityPreservation`, `ReconstructionSymmetric`.
  - DG files add interface-state, interface-gradient, modal-reconstruction and volume-integral consistency and linearity theorems.
  - No published Lean theorem states TVD, L² stability of the update, positivity, a CFL condition, convergence, order of accuracy, energy conservation or divergence preservation.

  Method: theorem names extracted from all files in the 7 solver repos on 2026-10-04. — [AdvectionDiffusion proofs](https://github.com/lanyonai/AdvectionDiffusion/tree/main/proofs); [CompressibleEuler proofs](https://github.com/lanyonai/CompressibleEuler/tree/main/proofs); [IdealMHD proofs](https://github.com/lanyonai/IdealMHD/tree/main/proofs)
- **[Fact]** The published "hyperbolicity" theorems are tautological as stated.
  - The 1D advection version is `theorem xHyperbolicity ... : (∃ r1 : Real, r1 = (xFluxJacobianEigenExprs C P U).lambda1) := by refine ⟨..., rfl⟩`.
  - The Euler and ideal-MHD versions are conjunctions of the same form over 3 and 8 eigenvalues.
  - No published Lean file contains a theorem relating these eigenvalue definitions to the flux Jacobian; there are zero uses of `Matrix.det`, `charpoly`, `mulVec` or eigen predicates.
  - The MHD eigenvalue definitions use `Real.sqrt`, which Mathlib defines as 0 for negative arguments. "Real eigenvalues" therefore holds by typing, even where the true eigenvalues would be complex.
  - The research notes link these exact lines as their hyperbolicity proofs.

  Sources: [linear_advection_1d.lean#L100](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean#L100); [advection_diffusion_full_2d.lean#L197](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/advection_diffusion_full_2d.lean#L197); [compressible_euler_3d.lean#L473](https://github.com/lanyonai/CompressibleEuler/blob/main/proofs/compressible_euler_3d.lean#L473); [hyperbolic_ideal_mhd_3d.lean#L1417](https://github.com/lanyonai/IdealMHD/blob/main/proofs/hyperbolic_ideal_mhd_3d.lean#L1417); cited as proofs in [Ideal MHD note](https://lanyon.ai/research/ideal-mhd/) and [Advection-Diffusion note](https://lanyon.ai/research/advection-diffusion/)
- **[Fact]** The substantive theorems are local, per-interface algebraic identities. Examples:
  - `xFluxConservative`: the sum of left and right fluctuations equals the flux jump, under non-degenerate wave-speed hypotheses.
  - `xWaveJumpCondition`.
  - Reconstruction consistency: equal neighbouring states reproduce the state.
  - Linearity preservation: linear data are reconstructed exactly.
  - Minmod reconstruction symmetry.

  Proofs use `simp`, `ring_nf`, `field_simp`, `split_ifs <;> linarith`. — [linear_advection_1d.lean](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean)
- **[Fact]** The C kernels contain runtime "consistent" checks with a floating-point tolerance, for example `return ((fabs(f_left_reconstruction - f) < 1.0e-8));`. The Lean theorems state exact real equalities instead. — [linear_advection_1d.c (≈L418)](https://github.com/lanyonai/AdvectionDiffusion/blob/main/implementations/linear_advection_1d.c)
- **[Company claim]** On what is not checked: "the /verify step does not check if the equation or proofs of correctness are semantically correct. Semantic correctness is ensured by our sophisticated RAG pipeline ... This meta-level checking is, in general, not 100% foolproof". — [Advection-Diffusion note](https://lanyon.ai/research/advection-diffusion/)
- **[Fact]** Admitted limitations in the notes:
  - "At present, it does not correct for positivity violations or thermodynamic inconsistencies in the diffusion terms." — [Advection-Diffusion](https://lanyon.ai/research/advection-diffusion/)
  - For the DG Burgers solver, "Lanyon applies limiters in this case, though the limiters are not verified. A consequence of this is that we can't prove the total-variation bounded property, leading to single-cell oscillations in the solution." — [Burgers note, Jul 2026](https://lanyon.ai/research/burgers-equation/)
  - For the Euler system, Lanyon's solutions "heuristically satisfy TVD due to their solution structure even in the abscence of formal proofs of TVD for vector nonlinear hyperbolic PDEs." — [Nonlinear benchmarking, Aug 2026](https://lanyon.ai/research/nonlinear-benchmarking/)
- **[Fact]** V&V framing: none of the 17 site pages or the two founders' papers mentions ASME V&V 10/20/40, NAFEMS, DO-178C, NQA-1, Roache, Oberkampf, manufactured solutions, Richardson extrapolation or grid-convergence studies. Method: full-text search of the downloaded pages and papers on 2026-10-04. Comparison with exact solutions is visual, for example the Burgers caption: "The exact and Lanyon-computed solutions are essentially identical". — [Burgers note](https://lanyon.ai/research/burgers-equation/); [lanyon.ai research index](https://lanyon.ai/research/)
- **[Company claim]** Long-term framing: "This remains a fundamental aspect of the Lanyon R&D roadmap: to build a type system for the physical universe." — [Our Vision, J. Gorard, Sep 2026](https://lanyon.ai/blog/vision/)

### Inferences
- The guarantee is "**the generator's Lean statements are proved, and the C was emitted from the same spec**". It is not "the C program is proved correct".
  - No published artefact links the C to the Lean: no refinement proof, no verified extraction, no VST/Frama-C proof, no verified compiler.
  - The generator itself (DSL → Lean and DSL → C emitters, plus the Lisp prover) is proprietary and unverified, so it belongs to the trusted computing base.
  - The misformalization risk Lanyon criticises in LLM autoformalization has therefore moved into the template generator rather than been eliminated. The hyperbolicity theorems are a concrete case: the generator emits a statement too weak to mean "hyperbolic", and Lean accepts it.
- Lanyon's own benchmark rubric counts "true-by-construction theorems presented as substantive" as a misformalization escape hatch. By that rubric, the published `Hyperbolicity` theorems (and trivial ones such as `|a| ≥ |a|` in 1D advection `WaveStability`) would arguably be flagged.
- The floating-point claim is not certified by Lean. The proofs hold over ℝ, where `ring_nf`/`field_simp` are always sound. Whatever IEEE-754 discipline exists lives in the unverified Racket rewriting rules, which the 2025 paper itself calls "ad hoc" and applied "wherever possible".
- In ASME V&V terms, the closest fit is a formal, partial form of **code verification**: consistency and conservation of the discrete operators. Lanyon offers nothing for solution verification (error estimates for a given run) or validation.

### Gaps
- The internals of the Lanyon-era Lisp prover. Lanyon promised a future research note "about the internal operation of Lanyon's symbolic theorem-prover, and its interaction with other proof assistant languages"; none had appeared by 2026-10-04. — [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/)
- Whether the Racket-side checks (for example eigenvalue/Jacobian consistency, CFL, TVD, as in the 2025 public pipeline) still run inside Lanyon before Lean emission. Not stated, and not inspectable.
- How BEACONS (neural solvers with L∞ certificates) relates to the Lanyon product. It is not mentioned on lanyon.ai.

## 3. Evidence: papers, patents, posts, benchmarks, repositories, demos, talks, customers, standards

### Takeaway
The evidence is entirely first-party:
- two arXiv preprints (2025, 2026), apparently not peer-reviewed as of 2026-10-04;
- nine published research notes (a tenth is "in preparation") and four blog posts (July to September 2026);
- two self-run benchmark posts that grade frontier LLMs, not Lanyon;
- seven public GitHub repositories of generated outputs, with no generator, no build setup and no drivers.

The published Lean files do type-check, and are free of `sorry`, `axiom` and `native_decide`; I confirmed this on two files. I found no patents, customers, case studies, conference talks, independent replications or independent technical critiques. Third-party coverage is press-release syndication plus a podcast.

### Cited Findings
- **[Fact]** Research notes and posts, with dates and authors (all lanyon.ai):
  - [Advection-Diffusion](https://lanyon.ai/research/advection-diffusion/) (Hakim, Jul 2026)
  - [Maxwell](https://lanyon.ai/research/maxwell-equations/) (Hakim, Jul 2026)
  - [Burgers](https://lanyon.ai/research/burgers-equation/) (Hakim, Jul 2026)
  - [Euler](https://lanyon.ai/research/euler-equations/) (Hakim, Jul 2026)
  - [Formulary](https://lanyon.ai/research/formulary/) (Hakim, Jul 2026)
  - [GR Maxwell](https://lanyon.ai/research/gr-maxwell/) (Gorard, Jul 2026)
  - [Linear PDE benchmarking](https://lanyon.ai/research/linear-benchmarking/) (Juno, Jul 2026)
  - [Nonlinear PDE benchmarking](https://lanyon.ai/research/nonlinear-benchmarking/) (Juno, Aug 2026)
  - [Ideal MHD](https://lanyon.ai/research/ideal-mhd/) (Hakim, Aug 2026)
  - Blog: [Announcing](https://lanyon.ai/blog/welcome/) (Jul 2026), [Partners](https://lanyon.ai/blog/fundraising/) (Aug 2026), [Dorland Prize](https://lanyon.ai/blog/dorland/) (Sep 2026), [Vision](https://lanyon.ai/blog/vision/) (Sep 2026)
- **[Fact]** Self-reported generation statistics (README / note):

  | Repo | Time | Lean lines | C lines | Defs / theorems |
  |---|---|---|---|---|
  | AdvectionDiffusion | ~506 s | 20,609 | 17,986 | 980 / 582 |
  | MaxwellEquations | ~130 s | 15,209 | 8,516 | 270 / 156 |
  | IdealMHD | ~434 s | 49,317 | 31,964 | — |
  | CompressibleEuler | "just under 3 minutes" | >10,000 | >6,500 | — |
  | GR Maxwell | "a little over 5 minutes" | ~25,000 | ~12,000 | — |

  The MHD Lean takes "just under four minutes to typecheck". Sources: [AdvectionDiffusion README](https://github.com/lanyonai/AdvectionDiffusion); [MaxwellEquations README](https://github.com/lanyonai/MaxwellEquations); [IdealMHD README](https://github.com/lanyonai/IdealMHD); [Euler note](https://lanyon.ai/research/euler-equations/); [GR Maxwell note](https://lanyon.ai/research/gr-maxwell/); [Ideal MHD note](https://lanyon.ai/research/ideal-mhd/)
- **[Fact]** My counts match these figures closely: 20,588 / 977 / 582 for AdvectionDiffusion; 49,311 Lean lines for IdealMHD; 15,204 for Maxwell. Identical theorem counts across repos (for example 156 for Euler, Maxwell, GR Maxwell and ideal MHD) follow from the fixed per-file template. Method: line and declaration counts on 2026-10-04. — [lanyonai repos](https://github.com/lanyonai)
- **[Fact]** Independent type-check by this researcher (2026-10-04), using Lean v4.30.0 and Mathlib commit c5ea00351c (2026-05-26):
  - [linear_advection_1d.lean](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean) (469 lines) compiled in 34 s with no errors or warnings.
  - [compressible_euler_1d.lean](https://github.com/lanyonai/CompressibleEuler/blob/main/proofs/compressible_euler_1d.lean) (871 lines) compiled in 17 s.
  - `#print axioms` on `xHyperbolicity`, `xFluxConservative` and `xReconstructionSymmetric` shows only `[propext, Classical.choice, Quot.sound]`.
  - Across all 72 Lean files there are zero occurrences of `sorry`, `admit`, `native_decide` or `axiom` declarations.
- **[Fact]** The repos ship no `lakefile` or `lean-toolchain`, so Lean and Mathlib versions are unpinned. They ship no generator, no DSL compiler and no simulation drivers. Six of the seven solver repos are MIT-licensed; AdvectionDiffusion, the most-starred, has no license file (GitHub API, 2026-10-04). — GitHub API tree listings, 2026-10-04; [lanyonai](https://github.com/lanyonai)
- **[Fact]** GitHub org `lanyonai`: created 2026-04-04; 139 followers; 9 public repos. Solver repos were created 2026-07-17 to 2026-08-27; stars as of 2026-10-04:

  | Repo | Stars |
  |---|---|
  | AdvectionDiffusion | 73 |
  | MaxwellEquations | 44 |
  | GeneralRelativisticMaxwell | 20 |
  | BurgersEquation | 14 |
  | CompressibleEuler | 12 |
  | ElectrostaticVlasov | 12 |
  | IdealMHD | 9 |

  — [github.com/lanyonai](https://github.com/lanyonai)
- **[Fact]** Benchmark design (Lanyon-run, July and August 2026):
  - Models: Kimi K3, Claude Opus 4.8, Claude Fable 5, GPT-5.5 and GPT-5.6 Sol, at "max" or "xhigh" effort.
  - Prompts: "detailed" and "terse" variants, 3 trials each.
  - Tasks: 1D linear advection, 2D Maxwell plane wave, Burgers, and a 2D Euler Riemann problem.
  - Grading: "every agent review every output, so that all frontier models all review each other", on anonymized snapshots, with verdicts faithful / partial / misformalized / disconnected.
  - The tables list only the five frontier models. Lanyon's own outputs are reported as single reference numbers and are not graded by the panel.

  [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/); [Nonlinear benchmarking](https://lanyon.ai/research/nonlinear-benchmarking/)
- **[Fact]** Benchmark results as published, Lanyon first:

  | Task | Lanyon | Frontier models (mean output tokens) | Lanyon's stated reduction |
  |---|---|---|---|
  | 1D advection | ~7 s, ~800 output tokens | ~18k–153k | "upwards of a factor of 20-100" |
  | 2D Maxwell | ~23 s, ~600 tokens | 25k–219k | "50 to 250" |
  | 1D Burgers | ~6 s, 663 tokens | 21k–142k | "30-200" |
  | 2D Euler | ~23 s, 2,940 tokens | — | — |

  - On Euler: "the only model which came even close to Lanyon in terms of accuracy ... was Fable 5, and only with the detailed prompt, at ~70 times as many output tokens and >130 times the wallclock time".
  - With detailed prompts, frontier models scored "faithful x3" on 1D advection and Maxwell (all five models), and on Burgers (all five). On Euler, three models scored faithful x3; Opus 4.8 and Kimi K3 scored partial x2, faithful x1.
  - Failures concentrated in terse prompts and in the 2D Euler C code: "X" (an incorrect solver) and "2nd / OSC FE" (an unstable forward-Euler update).

  [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/); [Nonlinear benchmarking](https://lanyon.ai/research/nonlinear-benchmarking/)
- **[Fact]** Benchmark caveats stated by Lanyon itself:
  - Lanyon's figure "is inflated once you include the cost of writing the simulation driver" (footnote 2).
  - Kimi's cost accounting was unreliable: Moonshot's CLI did not report it, and Moonshot had credited the account $20.
  - Opus 4.8 was used instead of Opus 5 for consistency.

  [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/); [Nonlinear benchmarking](https://lanyon.ai/research/nonlinear-benchmarking/)
- **[Fact]** Customers and partners: none are named on the website, in the blog or in the press release. Named relationships are investors only:
  - Dimension (lead), Industrious Ventures, and angel Siqi Chen.
  - Lanyon also contributed to the APS William D. Dorland Prize endowment.

  [Partners post](https://lanyon.ai/blog/fundraising/); [Dorland post](https://lanyon.ai/blog/dorland/); [PR Newswire, 2026-08-17](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- **[Third party]** Theories of Everything / Curt Jaimungal episode with Gorard, 2026-10-01. The page summary: "an AI states a physical model in a formal language, and the system then generates both the simulation code and a proof, checkable by computer, that the code does what the model says"; Gorard calls this "type systems for physics". The page also mentions both papers. — [Substack episode page](https://curtjaimungal.substack.com/p/the-physicist-revolutionizing-physics); listed in [The Neuron digest, 2026-10-01](https://theneuron.ai/digest/everything-that-happened-in-ai-today-thursday-october-1-2026)
- **[Third party]** Gorard's X posts, known only from search-result titles because X returned HTTP 402:
  - The stealth reveal: "We're building a radically new kind of formally verified AI for science, math, engineering, and everything else, over at @lanyon_ai." — [X post](https://x.com/getjonwithit/status/2078105967716082089)
  - The 2025 paper thread: "We developed the first automated theorem-proving framework for (hyperbolic) PDE solvers". — [X post](https://x.com/getjonwithit/status/1902158541839856071)
- **[Third party]** Press coverage found in search results (Yahoo Finance, AOL, AI Journal, TechEdgeAI, Aerospace Trends, The SaaS News) carries the PR Newswire headline or its claims. I found no independent technical analysis. — [Yahoo Finance](https://finance.yahoo.com/technology/ai/articles/lanyon-ai-emerges-stealth-build-070000418.html); [TechEdgeAI](https://techedgeai.com/lanyon-ai-emerges-with-10-6m-bet-on-provably-correct-scientific-ai/); [Aerospace Trends](https://www.aerospace-trends.com/lanyon-ai-announces-stealth-exit-and-advances-in-high-precision-scientific-computing/)

### Inferences
- The strongest independently checkable evidence is narrow but real. The published Lean files compile with no escape hatches, and their reported sizes are accurate.
- The weakest points are:
  - what the theorems actually assert (template, local, real-valued; some tautological);
  - the absence of any published link from proof to executable code;
  - that the benchmarks grade only competitors, using competitors as graders, with n=3.
- No standards are cited, and the vocabulary (formal methods and CFL/TVD theory) does not map to the ASME/NAFEMS processes regulated industries use. This matters for the stated aerospace and nuclear targets.

### Gaps
- Patents: a quick web search on Gorard/Hakim and formally verified code generation found none. Google Patents and Espacenet were not queried directly.
- Talks: none found at SC, NAFEMS, the ASME V&V Symposium, SIAM CSE or NeurIPS/ICML workshops. The search was not exhaustive.
- The Jaimungal episode audio/video was not reviewed; only the page text was.
- Peer-review status of arXiv 2503.13877 and 2602.14853: no journal reference is listed on arXiv.
- Hacker News: a domain-restricted search found no Lanyon thread. Reddit could not be searched; the tool blocks reddit.com.
- LinkedIn staff posts were not accessible.
- No job postings were found beyond the company page's generic "Research Scientists" and "Research Engineers" calls, which give mailto only. Job postings therefore reveal nothing about the stack. — [Company page](https://lanyon.ai/company/)

## 4. Claims against evidence

### Takeaway
The headline marketing claims run well ahead of the public evidence:
- "mathematically impossible ... to make a mistake"
- "100% reliable"
- "prove sophisticated mathematical theorems"
- "invent new state-of-the-art algorithms"
- IEEE-754-respecting proofs
- "formally verified C"
- "end-to-end" correctness

The company's own footnotes narrow the first claim to "syntactic correctness". The claims about generation speed, artifact size and proof hygiene hold up. The cost advantage over frontier LLMs is plausible but entirely self-measured.

### Cited Findings
- **Claim A, "mathematically impossible for Lanyon ever to make a mistake" / "provably correct by its very construction"** — [Welcome, Jul 2026](https://lanyon.ai/blog/welcome/); [PR, 2026-08-17](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
  - **[Fact]** The company qualifies it in a footnote: "It is perhaps more accurate to say that it is mathematically impossible for Lanyon ever to commit a misformalization ... Lanyon guarantees perfect syntactic correctness. The issue of semantic correctness ... remains an open area of research for us". — [Welcome](https://lanyon.ai/blog/welcome/)
  - Evidence: the published Lean type-checks (my check), but its statements concern Real-valued definitions and some are tautological; no proof links them to the C ([linear_advection_1d.lean](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean)).
  - **Status:** uncheckable as stated, because the generator is proprietary and unverified. It is supported only in the narrow sense that the C and Lean are claimed to come from one spec.
- **Claim B, "Lanyon is not only 100% reliable"** — [PR](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
  - **[Fact]** Lanyon's own notes concede that limiters for DG Burgers are "not verified", that the RAG semantic check is "not 100% foolproof", and that TVD for the Euler system holds only "heuristically". — [Burgers](https://lanyon.ai/research/burgers-equation/); [Advection-Diffusion](https://lanyon.ai/research/advection-diffusion/); [Nonlinear benchmarking](https://lanyon.ai/research/nonlinear-benchmarking/)
  - **Status:** contradicted by the company's own technical notes.
- **Claim C, the generated solvers are "second order and TVD by construction"** — [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/)
  - **[Fact]** No published Lean theorem states TVD or second-order accuracy (my theorem-name census). The 2025 paper reports executable Racket proofs of second-order TVD for the minmod and monotonized-centered limiters, with partial proofs for superbee and van Leer. — [arXiv:2503.13877](https://arxiv.org/abs/2503.13877)
  - **Status:** partly supported by pre-Lanyon Racket work; absent from Lanyon's public Lean artifacts.
- **Claim D, "Hyperbolicity" proved for Euler, MHD and others** — [Ideal MHD](https://lanyon.ai/research/ideal-mhd/); [Euler](https://lanyon.ai/research/euler-equations/)
  - **[Fact]** The linked theorems have the form `∃ r : Real, r = λᵢ`, proved by `rfl`. — [hyperbolic_ideal_mhd_3d.lean#L1417](https://github.com/lanyonai/IdealMHD/blob/main/proofs/hyperbolic_ideal_mhd_3d.lean#L1417)
  - **Status:** not supported by the published proof statements.
- **Claim E, the symbolic prover "respects the axioms of floating point arithmetic" (IEEE-754)** — [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/)
  - **[Fact]** The Lean is over ℝ, and the 2025 paper says IEEE-safe rewriting was "tried wherever possible" with "ad hoc" rules. — [arXiv:2503.13877 PDF](https://arxiv.org/pdf/2503.13877)
  - **Status:** not machine-checked; uncheckable without the generator.
- **Claim F, "formally verified C code" / "end-to-end formally verified solvers"** — [repo READMEs](https://github.com/lanyonai/AdvectionDiffusion)
  - **[Fact]** The C contains only spatial kernels with an empty `main` ("Insert simulation drivers here"). The time-stepper, boundary conditions and driver are generated at `/simulate` and are not published or verified. — [linear_advection_1d.c](https://github.com/lanyonai/AdvectionDiffusion/blob/main/implementations/linear_advection_1d.c)
  - **Status:** overstated.
- **Claim G, "a tiny fraction of the token and compute cost of frontier models like GPT-5.6 and Fable 5"** — [PR](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
  - **[Fact]** Self-run benchmarks show 600–2,940 Lanyon output tokens against roughly 18k–219k for frontier models, excluding the driver. Lanyon is not scored by the panel. Lanyon's own LLM, inference compute and symbolic-expansion compute are not costed. — [Linear](https://lanyon.ai/research/linear-benchmarking/); [Nonlinear](https://lanyon.ai/research/nonlinear-benchmarking/)
  - **Status:** plausible but unreplicated.
- **Claim H, "thousands of lines of type-correct mathematical proof every second" ("at top speed")** — [Welcome](https://lanyon.ai/blog/welcome/)
  - **[Fact]** The published end-to-end rates are about 40–120 Lean lines/s: 20,609 lines in ~506 s; 15,209 in ~130 s; 49,317 in ~434 s. These totals include LLM time. Type-checking takes minutes. — [AdvectionDiffusion](https://github.com/lanyonai/AdvectionDiffusion); [MaxwellEquations](https://github.com/lanyonai/MaxwellEquations); [IdealMHD](https://github.com/lanyonai/IdealMHD)
  - **Status:** peak symbolic-expansion speed is uncheckable; the end-to-end figures are well below it.
- **Claim I, Lanyon can "prove sophisticated mathematical theorems, and invent new state-of-the-art algorithms"** — [Welcome](https://lanyon.ai/blog/welcome/); [PR](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
  - **[Fact]** The public theorems are a fixed template of local algebraic identities. No new algorithm is claimed in any note; the schemes are standard (Lax/Rusanov, Roe, minmod, Dedner cleaning, modal DG). — [lanyonai repos](https://github.com/lanyonai); [Ideal MHD](https://lanyon.ai/research/ideal-mhd/)
  - **Status:** unsupported by public evidence.
- **Claim J, the proofs are free of `sorry`, `axiom` and `native_decide`** (implicit in the benchmark rubric)
  - **[Fact]** Verified for all 72 files by scan, and for two files by compile plus `#print axioms`. — [lanyonai repos](https://github.com/lanyonai)
  - **Status:** supported.
- **Claim K, "Lanyon can one-shot a provably-correct C implementation of a complex 3D nonlinear PDE solver (e.g. resistive MHD), and synthesize a ~10k-line Lean + Rocq proof ... in ~3 minutes"**
  - Seen only in search-engine snippets tied to a Digg aggregation of a Gorard post; the page was unreachable (404 / empty). — [Digg](https://digg.com/tech/31wa8sqc)
  - **[Fact]** It conflicts with the company's own July 2026 note that Rocq support is future work ([Advection-Diffusion](https://lanyon.ai/research/advection-diffusion/)). It also conflicts with the August 2026 MHD note, which covers *ideal* MHD only (~50k Lean lines, ~8 min) and says resistive and reconnection physics need "non-ideal effects or ... more complex models" ([Ideal MHD](https://lanyon.ai/research/ideal-mhd/)).
  - **Status:** unverified and inconsistent with primary sources.
- **Claim L, target sectors include "GPU kernel optimization, frontier AI inference"** — [PR](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
  - **[Fact]** All public outputs are CPU C kernels for PDEs; there is no GPU or inference artifact. The Formulary lists GPU use only as an "Ultimate Discrete Scheme" desideratum. — [Formulary](https://lanyon.ai/research/formulary/); [lanyonai](https://github.com/lanyonai)
  - **Status:** aspirational.
- **Claim M, near-term roadmap: "very soon ... a provably L² stable, energy conserving discontinuous Galerkin discretization of the Vlasov-Maxwell system"**, and "Lanyon supports such extended plasma models, eventually including the full Vlasov-Maxwell system"
  - Sources: [Dorland post, Juno, Sep 2026](https://lanyon.ai/blog/dorland/); [Ideal MHD](https://lanyon.ai/research/ideal-mhd/)
  - **Status:** forward-looking. Only the ElectrostaticVlasov template repo exists. — [ElectrostaticVlasov](https://github.com/lanyonai/ElectrostaticVlasov)

### Inferences
- A fair one-line summary: Lanyon demonstrably auto-generates large, checkable Lean files of *local algebraic identities* about standard FV/DG building blocks, together with matching C kernels, quickly and cheaply. "Provably correct simulation" is a much stronger claim than those artifacts support.
- Lanyon's benchmarks test a property Lanyon defines: whether a frontier model's Lean matches its own C. They do not compare against Lanyon's own artifacts on the same rubric, or against conventional verification (MMS/order-of-accuracy test suites, Frama-C, FPTaylor, interval methods).

### Gaps
- No independent replication of any benchmark exists.
- No third-party audit of the generator exists.
- The original text and date of the "Lean + Rocq / resistive MHD" claim could not be retrieved.

## 5. Business model and positioning

### Takeaway
As of 2026-10-04 Lanyon presents itself as a "fundamental R&D lab". It publishes no pricing, deployment model (SaaS or on-premises), API, waitlist, customers or case studies; contact is by email only. Its stated targets are aerospace, space and atmospheric propulsion, nuclear energy and defense. It positions itself against "autoformalization" by general-purpose frontier LLMs, not against named companies or traditional V&V vendors.

### Cited Findings
- **[Company claim]** "We think of Lanyon AI as a fundamental R&D lab, guided by a clear philosophical vision: to formally verify the physical universe". Investors' "initial questions were not about run rates or market share ('those things will come,' they said)". — [Introducing our Partners, Gorard, Aug 2026](https://lanyon.ai/blog/fundraising/)
- **[Company claim]** On customers: "As three nerdy theorists fresh out of academia, we knew that convincing people in critical industries like aerospace, nuclear, and defense that Lanyon would be the solution to all of their woes was going to be a tall order." — [Partners post](https://lanyon.ai/blog/fundraising/)
- **[Company claim]** "Our initial target areas will be critical industries: aerospace engineering, space and atmospheric propulsion, nuclear energy." — [Welcome](https://lanyon.ai/blog/welcome/)
- **[Company claim]** The press release adds "physics, engineering, GPU kernel optimization, frontier AI inference". It quotes investor Simon Barnett (Dimension): next-token prediction "still ends at 'mostly right', a standard that doesn't clear the bar for flight controls, nuclear systems, or simulating chip tape-outs." — [PR Newswire, 2026-08-17](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- **[Company claim]** Education is a secondary vision: "Lanyon provides the foundation for a variety of computational coursework in which verified solvers can be used to validate physical intuition about the world". — [Dorland post, Juno, Sep 2026](https://lanyon.ai/blog/dorland/)
- **[Fact]** The website's only calls to action are "Start a conversation" and "Come and join us", both `mailto:contact@lanyon.ai`. There is no pricing or product-access page. — [lanyon.ai](https://lanyon.ai/); [Company](https://lanyon.ai/company/)
- **[Fact]** Comparisons Lanyon itself draws:
  - The PR names GPT-5.6 and Fable 5 as cost comparators.
  - The benchmarks cover Kimi K3, Claude Opus 4.8, Claude Fable 5, GPT-5.5 and GPT-5.6 Sol.
  - It contrasts itself with an unnamed group: "Many AI companies have proposed a workflow of 'autoformalization'", and "Lanyon is much more than an LLM harness or a coding agent wrapper."

  [PR](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html); [Welcome](https://lanyon.ai/blog/welcome/); [Linear benchmarking](https://lanyon.ai/research/linear-benchmarking/)
- **[Fact]** Corporate details relevant to positioning (funding is covered elsewhere):
  - "$10.6 million initial fundraising round led by Dimension, with participation from Industrious Ventures".
  - Address: 100 Overlook Center, Suite 2145, Princeton, NJ.
  - Founders are "all formerly of Princeton University and/or the Princeton Plasma Physics Laboratory".

  [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

### Inferences
- Lanyon appears to have no customers or revenue yet, and no productized deployment: no pricing, no named pilots, an R&D-lab framing, and investors steering it towards industry contacts. Industrious Ventures is described as providing "contact with physical reality".
- For regulated aerospace and nuclear buyers, the missing pieces are a mapping to ASME V&V 10/20/40, NQA-1 or DO-178C evidence, plus solution verification and UQ. These would likely be gating; that is my inference.
- Lanyon's positioning is AI-native ("formal methods beat LLM autoformalization"). It does not position against existing V&V or floating-point analysis tooling, or against commercial CFD vendors.

### Gaps
- Pricing, licensing (SaaS, on-premises or air-gapped) and deployment are not published.
- Whether the Gkeyll fork signals a product integration path is not stated.
- The headcount beyond the three founders is not stated.

## 6. Has the company said anything about quantum computing or quantum simulation, or about the Wolfram Language or Mathematica?

### Takeaway
No product or technical statement concerns quantum computing, quantum simulation, the Wolfram Language or Mathematica. The only mentions are historical or philosophical asides in the September 2026 Vision essay and a biographical line in the press release about Gorard co-founding the Wolfram Physics Project. The technical stack uses Racket, Lean 4/Mathlib, C and Maxima, with no Wolfram components. Full-text search covered all 17 lanyon.ai pages and both founders' papers.

### Cited Findings
- **[Fact]** Wolfram, exact quote: "and Zuse, Fredkin, and Wolfram investigated what alternative foundations for physics might look like if they were to be based on discrete computational rules, as opposed to continuous analytical ones (resolving speculations first posed by Democritus and Epicurus)." — [Our Vision, J. Gorard, Sep 2026](https://lanyon.ai/blog/vision/)
- **[Fact]** Quantum, exact quote: "Almost simultaneously with Turing and Church, physicists such as Planck were coincidentally discovering that many facets of nature (such as the energy levels of hydrogen) also appear to be fundamentally discrete. Thanks to this fortuitous co-evolution of the computational and the quantum, many of our most fundamental models of reality, from black hole event horizons to condensed matter systems, are now discrete, computational, and information-theoretic in nature: an abrupt transition from the earlier Newtonian paradigm." — [Our Vision](https://lanyon.ai/blog/vision/)
- **[Fact]** Quantum, exact quote: "We built compasses long before we understood terrestrial magnetism, steam engines long before we understood thermodynamics, and light bulbs long before we understood quantum mechanics." — [Our Vision](https://lanyon.ai/blog/vision/)
- **[Fact]** Wolfram, exact quote: "Gorard is an award-winning applied mathematician, known previously for co-founding the Wolfram Physics Project with Stephen Wolfram." — [PR Newswire, 2026-08-17](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- **[Fact]** No other occurrences of "Wolfram", "Mathematica", "quantum", "qubit" or "Schrödinger" appear on lanyon.ai (17 pages: homepage, company page, blog and research indexes, 4 blog posts, 9 published research notes). — full-text search on 2026-10-04 of [lanyon.ai](https://lanyon.ai/), [blog](https://lanyon.ai/blog/), [research](https://lanyon.ai/research/)
- **[Fact]** Pre-company context, not a company statement: the founders' 2025 paper cites Gorard's earlier Wolfram Language work, "[9] Jonathan Gorard. 2024. Computational General Relativity in the Wolfram Language using Gravitas II: ADM Formalism and Numerical Relativity. arXiv:2401.14209". — [arXiv:2503.13877](https://arxiv.org/abs/2503.13877)
- **[Fact]** The stack has no Wolfram components:
  - Racket DSL — [maxwell_1d.rkt](https://github.com/lanyonai/MaxwellEquations/blob/main/specifications/maxwell_1d.rkt)
  - Lean 4 / Mathlib — [linear_advection_1d.lean](https://github.com/lanyonai/AdvectionDiffusion/blob/main/proofs/linear_advection_1d.lean)
  - C — [linear_advection_1d.c](https://github.com/lanyonai/AdvectionDiffusion/blob/main/implementations/linear_advection_1d.c)
  - Maxima for computer algebra — [LanyonScripts](https://github.com/lanyonai/LanyonScripts)
  - The upstream `gkylcas` repo holds 850 Maxima `.mac` files and 45 Racket files. — [gkeyllorg/gkylcas](https://github.com/gkeyllorg/gkylcas)

### Inferences
- In the published proof templates, every state variable is real-valued and every equation is a first-order hyperbolic (or advection-diffusion) conservation law, so the published scope does not cover Schrödinger, Lindblad, quantum circuits or tensor networks. The Vlasov and DG roadmap stays classical. This is inferred from the template and the domains published so far.
- The Wolfram link is biographical (Gorard's Wolfram Physics Project history, and his earlier Wolfram Language numerical-relativity package). It is not a technical dependency or an announced integration. Lanyon's CAS lineage is Maxima, via Gkeyll.

### Gaps
- The Jaimungal interview audio (2026-10-01) was not reviewed. Gorard may discuss quantum topics or the Wolfram Language there; the page text does not.
- Gorard's X posts could not be read in full (HTTP 402), so any informal quantum or Wolfram remarks there are unchecked.
