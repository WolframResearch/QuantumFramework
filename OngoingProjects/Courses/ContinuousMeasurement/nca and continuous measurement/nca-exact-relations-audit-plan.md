---
Template: Default
---

# Plan: Exact Relations for Continuously Monitored Systems from Operator Algebra

**The question.** A system observable $\hat A$ is continuously monitored; the conditional state ($|\psi\rangle$ or $\rho$) obeys a stochastic differential equation. Which of its properties are *exactly* computable by operator algebra, before any grid or trajectory, and how do the Wolfram Language 15.0 `NonCommutativeAlgebra` (NCA) features produce them as benchmarks for the numerics? The purpose is **validation, aimed squarely at the conditioning term**: the $\sqrt{\lambda\eta}\,\mathcal H[\hat A]\rho\,dW$ part is what a continuous-measurement simulation can actually get wrong (a wrong sign, prefactor, efficiency, or Itô/Stratonovich reading), so the load-bearing benchmarks are the conditioning-specific ones (classes 2 and 4, and their dependence on the efficiency $\eta$); the unconditional Lindblad-drift relations are kept only as necessary sanity checks that a correct drift does not certify a correct trajectory. The deterministic-variance result in *Watching a Quantum Particle* is one such conditioning benchmark (class 2 below), a hint, not the template.

**Deliverable.** A claims ledger: one row per exact relation, mapping it to a primary reference (resolved identifier), the reproducing kernel snippet, and an independent concrete carrier that confirms it. No relation enters on memory: it enters when it reduces to zero or evaluates to a c-number modulo a basis shown complete, a concrete representation agrees, and it is attributable to a paper actually read.

---

## 1. The two layers

Averaging over detector records (or discarding them) gives the Lindblad master equation; conditioning on one record gives the stochastic equation.

$$
\text{unconditional:}\quad \dot\rho=-\tfrac{i}{\hbar}[\hat H,\rho]+\lambda\,\mathcal D[\hat A]\rho,
$$

$$
\text{conditional:}\quad d\rho=\Big(-\tfrac{i}{\hbar}[\hat H,\rho]+\lambda\,\mathcal D[\hat A]\rho\Big)dt+\sqrt{\lambda\eta}\,\mathcal H[\hat A]\rho\,dW,
$$

with $\mathcal D[\hat A]\rho=\hat A\rho\hat A^\dagger-\tfrac12\{\hat A^\dagger\hat A,\rho\}$ and $\mathcal H[\hat A]\rho=\hat A\rho+\rho\hat A^\dagger-\langle\hat A+\hat A^\dagger\rangle\rho$. Efficient monitoring of a Hermitian $\hat A$ keeps the conditional state pure (the stochastic Schrödinger equation); efficiency $\eta<1$ or a non-Hermitian $\hat A$ forces the density-matrix form. The $dW$ here is the **diffusive** (Wiener) case, homodyne/heterodyne-like; the canonical *counting* measurement replaces it with a Poisson jump, and a single Kraus construction yields both (Rouchon, `arXiv:2208.07416`). This plan works the diffusive case only, a recorded narrowing (see the revision log), not a claim the algebra stops at jumps: the same $\mathcal L^\dagger$ engine covers them.

---

## 2. The engine: the adjoint Lindbladian, and closure

Every exact statement is a consequence of one superoperator acting on observables:

$$
\frac{d}{dt}\langle\hat O\rangle=\langle\mathcal L^\dagger\hat O\rangle,\qquad
\mathcal L^\dagger\hat O=\tfrac{i}{\hbar}[\hat H,\hat O]+\lambda\big(\hat A^\dagger\hat O\hat A-\tfrac12\{\hat A^\dagger\hat A,\hat O\}\big).
$$

$\mathcal L^\dagger$ is pure commutator and anticommutator algebra: the object NCA exists to normal-order. The one property that governs how much is solvable in closed form is **closure**: does $\mathcal L^\dagger$ map the span of a finite monomial set $\{\hat O_i\}$ into itself? If yes, the expectations obey a c-number linear system

$$
\dot{\vec m}=M\vec m+\vec b,
$$

and the exact transient $e^{Mt}$, the exact rates $\operatorname{spec}M$, and the exact steady state $-M^{-1}\vec b$ all follow by ordinary linear algebra. Whether such a finite set exists, and the smallest one, is an **elementary** question: it asks for a finite invariant monomial subspace of the *unconditional* adjoint Lindbladian, standard for Gaussian-preserving Lindbladians and for spin, decided by the moment computation itself. Do not confuse it with the quantum *finite-dimensional-filter* problem and its **estimation Lie algebra** (Amini and Gough, `arXiv:1804.10575`; classically Brockett-Mitter-Beneš): that is a different and harder *conditional* question, a necessary condition only, developed for homodyne observation, and it belongs to class 2, not here.

**What NCA adds is mechanization, not a new question.** The one-mode closures are pen-and-paper and were done by hand decades ago; the value is a symbolic engine that, for spin, multimode, or nonlinear $\hat H$, builds $M$, decides closure, and emits the benchmarks without a human tracking the operator ordering. The one regime where automation genuinely beats hand algebra is when the hierarchy does **not** close (generic nonlinear $\hat H$): the same computation is what detects non-closure and does the elimination. A caution on the machinery: a finite PBW basis for the *algebra* does not imply the *moment hierarchy* closes, and noncommutative Gröbner bases need not terminate, so closure is decided by the moment computation itself, not assumed from the encoding. And closure is searched over *monomial* spans, so a system that closes on a non-monomial set (projectors, a Lie-algebra or Casimir-reduced basis) is reported as non-closing, a false negative; "closure" here means monomial-span closure.

---

## 3. Three engines, two corollaries, and the failure case

This is **not** a clean five-way partition and should not be sold as one. Three independent engines do the work: (1) moment-hierarchy closure, (2) the Gaussian deterministic covariance, (3) multi-time correlations and output spectra. Two further items are corollaries, not separate results: (5) the steady-state identities are engine 1's fixed point plus a one-line saturation check; (4) the measurement-statistics martingale and unconditional-variance claims are engine 1's closure on $\{\hat A,\hat A^2\}$, with the Born-weight collapse a measurement-theory input rather than an algebra output. One genuine question sits outside all of them, the non-closing case, treated at the end. The numbered items below keep their labels for cross-reference; the ordering within is weakest premise first.

**The caveat that governs what these benchmarks are worth.** The engine is the *unconditional* adjoint Lindbladian, with no measurement-record term. So engines 1 and 3 and corollary 5 (the heating slope $\lambda\hbar^2$, the rates, the steady width and purity) are properties of the *dissipative* system: they hold with the detector unplugged, and because $\mathbb E[dW]=0$ they are independent of the conditioning coefficient $\sqrt{\lambda\eta}$. A simulator with a wrong sign, wrong prefactor, or wrong Itô/Stratonovich convention on the $dW$ term reproduces them anyway after averaging. Only class 2 (deterministic conditional covariance) and class 4 (martingale, collapse) test the conditioning that is actually in doubt in a continuous-measurement code, and the efficiency $\eta$ is a knob none of the engines yet track. A genuine *measurement*-validation suite must be aimed at class 2/4 and at $\eta$; the drift benchmarks catch Lindblad-drift bugs only. **That is now the plan's primary aim** (revision R7): classes 2 and 4, the efficiency $\eta$, and the convention-sensitivity check below are load-bearing, and engines 1, 3 and corollary 5 are demoted to drift sanity checks, necessary but not sufficient because they pass under a wrong conditioning term.

**The primary benchmark: a conditioning check that fails under a wrong $dW$ term.** Build an exact conditional quantity that depends on the sign, prefactor, efficiency, and Itô/Stratonovich reading of $\sqrt{\lambda\eta}\,\mathcal H[\hat A]\rho\,dW$. Two concrete ones: the class-2 conditional-covariance Riccati carries $\eta$ in its coefficients, so its steady covariance scales with $\eta$; and the short-time growth of the conditional purity, or of $\mathrm{Var}(\hat A)$, is linear in $\lambda\eta\,dt$ and flips sign under the wrong convention. Assert these against the trajectory simulator: a code with the drift right and the conditioning wrong reproduces the heating slope $\lambda\hbar^2$ and *fails* this check. That asymmetry is the whole reason the suite exists.

1. **Closed moment systems: transients, rates, steady states.** *Assumes a closing monomial set exists.* Build $M$ from $\langle\mathcal L^\dagger\hat O_i\rangle$; read the exact time course, the exact heating and relaxation rates, and the exact asymptote. The position heating slope $\frac{d}{dt}\mathbb E\langle\hat p^2\rangle=\lambda\hbar^2$ is one entry of one $M$. *Benchmark:* the whole curve, slope and asymptote, nothing to tune.

2. **Deterministic conditional variances.** *Adds quadratic $\hat H$, linear $\hat A$.* The $dW$-coefficient of the covariance cancels identically, so every trajectory shares one covariance obeying a noise-free Riccati ODE (the Kalman-Bucy filter; Jacobs 2014 §3.1.3). NCA proves the cancellation symbolically, and for the free particle the whole conditional evolution is solved in closed form (Belavkin and Staszewski, Phys. Rev. A 45, 1347, 1992): the posterior stays Gaussian with steady dispersion tending to $(\hbar/2\lambda m)^{1/2}$. Their $\lambda$ and the essay's differ by a factor of two, so the convention-independent content is the scaling $\sigma_x^2\propto(\hbar/\lambda m)^{1/2}$, which matches the essay's $\sigma_x=(\hbar/\lambda m)^{1/4}/\sqrt2$; do not claim the prefactor without fixing the strength convention. *Benchmark:* independent runs agreeing to machine precision, a test with no statistics in it.

3. **Output spectra and two-time correlations.** *Adds the regression theorem, or the characteristic operator.* The measurement record's power spectrum and any two-time correlator follow exactly from the same $\mathcal L$: through the quantum regression theorem, or through Barchielli's characteristic functional, an operator-valued generating function that yields every output moment by differentiation. No trajectory averaging. *Benchmark:* the spectrum of the simulated record.

4. **Measurement statistics.** $\langle\hat A\rangle$ is a martingale when $[\hat H,\hat A]=0$; $\operatorname{Var}(\hat A)$ decays deterministically; the trajectory ensemble of $\langle\hat A\rangle$ collapses onto $\operatorname{spec}\hat A$ with Born weights. *Benchmark:* ensemble means, variances, and the collapse histogram.

5. **Steady-state identities.** Purity and the uncertainty product at the fixed point are exact algebraic identities of the closed system, e.g. minimum-uncertainty saturation $V_xV_p-C^2=\hbar^2/4$. *Benchmark:* the steady width and purity the numerics must hit.

**The failure case: the question the engines do not answer.** For generic nonlinear $\hat H$ the moment hierarchy does not close, and then the exactly-computable content is a different set: conservation laws and symmetry-fixed quantities, short-time (Taylor) coefficients to any finite order, moment *inequalities* and bounds, and the output counting statistics from the characteristic functional (class 3), which stays exact when no finite moment set closes. This is where the tool earns more than pen and paper, so a non-closing case belongs in the plan (step 5 of the sequencing), not only the closing ones. Note that class 5 is engine 1's asymptote, not an independent result: the steady state is a symbolic fixed-point `Solve` on $M$, and only the saturation identity ($V_xV_p-C^2=\hbar^2/4$) is genuinely new content.

---

## 4. `nca-route`: the encoding is chosen by the monitored system

The whole routing decision is how the `NonCommutativeAlgebra` object is built, and it is dictated by the algebra $\hat A$ lives in. Probes are run, not recalled.

| System monitored | Fundamental relations | Route | Concrete carrier |
|---|---|---|---|
| Position / quadrature (1 boson) | $[\hat x,\hat p]=i\hbar$ or $[\hat a,\hat a^\dagger]=1$ | **E3** `"CommutationRelations"` (PBW; proven working) | truncated Fock / `QuantumOperator` |
| $N$ bosonic modes | Heisenberg-Weyl, $2N$ generators | **E3** | multimode Fock truncation |
| Spin-$j$ / qubit | $\mathfrak{su}(2)$; Casimir $\hat J^2=\hbar^2 j(j{+}1)$ a quotient | **E3 (kernel-verified)** | spin-$j$ matrices / QF Pauli |
| Fermion / Majorana | CAR $\{\hat c_i,\hat c_j^\dagger\}=\delta_{ij}$ | **E4** `"Structure"->{"Clifford",{pvars,nvars,zvars}}` (doc-grounded) | gamma matrices, Jordan-Wigner |
| Generic finite-$d$ | none exploitable | **E5** matrix carrier `{Dot,n}` | the matrices / QF |

The decision that matters is **E3 versus E2**. E3 (commutation shape) buys a finite Gröbner basis, left ideals, and a PBW normal form where `NonCommutativeExpand` alone settles equality; E2 (`"Relations"`) is general but its basis may be infinite and silently truncated, and a truncated basis gives false non-membership verdicts above its degree cap. Reaching for E3/E4 first is free: a wrong shape refuses to build and says so, rather than truncating.

Gate probes to run under `wolframscript -file`, no `Quiet`, no `Print`:

```wl
(* Baseline, proven in Part II: the CCR builds as E3. *)
xp = NonCommutativeAlgebra[<|
    "Generators" -> {{xo, po}},
    "CommutationRelations" -> {po ** xo - xo ** po + I \[HBar]}|>];
ListQ[xp["Generators"]]                (* G3 liveness: True => built; unevaluated => ::ncalg *)

(* Kernel-verified (WL 15.0.0): su(2) builds directly as E3. *)
su2 = NonCommutativeAlgebra[<|
    "Generators" -> {{jx, jy, jz}},
    "CommutationRelations" -> {
      jy ** jx - jx ** jy + I jz,        (* [jx,jy]=i jz, leading-monomial-positive *)
      jz ** jy - jy ** jz + I jx,        (* [jy,jz]=i jx *)
      jz ** jx - jx ** jz - I jy}|>];    (* [jz,jx]=i jy *)
ListQ[su2["Generators"]]
(* Result True: the three cyclic relations are already a Groebner basis.
   The proof is that NonCommutativeGroebnerBasis EQUALS the input set, not that
   its length is 3 (a non-Groebner set can also return length 3). No completion needed.
   Pauli/2 carrier anchored: correct relations -> 0, wrong sign -> nonzero. *)

(* Infinite-basis detector whenever E2 is used: length tracking the cap = truncation. *)
Length[NonCommutativeGroebnerBasis[ideal, alg, MaxIterations -> #]] & /@ {20, 40, 60}
```

Criterion, stated before scoring (unchanged from `nca-route`): **finiteness > canonicity > independence > inspectability > cost.** Keep $\hbar$, $\lambda$, mass, $j$, $n$ symbolic; `{Dot,n}` with symbolic $n$ proves a matrix identity for every truncation dimension at once.

Verify by leaving the algebra: the abstract result and any cross-check inside the same object share its monomial order, so they agree while both are wrong. Use a concrete carrier anchored to a known value, and assert the *wrong* value fails (the su(2) Casimir, the Clifford metric signature).

---

## 5. Literature audit (ground every class; resolve identifiers via the fetcher)

Protocol: build the complete fetch list first, then read it (no sampling; for the two books, the relevant chapters below in full); run the `arxiv-fetcher` pipeline to resolve identifiers and pull PDF plus TeX; cite DOI/ISBN for books and pre-1994 papers. Identifiers were confirmed to exist by web search on 2026-07-28 unless marked *resolve*.

**Primary anchors (read first; between them they carry all five classes).**

- **Jacobs, *Quantum Measurement Theory and its Applications* (Cambridge, 2014), ISBN 978-1-107-02548-6.** The definitive pedagogical source, superseding the 2006 review for depth. Chapter 3 is the spine: Gaussian continuous measurement (3.1); the SME **as the classical Kalman-Bucy filter** (3.1.3) = class 2; the **power spectrum of the record** (3.1.4, 3.7.2, App. F.8) = class 3; **equations of motion for Gaussian states** (3.7.1) = class 1; **steady states for linear open systems** (App. G.10) = class 5; and **operator reordering/splitting relations** (App. G.7) = the normal-ordering the NCA route automates. Read Chapter 3 and Appendices D, F, G.
- **Barchielli and Gregoratti, *Quantum Trajectories and Measurements in Continuous Time: The Diffusive Case*, LNP 782 (Springer, 2009), DOI 10.1007/978-3-642-01298-3.** The rigorous monograph, finite-dimensional and diffusive (the spin side). Linear and nonlinear SSE/SME (Ch. 2, 3, 5); the **characteristic operator/functional** (Ch. 4.2), an exact operator-valued generating function for every output moment, a stronger class-1/3/4 tool than moment-by-moment; the **spectrum of the output process** (Ch. 4.5); worked **two-level-atom homodyne/heterodyne spectra and squeezing** (Ch. 8, 9) = class 3 done rigorously on spin. Read Ch. 1-5, 8-9.
- **Belavkin and Staszewski, "Nondemolition observation of a free quantum particle," Phys. Rev. A 45, 1347 (1992).** The exact closed-form solution of the essay's own model: the posterior stays Gaussian and its dispersion tends to $(\hbar/2\lambda m)^{1/2}$ (prefactor convention-dependent, see class 2). The single best external check for classes 1, 2 and 5 in the position case, and the historical origin of the deterministic-variance result.

**Thematic clusters (fill in and cross-check).**

- **A. General SME and moment equations for observable $\hat A$** (classes 1, 3, 4). Jacobs and Steck, Contemp. Phys. 47, 279 (2006), `arXiv:quant-ph/0611067`; Wiseman and Milburn, *Quantum Measurement and Control* (Cambridge, 2010), including the quantum regression theorem for class 3; Gisin and Percival, J. Phys. A 25, 5677 (1992); Wiseman, "Linear quantum trajectories," `arXiv:quant-ph/9801042`.
- **B. Gaussian conditional dynamics, filtering, and the covariance-matrix toolbox** (classes 2, 5). Doherty and Jacobs, Phys. Rev. A 60, 2700 (1999), `arXiv:quant-ph/9812004`; Genoni et al., Phys. Lett. A 494, 129260 (2024), `arXiv:2312.13214` (verified in source: general-dyne monitoring gives a stochastic first-moment evolution and a **deterministic** conditional covariance obeying a Riccati equation, class 2 in full generality); Bouten, van Handel and James, SIAM J. Control Optim. 46, 2199 (2007), `arXiv:math/0601741`; Belavkin, J. Multivariate Anal. 42, 171 (1992) and Comm. Math. Phys. 146, 611 (1992); Weedbrook et al., *Gaussian quantum information*, Rev. Mod. Phys. 84, 621 (2012), `arXiv:1110.3234`.
- **C. Symbolic prior art, and why this plan makes no novelty claim.** Symbolic engines already joined to continuous measurement: QNET (Tezak, Goerz, Mabuchi; SymPy operator and superoperator algebra over the Gough-James SLH formalism), SymPsi, NCAlgebra (Helton et al., UCSD), QuAlg (`arXiv:2008.06467`), QuantumCumulants.jl (`arXiv:2105.01657`, unconditional mean-field). Two framings to keep straight: **engine-1 closure is elementary** (an invariant monomial subspace of the *unconditional* adjoint Lindbladian, standard for Gaussian-preserving Lindbladians), whereas the *conditional* finite-dimensional-filter question is the harder one, with the **estimation Lie algebra** as its apparatus (Amini and Gough, `arXiv:1804.10575`, a necessary condition only, homodyne) and it belongs to class 2, not engine 1. The deterministic-covariance cancellation is a known by-hand result (Genoni, Lami and Serafini, `arXiv:1607.02619`; Jacobs §3.1.3). **The novelty check was run, and the residual novelty does not survive it.** QNET's `symbolic_heisenberg_eom` already emits both halves the plan called new: the ensemble-averaged (unconditional) equation of motion, from which $M$ is read, and the conditional QSDE when noises are supplied. SymPsi generates operator equations of motion from a Hamiltonian plus collapse operators, and a *conditional* stochastic cumulant mean-field for homodyne monitoring is published (Zhang et al., `arXiv:2306.00868`), so even the conditional gap is not clean. The only piece not clearly pre-existing is an exact closure *decision* plus non-closing *elimination* by noncommutative Gröbner/PBW, and the draft concedes that need not terminate. **So this plan makes no novelty claim.** Its worth is validation: an independent symbolic route to benchmark a simulation against, worth building even when the method is not new, and to be run head-to-head against `symbolic_heisenberg_eom` on the same spin and boson examples to confirm agreement.
- **D. Output spectra and quantum noise** (class 3). Clerk, Devoret, Girvin, Marquardt, Schoelkopf, *Introduction to quantum noise, measurement, and amplification*, Rev. Mod. Phys. 82, 1155 (2010), `arXiv:0810.4729`; plus the spectra chapters of the two primary anchors.
- **E. External benchmark values** (classes 1, 5). Standard-quantum-limit and backaction bounds: Caves (1980-1982); Braginsky and Khalili, *Quantum Measurement* (Cambridge, 1995). Diósi continuous-measurement and CSL heating/diffusion rates and smeared-position operators for the smeared-$\hat A$ continuation. Scan recent "exact under continuous monitoring" results for priority, e.g. `arXiv:2606.29042`, `arXiv:1611.02077`, judged on reading not title.

---

## Verification log (2026-07-28)

Line-by-line checks run so far. Method: K = kernel (WL 15.0.0), D = doc page, S = source paper, T = book table of contents.

| Claim | Method | Result |
|---|---|---|
| CCR builds as E3; $[\hat x,\hat p]=i\hbar$, $[\hat x,[\hat x,\hat p^2]]=-2\hbar^2$, $\hat x\hat p\hat x=\hat x^2\hat p-i\hbar\hat x$ | K | confirmed |
| su(2) builds directly as E3, no Gröbner completion (the Gröbner basis *equals* the three input relations, verified set-equal, not merely length 3); Pauli/2 carrier anchored, wrong sign rejected | K | confirmed; resolves the plan's open question |
| Fermion/CAR encoding is E4 `"Structure"->{"Clifford",{pvars,nvars,zvars}}` | D | confirmed (worked examples in the doc pages) |
| Class 2: general-dyne monitoring gives a deterministic conditional covariance obeying a Riccati equation | S (Genoni `2312.13214`, l. 912, 1010) | confirmed verbatim |
| QuantumCumulants is unconditional mean-field / cumulant only | S (title + abstract) | confirmed; supports the conditional-tooling gap |
| Belavkin-Staszewski free-particle steady dispersion $\propto(\hbar/\lambda m)^{1/2}$ | S (PRA 45, 1347 abstract) | scaling confirmed; prefactor convention-dependent |
| Barchielli characteristic operator (Ch. 4.2) and output spectrum (Ch. 4.5); Jacobs Kalman-Bucy (§3.1.3) | T | confirmed |
| Identifiers: Doherty-Jacobs `quant-ph/9812004`, Clerk et al. `0810.4729` | S | resolved and fetched |
| Paper library: 11 arXiv items fetched (TeX or PDF) | - | complete |

Still to verify: the Belavkin filtering references (J. Multivariate Anal. 42, 171; CMP 146, 611) by full read, the smeared-$\hat A$ (CSL) heating rates, and the class-3 output-spectrum formulas against a concrete carrier.

## 6. Verification and cross-carriers

Each class closes when two routes sharing no assumptions agree on an anchored value.

1. **Kernel NCA vs QF matrices.** Build $\hat x,\hat p,\hat A,\hat H$ as truncated matrices or `QuantumOperator`s; evolve the Lindblad equation with `QuantumEvolve` (Lindblad rule-form `{L}->{g}`; a zero Hamiltonian passed as an explicit zero matrix) and confirm $M$. Watch the truncation failure mode: heating drives $\langle\hat p^2\rangle$ off any fixed cutoff, and a coarse cutoff saturates silently at a plausible wrong number, so the carrier dimension is anchored to $\sqrt{\lambda\hbar^2 t}$ from above.
2. **Kernel NCA vs QuantumCumulants.jl** (cluster C): two independent symbolic derivations agreeing on a symbolic-parameter law.
3. **Structural anchors**: assert the Casimir / metric value and that the wrong value fails.

---

## 7. Sequencing

1. **Gate probes (de-risks everything, one short kernel session).** Run the CCR baseline, the su(2) build, and a Clifford/CAR build; record which encoding each system admits. Until this runs, the spin and fermion picks are hypotheses.
2. **Build the conditioning benchmarks first, for position.** The class-2 conditional-covariance Riccati with the efficiency $\eta$ carried in its coefficients, and the convention-sensitivity check (the short-time $\lambda\eta\,dt$ growth of conditional purity / $\mathrm{Var}(\hat A)$ that flips sign under a wrong $dW$ term); assert both against the trajectory simulator, since these are what the suite exists to catch. Then, as necessary-not-sufficient regression against known-good physics, reproduce the drift sanity checks (the class-1 $M$ matrix, the heating slope $\lambda\hbar^2$).
3. **Fetch and read clusters A-D; fill the ledger.**
4. **Extend to spin-$j$** (dephasing rates, Bloch-vector contraction, collapse statistics) **and a linear bosonic mode with nontrivial quadratic $\hat H$** (measured quadrature, squeezed steady state), each with its carrier cross-check and its output spectrum (class 3).
5. **Do at least one non-closing case** (a monitored anharmonic or Kerr oscillator, or a spin with nonlinear coupling), where the hierarchy does not close: show the tool detecting non-closure, running the elimination, and delivering the short-time coefficients and moment bounds that survive. Without this the plan only reproduces textbook benchmarks a practitioner already owns; this is the step that makes automation worth more than pen and paper.

**Ledger schema** (one row per relation): relation | engine/corollary | assumptions | NCA encoding + probe | reproducing snippet | concrete carrier + anchor | primary reference (resolved id) | status.

**Spawn boundary.** The fetch-and-read (step 3) and the extensions (steps 4-5) are each self-contained and heavy; either is a separate session once steps 1-2 fix the method.

**Scope guard.** This plan covers exact relations obtainable by algebra as benchmarks for numerics, in the **diffusive** case. Set aside deliberately, each a recorded narrowing rather than a claim the method fails: (a) jump/counting (Poisson) monitoring, which the same $\mathcal L^\dagger$/Kraus construction covers (Rouchon `arXiv:2208.07416`); (b) the filtering-*optimality* question (is the conditional state the least-squares estimator, and its innovations structure), the estimation-theory half of the same object; (c) propagators, feedback-control design, and collapse-model phenomenology beyond borrowing their rates as check-points.

## Revision log

Append-only; records deviations from the original frame so they are revision, not drift.

- **2026-07-28, R1.** su(2) encoding resolved from open question to E3, by kernel evidence: the Gröbner basis of the three cyclic relations *equals* the input set (so no completion; length 3 alone would not prove it).
- **2026-07-28, R2.** Added the Belavkin-Staszewski exact-solution anchor; elevated the Jacobs and Barchielli-Gregoratti books to primary anchors.
- **2026-07-28, R3.** Narrowed scope to the diffusive case explicitly (the draft used $dW$ throughout without saying so); jump/counting set aside as (a) above. An unprincipled narrowing, recorded here rather than left silent.
- **2026-07-28, R4.** Reframed the closure question as the quantum finite-dimensional-filter problem (estimation Lie algebra, Amini-Gough), and tempered the novelty from "answers when closure holds" to "mechanizes a conceptually-solved decision and emits the benchmark generator," pending a check against QNET and SymPsi.
- **2026-07-28, R5.** Relabeled the "five classes" as three engines, two corollaries, and the non-closing failure case, since class 5 is engine 1's asymptote and class 4 is closure on $\{\hat A,\hat A^2\}$ plus a Born-rule input.
- **2026-07-28, R6 (corrects R4, from the second adversarial pass).** The estimation Lie algebra (Amini-Gough) was mis-applied in R4: it decides the *conditional* finite-dimensional-filter problem (a necessary condition, homodyne), not engine 1's *unconditional* moment closure, which is elementary. Amini-Gough is now confined to class 2. The novelty claim is **retired**: the deferred QNET/SymPsi check was run, and QNET's `symbolic_heisenberg_eom` plus the published conditional cumulant mean-field (`arXiv:2306.00868`) already emit the benchmark generator; the plan is repositioned as a validation tool. Recorded, as an unlogged scope slippage now caught: the opening question is *conditional* ("properties of the conditional state"), but engines 1/3 and corollary 5 compute *unconditional* ensemble averages and only class 2/4 are conditioning-dependent, so the flagship benchmarks certify the Lindblad drift, not the conditioning term. Re-aiming the suite at class 2/4 and at $\eta$ is an open design decision for the author, not yet taken.
- **2026-07-28, R7 (design decision taken).** Re-aimed the plan at the conditioning term: classes 2 and 4, the efficiency $\eta$, and a $dW$-convention-sensitivity benchmark are now primary; engines 1/3 and corollary 5 are demoted to drift sanity checks (necessary, not sufficient). The opening, the Section 3 caveat, and the sequencing (step 2 builds the conditioning benchmarks first) are reframed accordingly. Validation, not novelty, remains the stated worth (R6).
