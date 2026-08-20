# Phase 2 handoff: Korotkov quantum-Bayesian sibling essays

Phase 2 wrote the two learning-by-computing sibling essays the core essay deliberately left
out (its close names them), each grounded in papers read in full first. This is the durable
handoff; an unattended session produced it.

## Papers read in full (Mads's read-all policy; TeX by arXiv ID, no sampling)

All 7 read completely from `scratchpad/refs/<id>/`. Every essay/sheet equation is anchored to a
`<file>:<line>` in these sources.

Non-commuting (priority):
- **1702.08077** (`XZcorr-arXiv.tex`, Atalaya-Korotkov-Hacohen-Gourgy-Martin-Siddiqi) - the theory:
  two output signals, the Ito Bayesian SDEs for simultaneous sigma_z + sigma_phi (:80-90), the
  closed-form self/cross correlators K_zz, K_zphi and Gamma_pm (:135-150), the collapse recipe
  (:117, :631), the small-time K_zphi(0)=cos phi (:160), the antisymmetrized cross-correlator that
  reads a residual Rabi (:204). Supplement: single-observable Ito form (:497), rotate-and-add (:507-526).
- **1710.05249** (`MTC-Arxiv.tex`, multi-time correlators) - the cleanest general update: Bloch form
  r-dot = Lambda_ens(r-r_st) + sum_l [n_l-(n_l.r)r]/sqrt(tau_l) xi_l (:96); density-matrix general form
  (:409, :414); even-N factorization into 2-time correlators (:168), odd-N (:173).
- **1608.06652** (Hacohen-Gourgy experiment) - the SQM scheme, two-channel SME (:361), the Gaussian
  measurement operator Omega(V) (:369), the disturbance/uncertainty bound (:169-173), steady radius
  r=sqrt(eta) and isotropic diffusion for orthogonal axes (:164).

Entanglement:
- **1402.1868** (`measurement-based-entanglement-V2.tex`, Roch) - remote cascaded qubits, half-parity
  meter, X-state concurrence (:344), diagonal classical Bayes (:406) + odd-parity coherence ride-along
  (:411), Gamma_meas (:181), experimental peak concurrence ~0.35 (:194).
- **1603.09623** (`Twoqubitpaper-Mar30.tex`, Chantasri) - the two-qubit Bayes update (:136), X-state
  concurrence (:153), the concurrence-readout closed form (:177) and the sharp max-concurrence bound
  C_max(t) (:207), the three most-likely-path branches (:271-395).

QEC capstone (Mads's own line, read for framing):
- **1612.02096** (Atalaya-Bahrami, 4-qubit Bacon-Shor) - the gauge qubit under simultaneous X-and-Z
  continuous measurement diffusing on the great circle (:576, :616-626), error syndromes as
  cross-correlators <I_X12 I_X34>, <I_Z13 I_Z24> = exp(-2 Gamma_m|t1-t2|) (:771).
- **1910.08272** (Atalaya-Korotkov-Whaley, 9-qubit error-correcting) - 12 gauge ops, 4 gauge qubits,
  gauge SME (:471), triple cross-correlators (:618-635), logical error rate ~ Gamma_d^1.88, crossover
  Gamma_d ~ 10^-3/tau_coll (:1236).

## New physics-sheet sections (appended to `korotkov_physics_sheet.md`)

- **Section L**: simultaneous non-commuting z-and-x Bayesian update + output correlators. The two
  readouts and Ito SDEs, the general Bloch/density-matrix update, the discrete Gaussian-Kraus update
  M(Ibar,A)=cosh(Ibar/2D)I+sinh(Ibar/2D)A and why two of them do not commute, the collapse-into-diffusion
  dynamics, r=sqrt(eta), the disturbance bound, the closed-form K_zz/K_zphi and cos-phi law, and the
  Bacon-Shor capstone link. Every line anchored.
- **Section M**: two-qubit measurement-induced entanglement trajectories. The two-qubit Bayes update
  with geometric-mean coherence, half-parity likelihoods, X-state concurrence, the concurrence-readout
  closed form and the sharp C_max(t) bound, the three most-likely branches, the parity-meter exact
  solutions, and the Roch cascaded-SME anchors. Every line anchored.

## Sibling essays written (both Template: Default, this folder)

- **`Watching-Two-Axes.md`** (priority) + built `Watching-Two-Axes.nb`. Simultaneous continuous
  measurement of non-commuting z and x on one qubit. Through-line, delivered: the single-axis
  Gaussian-Bayes update is a diagonal closed form, and it BREAKS for two axes because the two Gaussian
  Kraus operators fail to commute by 2i sinh sinh sigma_y (shown symbolically); the general matrix
  update takes over. Then: orthogonal measurement gives isotropic great-circle diffusion (no collapse);
  the angle phi tunes collapse -> localized -> isotropic (mean|z|_end 1 -> 0.83 -> 0.64=2/pi); the
  uncertainty floor forbids the diffusion from stopping (sum-of-variances >= |commutator|, saturating
  only along y); efficiency caps the steady radius at sqrt(eta); and the two output correlators K_zz,
  K_zphi are recovered from simulated records, with the zero-lag cross-correlator reading cos(phi).
  Closes to the Bacon-Shor gauge qubit. 19 wl cells.
- **`Watching-Two-Qubits.md`** + built `Watching-Two-Qubits.nb`. Two-qubit measurement-induced
  entanglement by a half-parity meter. Delivered: the two-qubit Bayes update is the single-qubit one
  one level up, the odd-parity coherence rides the geometric mean of its two populations (ride-along
  ratio conserved, shown symbolically); a half-parity meter sorts records into |00>, |11>, or the odd
  subspace, steering odd-bin records into a Bell state; the concurrence has a sharp closed-form ceiling
  C_max(t) recovered from the ensemble envelope (MC 0.99-quantile matches C_max to <0.001); the
  final-concurrence distribution is bimodal; and the ceiling height is a measurement-rate-vs-dephasing
  trade-off (peak ~0.37 at the demo gamma, matching the Roch ~1/3 experiment). 9 wl cells.

## Verification status (IMPORTANT)

Both essays are **wolframscript-and-lint-checked but NOT adversarially verified**:
- Every wl cell runs cleanly via `wolframscript -file` (extracted and evaluated in order; no Quiet, no
  Print; symbolic checks return True, Monte-Carlo checks match the closed forms).
- `outcell_lint.py` clean on both (0 hard, 0 soft).
- md2nb built both `.nb` input-only (`"Evaluate"->False`); all section headings render (none silently
  dropped): Two-Axes 11 headings / 19 input cells, Two-Qubits 9 headings / 9 input cells.

**DEFERRED for Mads (this was an unattended session, so the ask-first loops were not run):**
- `/wl-verify` (fresh-context WL correctness reviewer) on both essays.
- `/essay-verify` (fresh-context pedagogy-DNA reviewer) on both essays.
Run both when you pick this up. The essays are self-validating in spirit (symbolic checks + MC recovery
of the closed forms), but no independent adversarial pass has graded them.

## What remains / not done

- The two adversarial loops above (deferred, ask-first).
- A possible third sibling was not in scope; the two requested siblings (non-commuting priority + two-qubit
  entanglement) are both done.
- The core `Korotkov-Quantum-Bayesian.md` and its commit `f2098851` were NOT touched.
- Heavy cells flagged in-text: the Two-Axes correlator-recovery cell (~2 min ensemble) and its two
  steady-radius/transition ensembles; the Two-Qubits ensembles are a few seconds each.

## Provenance note

This session initially wrote into the main working tree by absolute path before relocating the
deliverables into the worktree `brave-kapitsa-34b6d1` and committing there; the main tree's
pre-existing uncommitted changes (core essay, OUTLINE, core .nb) were left exactly as found.
