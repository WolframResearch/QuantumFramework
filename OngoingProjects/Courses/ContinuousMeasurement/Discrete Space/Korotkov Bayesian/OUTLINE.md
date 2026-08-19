# Outline: Watching a Qubit, Korotkov's Quantum-Bayesian Formalism

Genre: "learning by computing" (computation-first .md, built to .nb via md2nb).
Ground truth: `korotkov_physics_sheet.md` (this folder). Every equation below cites a sheet section.
Status: N2 (spine). Cells are PLANNED here, not written; write them in N3+ against the exemplar.

Template: Default (per Mads). Exact front-matter keys follow the MarkdownToNotebook docs;
confirm them when the essay .md is created (N3), do not recall.

--------------------------------------------------------------------------------
## Scope
Single qubit, one symmetric broadband detector (QPC-like), "bad cavity"/Markov limit.
Symmetric detector => correlation K = 0, so ONLY the informational ("spooky") back-action;
no phase/"realistic" back-action. Defer to a later pass (N7): non-ideal-with-correlation,
finite detector bandwidth (1606.07162), feedback, two-qubit entanglement.

## Organizing question
A weak continuous measurement never collapses the qubit outright. What does the noisy
record actually contain, and why can the qubit's coherent oscillation never rise more
than 4x above the detector noise?

## Notation & conventions (preamble cell + text; shared by every section)
- units hbar = 1; measured ("localized") basis {|1>,|0>}; z = rho_11 - rho_00.
- output   I(t) = I_c + (Delta I / 2) z(t) + xi(t),  single-sided noise density S.  [sheet K]
- step avg Ibar = (1/tau) \int I dt',  variance D = S/(2 tau).                        [sheet E]
- rates    Gamma_m = (Delta I)^2 / (4 S)  (information rate),                          [sheet B]
           Gamma_d = total ensemble dephasing,  gamma = Gamma_d - Gamma_m (>=0, extra),
           ideality eta = Gamma_m / Gamma_d in [0,1]  (ideal <=> eta = 1, gamma = 0).
- qubit    symmetric (eps = 0):  H_qb = (Omega_R/2) sigma_x,  Rabi frequency Omega_R.  [sheet A]

## Cross-links (state once, up front)
Sibling essays under `Discrete Space/`: `NDSolve vs Ito/` (the Stratonovich Bloch SDE, x,y,z,Q)
and `Manual StepLike Ito/` (Rouchon trajectories). This essay is the Bayesian-update door into
the same physics; Section 7 closes the loop by matching the SDE.

================================================================================
## Section ladder (weakest assumption first)
================================================================================

### 1. A noisy detector cannot separate the two states at once
- establishes: distinguishability is gradual, not instantaneous. [sheet A, B, K]
- bridge: For a state |1> or |0> the detector current is a Gaussian about I_1 or I_0; over a
  window tau their spread D = S/(2 tau) shrinks, so the two only pull apart with time:
- cell: plot P_j(Ibar) = Normal(I_j, sqrt D) for a few tau. Out = the overlapping curves.
- interpretation: a single short look is ambiguous; information is the shrinking overlap,
  set by the measurement time tau_m = 2 S/(Delta I)^2.  [sheet B]

### 2. The populations follow Bayes' rule (frozen qubit, H = 0)
- establishes: diagonal update = classical Gaussian Bayes; gradual localization. [sheet D, E]
- bridge: With no Hamiltonian the record only tells us which state we are in, so the
  populations update by the likelihood ratio and drift toward one pole:
- cell: one trajectory. Start rho = diag(1/2,1/2); each step draw Ibar from the mixture
  rho_11 N(I_1,D) + rho_00 N(I_0,D), then rho_jj -> rho_jj P_j(Ibar)/Norm. Out = z(t).
- interpretation: "collapse" here is accumulation of log-likelihood, a noisy drift to z = +/-1,
  not a jump; the final pole is random with probability rho_jj(0).  [sheet E]

### 3. A pure state stays pure in a single run (the purity ride-along)
- establishes: ideal detector conserves |rho_10|/sqrt(rho_11 rho_00). [sheet D]
- bridge: The coherence is pinned to the populations it rides on, so for an ideal detector
  a state that starts pure never leaves the Bloch sphere:
- cell: rerun Section 2 from a pure state, tracking purity |rho_10|^2/(rho_11 rho_00) each
  step (coherence update rho_10 -> rho_10 sqrt(rho_11' rho_00'/(rho_11 rho_00))). Out = purity(t).
- interpretation: information gain alone does not decohere a single run; the ensemble
  decoheres only because different runs localize to different poles.  [sheet D, C]

### 4. Rabi under continuous measurement (add the Hamiltonian)
- establishes: the full transparent integrator = Bayes kick between unitary half-steps. [sheet C, E]
- bridge: Switch on H_qb and interleave it with the measurement by a symmetric (Strang) split,
  a half rotation, the Bayesian update, another half rotation:
- cell: the ~20-line integrator; one trajectory of (x,y,z). Second run: start fully mixed and
  watch it purify into phase-locked-but-diffusing oscillation. Out = x,y,z(t).
- interpretation: measurement localizes toward the z-poles while H rotates; the trajectory is
  noisy Rabi whose phase slowly diffuses, and a mixed start sharpens as the record accrues.

### 5. A non-ideal detector caps the purity (eta < 1)
- establishes: non-ideality = a pure-dephasing factor e^{-gamma tau}; eta sets the ceiling. [sheet J]
- bridge: A realistic detector loses some information internally, which adds a dephasing
  factor to the coherence and stops the state short of pure:
- cell: same trajectory with rho_10 *= e^{-gamma tau}, gamma = Gamma_d(1-eta); overlay
  purity(t) for eta = 1, 0.7, 0.4. Out = purity(t) family.
- interpretation: gamma is potential information dissipated before readout; the steady purity
  is set by eta, and eta = 1 is the only value that lets you monitor the wavefunction.  [sheet J, B]

### 6. PAYOFF: the output spectrum and the factor of 4
- establishes: the coherent line is bounded, peak-to-pedestal = 4 eta. [sheet I, K, 0003225]
- bridge: Rabi oscillation modulates the record, so its spectrum carries a peak at Omega_R;
  simulate many long records, average the periodogram, and compare with the closed form:
- cell: ensemble of long r(t); ensemble-averaged PowerSpectralDensity; overlay
  S_I(omega) = S + 4 eta S * Omega_R^2 Gamma^2 / [(omega^2 - Omega_R^2)^2 + Gamma^2 omega^2],
  Gamma = Gamma_d. Sweep eta; read the peak height. Out = PSD with peak + closed-form overlay.
- interpretation: the line cannot exceed 4x the noise floor, and it reaches 4 only for an
  ideal detector; half the peak is the correlation between the detector noise and the
  measurement back-action, which no classical harmonic signal can produce.  [sheet I, 0003225]

### 7. Cross-check: the Stratonovich SDE and the Ito trap
- establishes: the Bayesian update and the Stratonovich SDE agree; a naive Ito stepper does not. [sheet C]
- bridge: The same physics is written as a Bloch-Langevin SDE in `NDSolve vs Ito/`; run both
  on matched ensembles and let tau -> 0 to see them converge:
- cell: overlay the coherence decay from (a) the Bayesian update, (b) the Stratonovich SDE with
  its drift correction, (c) a plain Ito-Euler step on the same RHS WITHOUT the +(Delta I)^2/4S
  drift term. Out = three coherence-vs-time curves; (c) peels off at the wrong rate.
- interpretation: the Bayesian form has no calculus ambiguity at all (discrete Gaussian Bayes +
  unitary split); the SDE reproduces it only when the Stratonovich->Ito drift is kept.  [sheet C]

### 8. The other back-action: the phase kick and the second quadrature
- establishes: a general detector adds a unitary, record-driven z-rotation (phase/"realistic"
  back-action); the symmetric K=0 detector used so far was the special case U_r = I. [sheet J, K]
- bridge: A general reading rotates the coherence about z as well as squeezing the populations;
  turn the correlation K on and watch the same record wind the coherence phase:
- cell: one record, frozen qubit from |+>, coherence path in the (x,y) plane with phase
  back-action (K != 0, factor e^{-iK Ibar tau} on rho_10) beside the informational-only run
  (K = 0). Out = the radial K=0 meridian track vs the winding K!=0 parallel spiral.
- interpretation: the phase kick is unitary (single run stays pure) but scrambles the ensemble
  phase, adding K^2 S/4 to Gamma = (Delta I)^2/4S + K^2 S/4 + gamma, so the Section-6 peak is
  4 eta-tilde, at most 2 for a phase-preserving cQED amp (two quadratures, split 50/50). Lands on
  the Kraus polar form M_r = U_r sqrt(M_r^dag M_r): U_r = e^{-iK Ibar tau sigma_z/2} (phase),
  sqrt(M_r^dag M_r) = the informational Bayes factor.  [sheet J, K; POVM bridge]

================================================================================
## Validation & close
================================================================================
The essay is self-validating: Section 6 recovers the closed-form Korotkov-Averin spectrum from
simulated trajectories, and Section 7 matches the independent Stratonovich SDE. Both are the
"is it right" gates before the ask-first /wl-verify and /essay-verify loops.

================================================================================
## Sibling roadmap (Phase 2)
================================================================================
This essay is the K=0 core (informational back-action, one symmetric broadband detector).
Section 8 opens the door to the phase back-action at the abstract-detector level; each Phase-2
sibling turns one deferred piece of the SAME Bayesian update into its own essay:
- DONE, `Watching-Two-Axes.md` (+ .nb): simultaneous non-commuting measurement, two detectors
  reading sigma_z and sigma_phi at once; the two Gaussian kicks fail to commute by a rotation
  about the third axis, no measured axis wins, the state diffuses on the sphere and the purity
  ride-along of Section 3 is lost (steady Bloch radius sqrt(eta)); zero-lag cross-correlator = cos phi;
  Bacon-Shor gauge-qubit connection. /wl-verify + /essay-verify both OPEN 0.
- DONE, `Watching-Two-Qubits.md` (+ .nb): measurement-induced entanglement, one half-parity probe
  reading a joint operator; equal coupling of |01>,|10> (dv_2 = dv_3) => the odd-parity coherence
  rides the geometric mean sqrt(P2 P3) (the single-qubit purity ride-along one level up) and
  survives localization as a Bell state; closed-form maximum concurrence recovered from the
  ensemble envelope. /wl-verify + /essay-verify both OPEN 0.
- Continuous quantum error correction (still future): continuously monitored stabilizers with
  Bayesian tracking of the syndrome, the same update run on the code space.
Also deferred (named in the essay's own close): finite detector bandwidth (1606.07162, track the
qubit and the cavity field together) and measurement feedback to lock the Rabi oscillation.

## Build/verify checklist (for N3+)
- write each subsection against the exemplar (introduction-to-QIS-revised.md) + learning-by-computing skill.
- one bridge sentence ending ":" above each ```wl fence; one interpretation sentence below; Out cell = bare value only.
- no coined nouns (density matrix, purity, Bloch vector, likelihood, dephasing rate, ideality).
- run every wl cell through wolframscript before building.
- Template: Default (set); confirm exact front-matter keys from MarkdownToNotebook docs; build .nb only at the end.
- outcell_lint.py (deterministic) BEFORE the /essay-verify loop; /wl-verify + /essay-verify are ask-first.
