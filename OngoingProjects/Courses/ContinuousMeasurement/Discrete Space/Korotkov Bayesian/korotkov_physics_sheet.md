# Korotkov quantum-Bayesian formalism: grounded physics sheet

Source of truth = the TeX under `scratchpad/korotkov_refs/`. Every equation below carries a
`file:line` anchor. Notation kept in Korotkov's own symbols. Built for the essay-only pass.

Status of the 8-paper read:
- [x] cond-mat/0008461  Meas2-3.tex          (PRIMARY: full Bayesian formalism) — READ FULL
- [x] cond-mat/0002203  short.tex            (Korotkov-Averin spectrum, S/N<=4) — READ FULL
- [x] cond-mat/0211647  NONIDEA2.TEX         (nonideal detectors in Bayesian formalism) — READ FULL
- [x] cond-mat/0209629  ARW-XXX.TEX          (Bayesian-approach review chapter) — READ FULL
      Kor-rev: pedagogical synthesis (Bloch+pointer AND correspondence-principle derivations both given).
      Clean NORMALIZED spectrum (ARW-XXX.TEX:1060-1062): S_I(omega)/S_0 = 1 + 4 eta/[(omega/Omega)^2
        + (omega^2-Omega^2)^2/(Omega^2 Gamma_d^2)]  -- nicest form for the essay's payoff plot.
      Fundamental bound Gamma_d tau_m >= 1/2 (ARW-XXX.TEX:778); eta = 1/(2 Gamma_d tau_m).
      Feedback (ARW-XXX.TEX:1131-1210): linear H_fb=H(1-F*Dphi), synchronization D=exp(-C/32F)->1 at F>>C,
        C = hbar(Delta I)^2/(S_0 H); direct feedback H_fb/H-1 = F{2[I(t)-I_0]/Delta I - cos(Omega t)} sin(Omega t),
        best near F/C=1/4 (no real-time Bayesian solve needed). Dephasing suppressed if d_e=gamma_e/[(Delta I)^2/4S_0] small.
      N-qubit Bayesian eq (ARW-XXX.TEX:1225-1244): gamma_ij=(eta^{-1}-1)(I_i-I_j)^2/4S_0; NO decoherence
        between equally-coupled states (I_i=I_j). Two-qubit: spontaneous Bell-singlet entanglement
        (prob 1/4) vs oscillating subspace (prob 3/4, peak-to-pedestal 32 eta/3).
- [x] cond-mat/9909039  Meas-s3m.tex         (origin: double dot) — READ FULL
      Origin/confirm: nonlinear Langevin sigma_dot_11 = -sigma_11 sigma_22 (2 Delta I/S_I)[I(t)-I_0]
      (Meas-s3m.tex:319-326); Gamma_d=(Delta I)^2/4S_I, S_I=2eI_0; tau_loc=tau_dis=2S_I/(Delta I)^2,
      tau_loc=tau_d/2; coupling C = hbar(Delta I)^2/(S_I H); nonideal gamma_d=Gamma_d-(Delta I)^2/4S_I
      as "two ideal detectors in parallel, one output inaccessible". Purification + Zeno first shown.
      (Sign note: here Delta I = I_2 - I_1, opposite to 0008461. Density matrix called sigma.)
- [x] cond-mat/0003225  Meas-sp3.tex         (output spectrum of oscillations) — READ FULL
      Kor-spec: THE spectrum derivation. Symmetric-qubit output spectrum (Meas-sp3.tex:470-474):
        S_I(omega) = S_0 + Omega^2 (Delta I)^2 Gamma / [ (omega^2-Omega^2)^2 + Gamma^2 omega^2 ]
        Gamma = eta^{-1} (Delta I)^2/4 S_0 = alpha eta^{-1} Omega,  coupling alpha = hbar(Delta I)^2/(8 S_0 H)
        peak height = 4 eta S_0  (peak-to-pedestal = 4 eta); width = alpha eta^{-1} Omega.
      Weak-coupling Lorentzian (Meas-sp3.tex:369): S_I = S_0 + 4 S_0/[1+[8 S_0(omega-Omega)/(Delta I)^2]^2].
      WHY 4 (Meas-sp3.tex:371-384, 793-807): half the peak = z(t) evolution, half = nonclassical
        correlation K_{xi z} between detector noise xi and the qubit backaction z. Conventional view:
        discrete eigenvalues of z-hat (collapse) give the extra factor 2. Integral under peak=(Delta I)^2/4,
        twice the classical harmonic value.
      Zeno (strong coupling): overdamped at alpha/eta > 2, peak gone at alpha/eta > sqrt2; Lorentzian
        omega_1 = Omega^2/Gamma = Omega eta/alpha. Bayesian == conventional spectrum proven (any alpha,eta,eps).
      Finite-T env (Meas-sp3.tex:592-748) is NOT generally the same as detector nonideality eta<1
        (coincide only for |eps/H|<<1, or high-T with tuned gamma_1,gamma_2). => keep env and eta distinct.
- [x] 1111.4016         LesHouches2.tex      (pedagogical, circuit QED) — READ FULL (+ appended draft)
- [x] 1606.07162        finite-kappa-4.tex   (moderate bandwidth) — READ FULL
      Finite-bandwidth (kappa ~ Gamma, NOT bad-cavity): the SAME Bayesian equations (rho-diag/rho-off,
      Stratonovich forms identical to 0008461) but applied to the ENTANGLED qubit-resonator matrix
      elements rho_jj', with the two resonator fields evolving classically
        alpha_dot_j = -i(omega_r +/- chi - omega_d) alpha_j - (kappa/2) alpha_j - i eps   (finite-kappa-4.tex:369,422)
      and time-dependent dephasing Gamma_d(t) = (kappa/2)|alpha_1(t)-alpha_0(t)|^2 (finite-kappa-4.tex:509,960).
      State = 4 numbers (rho_00,rho_11,rho_10,rho_01) + 2 fields (alpha_0,alpha_1). Main result eqs
      finite-kappa-4.tex:942-964 (phase-sensitive), :1004-1019 (phase-preserving); Strat/Ito diff forms :1055-1112.
      Numerical WIN: exact for arbitrary Delta t when params constant (unlike Wiener-process quantum
      trajectories) -> larger steps. NOT valid for an evolving (Rabi) qubit when Rabi ~ kappa (cat states form).
      Key coherent-state fact (Appendix A, finite-kappa-4.tex:1472,1476-1483): a coherent state stays
      COHERENT (pure!) under damping, alpha(t)=alpha(0) e^{-i omega_r t} e^{-kappa t/2}; no decoherence.
      This is why the whole cQED-measurement analysis stays simple.
      => THIS is the concrete "improve the formalism" candidate (N7): finite-bandwidth entangled
         qubit-resonator Bayesian update, the extension past the essay-only broadband toy.

===============================================================================
## READING COMPLETE: all 8/8 read in full. Every essay claim is anchored above. ##
===============================================================================
Essay-only build rests on: sections A-E + I + J + K (the discrete Bayesian update, the ideal-detector
purity ride-along, non-ideality eta, the Korotkov-Averin spectrum, peak-to-pedestal = 4 eta_tilde).
Use 1111.4016's z=rho_11-rho_00 notation. Validate against the sibling Stratonovich SDE (NDSolve vs Ito/).
Deferred to N7 (formalism extension, NOT essay-only): finite bandwidth (1606.07162), feedback, entanglement.

===============================================================================
## A. Setup and notation  (0008461)
===============================================================================

Qubit density matrix rho_ij in the LOCALIZED basis |1>,|2> (the basis the detector couples to).
rho_11 + rho_22 = 1; rho_21 = rho_12*; pure state <=> |rho_12|^2 = rho_11 rho_22.

Qubit Hamiltonian (Meas2-3.tex:154-157):
  H_QB = (eps/2)(c1+ c1 - c2+ c2) + H (c1+ c2 + c2+ c1)
  eps = energy asymmetry between levels; H = tunneling/mixing (real WLOG).
  Diagonal-basis level splitting: hbar*Omega = (4 H^2 + eps^2)^{1/2}   (Meas2-3.tex:166)
  Symmetric qubit eps=0 => Omega = 2H/hbar   (Meas2-3.tex:1493)

Detector: output current I(t). Averages I_1 (if state |1>) and I_2 (if |2>).
  Delta I = I_1 - I_2                                   (Meas2-3.tex:217-219)
  I_0 = (I_1 + I_2)/2                                   (Meas2-3.tex:256)
  S_1, S_2 = low-freq spectral densities of detector noise for the two states.
  S_0 = (S_1 + S_2)/2  (symmetrized)                    (Meas2-3.tex:258)
  Weakly-responding (linear) regime: |Delta I| << I_0, |S_1 - S_2| << S_0.

Shot noise forms (context, not needed for the toy essay):
  tunnel junction: S_{1,2} = 2 e I_{1,2}                (Meas2-3.tex:244, Schottky)
  QPC:             S_{1,2} = 2 e I_{1,2} (1 - T_{1,2})  (Meas2-3.tex:283)

===============================================================================
## B. Rates and the ideality eta   (0008461)  -- LOAD-BEARING
===============================================================================

Typical measurement time (weakly responding)   (Meas2-3.tex:262):
  tau_m = 2 S_0 / (Delta I)^2

Ensemble decoherence rate (the "conventional"/deterministic dephasing):
  Gamma_d  (general),  with lower bound Gamma_d >= (2 tau_m)^{-1}
  ideal:   Gamma_d = (Delta I)^2 / (4 S_0)             (Meas2-3.tex:425, GdS0)
  relation Gamma_d = 1/(2 tau_m) for ideal              (Meas2-3.tex:420, GdTm)

Measurement ("information acquisition") rate  -- Korotkov writes it as (Delta I)^2/4S_0:
  Gamma_m := (Delta I)^2 / (4 S_0)         [= 1/(2 tau_m)]

Pure-environment (extra) decoherence   (Meas2-3.tex:581-583, gamma_d):
  gamma_d = Gamma_d - (Delta I)^2/(4 S_0) = Gamma_d - Gamma_m >= 0

Ideality factor  (Meas2-3.tex:589-591):
  eta := 1 - gamma_d/Gamma_d = 1/(2 Gamma_d tau_m)  =  Gamma_m/Gamma_d,   0 < eta <= 1
  eta = 1  <=> ideal (quantum-limited) detector: tunnel junction, symmetric QPC,
              quantum-limited SQUID.  eta<1 => extra dephasing gamma_d>0.

Finite-temperature tunnel junction ideality  (Meas2-3.tex:1317):
  eta = [tanh(beta eV / 2)]^2

Correlated output+backaction (SET) generalization (Meas2-3.tex:1400-1422): a
quantum-limited TOTAL energy sensitivity (eps_I eps_phi - eps_{I phi}^2)^{1/2} = hbar/2
is equivalent to ideality eta~=1. Uncorrelated case reduces to eta = (hbar/2)^2/(eps_I eps_phi)
= 1/(2 Gamma_d tau_m).  -- context for "why eta<=1", not needed for the toy build.

===============================================================================
## C. The Bayesian equations (Stratonovich)   (0008461)  -- THE FORMALISM
===============================================================================

Stratonovich SDEs for a qubit continuously measured (Meas2-3.tex:566-578, Bayes1/Bayes2):
  d rho_11/dt = - d rho_22/dt
              = -2 (H/hbar) Im rho_12  +  rho_11 rho_22 (2 Delta I / S_0) [I(t) - I_0]
  d rho_12/dt = i (eps/hbar) rho_12  +  i (H/hbar)(rho_11 - rho_22)
                - (rho_11 - rho_22)(Delta I / S_0)[I(t) - I_0] rho_12
                - gamma_d rho_12
  (interpreted in Stratonovich sense; explicitly stated Meas2-3.tex:578,629-635)

Detector output that closes the system (Meas2-3.tex:607-609, I(t)):
  I(t) - I_0 = (Delta I / 2)(rho_11 - rho_22) + xi(t),   xi white, S_xi = S_0

Ito form of the SAME system (Meas2-3.tex:667-679, Ito1/Ito2):
  d rho_11/dt = -2 (H/hbar) Im rho_12 + rho_11 rho_22 (2 Delta I/S_0) xi(t)
  d rho_12/dt = i(eps/hbar) rho_12 + i(H/hbar)(rho_11-rho_22)
                - (rho_11-rho_22)(Delta I/S_0) rho_12 xi(t)
                - [ gamma_d + (Delta I)^2/(4 S_0) ] rho_12
  NOTE the extra +(Delta I)^2/4S_0 in the Ito coherence decay: it is the Strat->Ito drift
  correction, NOT physical decoherence (Meas2-3.tex:694-698). Get this wrong (feed Strat
  eqs to an Ito stepper without it) and rho_12 decays at the wrong rate. THE trap the essay
  demonstrates.

Strat<->Ito translation rule used (Meas2-3.tex:649-664):
  Strat  dx_i = G_i + F_i xi   <=>  Ito  dx_i = G_i + (S_xi/4) sum_k (dF_i/dx_k) F_k + F_i xi

Ensemble average of the Ito eqs over xi  ->  the conventional eqs (Meas2-3.tex:397-404):
  d rho_11/dt = -2(H/hbar) Im rho_12
  d rho_12/dt = i(eps/hbar) rho_12 + i(H/hbar)(rho_11-rho_22) - Gamma_d rho_12

===============================================================================
## D. Exact H=0 solution: Bayes update + purity ride-along  (0008461) -- CONFIRMS MEMORY
===============================================================================

For H=0 the SDEs integrate exactly over a step tau (Meas2-3.tex:708-723, sol1/sol2):
  rho_11(t+tau)/rho_22(t+tau) = [rho_11(t)/rho_22(t)]
                                * exp[-(Ibar - I_1)^2 tau / S_0] / exp[-(Ibar - I_2)^2 tau / S_0]
  rho_12(t+tau)/sqrt(rho_11(t+tau) rho_22(t+tau))
       = [rho_12(t)/sqrt(rho_11(t) rho_22(t))] * exp(i eps tau/hbar) * exp(-gamma_d tau)
  where  Ibar(tau) = (1/tau) \int_t^{t+tau} I(t') dt'   (Meas2-3.tex:722)

=> sol1 IS the Bayes formula on populations; sol2 says |rho_12|/sqrt(rho_11 rho_22)
   is CONSERVED for an ideal detector (gamma_d=0), decaying only via exp(-gamma_d tau).
   [memory claim 3 CONFIRMED]

===============================================================================
## E. Monte-Carlo algorithm (the Bayesian update to implement) (0008461) -- CONFIRMS MEMORY
===============================================================================

Bayesian algorithm, step Delta t (Meas2-3.tex:730-754):
  1. draw the averaged current Ibar over Delta t from the Gaussian MIXTURE
       P(Ibar) = rho_11 N(I_1, D) + rho_22 N(I_2, D),   D = S_0/(2 Delta t)   (Meas2-3.tex:740-747, Gauss1)
     [memory claim 1 CONFIRMED: averaged-current variance D = S_0/(2 Delta t)]
  2. substitute Ibar into sol1/sol2 -> rho_ij(t+Delta t)
  3. apply the finite-H unitary rotation for Delta t (rotation in rho_11 - rho_12 plane)
  Requires Delta t << hbar/H.  "Allows longer timesteps" than the Ito stepper.

Ito algorithm alternative (Meas2-3.tex:756-779): draw xibar ~ N(0,D), D=S_0/(2 Delta t),
  Euler-step Ito1/Ito2. Requires Delta t << all of {S_0/(Delta I)^2, hbar/H, hbar/eps, 1/gamma_d}.
  The two algorithms are equivalent as Delta t->0 (Meas2-3.tex:781-789). <-- THIS is the
  Bayesian-vs-SDE equivalence the essay makes felt, and the natural place to expose the Ito trap.

Validity window (from the pointer derivation, Meas2-3.tex:1084-1099):
  e/I_0 << Delta t << min[ e I_0/(Delta I)^2 , hbar/H ]
  i.e. tau_detector << Delta t << (Gamma_m^{-1}, Omega^{-1})   [memory claim 5 CONFIRMED in spirit]

===============================================================================
## F. Quantum jumps (finite Delta I), pointer/collapse derivation (0008461) -- context
===============================================================================

Frequent one-electron readout gives quantum-jump eqs (Meas2-3.tex:1147-1169, QJ1-QJ5):
  between jumps: drift with -(Delta I/e) rho_11 rho_22 etc.;
  at an electron passage: rho_11 -> I_1 rho_11 / (I_1 rho_11 + I_2 rho_22)  (Bayesian update),
  rho_12 rides purity. Weakly-responding limit |Delta I|<<I_0 -> Bayes1/Bayes2 (diffusion).
Microscopic justification of the whole Bayesian picture via a "pointer" that reads electron
number n and collapses (Meas2-3.tex:981-1000). Not needed for the toy essay; good for a
"where does this come from" aside.

===============================================================================
## G. Feedback and the S/N<=4 statement (0008461)  -- payoff anchor (defer detail to 0002203)
===============================================================================

Feedback to lock symmetric-qubit oscillations (Meas2-3.tex:1490-1515):
  desired: rho_11 = [1+cos(Omega t)]/2, rho_12 = i sin(Omega t)/2, Omega = 2H/hbar
  control: H_fb(t) = H[1 - F * Dphi(t - tau_d)],  Dphi = desired phase - measured phase
  measured phase phi(t) = arctan[ 2 Im rho_12^a / (rho_11^a - rho_22^a) ] on a window-averaged rho^a.
  d := 4 S_0 gamma_d/(Delta I)^2 = ratio of extra-env decoherence to measurement rate.

KEY (Meas2-3.tex:1563-1567): WITHOUT feedback the spectral density of the oscillations has a
peak-to-pedestal ratio with MAXIMUM VALUE 4; with good feedback it sharpens to a near
delta-function. Cited to Kor-osc (0003225) and Kor-Av (0002203). [S/N<=4 CONFIRMED as
peak-to-pedestal of the detector-output spectrum; exact conditions -> read 0002203 next.]

===============================================================================
## H. Memory-reconciliation (my first-message claims vs the TeX)
===============================================================================
1. averaged-current variance S/(2 tau)      -> CONFIRMED  D = S_0/(2 Delta t)   (0008461:747)
2. eta = Gamma_m/Gamma_d, gamma=Gamma_d-Gamma_m -> CONFIRMED  (0008461:582-591)
3. |rho_12|/sqrt(rho_11 rho_22) invariant (ideal) -> CONFIRMED  (0008461:714-719)
4. S/N <= 4 peak-to-pedestal                 -> CONFIRMED as stated (0008461:1566); exact
                                                 conditions pending 0002203.
5. validity window tau_det << tau << (Gamma_m^-1, Omega^-1) -> CONFIRMED (0008461:1084-1099)
All five recalled relations survive contact with the primary source.

===============================================================================
## I. The Korotkov-Averin spectrum and the S/N<=4 bound  (0002203) -- PAYOFF, grounded
===============================================================================

Notation in THIS paper: Hamiltonian H = -(1/2)(eps sigma_z + Delta sigma_x + sigma_z U).
So Delta = tunnel coupling (maps to 2H of 0008461), eps = bias.
  Oscillation frequency: Omega = (Delta^2 + eps^2)^{1/2}   (short.tex:293-294)   [eps=0 => Omega=Delta]
  delta I = current response to the oscillation (== Delta I of 0008461, linear regime)
  S_0 = shot-noise pedestal = 2 e <I> R                     (short.tex:268-269)
  Backaction dephasing rate (short.tex:252, Eq.8):
    Gamma = eV[(delta D)^2 + u^2]/(8 pi D R)
    symmetric coupling u=0 (ideal): Gamma = (delta I)^2/(4 S_0)  == Gamma_m   (short.tex:266-268)
    u =/= 0 (asymmetric = nonideal) makes Gamma larger at fixed response -> extra dephasing.

Output spectral density, eps=0 (short.tex:283-286, Eq.9):
    S_I(omega) = S_0 + Gamma Omega^2 (delta I)^2 / [ (omega^2 - Omega^2)^2 + Gamma^2 omega^2 ]
  Peak sits at omega = Omega, width ~ Gamma.
  Peak height at eps=0 (max amplitude): S_max = (delta I)^2 / Gamma      (short.tex:298-299)

THE BOUND (short.tex:301-304, Eq.10):
    S_max/S_0 = 4 (delta D)^2 / [ (delta D)^2 + u^2 ]  <=  4
  Universal: independent of coupling strength. Saturates at 4 for the symmetric/ideal
  detector (u=0); degrades with asymmetry u (i.e. with non-ideality / eta<1).
  One-line origin: ideal Gamma = (delta I)^2/4S_0 => S_max = (delta I)^2/Gamma = 4 S_0. Done.

Line intensity (short.tex:339-342, Eq.11):  int_0^inf [S_I - S_0] domega/2pi = (delta I)^2/4.
  Depends on coupling; classical oscillations of amplitude delta I/2 would give HALF this
  (a quantum signature; nice aside for the essay).

Strong-measurement / Zeno limit (short.tex:356-362): Gamma >> Omega -> coherent peak turns
into a zero-frequency Lorentzian, incoherent tunneling rate gamma = Delta^2/2Gamma DECREASES
with Gamma (Zeno). Energy relaxation Gamma_e broadens the peak and pushes S_max/S_0 below 4
(short.tex:419-432). => in the essay, S/N approaches 4 for ideal + weak coupling, and any extra
relaxation / asymmetry / non-ideality drops it below 4. This is exactly the eta-dependence to show.

Essay note: this paper derives the spectrum NON-selectively (ensemble/Bloch), and states all
results "can be reproduced within the selective description" (short.tex:88-90). The essay does
the reverse-and-forward: simulate the SELECTIVE (Bayesian) trajectories, ensemble-average the
output periodogram, and RECOVER Eq.9 + the bound. That is the felt cross-check.

===============================================================================
## J. Non-ideality, the clean single-qubit picture  (0211647) -- feeds N4
===============================================================================

Model = ideal detector + classical output noise xi_1 (density S_1) + backaction noise.
For the ESSAY (toy, single qubit, no correlation K, no asymmetry theta) only this is needed:

Ideality (general)   (NONIDEA2.TEX:519-528, eta-gen):
  eta = [ (Delta I)^2 / 4 S_Sigma ] / Gamma_Sigma ,   S_Sigma = S_0 + S_1
  simplest case (extra output noise only, xi_3=0):  eta = S_0/(S_0 + S_1).

Coherence gets a pure-dephasing factor. Extra output noise adds  exp(-gamma_1 tau) with
  gamma_1 = (Delta I)^2 S_1 / (4 S_0 S_Sigma)                 (NONIDEA2.TEX:471)
Uncorrelated backaction xi_3 adds  -gamma_3 rho_12,  gamma_3 = S_3/(4 hbar^2)  (NONIDEA2.TEX:337)
Single-qubit vs ensemble dephasing:  gamma_single = Gamma_Sigma - (Delta I)^2/4 S_Sigma
  => in toy terms:  gamma = Gamma_Sigma (1 - eta).  Ideal eta=1 => gamma=0 => purity preserved.

Frozen-qubit exact Bayesian solution restated (NONIDEA2.TEX:368-390, "Quantum Bayes theorem"):
  rho_11(tau) = [ 1 + (rho_22(0)/rho_11(0)) exp{-(Ibar-I2)^2 tau/S_0}/exp{-(Ibar-I1)^2 tau/S_0} ]^{-1}
  rho_12(tau) = rho_12(0) sqrt(rho_11(tau) rho_22(tau)/(rho_11(0) rho_22(0)))  [times e^{-gamma tau}]
  -- same content as 0008461 sol1/sol2; cleaner algebraic form to code the update.

General N-qubit "Quantum Bayes" update (NONIDEA2.TEX:924-935, QBayes-ent) -- compact statement:
  rho_ij(tau) = rho_ij(0) sqrt(P_i P_j)/sum_k rho_kk P_k ,  P_i = Gaussian(Ibar; I_i, D0), D0=S_0/2tau
  populations Bayes-update; each coherence carries the geometric mean of the two likelihoods.
  THIS one line + a unitary step IS the whole implementable formalism.

Beyond the toy (note only, NOT essay-only): correlated K adds iK[I-I0]rho_12 with bound
Gamma_Sigma >= (Delta I)^2/4S_Sigma + K^2 S_Sigma/4 (NONIDEA2.TEX:685); asymmetric coupling shifts
eps -> eps + theta I_0. Defs eta, eta-tilde, eta-tilde_2 coincide at K=0.

===============================================================================
## K. Modern framing + clean z-update + Rabi spectrum  (1111.4016) -- USE THIS in the essay
===============================================================================

Bloch coordinate:  z = rho_11 - rho_00  (z-basis = measured basis). Output (LesHouches2.tex:302-304):
  I(t) = I_c + (Delta I/2) z(t) + xi(t),   S_xi = S = 2 e I_c (single-sided),  I_c=(I_0+I_1)/2

THE update, cleanest form (LesHouches2.tex:356-406, rho-diag / rho-off-3):
  Ibar_m(t) = (1/t) int_0^t I dt'
  rho_11/rho_00 -> (rho_11/rho_00) * exp[-(Ibar_m - I_1)^2/2D] / exp[-(Ibar_m - I_0)^2/2D],  D=S/2t
  rho_01(t) = rho_01(0) sqrt(rho_00(t) rho_11(t)/(rho_00(0) rho_11(0))) * exp[iK Ibar_m t] * exp[-gamma t]
  Gaussian likelihood P_{|j>}(Ibar_m) = (1/sqrt(2 pi D)) exp[-(Ibar_m - I_j)^2/2D].

TWO back-actions (the pedagogical spine, LesHouches2.tex:186-251):
  "spooky" (informational, non-unitary): the Bayes factor; moves state along MERIDIANS (changes z).
      No physical mechanism, pure information. THIS is what the toy essay is about.
  "realistic" (unitary): the exp[iK Ibar_m t] phase; moves along PARALLELS (rotates phase).
      Physical mechanism (photon-number / asymmetric-detector backaction). K=0 for a SYMMETRIC QPC.
  => the toy essay uses a symmetric detector (K=0): only spooky back-action, purity ride-along exact.

Ensemble dephasing decomposition (LesHouches2.tex:433-437, Gamma):
  Gamma = (Delta I)^2/4S  +  K^2 S/4  +  gamma
          [spooky]           [realistic]  [extra env]
  eta = 1 - gamma/Gamma  (quantum efficiency; purity preserved iff eta=1, i.e. gamma=0)
  eta_tilde = (Delta I)^2/4 S Gamma  = fraction of Gamma from SPOOKY only.

PEAK-TO-PEDESTAL of the Rabi spectral peak = 4 * eta_tilde   (LesHouches2.tex:452-453, 1010-1012)
  => the "S/N <= 4" bound is really "= 4 eta_tilde". It hits 4 only when ALL dephasing is spooky
     (symmetric QPC / broadband, K=0, gamma=0 => eta_tilde=1). Phase-preserving cQED amp: eta_tilde<=1/2
     => ratio <= 2. This SHARPENS the 0002203 statement and gives the eta-knob for N6.

Boxed Rabi output spectrum, exact-resonance, from the appended draft (LesHouches2.tex:1900-1916):
  S_I(omega) = S + 4 eta S * Omega_R^2 Gamma^2 / [ (omega^2 - Omega_R^2)^2 + Gamma^2 omega^2 ]   (peak-1)
  Good-oscillation limit Omega_R >> Gamma (Lorentzian of height 4 eta S at omega=Omega_R):
  S_I(omega) = S + 4 eta S / [ 1 + 4 (omega - Omega_R)^2 / Gamma^2 ]                              (peak-2)
  Consistency check with 0002203 Eq.9: peak height = 4 eta S; for ideal broadband eta=1 => ratio 4. OK.
  (eta here plays the role of eta_tilde for K=0; matches section I.)

POVM/Kraus bridge (LesHouches2.tex:494-517): the measurement (Kraus) operator factorizes
  M(Ibar_m) = exp(-iK Ibar_m t sigma_z/2) [ sqrt(P_{|0>}) |0><0| + sqrt(P_{|1>}) |1><1| ]
  = (unitary "realistic") x (positive "spooky" Bayes sqrt). Nice one-line link to POVM language.

cQED context (not needed for toy, note only): dispersive H = (hbar w_qb/2) sigma_z + hbar w_r a+a
+ hbar chi a+a sigma_z; ensemble dephasing Gamma = 8 chi^2 nbar/kappa (LesHouches2.tex:699);
phase-sensitive amp (info quadrature) = pure spooky, phase-preserving = spooky+realistic split 50/50.
Feedback controller (2nd channel) Delta(w_qb - w_R) = -K[Q(t)-<Q>] (LesHouches2.tex:1269).

WHICH NOTATION FOR THE ESSAY: use 1111.4016's z=rho_11-rho_00, I(t)=I_c+(Delta I/2)z+xi, S single-sided,
the rho-diag/rho-off update, and Gamma/eta/eta_tilde. Symmetric detector K=0 keeps it pure-spooky.
Map to 0008461 (rho_11 rho_22, S_0) only when citing the exact-solution anchors.

================================================================================
### PHASE 2 EXTENSION (sibling essays). New refs under scratchpad/refs/<id>/.  ###
================================================================================
Read-all done for 7 more papers (non-commuting: 1702.08077, 1710.05249, 1608.06652;
entanglement: 1402.1868, 1603.09623; QEC capstone: 1612.02096, 1910.08272). Every
equation below carries a <file>:<line> anchor into those TeX sources.

===============================================================================
## L. Simultaneous non-commuting z-and-x Bayesian update + output correlators ##
    (1702.08077 = XZcorr-arXiv.tex; 1710.05249 = MTC-Arxiv.tex; 1608.06652 = Hacohen-Gourgy)
===============================================================================

L0. THE SCHEME. One qubit, TWO linear detectors measuring simultaneously sigma_z and
sigma_phi = sigma_z cos phi + sigma_x sin phi (angle phi between the two Bloch axes).
Rabi-rotated effective qubit; phase-sensitive amp on the optimal quadrature => NO phase
backaction (K=0), only informational backaction. Bloch rho = (1 + x sigma_x + y sigma_y + z sigma_z)/2.

L1. Two output signals (XZcorr:64-69, eq:outputs-Iz/Iphi):
  I_z(t)     = Tr[sigma_z rho]     + sqrt(tau_z) xi_z(t)
  I_phi(t)   = Tr[sigma_phi rho]   + sqrt(tau_phi) xi_phi(t)
  white, uncorrelated: <xi_z xi_z>=<xi_phi xi_phi>=delta(t-t'), <xi_z xi_phi>=0  (XZcorr:75)
  tau_z, tau_phi = "measurement" (collapse) times for informational SNR=1 per channel.
  quantum efficiencies eta_z = 1/(2 tau_z Gamma_z), eta_phi = 1/(2 tau_phi Gamma_phi) (XZcorr:92-93).
  Experiment (Hacohen-Gourgy): eta_z=0.49, eta_phi=0.41.

L2. Ito SDEs, measurement only (XZcorr:80-90, eq:Ito-x/y/z) -- THE non-commuting update:
  x. = -Gamma_z x - Gamma_phi cos phi (x cos phi - z sin phi)
        - tau_z^{-1/2} xz xi_z - tau_phi^{-1/2}[xz cos phi - (1-x^2) sin phi] xi_phi
  y. = -(Gamma_z+Gamma_phi) y - tau_z^{-1/2} yz xi_z - tau_phi^{-1/2} y[z cos phi + x sin phi] xi_phi
  z. =  Gamma_phi sin phi (x cos phi - z sin phi)
        + tau_z^{-1/2}(1-z^2) xi_z + tau_phi^{-1/2}[(1-z^2) cos phi - xz sin phi] xi_phi
  Derivation: single-sigma_z Ito (XZcorr:497-504), rotate to sigma_phi basis (XZcorr:507-519),
  ADD the two channels' terms (XZcorr:526, uncorrelated noises).
  Extra (non-measurement) evolution: H = hbar Omega_tilde_R sigma_y/2, T1/T2 decoherence (XZcorr:97-101):
    x. += Omega_tilde_R z - gamma x, y. += -T2^{-1} y, z. += -Omega_tilde_R x - gamma z,
    gamma = (T1^{-1}+T2^{-1})/2.

L3. CLEANEST GENERAL FORM (Bloch, unital, K=0)  (MTC:94-97, eq:Bayesian-eq):
    r. = Lambda_ens (r - r_st) + sum_ell [ n_ell - (n_ell . r) r ] / sqrt(tau_ell) * xi_ell(t)
  n_ell = ell-th measurement axis. The back-action term [n_ell - (n_ell.r) r] pulls r toward n_ell.
  Ensemble (Lindblad) part: L_m[rho] = sum_ell Gamma_ell [sigma_ell rho sigma_ell - rho]/2,
  Gamma_ell = 1/(2 eta_ell tau_ell)  (MTC:104). Unital <=> r_st = 0 (MTC:106).

L3b. GENERAL DENSITY-MATRIX form, arbitrary Hermitian A_ell (MTC:409,414, eq:Bayes-general/-2):
  rho. = L[rho] + sum_ell [ A_ell rho + rho A_ell - 2 rho Tr(A_ell rho) ] / sqrt(2 S_ell) * xi_ell(t),
  I_ell = Tr(A_ell rho) + sqrt(S_ell/2) xi_ell,  L_m = sum_ell [A_ell rho A_ell -(A_ell^2 rho+rho A_ell^2)/2]/(2 eta_ell S_ell).

L4. DISCRETE KRAUS UPDATE (the implementable form; reduces to core essay's single-basis update).
  For a measured Pauli A (eigenvalues +/-1), reading Ibar over step Dt, window variance D:
    M(Ibar, A) = cosh(Ibar/2D) I + sinh(Ibar/2D) A        [= exp(Ibar A /2D), since A^2=I]
    rho -> M rho M / Tr(M rho M).
  For A=sigma_z: M diagonal, populations Bayes-update + coherence rides -> IS the core essay's bayesUpdate.
  For TWO non-commuting A's (sigma_z, sigma_x): M_z, M_x DO NOT COMMUTE, so M_x M_z != M_z M_x and
  there is NO closed-form population update. Carry the full 2x2 rho. Symmetric (Strang) split
  M_z^{1/2} M_x M_z^{1/2} is 2nd-order in Dt; agrees with the Ito SDE as Dt->0.
  (Hacohen-Gourgy measurement operator, HG:369-373, eq:XYMeasOpp: Omega(V)=exp[sum_i -(Gamma_i eta_i/2)(V_i - sigma_delta,i)^2 dt],
   rho(t+dt)=E_{1-eta}[Omega rho Omega^dag/Tr]. Same content; the (V-A)^2 Gaussian Kraus.)
  Honest draw: Ibar_A = RandomChoice[{(1+<A>)/2,(1-<A>)/2}->{+1,-1}] + N(0, sqrt D). D = 1/(2 Gamma Dt) ideal.
  Two-channel SME (HG:361-365, eq:SMEFinal): drho = sum_i (Gamma_i/2) D[sigma_delta_i] rho dt + sqrt(Gamma_i eta_i/2) H[sigma_delta_i] rho dW_i.

L5. QUALITATIVE DYNAMICS (Hacohen-Gourgy experiment):
  phi=0 (commuting): standard collapse to the two z-poles.
  0<phi<pi/2: state localizes to finite regions near +/- axes, no point collapse.
  phi=pi/2 (orthogonal, maximally incompatible): ISOTROPIC persistent diffusion, uniform random walk
    on the Bloch sphere, NO imprint of the axes, no collapse. Azimuthal diffusion coeff 2Gamma (great circle).
  Steady-state radius r = sqrt(eta) for orthogonal ideal-ish detectors (HG:164); most-likely steady purity 0.89.
  Measurement-induced disturbance / uncertainty floor (HG:169-173, eq:distmap):
    Tr[drho^dag drho] = (Var(sigma_delta,1) Gamma_1 eta_1 + Var(sigma_delta,2) Gamma_2 eta_2) dt
                       >= |<[sigma_delta1, sigma_delta2]>| sqrt(eta_1 eta_2 Gamma_1 Gamma_2) dt.
    Sum-of-variances form (MacCone-Pati), never trivial for non-commuting => state must diffuse forever.

L6. OUTPUT SIGNAL CORRELATORS K_ij(tau) = <I_j(t1+tau) I_i(t1)>, tau>0 (XZcorr:107-115).
  Collapse recipe (XZcorr:117-123, 631-637): replace continuous meas at t1 by PROJECTIVE meas of sigma_i
  (result +/-1 w.p. {1 +/- Tr[sigma_i rho(t1)]}/2, collapse to |1_i> or |0_i>), evolve rho_av (noiseless)
  by the ensemble eqs, then measure sigma_j:
    K_ij(tau) = Tr[sigma_j rho_av(t1+tau|1_i)](1+Tr[sigma_i rho(t1)])/2 - Tr[sigma_j rho_av(t1+tau|0_i)](1-...)/2.
  For unital evolution simplifies to K_ij(tau)=Tr[sigma_j rho_av(tau|1_i)], INDEPENDENT of the state (XZcorr:676).
  Closed forms (XZcorr:135-150, eq:K-zz/K-zphi/Gamma_pm):
    K_zz(tau)  = (1/2)[1 + (Gamma_z+cos2phi Gamma_phi)/(Gp-Gm)] e^{-Gm tau}
               + (1/2)[1 - (Gamma_z+cos2phi Gamma_phi)/(Gp-Gm)] e^{-Gp tau}
    K_zphi(tau)= [(Gamma_z+Gamma_phi)cos phi + 2 Omega_tilde_R sin phi]/[2(Gp-Gm)] (e^{-Gm tau} - e^{-Gp tau})
               + (cos phi/2)(e^{-Gm tau} + e^{-Gp tau})
    Gamma_pm = [Gamma_z+Gamma_phi +/- sqrt(Gamma_z^2+Gamma_phi^2+2 Gamma_z Gamma_phi cos2phi - 4 Omega_tilde_R^2)]/2
               + (T1^{-1}+T2^{-1})/2.
  Rotational symmetry: K_phiphi, K_phiz from K_zz, K_zphi via Gamma_z<->Gamma_phi, phi->-phi (XZcorr:152).
  Do NOT depend on tau_z,tau_phi (=> not on eta): non-ideal = ideal + extra output noise, only hits K_ii(0) (XZcorr:155).

L7. SPECIAL CASES (XZcorr:157-166):
  (i) small time: K_zz(+0)=1, K_zphi(0)=K_phiz(0)=cos phi. [the felt cross-correlator-at-zero payoff]
  (iii) Omega_tilde_R=T1^{-1}=T2^{-1}=0: phi=0 full corr (K_zphi=K_zz=1); phi=pi (anti, =-1);
        phi=pi/2 no corr (K_zphi=0), K_zz=e^{-Gamma_phi tau}, K_phiphi=e^{-Gamma_z tau}.
  (iv) Omega_tilde_R=0: cross-correlator symmetric K_zphi=K_phiz.
  Zeno pinning |phi|<<1: K_zphi(tau) ~ exp(-2 Gamma_jump tau), Gamma_jump=(phi^2 Gamma_z Gamma_phi + Omega_tilde_R^2)/[2(Gamma_z+Gamma_phi)]+... (XZcorr:162-165).
  ANTISYMMETRIZED cross-correlator estimates a small residual Rabi (XZcorr:204, eq:Kzphi-antisym):
    K_zphi(tau) - K_phiz(tau) = [2 Omega_tilde_R sin phi/(Gp-Gm)](e^{-Gm tau} - e^{-Gp tau}).  Sensitive Omega_tilde_R probe.
  Multi-time (MTC): even-N factorizes into product of 2-time correlators (MTC:168), odd-N adds <I(t1)> (MTC:173);
    2-time K_{ij}(t_i,t_k)=n_k[exp(int Lambda_ens) n_i] (MTC:178-180). Verified vs experiment (3-time K_phizphi, 4-time K_zphizphi).

L8. CAPSTONE LINK (Mads's own line). 4-qubit Bacon-Shor code (1612.02096): the four gauge operators
  X12,X34,Z13,Z24 continuously measured reduce to a GAUGE QUBIT under simultaneous X-and-Z measurement,
  diffusing on the great circle (4qubit:576, 616-626), errors read from cross-correlators
  <I_X12 I_X34>, <I_Z13 I_Z24> = exp(-2 Gamma_m|t1-t2|) (4qubit:771; same-time = +1, NOT x_g^2, XZcorr:155 physics).
  9-qubit code (1910.08272): 12 gauge ops, 4 gauge qubits, gauge SME (main:471),
  triple cross-correlators (main:618-635), logical error rate ~ Gamma_d^1.88, crossover Gamma_d~10^-3/tau_coll (main:1236).

===============================================================================
## M. Two-qubit measurement-induced entanglement trajectories  (1603.09623 = Chantasri; 1402.1868 = Roch) ##
===============================================================================

M0. THE SCHEME. Two remote qubits, each in its own cavity, jointly measured in a "bounce-bounce" geometry
  so the probe reflects off cavity 1 then cavity 2. Tuned so that the phase shifts for |01> and |10> are
  EQUAL: the joint measurement is a HALF-PARITY meter distinguishing three groups {|00>}, {|11>}, {odd = |01>,|10>}
  but NOT within the odd subspace. Post-select the odd outcome -> Bell state (|01>+|10>)/sqrt2. No local coupling;
  entanglement is purely measurement-induced.

M1. THE TWO-QUBIT BAYESIAN UPDATE (Chantasri, Twoqubit:136, eq-bayes) -- clean generalization of the
  single-qubit purity ride-along:
    rho_ij(t) = rho_ij(0) sqrt( p(V_t|i) p(V_t|j) ) e^{-gamma_ij t} / sum_k rho_kk(0) p(V_t|k)
  i,j,k in {1,2,3,4} = {|00>,|01>,|10>,|11>}. POPULATIONS Bayes-update; each COHERENCE carries the
  GEOMETRIC MEAN sqrt(P_i P_j) of the two likelihoods, decaying by e^{-gamma_ij t}. (Same one-line N-qubit
  "quantum Bayes" as sheet J, NONIDEA2.TEX:924-935.)
  Gaussian likelihoods (Twoqubit:143): p(V_t|i) = (t/pi s)^{-1/2} exp{-(V_t - dv_i)^2 t/s}, s = 1/(2 eta_m),
  eta_m ~ 0.22 (Chantasri exp). Half-parity: dv_2 ~ dv_3 ~ 0, -dv_1 ~ dv_4 ~ dv (Twoqubit:143).
  Extra dephasing gamma_ij ~ (eta_m^{-1}-1)(dv_i-dv_j)^2/(4s) (Twoqubit:143); rho_23 = odd-parity coherence
  survives (dv_2=dv_3 => gamma_23 minimal), all others damp fast.

M2. CONCURRENCE (X-state; only rho_23 and diagonals survive) (Twoqubit:153, eq-conc0):
    C(t) = 2 max{ 0, x_5 - sqrt(x_1 x_4) },  x_1=rho_11(=|00>), x_4=rho_44(=|11>), x_5=|rho_23|(=|rho_{01,10}|).
  (Roch simplified concurrence, Roch:344: C = 2 max(0, |rho_{01,10}| - sqrt(rho_00,00 rho_11,11)).)

M3. CONCURRENCE-READOUT CLOSED FORM, perfectly-symmetric half-parity, product-of-x-states init (Twoqubit:177, eq-conc3):
    C_ps,x(V_t, t) = [ e^{-gamma t} - e^{-dv^2 t/s} ] / [ 1 + cosh(2 V_t dv t/s) e^{-dv^2 t/s} ].
  MAX over readout (at V_t=0) = the sharp upper concurrence BOUND (Twoqubit:207, eq-maxconc):
    C_max,ps,x(t) = [ e^{-gamma t} - e^{-dv^2 t/s} ] / [ 1 + e^{-dv^2 t/s} ].
  Two competing rates: extra dephasing gamma (kills the numerator) vs measurement rate dv^2/s (grows it,
  then the -e^{-dv^2 t/s} decays it). Rises from 0, peaks, decays. The concurrence distribution p_C,t(c) has a
  sharp cutoff at C_max(t): entanglement cannot be created faster than this even in rare records.

M4. MOST-LIKELY-PATH (Twoqubit:231-395). Log-likelihood log P = S0 - int (1/s) sum_k (v_t - dv_k)^2 x_k dt' (Twoqubit:244).
  Optimal readout v_t is CONSTANT in time for no-drive (Twoqubit:343). THREE branches: v_t=dv_{2,3}=0 -> high-concurrence
  (odd subspace), v_t=+/-dv -> low-concurrence (collapse to |00> or |11>). High branch analytic (Twoqubit:271-277):
    x_{1,4} ~ x0_{1,4} e^{-dv^2 t/4s}, x_{2,3} ~ x0_{2,3}, x_5 ~ x0_5 e^{-gamma t}, norm N=(1-x0_1-x0_4)+(x0_1+x0_4)e^{-dv^2 t/4s};
    its concurrence coincides with C_max(t). Bimodal time-to-max-concurrence.
  PARITY-METER (full parity, distinguishes even from odd but not within) exact solutions (Twoqubit:369-395):
    x_o(t)=x0_o e^{-lambda t}/[1 - x0_o(1-e^{-lambda t})], lambda = 2 v_t dv/s - dv^2/s; x_e=1-x_o;
    C(t)=2 max{0, [x0_5 |e^{-(gamma+lambda)t}| - sqrt(x0_1 x0_4)]/|1-(x0_2+x0_3)(1-e^{-lambda t})|}.

M5. ROCH cascaded-SME anchors (1402.1868): measurement rate Gamma_meas = (1/2) eta_meas eta_loss |alpha_in|^2 sin(2 Dphi)^2 (Roch:181).
  Diagonal classical Bayes rho^fin_ij,ij = rho^in_ij,ij p_sel(i,j)/p_ent (Roch:406). Odd-parity coherence ride-along
  |rho^fin_01,10| = |rho^in_01,10| sqrt(rho^fin_01,01 rho^fin_10,10)/sqrt(rho^in_01,01 rho^in_10,10) * exp[-distinguishability dephasing] (Roch:411).
  Concurrence peaks ~0.35 in experiment (SNR-limited vs decoherence trade-off, Roch:194). eta_meas=0.4, eta_loss=0.81.
  Cascaded (one-way) SME + polaron transformation (Roch:452-549); X-state so only rho_00,00, rho_11,11, rho_01,10 matter for C.

M6. ESSAY BUILD NOTE. The two-qubit essay implements M1 directly: 4x4 rho, populations Bayes + geometric-mean
  coherence, gamma_23 the only surviving off-diagonal decay. Recover C_max(t) (M3) from an ensemble of records;
  show the three most-likely branches (M4); the bounded, bimodal concurrence distribution. Symmetric half-parity,
  product-of-x init x0_i=1/4. Ideal-ish: gamma=extra dephasing knob (eta_m). Do NOT need the cavity/polaron layer.

===============================================================================
## Experimental Anchors (grounded from full-text arXiv reads)
===============================================================================
Each subsection below was extracted from the TeX source under scratchpad/korotkov_refs/<id>/,
read end to end (main text AND supplement). Anchors are <file>:<line> into those sources plus the
paper's own equation/section labels. Numbers are quoted as the paper states them (with units). All
five e-prints came down as TeX tarballs (no PDF fallback). NOTE: 1612.02096 is a THEORY paper
(analytics + Monte Carlo), not a lab experiment; it is anchored here for its formal reduction and
its stated results, not for a measured number.

-------------------------------------------------------------------------------
### N1. 1305.7270  Murch-Weber-Macklin-Siddiqi (single-qubit trajectories)
-------------------------------------------------------------------------------
K. W. Murch, S. J. Weber, C. Macklin, I. Siddiqi, "Observing single quantum trajectories of a
superconducting qubit", Nature 502, 211 (2013). arXiv:1305.7270.
(main = sqwm_v5_arxiv.tex; supplement = z_sqwm_sup.tex)

SETUP. 3D transmon dispersively coupled to a copper waveguide cavity, H_int = -hbar chi a+a sigma_z
(sqwm:75). A near-quantum-limited lumped-element Josephson parametric amplifier (LJPA) run
phase-sensitively amplifies ONE quadrature of the reflected tone and de-amplifies the other, so the
experimenter chooses WHICH backaction the qubit feels. Two choices:
  "Z-measurement": amplify the quadrature carrying qubit-state info (the state-dependent cavity phase
     shift). State is driven toward the poles/eigenstates; trajectory rides a MERIDIAN of the Bloch
     sphere (sqwm:78,107; z^Z=tanh, x^Z=sqrt(1-z^2), Fig.3a,b).
  "phi-measurement": amplify the quadrature carrying intracavity photon-number info (amplitude). The
     superposition phase diffuses (AC-Stark), no projection; trajectory confined to the EQUATOR
     (sqwm:78,170). Physical-quadrature identity: qubit-state info sits in the quadrature IN
     QUADRATURE with the input tone, photon-number info sits IN PHASE with it (sqwm:75). Abstract's
     compressed wording: "selectively measure either the phase or amplitude of the cavity field, and
     thereby confine trajectories to either the equator or a meridian" (sqwm:57).
Bayesian update: z^Z = tanh(V_m S / 2 Delta V), x^Z = sqrt(1-(z^Z)^2) e^{-gamma tau}  (sqwm:117, eq:zz);
  x^phi = cos(S V_m/(2 Delta V)) e^{-gamma tau}, y^phi = -sin(S V_m/(2 Delta V)) e^{-gamma tau} (sqwm:130-131).

NUMBERS (with anchors):
- Quantum efficiency eta = 0.49, from linear fit of S vs nbar (sqwm:119, also 145,176; Fig.1e inset).
  Decomposed eta = eta_col * eta_amp, eta_col = 0.72 (collection), eta_amp = 0.68 (amplifier) [z_sqwm_sup:140].
  "among the highest reported for a continuous variable" (sqwm:176).
- Dimensionless measurement strength S = 64 tau chi^2 nbar eta / kappa = Delta V^2 / sigma^2 (SNR) (sqwm:119).
- Dephasing rate gamma = 8 chi^2 nbar (1-eta)/kappa + 1/T_2* ; first term = measurement-induced dephasing
  from the 1-eta undetected fraction (sqwm:122).
- Dispersive coupling chi/2pi = -0.49 MHz; cavity decay kappa/2pi = 10.8 MHz (Fig.1 caption, sqwm:67).
  (A commented Methods draft line lists -0.52 MHz; the live Fig.1 value is -0.49 MHz.)
- LJPA: 10 dB gain, 20 MHz instantaneous bandwidth; two-junction SQUID from 2 uA junctions, 3 pF shunt;
  pumped by sidebands +/-300 MHz about cavity (sqwm:67,182).
- Qubit: E_c/h = 200 MHz, E_J/h = 11 GHz, omega_q/2pi = 3.999 GHz; cavity 6.8316 GHz; T_2* = 20 us (sqwm:122,180).
- Weak-measurement window: integrate V_m for 1.8 us; trajectory time step tau_{i+1}-tau_i = 16 ns;
  ~10^5 repetitions per point (sqwm:107,144,164).
- Herald/tomography readout S = 42, 800 ns; pi/2 rotation 16 ns; readout fidelity 95%; ~4% of shots
  outside {|0>,|1>} discarded (sqwm:184).
- TOMOGRAPHIC VALIDATION: tomographically reconstructed <sigma_i> vs V_m match theory with eta=0.49 for
  both Z (nbar=0.4, S=3.15) and phi (nbar=0.46, S=3.62) measurements, "excellent agreement" (sqwm:144-145,
  Fig.2c,d); single-shot reconstructed trajectories track the tomographic ensemble (Fig.3).
- EPR-steering witness S^{Z,phi} = <z^Z>^2 + <x^phi>^2 + <y^phi>^2 <= 1; experiment gives ~0.8; 1.4 dB
  loss between cavity port and amp; eta_amp>0.5 => most backaction is quadrature-dependent (z_sqwm_sup:136-140).
grounds: single-qubit core (essay Sections 3-6) and the Section-8 informational-vs-phase quadrature split
  (Z-measurement = pure informational backaction on a meridian; phi-measurement = phase backaction on the equator).

-------------------------------------------------------------------------------
### N2. 1402.1868  Roch et al. (remote measurement-induced entanglement)
-------------------------------------------------------------------------------
N. Roch, M. E. Schwartz, F. Motzoi, C. Macklin, R. Vijay, A. W. Eddins, A. N. Korotkov, K. B. Whaley,
M. Sarovar, I. Siddiqi, "Observation of measurement-induced entanglement and quantum trajectories of
remote superconducting qubits", Phys. Rev. Lett. 112, 170501 (2014). arXiv:1402.1868.
(measurement-based-entanglement-V2.tex; main + Supplemental Material read in full)

SETUP. Two 3D-transmon qubits, each in its own copper cavity, connected via two microwave circulators
and 1.3 meters of ordinary coaxial cable so a probe reflects off cavity 1 THEN cavity 2 (cascaded /
"bounce-bounce", meas:148,150,565). Single-qubit reflection coefficient r^+/- (meas:158). Engineer the
cavities/dispersive shifts so the RELATIVE phase shifts match, Delta_phi_1 = Delta_phi_2: then |01> and
|10> acquire the SAME phase shift delta_1+delta_2 and become indistinguishable, while |00>,|11| stay
separable = a HALF-PARITY meter (meas:161,163). Post-select the odd-parity outcome -> Bell state
(|01>+|10>)/sqrt2. Homodyne via a phase-sensitive LJPA. Validated Fig.2b: |00>,|11> histograms well
separated, |01>,|10> fully overlapping.

NUMBERS (with anchors):
- Physical separation: 1.3 meters of coaxial cable between the two cavities (meas:148,150,565).
- PEAK CONCURRENCE = 0.35, reached at intermediate t_m where the SNR-improvement rate ~ Gamma_loss
  (meas:194). "comparable to optical remote-entanglement; but the creation RATE is orders of magnitude
  higher, Gamma_creation/2pi = 1 kHz" (meas:194). Models predict ~70% is reachable with better hardware (meas:196).
- Measurement efficiency eta_meas = 0.4 +/- 0.10 (meas:161,394; Table S1, meas:634).
- Inter-cavity power-transfer efficiency eta_loss ~ 0.81 +/- 0.05 (main text meas:161; Table S1 meas:633).
  [The Simplified-Theory section rounds it to 0.75 in prose (meas:394); the calibrated value is 0.81 +/- 0.05.]
- Entanglement-generation rate Gamma_meas = (1/2) eta_meas eta_loss |alpha_in|^2 sin(2 Delta_phi)^2 (meas:181, Eq.3).
  nbar_1 = 1.2 => Gamma_meas/2pi ~ 210 kHz, tau_meas = 1/Gamma_meas ~ 750 ns (meas:184).
- SNR ~ 2|alpha_in| sin(2 Delta_phi) sqrt(eta_loss eta_meas t_m) (meas:189, Eq.5).
- Loss-induced dephasing of qubit 1: Gamma_loss ~ 2(1-eta_loss)|alpha_in|^2 sin(Delta_phi)^2 (meas:193, Eq.6).
- Three regimes vs t_m: SNR-dominated t_m < 0.75 tau_meas; stabilization 0.75-1.25 tau_meas; decoherence
  decay t_m > 1.25 tau_meas (meas:192-194).
- Concurrence (X-state) C = 2 max(0, |rho_01,10| - sqrt(rho_00,00 rho_11,11)) (meas:344, Eq.concurrence-simple;
  main-text simplified form meas:188).
- HALF-PARITY / joint dispersive readout: choose omega_m where B_out^{(01)} = B_out^{(10)}; found by tuning
  omega_m until the |01>,|10> single-shot histograms fully overlap (meas:384,581); omega_m/2pi = 7.19326 GHz (meas:585).
- Post-selection p_ent = 10% (compensates finite eta; 50% for perfectly separated histograms) (meas:186,401).
- 8,000 reps per tomographic rotation per t_m; 30 tomography rotations; 17 data sets for error bars;
  tomographic reconstruction at t_m = 0.65 us; state-prep+tomo fidelity 98.8% (meas:186,163,651).
- Projective post-selection of |00>: nbar_1 = 6.2, 1 us readout (meas:186).
- QUANTUM-BAYESIAN VALIDATION FOR A CASCADED SYSTEM (explicit): abstract "confirming the validity of the
  quantum Bayesian formalism for a cascaded system" (meas:135); body "demonstrates the validity of quantum
  trajectory theories for cascaded quantum systems" (meas:207). Diagonal classical Bayes update
  rho^fin_ij,ij = rho^in_ij,ij p_sel(i,j)/p_ent (meas:406, Eq.rho_fin); odd-parity coherence rides the
  geometric mean |rho^fin_01,10| = |rho^in_01,10| sqrt(rho^fin_01,01 rho^fin_10,10)/sqrt(rho^in_01,01 rho^in_10,10)
  x exp[-distinguishability dephasing] (meas:411, Eq.rho-01,10-num).
- THREE-BRANCH / most-likely-path structure: single trajectories are projected onto the Bell state OR onto
  the non-entangled |00> or |11> (meas:207, Fig.4b). Cascaded SME + polaron transform (meas:452-549); SME
  simulated with 10^5 Wiener instances at 1 ns step (meas:202,704).
- Table S1 (meas:626-635): omega_q/2pi = 4.31143 / 4.46143 GHz; omega_r/2pi = 7.1864 / 7.1984 GHz;
  kappa/2pi = 18.5 / 21 MHz; chi/2pi = 1.275+/-0.025 / 1.085+/-0.035 MHz; T_1 = 27+/-5 / 20+/-3 us;
  T_2* = 16+/-3 / 12+/-2 us (qubits 1 / 2).
grounds: two-qubit generalization (essay Section M / two-qubit essay): 4x4 rho, populations Bayes-update,
  odd-parity coherence ride-along, half-parity meter, three most-likely branches, cascaded-system validation.

-------------------------------------------------------------------------------
### N3. 1608.06652  Hacohen-Gourgy et al. (simultaneous non-commuting observables)
-------------------------------------------------------------------------------
S. Hacohen-Gourgy, L. S. Martin, E. Flurin, V. V. Ramasesh, K. B. Whaley, I. Siddiqi, "Dynamics of
simultaneously measured non-commuting observables", Nature 538, 491 (2016). arXiv:1608.06652.
(Hacohen-Gourgy_Non-Commuting_Measurements-ArxivFinal.tex; main + Methods read in full)

SETUP. ONE transmon dispersively coupled to a MULTIMODE aluminum cavity; each of the two lowest cavity
modes is a separate readout channel monitored by its own phase-sensitive LJPA (~90% of each signal routed
to its amp) (HG:91,94,101,201). "Single-quadrature measurement" (SQM): drive Rabi oscillations at
Omega_R/2pi = 40 MHz so the qubit becomes an effective low-frequency qubit; apply a pair of sidebands at
+/-Omega_R to each mode; the RELATIVE sideband phase delta sets the measured axis
sigma_delta = sigma_x cos delta + sigma_y sin delta (HG:107,122,125,317). Two modes therefore measure two
independently-chosen axes sigma_{delta_1}, sigma_{delta_2}; the ANGLE between them is delta_2 - delta_1.
Effective-qubit rate Gamma eta = 2 chi^2 abar_0^2 eta/kappa (HG:125); Gamma = 8 g_tilde^2/kappa =
2 chi^2 nbar_0/kappa with g_tilde = chi abar_0/2 (HG:356,318).

COLLAPSE-TO-DIFFUSION TRANSITION as the axes go parallel -> orthogonal (HG:160-165, Fig.3a):
  delta_1 = delta_2 (angle 0, commuting): standard collapse to the two poles.
  0 < angle < 90 deg: state localizes to finite regions near the +/- axes, NO point collapse.
  angle = 90 deg (sigma_x and sigma_y, maximally incompatible): the axes leave NO imprint; ISOTROPIC
    PERSISTENT DIFFUSION = uniform random walk on the Bloch sphere, no collapse.

NUMBERS (with anchors):
- Rabi drive Omega_R/2pi = 40 MHz, stabilized to within 10 kHz by feedback; LO leakage suppressed to 10^-4
  photon level (HG:107,151,220).
- Quantum efficiencies eta_1 = 0.41 (mode 1), eta_2 = 0.49 (mode 2) (HG:101,201). Mode-1 dephasing rate
  Gamma_1/2pi = 122 kHz (HG:101, Fig.1c). eta_i = (mu_up-mu_down)^2/(8 tau sigma^2 Gamma_i) (HG:239, from Korotkov 1111.4016).
- Steady-state radius (orthogonal, ideal-ish) r = sqrt(eta) (HG:164). Most-likely steady-state purity P = 0.89
  (HG:135,164-165).
- Azimuthal diffusion on a ring (inner radius 0.86, outer 0.92): measured variance slope 1.4 us^-1 vs
  expected 1.5 us^-1 for a perfect random walk; ~10% measurement-rate uncertainty (HG:136-137).
- Measurement-induced disturbance / uncertainty floor (HG:169-173, eq:distmap):
    Tr[drho^dag drho] = (Var(sigma_delta,1) Gamma_1 eta_1 + Var(sigma_delta,2) Gamma_2 eta_2) dt
                       >= |<[sigma_delta1, sigma_delta2]>| sqrt(eta_1 eta_2 Gamma_1 Gamma_2) dt.
  SUM-of-variances form (MacCone-Pati), never trivial for non-commuting => state must diffuse forever (HG:175).
  disturbance mapped over Bloch sphere for angles delta_2-delta_1; no zero-disturbance point once non-commuting
  (measured over a 64 ns interval, HG:167, Fig.4).
- Two-channel SME (HG:361-365, eq:SMEFinal): drho = sum_i (Gamma_i/2) D[sigma_delta_i] rho dt
    + sqrt(Gamma_i eta_i/2) H[sigma_delta_i] rho dW_i;  V_i dt = <sigma_delta_i> dt + dW_i/sqrt(2 eta_i Gamma_i).
  Discrete Kraus / measurement operator (HG:371, eq:XYMeasOpp): Omega(V) = exp[ sum_i -(Gamma_i eta_i/2)(V_i - sigma_delta,i)^2 dt ],
    rho(t+dt) = E_{1-eta_i}[ Omega rho Omega^dag / Tr ].
- Trajectory update interval Delta t = 16 ns; tomography every 200 ns; 16 initial states reconstructed from
  ~10,000 trajectories each; 57% of tomographic points within error bars (HG:153,129,158,244).
- Methods hardware (HG:198): E_c/h = 220 MHz, omega_q/2pi = 4.262 GHz, T_1 = 60 us, T_2echo = 40 us, Rabi
  decay 25 us; cavity 81x51x20 mm Al; modes omega_1/2pi = 6.666 GHz, omega_2/2pi = 7.391 GHz;
  kappa_1/2pi = 7.2 MHz, kappa_2/2pi = 4.3 MHz; chi_1/2pi = 0.18 MHz, chi_2/2pi = 0.23 MHz; LJPA gains 15/18 dB.
grounds: two-axes generalization (essay Section L). This is the experimental collapse->localized->isotropic-
  diffusion transition and the r = sqrt(eta), purity 0.89 steady state that the two-axis essay reproduces.

-------------------------------------------------------------------------------
### N4. 1702.08077  Atalaya-Hacohen-Gourgy-Martin-Siddiqi-Korotkov (output correlators)
-------------------------------------------------------------------------------
J. Atalaya, S. Hacohen-Gourgy, L. S. Martin, I. Siddiqi, A. N. Korotkov, "Correlators in simultaneous
measurement of non-commuting qubit observables", arXiv:1702.08077 (Phys. Rev. Lett. format; own journal
line not present in fetched source). (XZcorr-arXiv.tex; main + Supplemental Material read in full)

SETUP. Same apparatus as N3 (Rabi-rotated effective qubit, stroboscopic sideband measurement, two lowest
cavity modes). Two linear detectors continuously and simultaneously measure sigma_z and
sigma_phi = sigma_z cos phi + sigma_x sin phi, phi = angle between the two Bloch axes (XZcorr:36,53).
Phase-sensitive amps on the optimal (informational) quadrature => NO phase backaction, only informational
backaction (XZcorr:78). Output signals (XZcorr:64-68):
  I_z(t)   = Tr[sigma_z rho]   + sqrt(tau_z) xi_z(t),
  I_phi(t) = Tr[sigma_phi rho] + sqrt(tau_phi) xi_phi(t),  white, uncorrelated (XZcorr:75).
Quantum efficiencies eta_z = 1/(2 tau_z Gamma_z), eta_phi = 1/(2 tau_phi Gamma_phi) (XZcorr:92).

CLOSED-FORM CORRELATORS (this is the deliverable for the essay's kzzCF / kzpCF).
Definition K_ij(tau) = <I_j(t1+tau) I_i(t1)>, tau>0 (XZcorr:111, eq:corr-def). Results (XZcorr:136-150):

  K_zz(tau) = (1/2)[ 1 + (Gamma_z + cos(2 phi) Gamma_phi)/(Gamma_+ - Gamma_-) ] exp(-Gamma_- tau)
            + (1/2)[ 1 - (Gamma_z + cos(2 phi) Gamma_phi)/(Gamma_+ - Gamma_-) ] exp(-Gamma_+ tau)

  K_zphi(tau) = [ (Gamma_z + Gamma_phi) cos(phi) + 2 OmegaR_t sin(phi) ] / [ 2 (Gamma_+ - Gamma_-) ]
                  * ( exp(-Gamma_- tau) - exp(-Gamma_+ tau) )
              + (cos(phi)/2) * ( exp(-Gamma_- tau) + exp(-Gamma_+ tau) )

  Gamma_pm = { Gamma_z + Gamma_phi
               +/- sqrt[ Gamma_z^2 + Gamma_phi^2 + 2 Gamma_z Gamma_phi cos(2 phi) - 4 OmegaR_t^2 ] } / 2
             + (T1^-1 + T2^-1)/2

  where OmegaR_t = tilde-Omega_R = Omega_R - Omega_rf is the small residual Rabi mismatch (H = hbar OmegaR_t sigma_y/2).
  K_phiphi, K_phiz follow by Gamma_z <-> Gamma_phi and phi -> -phi (XZcorr:152). The correlators do NOT depend
  on tau_z, tau_phi (hence not on eta): non-ideality only hits the zero-lag self-correlator K_ii(0) (XZcorr:155).

KEY SPECIAL CASES:
- ZERO-LAG CROSS-CORRELATOR = cosine of the angle between the axes:
    K_zz(+0) = 1,  K_zphi(0) = K_phiz(0) = cos(phi)   (XZcorr:160, eq:small_time_limit). Verified vs experiment,
    Fig.2b (crosses = data, dashed = cos phi).
- OmegaR_t = T1^-1 = T2^-1 = 0: phi=0 full correlation (K_zphi=K_zz=1); phi=pi full anticorrelation (-1);
    phi=pi/2 zero cross-correlation, K_zz = exp(-Gamma_phi tau), K_phiphi = exp(-Gamma_z tau) (XZcorr:166).
- OmegaR_t = 0: cross-correlator symmetric, K_zphi = K_phiz (XZcorr:166).
- Zeno pinning |phi|<<1: K_zphi(tau) ~ exp(-2 Gamma_jump tau),
    Gamma_jump = (phi^2 Gamma_z Gamma_phi + OmegaR_t^2)/[2(Gamma_z+Gamma_phi)] + (T1^-1+T2^-1)/4 (XZcorr:162-164).
- ANTISYMMETRIZED cross-correlator PROBES A WEAK COHERENT ROTATION (XZcorr:204, eq:Kzphi-antisym):
    K_zphi(tau) - K_phiz(tau) = [ 2 OmegaR_t sin(phi) / (Gamma_+ - Gamma_-) ] ( exp(-Gamma_- tau) - exp(-Gamma_+ tau) ).
  Fitting the phi=pi/2 experimental antisymmetrized correlator gives OmegaR_t/2pi ~ 12 kHz (XZcorr:207,212),
  a residual Rabi mismatch otherwise very hard to see under a 40 MHz drive.

NUMBERS (with anchors), from the experiment used for the comparison (XZcorr:189, 611-612):
- eta_z = 0.49, eta_phi = 0.41 (XZcorr:93). (z-channel = kappa 4.3 MHz mode; phi-channel = kappa 7.2 MHz mode.)
- T1 = 60 us, T2 = 30 us; Gamma_z^-1 = Gamma_phi^-1 = 1.31 us; Omega_R ~ Omega_rf = 2pi x 40 MHz (XZcorr:189,611).
- kappa_z/2pi = 4.3 MHz, kappa_phi/2pi = 7.2 MHz; kappa_z^-1 = 37 ns, kappa_phi^-1 = 22.1 ns; Omega_R^-1 = 4 ns;
  resonator modes omega_r,z/2pi = 7.4 GHz, omega_r,phi/2pi = 6.7 GHz (XZcorr:189,612,983).
- 11 angles phi_n = n pi/10, n=0..10, plus small correction delta_phi = (kappa_phi - kappa_z)/2 Omega_R ~ 0.036
  (~2 deg), so phi = phi_n + delta_phi (XZcorr:180,190,557).
- ~200,000 traces per angle, 5 us each, 4 ns sampling; correlators averaged over t1 in [1, 1.5] us (XZcorr:190-191).
- Detector responses Delta_I_z = 4.0, Delta_I_phi = 4.4 (arb units); JPA half-bandwidths 3.6 MHz (z) and 10 MHz
  (phi); analog filter cutoff ~25 MHz (~40 ns period ripple) (XZcorr:191,199). Markovian theory valid tau >~ 30 ns
  but agrees even below (XZcorr:194).
- Cooled to 30 mK (XZcorr:388). Weak self-correlator qubit tail Gamma/kappa_z = 0.028, Gamma/kappa_phi = 0.017 (XZcorr:1106).
grounds: two-axes Sections 6-7 (essay kzzCF, kzpCF). K_zz, K_zphi above ARE the closed forms to paste into a
  WL cell and check against the essay's correlators; zero-lag K_zphi(0)=cos(phi) is the felt payoff; the
  antisymmetrized combination is the OmegaR_t probe.

-------------------------------------------------------------------------------
### N5. 1612.02096  Atalaya-Bahrami-Pryadko-Korotkov (Bacon-Shor, THEORY)
-------------------------------------------------------------------------------
J. Atalaya, M. Bahrami, L. P. Pryadko, A. N. Korotkov, "Bacon-Shor code with continuous measurement of
non-commuting operators", arXiv:1612.02096 (Phys. Rev. A format; own journal line not present in fetched
source). (4qubit-paper_arxiv.tex; read in full). THEORY paper: analytics + Monte Carlo, NO lab experiment.
(Bahrami is the user; statements below are quoted from his own paper.)

SETUP. Four-qubit Bacon-Shor code, four two-qubit GAUGE operators measured simultaneously and continuously:
  X1X2 = X12 = G1,  X3X4 = X34 = G2,  Z1Z3 = Z13 = G3,  Z2Z4 = Z24 = G4   (4qubit:113-116, Eq.4-operators).
Out of the six pairs, four are non-commuting; conventional operation instead measures them projectively in
two steps (Z-pair then X-pair). Stabilizers X_all = X12 X34, Z_all = Z13 Z24 (4qubit:130-133). The 16-dim
space splits into code space Q_0 (X_all=+1, Z_all=+1) and three error subspaces Q_X, Q_Y, Q_Z (4qubit:138-140).

THE GAUGE-QUBIT REDUCTION (the explicit link to the single-qubit two-axis picture):
- Under continuous measurement in Q_0 the four-qubit state stays |psi(t)> = a(t)|z+> + b(t)|z->, i.e. the
  measurement drives a "GAUGE QUBIT" |a,b>_g while leaving the logical qubit |alpha,beta>_L untouched
  (4qubit:446, Eq.psi-a-b; 4qubit:574).
- In that 2-dim subspace, G3=Z13 and G4=Z24 are each a Z-MEASUREMENT of the gauge qubit, while G1=X12 and
  G2=X34 are each an X-MEASUREMENT of it. So the four gauge operators = simultaneous X and Z measurement of a
  single (gauge) qubit -- exactly the Ruskov-Korotkov-Molmer / Hacohen-Gourgy single-qubit two-axis problem
  (4qubit:453-455,576). Ideal case Gamma_m = 1/(2 tau_m): the gauge qubit diffuses uniformly on the GREAT
  CIRCLE of the Bloch sphere, a,b real (4qubit:455,627). Uniform-case gauge-qubit Ito SDE at 4qubit:616-626:
    dx_g = (1-x_g^2)(xi1+xi2)/sqrt(tau_m) - x_g z_g (xi3+xi4)/sqrt(tau_m) - 2 Gamma_m x_g   (and z_g symmetric;
    y_g decays at 4 Gamma_m). Output signals I1=I_X12=x_g+sqrt(tau1)xi1, ..., I4=I_Z24=z_g+sqrt(tau4)xi4 (4qubit:599-602).

ERROR SYNDROME FROM TIME-AVERAGED CROSS-CORRELATORS (the Bacon-Shor close):
- Projective parity X12 X34 = +1 and Z13 Z24 = +1 is replaced by POSITIVE cross-correlators of the noisy
  outputs, <I_X12 I_X34> = +1 and <I_Z13 I_Z24> = +1 in Q_0; a single-qubit error flips one or both signs
  (X_i -> Q_X flips the Z-correlator; Z_i -> Q_Z flips the X-correlator; Y_i flips both) (4qubit:457-458).
- SAME-TIME cross-correlator = +1, NOT the naively expected x_g^2: the output noise xi_1 feeds back on the
  state via quantum backaction, <sqrt(tau_1) xi_1(t) x_g(t+0)> = 1 - x_g^2, so the two terms sum to 1
  (4qubit:769, following Korotkov PRB 63, 085312).
- TWO-TIME cross-correlator (uniform case): <I_1(t1) I_2(t2)> = <I_3(t1) I_4(t2)> = exp(-2 Gamma_m |t1-t2|)
  (4qubit:771, Eq.corr-tau). In error subspaces the sign flips to -exp(...). Cross-correlators of ORTHOGONAL
  gauge components (I_1 with I_3, etc.) vanish in every subspace (4qubit:779).
- Monitored (smoothed) correlators C_12(t), C_34(t) built by double time-integration: inner exponential kernel
  time constant tau_c ~ Gamma_m^-1, outer kernel (rectangular T_c^r or exponential T_c^e) >= an order of
  magnitude longer than Gamma_m^-1 (4qubit:789-799). An error is flagged when C crosses (1-Theta)<C_tilde>.

THEORETICAL NUMBERS (with anchors):
- Quantum efficiency of each detector eta_k = 1/(2 Gamma_k tau_k) (4qubit:553).
- Optimal inner integration: eta=1 -> tau_c,opt = 0.342 tau_m, <C_tilde> = 0.745, A^2 = 2.13 tau_m;
  eta=0.5 -> tau_c,opt = 0.247 tau_m, <C_tilde> = 0.670, A^2 = 2.20 tau_m (4qubit:844).
- Response time T_R^r = (Theta/2) T_c^r, T_R^e = T_c^e ln[2/(2-Theta)] (4qubit:904). Optimal threshold for
  exponential integration Theta_opt = 1.43, but symmetric Theta = 1 is the recommended practical choice
  (4qubit:918,953). Exponential integration gives ~31% smaller logical error rate than rectangular at Theta=1
  for the same false-alarm rate (4qubit:916).
- Logical error rate scales with T_R. Depolarizing channel: projective gamma_L = (22/9) Gamma_d^2 Delta t
  (4qubit:354); continuous gamma_L = (28/9) Gamma_d^2 T_R (4qubit:712). Ratio gamma_L,cont/gamma_L,proj
  = (14 ln2 / 11)(T_c^e / Delta t) ~ 0.9 T_c^e / Delta t (4qubit:1053).
- For a target false-alarm rate 10^-5 tau_m^-1: T_R^e ~ 22.6 tau_m (eta=1) or 28.3 tau_m (eta=0.5) at symmetric
  threshold (4qubit:940).
- MAIN RESULT: continuous and projective operation are comparable when the collapse timescale tau_m is about an
  ORDER OF MAGNITUDE less than the projective inter-measurement period, Delta t ~ 20 tau_m (non-Gaussian
  correction ~ Delta t/20; abstract "an order of magnitude less") (4qubit:48,1064,1076,1100).
- Monte Carlo: full 4-qubit density matrix, 10^4-10^5 trajectories, time step delta t = 5e-3 Gamma_m^-1,
  ideal Gamma_m = 1/(2 tau_m), phase backaction neglected (4qubit:969,971).
grounds: the Bacon-Shor close of the essay. The gauge qubit under two Z- and two X- gauge measurements IS the
  single-qubit two-axis picture (great-circle diffusion), and the error syndrome is read exactly from the
  time-averaged cross-correlators whose two-time form is exp(-2 Gamma_m |t1-t2|) with same-time value +1.

===============================================================================
## READING LOG: all 5/5 Experimental-Anchor papers read END TO END. ##
##   1305.7270 (main+supp), 1402.1868 (main+supp), 1608.06652 (main+methods),
##   1702.08077 (main+supp), 1612.02096 (full). All fetched as TeX (no PDF fallback). ##
===============================================================================
