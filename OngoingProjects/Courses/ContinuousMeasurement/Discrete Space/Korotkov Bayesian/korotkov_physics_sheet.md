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
