# Caldeira-Leggett High-T Lindblad: Exploration Plan

Purpose: map and prioritize how to explore the new master equation of arXiv:2508.14262 (Pleasance-Aurell-Petruccione), both symbolically (closed form) and by simulation. This is a living document. Add to the status log at the bottom as we develop each item.

Companion documents in this folder:
- `caldeira-leggett-high-T-lindblad.md` (+ `.nb`): the finished essay, Layers 0-1 (kernel algebra, positivity table, harmonic covariance flow).
- `symbolic-exploration.md`: the working development of the symbolic track below.
- `related-papers/tex/2508.14262/`: the paper source.

## The master equation (reference)

    L rho = -(i/hbar)[H_S, rho]
            - i(gamma/hbar)[X, {P, rho}]          (friction)
            - D_PP [X,[X,rho]]                     (position-position diffusion)
            - D_XX [P,[P,rho]]                     (momentum-momentum diffusion, NEW)
            - 2 D_XP [X,[P,rho]]                   (cross diffusion, = 0 here)

with H_S = P^2/2M + V(X), and Dekker coefficients

    D_XX = gamma/(6 M kB T),   D_PP = 2 gamma M kB T/hbar^2,   D_XP = 0,   gamma = eta/M.

Lindblad operators L1 = X, L2 = P; Kossakowski matrix a = {{2 D_PP, -i gamma/hbar}, {i gamma/hbar, 2 D_XX}}; det a = gamma^2/(3 hbar^2) > 0 (guaranteed Lindblad). Holds for ANY V(X); the whole generator carries no cutoff Omega.

## The one fact that organizes everything

Every term is at most quadratic in the canonical pair (X, P) EXCEPT the Hamiltonian commutator [V(X), .]. So V(X) alone draws the line between what is symbolic and what is numeric:

- Quadratic V (free particle, harmonic, inverted oscillator, constant force): the dynamics closes on the first and second moments; Gaussian states stay Gaussian. FULLY SYMBOLIC.
- Anharmonic V (quartic, double-well, washboard/periodic, Morse): [V(X), rho] feeds new moments (<V'(X)>, <X V'(X)>, <V''(X)>) into the equations; the moment hierarchy never closes; states go non-Gaussian. NUMERIC.

Everything below follows from this split.

## Symbolic track (closed form)

Ranked by value. Tags: [done] essay already covers it; [next] recommended first; [dev] to develop.

1. [next] **Observable equations of motion / the moment hierarchy, for general V.** Tool: `/nca-route` (encode X, P in the Weyl algebra, [X,P] = i hbar; compute the adjoint generator L-dagger[O]; reduce to normal form; cross-check against a Fock-matrix carrier). Output: exact d<O>/dt for O in {X, P, X^2, P^2, {X,P}, H_S, ...} without assuming Gaussianity; the exact closure for quadratic V; the precise closure-breakers (<V'(X)>, <X V'(X)>) for anharmonic V; the energy balance and fluctuation-dissipation content. This maps the whole problem and is the highest-value symbolic step.

2. [next] **The Wigner phase-space PDE.** Tool: Wigner transform of the operator equation; `/pde-route` for solving it. Output: for quadratic V, an exact Klein-Kramers / Fokker-Planck equation with a momentum diffusion hbar^2 D_PP = 2 gamma M kB T (classical Einstein term) and, from the new (P,P) piece, a POSITION diffusion hbar^2 D_XX = hbar^2 gamma/(6 M kB T) that is O(hbar^2). That last term is a purely quantum position-space diffusion that smooths fine (sub-Planck) structure: "penalizing high-frequency fluctuations" made concrete. The single cleanest way to SEE the classicalization mechanism, and the bridge to simulation.

3. [done]/[dev] **Harmonic covariance, in full.** Essay has the steady covariance and the floor-piercing. To extend: closed-form time-dependent sigma(t) (matrix exponential of the drift); decoherence rate of a spatial superposition (~ D_PP dx^2) and of a momentum superposition (~ D_XX dp^2, the new channel); Gaussian purity hbar/(2 sqrt(det sigma)) and von Neumann entropy; the Robertson-Schrodinger boundary vs D_XX. Tool: `/pde-route` (DSolve / LyapunovSolve).

4. [dev] **Canonical Lindblad operators + CP parameter map.** Diagonalize the Kossakowski matrix symbolically to get the two canonical jump operators (complex combinations of X and P) and rates; map det a over parameter space and against Diosi / Breuer-Petruccione (their Omega dependence vs this generator's Omega freedom). Feeds the SSE.

5. [dev] **Classical limit hbar -> 0.** Reduce symbolically to the classical Kramers equation; confirm the new position diffusion is the O(hbar^2) correction that vanishes, i.e. this term is where the quantumness lives.

## Simulation track (numeric)

Simulation pays off ONLY for anharmonic V and non-Gaussian states (harmonic is closed-form). Do not spend simulation budget on harmonic V except the HPZ benchmark.

Representations, most useful first here:
- Position grid / DVR for rho(x, x'): X diagonal, P by FFT, the three double commutators as sparse actions. Best for arbitrary V.
- Wigner grid: solve the Fokker-Planck + Moyal PDE. Best for SEEING classicalization.
- Fock truncation + vectorized Liouvillian (MatrixExp for time-independent L, or QuantumEvolve). Fast, but D_PP pumps <P^2> so a fixed truncation saturates: watch this, use large N or a co-moving basis.
- SSE stochastic unraveling with the P-noise coupling. Cross-check and the paper's own object.

Physics to simulate:
1. [dev] **Double-well / washboard**: tunneling and coherence suppression; Kramers escape rate vs T and gamma; how D_XX shifts it.
2. [dev] **Classicalization demo (the title claim)**: decohere a position cat (via D_PP) and, decisively, a momentum cat (via D_XX), showing the new term specifically kills high-momentum interference.
3. [dev] **Wigner negativity**: generated by anharmonic V, erased by the diffusions; the quantum-to-classical picture in one movie.
4. [dev] **Full-rho positivity**: Min[Eigenvalues[rho]] >= 0 for this generator vs dipping negative for CL (D_XX = 0) on a non-Gaussian state; the Hilbert-space version of the essay's covariance result.
5. [dev] **HPZ benchmark**: for harmonic V, propagate the exact Hu-Paz-Zhang equation (time-dependent Lyapunov) against this Markovian generator. The direct test of the "correct Markovian limit" claim. Both are Gaussian, so semi-analytic and cheap; the most important CONCEPTUAL check in the program.

## Priority and recommended order (opinion)

1. Start symbolic with `/nca-route` on the moment hierarchy for general V. Decisive, general, cheap; tells you exactly where symbolic stops and simulation must start, with no Gaussian assumption. Verify against a Fock carrier.
2. Then derive the Wigner PDE symbolically. Turns the paper's central but abstract claim into a concrete O(hbar^2) position diffusion; the natural bridge into simulation. If only two things get done, these two.
3. Do not simulate harmonic V (fully closed-form) except the HPZ benchmark, which is the highest-value conceptual test and nearly analytic.
4. Build one anharmonic simulation (double-well or a momentum cat) on a grid, showing complete positivity holds AND the classicalization a Gaussian analysis cannot exhibit.
5. Lower priority: the SSE (a cross-check and trajectory picture, not a discovery tool) unless the momentum-noise coupling is the point.

## Caveats and open questions

- Transcendental V (cos, exp) is not a Weyl-algebra polynomial, so `/nca-route` gives the hierarchy STRUCTURE (the coupling to <V'(X)>) cleanly, but the closure itself depends on V; polynomial V is where the symbolic derivation is exact.
- Fock truncation vs D_PP heating: the momentum diffusion pumps <P^2>, so a fixed N saturates at long times; grid/DVR or a co-moving basis is safer.
- The classical limit of D_XX: confirm whether the phase-space position diffusion is genuinely O(hbar^2) (a quantum correction) as expected; this is the subtle bit and the essay already saw the O(hbar^2) correction in the steady covariance.
- Open: does the paper's "correct Markovian limit" claim survive the HPZ benchmark quantitatively, and in what regime does the Markovian generator track HPZ best/worst?
- Open: a cumulant / Gaussian-closure approximation for anharmonic V, informed by the exact closure-breakers from step 1, as a cheap semi-analytic route between symbolic and full simulation.

## Status log

- 2026-08-18: plan created from the audit; essay Layers 0-1 shipped. Symbolic track starting with the moment hierarchy (`symbolic-exploration.md`).
- 2026-08-18: explored and verified two high-impact symbolic features (`symbolic-exploration.md`, features A and B).
  - THERMALIZATION (top result): the steady momentum variance is sigma_PP = M kB T + D_XX hbar^2 M^2 w0^2/(2 gamma); the unique D_XX matching the quantum Gibbs value to O(hbar^2) is gamma/(6 M kB T), exactly the paper's coefficient. Minimal D_XX for positivity is gamma/(8 M kB T), so the thermal value exceeds it by 4/3: positivity is a corollary of correct thermalization, not an independent input. This-work matches Gibbs sigma_PP exactly and sigma_XX/sigma_XP to O(hbar^2 gamma); standard CL misses all O(hbar^2) quantum corrections; Breuer-Petruccione sits at the positivity boundary but undershoots the thermal correction (1/16 vs 1/12). Reframes the (P,P) term as an equilibrium fix, not just a spectrum fix.
  - DECOHERENCE/CLASSICALIZATION: position-superposition rate Gamma_x = D_PP dx^2 = 2 gamma M kB T dx^2/hbar^2 (grows with T); momentum-superposition rate Gamma_p = D_XX dp^2 = gamma dp^2/(6 M kB T) (shrinks with T). The new (P,P) channel is the title's "penalizing high-frequency fluctuations"; opposite T-signatures are a testable fingerprint CL cannot produce.
  - Still to develop: Wigner Fokker-Planck (feature C), the /nca-route moment hierarchy for general V (engine), relaxation spectrum and entropy production.
- 2026-08-18 (correction, after an inquiry-critic pass): `symbolic-exploration.md`'s "features" section was pushed toward research-grade, then corrected. Two attempted new results were WITHDRAWN as wrong and/or already carried out in the companion essay's Parts VI-VII:
  - A claimed "asymptotic hierarchy" was wrong. The fluctuation-path expansion CONVERGES, radius = first Matsubara frequency z = hbar w/kB T = 2 pi, because the delta^(2n) coefficient is the moment/(2n)!, not the moment. My claim conflated the two: the bare moments pi^(2n)|B_2n| do grow factorially, but the smooth-path 1/(2n)! rescues them. The correct convergent treatment is essay Parts VI-VII.
  - A claimed "thermalization ceiling" (a permanent O(hbar^4) Gibbs miss) was wrong: the delta'''' term cancels that miss exactly (essay Part VII). True but minor: the 3-coefficient constant family is affine in hbar^2; the higher-order generator is not confined to it.
  - The broken-detailed-balance reading of D_XX (sigma_XP = -M hbar^2 D_XX, nonzero entropy production) is Artini et al. (arXiv:2507.23322) and Bernad-Homa-Csirik (Eur. Phys. J. D 72, 212 (2018)), cited not claimed.
  - Root cause: novelty claimed on a stale read of the essay (it was extended after the first read) plus an unverified load-bearing premise (moments != coefficients). The section now honestly reports both directions as done/known and foregrounds the genuinely-open frontier: the mean-force Gibbs comparison (O(gamma) equilibrium, cutoff-dependent for Ohmic) and the HPZ steady-covariance benchmark (essay Part VIII). Kept as a correct illustration of the known tension: the exact-thermalizing generator is D_XX = 0 and not completely positive (det a = -gamma^2/hbar^2).
- 2026-08-18 (pivot, pursued and landed): the mean-force comparison is now WORKED in `symbolic-exploration.md` (section "Mean force, not bare Gibbs"), all closed form, 18/18 cells kernel-verified. Answer to the open question: the CP steady state's O(gamma) deviations from bare Gibbs are NOT mean-force corrections; each covariance entry is obstructed by a different mechanism.
  - sigma_XP: true equilibrium value 0 exactly (time reversal); CP forces |sigma_XP^ss| >= gamma hbar^2/(8 kB T) (the positivity floor read as an equilibrium error bound; Artini's constraint).
  - sigma_PP: true O(gamma) correction grows as (2 M hbar gamma/pi) ln(hbar Omega/2pi kB T) (slope exact); diverges in the paper's own scaling limit, so bare-Gibbs matching at O(gamma^0) is the strongest cutoff-free equilibrium statement (completes the paper's defense at its order).
  - sigma_XX: true O(gamma) correction is -zeta(3) gamma hbar^3/(2 pi^3 M (kB T)^2) at leading order, ODD in hbar; the family's steady covariance (and every finite kernel-tower truncation) is EVEN in hbar: parity obstruction.
  - Known vs new: MFG = exact late-time state (Subasi et al PRE 86, 061132), MFG covariances (Grabert-Schramm-Ingold), general Lindblad-vs-MFG mismatch (Lee-Yeo PRE 106, 054145; spin models), Redfield recovers MFG at 2nd order (Thingna-Wang-Hänggi JCP 136, 194110): all cited. New = the closed-form trichotomy for the cutoff-free QBM family. Remaining open: O(gamma^2) anatomy, finite part of the log, and the HPZ transient (essay Part VIII).
- 2026-08-18 (verification + quality loops run to their caps; final state): /wl-verify converged (round 1 OPEN ISSUES: 0 on the mean-force version; a second invocation on the revised document found 4 issues, fixed, and round 2 verified everything with the one remaining item being a stale .nb, rebuilt; round 3 verified all claims with the only refuted item being a miscount in the orchestrator's own brief, not the document). /wl-quality ran 3 rounds (15 -> 13 -> 12 findings, every physics result confirmed correct at each round; the findings were self-check gaps and idiom). All quality fixes through round 3 are applied: LyapunovSolve with residual check, general three-coefficient covariance displayed and specialized, kossakowski/diffusion/lyapInv as single objects, GKSL=Dekker+Lamb operator identity (symbolic hbar), full PSD spectrum and family-wide Robertson-Schrodinger Resolve proofs, 40-level Liouvillian numeric reference with legality certificates, both-branch three-route sxxMF check, momentum representation checked against the Drude FDT integral, critical-damping and free-particle limit cells, parity theorem cell, CP window solved by Reduce, five-equation ladder identity with symbolic hbar and two anharmonicities, Wigner PDE cross term restored with specialization label. 40 cells, all reproducing in one kernel pass; .nb rebuilt in parity. Because the quality loop capped without OPEN ISSUES: 0 and the post-cap fixes came after the last verifier stamp, no ledger pass is recorded and both gate deferral markers are set: the final bytes are self-verified (full kernel run) but carry no fresh agent stamp.
- 2026-08-18 (inquiry-critic round 3, final; deflated to its honest core): round 3 reassigned two more "additions" to the literature and the fixes are applied. (a) The renormalization-resistant det a <= -gamma^2/hbar^2 floor is contained in Isar's review (quant-ph/0602149, below Eq. 3.26: the thermal-stationarity constraint at mu=lam reads 0 >= lam^2, at every frequency and mass); now cited, the cells kept as packaging. (b) The slice/off-slice dichotomy IS the published theorem of Artini et al. (translation-covariant CPTP quadratic MEs + confining potential -> NESS; breaking TC by fine-tuning -> equilibrium); the mu=0 cell is now presented as their theorem in one closed form. (c) zeta(3) demoted: it is the next (first gamma-dependent) term of the standard GSI Matsubara expansion, not a discovery; the PARITY is what carries the argument. (d) New from the critique, verified here: the CP-equilibrating window is EXACTLY |mu| <= lam sech(hbar w0/2kBT) (knife-edge at high T, width ~ lam hbar^2 w0^2/8(kBT)^2); the free particle thermalizes on the slice (tension needs the potential); and the ledger's one falsifiable prediction, stated in-doc: the asymptotic HPZ generator obeys det a <= -gamma_inf^2/hbar^2 (exact QBM is asymptotically non-CP). (e) Next steps 3 re-posed: arXiv:2211.15722 is T=0-only, so it supports a T->0 consistency check, not the finite-T comparison; Gaussian relative entropy S(rho_ss||rho_MF) added as the operational error measure (open). (f) Practical sentence now aimed first at this document's own feature A (guarding against over-reading its thermalization motivation) and qualified to quadratic V; the anharmonic O(gamma^0) verdict is open (Engine + Artini entropy production). Surviving as this document's own: the parity scoping of the delta-tower, the mean-force D_XP in closed form, the mu* boundary, the HPZ prediction. 23/23 cells reproduce; loop capped at 3 rounds per protocol, residual issues folded in as above.
- 2026-08-18 (inquiry-critic round 2; rescoped): round 2 killed the "theorem" framing with two kernel-verified counterexamples, both reproduced here before applying. (a) The two-parameter friction family (Isar-Sandulescu-Scheid) at mu=0 (quantum-optical damped oscillator) is completely positive WITH an exact Gibbs steady state (Dekker-Valsakumar slack +hbar^2 lam^2 csch^2/4 > 0): "exact thermalization vs CP" is a property of the Caldeira-Leggett friction SLICE (mu=lam), not of Markovian generators. (b) Under a renormalized drift (M*,w*) every diagonal target above the vacuum floor is effective-thermal and det a returns to exactly -gamma^2/hbar^2: the quadratic "price" was the bare-Hamiltonian choice, absorbable by the mean-force renormalization (Timofeev-Trushechkin). SURVIVING sharp fact: on the CL friction slice, exact stationarity of ANY diagonal covariance costs det a <= -gamma^2/hbar^2 (Dekker-Valsakumar floor, PLA 104,67 (1984), now cited), renormalization or not. Also fixed: the momentum log is ordering-dependent (hbar Om >> kBT gives the log, odd in hbar; hbar Om << kBT gives gam M hbar^2 Om/6kBT, even, linear in cutoff, = HTL's high-T appendix form; HTL mis-credit for the log removed); the inversion formulas cited (ISS 1994; Toscano-Nicacio PRA 106,062207 with detailed balance noted as the stronger criterion); zeta(3) relabeled "believed new, search documented"; HPZ downgraded to structural signature with the quantitative check re-posed conventions-first (Next steps 3); the stated question demoted (answerable from feature A in two lines: the CP diagonal has NO O(gamma) deviation at all). Title now "where the equilibrium statement ends". 23/23 cells reproduce.
- 2026-08-18 (inquiry-critic round 1 on the mean-force section; restructured): the critic reproduced every number but found the framing misassembled. Fixes applied: (a) the three obstructions are now SYMPTOMS of one theorem, stated and verified: the unique constant-coefficient generator with steady covariance diag(a,b) has D_XX = 0 forced, D_PP = 2 gamma b/hbar^2, D_XP = (b - a M^2 w0^2)/(2 hbar^2 M), and det a = -gamma^2/hbar^2 - (b - a M^2 w0^2)^2/(hbar^4 M^2), strictly below Caldeira-Leggett (quadratic positivity price); bare Gibbs is isotropic so the cross coefficient vanishes at O(gamma^0). (b) A false claim was removed: "no FV term supplies D_XP" contradicted HPZ, which has no (P,P) term, carries a time-dependent cross coefficient (Halliwell-Yu quant-ph/9508004), and lands on the MFG (Subasi et al): the exact dynamics equilibrates by exactly the theorem's route. (c) Credits added: the momentum-log slope is Grabert-Weiss-Talkner Z. Phys. B 55, 87 (1984) and Hilt-Thomas-Lutz PRE 84, 031110 (2011); the position-squeezing physics is HTL 2011; only the high-T zeta(3) closed form and the theorem are claimed. (d) The essay contradiction fixed (O(hbar^2) -> O(gamma^0); the tower refines the O(gamma^0) covariance and says nothing about O(gamma)); "completes the paper's defense" dropped (the paper makes no equilibrium claim). (e) The floor conditioned on Einstein-pinned D_PP. wl-verifier: round 1 OPEN ISSUES: 0 on the pre-restructure version; all 21 cells reproduce after. New sharp open question (Next steps 3): asymptotic HPZ Gamma f vs the theorem's D_XP at O(gamma) via arXiv:2211.15722.
