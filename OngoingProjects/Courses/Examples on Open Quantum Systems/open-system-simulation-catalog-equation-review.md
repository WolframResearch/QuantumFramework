# Equation-Led Review of *Watching Quantum Things*

## Scope and method

This is a fresh review of `open-system-simulation-catalog.md`, organized around the governing master equation or stochastic differential equation (SDE) in each of its twenty cases. For every case I ask two questions:

1. What physical consequences does the catalog already extract, and are those demonstrations correct and well shaped?
2. What are the most important consequences of the equation itself, whether or not the catalog currently shows them?

The formulas below were independently derived from the stated generators and checked against the code conventions used in the catalog. Specialized interpretation was checked against primary literature where useful. No arXiv TeX archive was needed: none of the decisive results depended on an unresolved convention or a derivation unavailable from the equations themselves. Consequently, no source archive was downloaded.

Severity labels have the following meaning:

- **[FATAL]** The present calculation does not establish the physical claim it is meant to establish.
- **[MISLEADING]** The mathematics or code can run, but the stated physical interpretation is materially wrong.
- **[IMPRECISE]** The main result is sound, but its domain, convention, or numerical status needs qualification.
- **[MISSING PAYOFF]** The case is basically sound but omits a particularly important result that follows simply from its equation.

## Overall assessment

The catalog is strongest when it follows one equation through several levels of description: exact evolution, geometry, trajectories, ensemble recovery, and an observable. The cavity-loss case, photon counting, localization, the Kalman filter, thermalization, rapid purification, and feedback cooling mostly achieve that standard.

The main repairs are concentrated rather than diffuse:

1. **Case 4 does not actually simulate atom-mediated cat preparation or the two-atom Ramsey correlation.** It inserts a pure cavity cat by hand. The unmeasured atom instead leaves the reduced cavity in a mixture.
2. **The shared `regressionSpectrum` helper diagonalizes the transpose of the state-space Liouvillian.** It therefore uses left modes as though they were propagated density-matrix modes. This makes the Mollow and homodyne spectra quantitatively unreliable and is especially unsafe near exceptional points.
3. **The QPC spectrum merges before the Liouvillian exceptional point.** Spectral maxima, pole frequencies, and eigenvector coalescence are three different notions.
4. **Case 19 falsely says that restoring detuning makes an ideally monitored conditional state mixed.** Detuning is Hamiltonian and preserves purity.
5. Several otherwise good cases omit the cleanest analytic consequences of their equations: decay constants, steady occupations, fluxes, exact covariance fixed points, or measurement-rate formulas.

## Toolkit-level findings

### T1. Steady states in a degenerate kernel

The function `steadyState` reshapes every null vector and divides it by its trace. This is safe when the steady state is unique and the numerical null vector happens to be a physical density matrix. It is not a generally valid prescription when the kernel is multidimensional: arbitrary null-space basis vectors need not be Hermitian, positive, or even have nonzero trace.

For pure dephasing, the physical steady set is the convex set

$$
\rho_{\mathrm{ss}}=p\,|e\rangle\langle e|+(1-p)\,|g\rangle\langle g|,
\qquad 0\leq p\leq 1,
$$

not merely a list of normalized algebraic basis elements. The review should distinguish the linear kernel of $\mathcal L$ from its intersection with the density-operator cone.

**Assessment:** **[IMPRECISE]** The present examples usually land on sensible matrices, but the helper is not a general physical steady-state solver.

### T2. Positivity is not accuracy

The Rouchon-style update is valuable because each finite step remains positive after normalization. But finite-step positivity does not mean the result is accurate for arbitrary $dt$. The statement that averaging trajectories recovers the Lindblad equation is an order-in-$dt$ statement, becoming exact in the continuum limit. Large steps can remain physical while having substantial weak and strong numerical error.

**Recommended addition:** accompany positivity checks with step-doubling and ensemble-convergence checks. Phrase the result as “recovers the Lindblad evolution to the scheme's integration order and in the $dt\to0$ limit.”

### T3. “Unwatched” and “record discarded” are not identical physical situations

The unconditional master equation applies both when no detector exists and when a detector exists but its record is discarded. “The system leaks but is not watched” is therefore too narrow. The physically decisive operation is conditioning on versus averaging over the record.

### T4. The regression-spectrum helper is not generally correct

The catalog advances density matrices with

$$
\operatorname{vec}\rho(t)=e^{\mathcal L t}\operatorname{vec}\rho(0).
$$

The modes used to expand a regression seed must therefore be right eigenvectors of $\mathcal L$. The helper instead evaluates `Eigensystem[Transpose[big]]`, obtaining right eigenvectors of $\mathcal L^T$, i.e. left modes of $\mathcal L$, and reshapes them as if they were propagated states.

The desired expression is

$$
S_{AB}(\omega)=2\operatorname{Re}\!\left[
\operatorname{Tr}\!\left(A(i\omega-\mathcal L)^{-1}B\right)
\right],
$$

with the stationary mode projected out when $\omega=0$. A minimal modal repair is to diagonalize `big`, not `Transpose[big]`. A better repair is a trace-constrained resolvent, because an eigenvector expansion becomes ill-conditioned at an exceptional point. Also, a pole with a complex residue is not by itself a positive Lorentzian; only the complete physical sum is guaranteed to have the appropriate reality and positivity properties.

Direct resolvent checks show that the transpose error is not cosmetic. For the weak-drive squeezed homodyne example, the correct minimum is about $0.737$ in the catalog's normalization, whereas the current helper gives about $0.854$. Symmetry happens to mask much of the error in the special QPC seed, but it does not rescue the general routine.

**Assessment:** **[FATAL]** for quantitative claims in Cases 20 and the generic helper; Case 11 additionally has an independent analytic interpretation error discussed below.

### T5. Conditions for pure conditional trajectories

A pure initial state stays pure when every decohering channel is represented in an efficient monitored unravelling. A mixed initial state may purify. A Hamiltonian term—including detuning—does not by itself make a conditional trajectory mixed. The global toolkit states most of these conditions correctly; later prose should preserve them exactly.

---

## Case-by-case review

### 1. Pure dephasing

The equation is

$$
\dot\rho=\gamma\,\mathcal D[\sigma_z]\rho.
$$

**Already explored.** The catalog derives fixed populations, exponential decay of the off-diagonal element, contraction of the Bloch disk toward the $z$ axis, and the degenerate steady-state manifold. This is in good shape and is one of the clearest symbolic cases. The only caution is the toolkit's algebraic treatment of a degenerate null space.

**Most important results and consequences.** With $\rho_{eg}(0)$ arbitrary,

$$
\rho_{eg}(t)=e^{-2\gamma t}\rho_{eg}(0),
\qquad
\rho_{ee}(t)=\rho_{ee}(0),
\qquad
\rho_{gg}(t)=\rho_{gg}(0).
$$

Thus

$$
x(t)=e^{-2\gamma t}x_0,
\qquad
y(t)=e^{-2\gamma t}y_0,
\qquad
z(t)=z_0,
$$

and the dephasing time in this convention is

$$
T_2=\frac{1}{2\gamma},
\qquad T_1=\infty.
$$

The purity is

$$
P(t)=\frac12\left[1+z_0^2+e^{-4\gamma t}(x_0^2+y_0^2)\right].
$$

**What to add.** **[MISSING PAYOFF]** State $T_2$, the purity law, and the absence of energy relaxation explicitly. Say that the environment acquires which-$z$ information; calling this “continuous measurement” should not imply that a record is actually collected.

### 2. Driven, damped two-level atom

The equation is

$$
\dot\rho=-i\left[\frac{\Omega}{2}\sigma_x,\rho\right]
+\gamma\mathcal D[\sigma_-]\rho.
$$

**Already explored.** The optical Bloch equations, fixed point, Liouvillian modes, Bloch-space spiral, invariant line, and a physically realizable ensemble are developed in unusual depth. The Bloch equations and steady state are correct. The physically realizable two-state ensemble is plausible, but the catalog should show the monitored unravelling and transition rates that realize it rather than infer realizability from an invariant chord alone.

**Most important results and consequences.** The Bloch equations are

$$
\dot x=-\frac\gamma2x,
\qquad
\dot y=-\frac\gamma2y-\Omega z,
\qquad
\dot z=\Omega y-\gamma(z+1).
$$

The steady state is

$$
x_{\mathrm{ss}}=0,
\qquad
y_{\mathrm{ss}}=\frac{2\Omega\gamma}{\gamma^2+2\Omega^2},
\qquad
z_{\mathrm{ss}}=-\frac{\gamma^2}{\gamma^2+2\Omega^2}.
$$

Therefore

$$
\rho_{ee}^{\mathrm{ss}}=\frac{\Omega^2}{\gamma^2+2\Omega^2},
\qquad
R_{\mathrm{fl}}=\gamma\rho_{ee}^{\mathrm{ss}}.
$$

This displays the central physical crossover: $\rho_{ee}^{\mathrm{ss}}\sim\Omega^2/\gamma^2$ for weak drive and saturates at $1/2$ for strong drive. The nonzero modes are

$$
-\frac\gamma2,
\qquad
-\frac{3\gamma}{4}\pm i\sqrt{\Omega^2-\frac{\gamma^2}{16}},
$$

so damped Rabi oscillations begin when $\Omega>\gamma/4$.

**What to add.** **[MISSING PAYOFF]** Put saturation, fluorescence flux, weak/strong limits, and the oscillation threshold before the more specialized invariant-line discussion.

### 3. Cavity cat under photon loss

The equation is

$$
\dot\rho=\gamma\mathcal D[a]\rho.
$$

**Already explored.** The Kraus solution, Wigner-fringe decay, exact numerical extraction of the coherence envelope, separation-squared scaling, and Fock-cutoff convergence are all strong. This is close to a model case.

**Most important results and consequences.** Let

$$
\eta=e^{-\gamma t}.
$$

Then a coherent state remains coherent,

$$
|\alpha\rangle\longmapsto|\sqrt\eta\,\alpha\rangle,
\qquad
\langle n(t)\rangle=\eta\langle n(0)\rangle.
$$

For a coherence operator,

$$
|\alpha\rangle\langle\beta|
\longmapsto
\exp\!\left[-\frac{1-\eta}{2}
\left(|\alpha|^2+|\beta|^2-2\alpha\beta^*\right)\right]
|\sqrt\eta\alpha\rangle\langle\sqrt\eta\beta|.
$$

Its magnitude is reduced by

$$
|C_{\alpha\beta}(t)|
=\exp\!\left[-\frac{1-e^{-\gamma t}}{2}|\alpha-\beta|^2\right].
$$

At short times,

$$
\Gamma_{\mathrm{dec}}=\frac\gamma2|\alpha-\beta|^2,
\qquad
t_{\mathrm{dec}}\simeq\frac{2}{\gamma|\alpha-\beta|^2}.
$$

For $|\alpha\rangle\pm|-\alpha\rangle$, $\Gamma_{\mathrm{dec}}=2\gamma|\alpha|^2$.

**What to add.** Explicitly compare the energy-relaxation time $1/\gamma$ with the much shorter cat-decoherence time. Also correct any suggestion that the asymptotic state is a persistent two-component mixture: both branches eventually reach the same vacuum, so the final state is pure vacuum. The mixture is an intermediate-time description.

### 4. Atom-made cat and Ramsey readout

The dispersive interaction is

$$
H_C=\chi a^\dagger a\,\sigma_z,
\qquad \phi=\chi\tau.
$$

**Already explored.** The catalog compares Wigner-fringe decay for two phase separations and uses trace distance from a shrinking incoherent mixture as a coherence proxy. Those calculations are meaningful for a cavity cat that has already been prepared.

**Most important results and consequences.** Starting with an atomic superposition gives the joint state

$$
|\Psi\rangle
=\frac{1}{\sqrt2}\left(
|e\rangle|\alpha e^{-i\phi}\rangle
+|g\rangle|\alpha e^{i\phi}\rangle
\right).
$$

If the atom is ignored, the cavity state is

$$
\rho_{mathrm{cav}}
=\frac12\left(
|\alpha e^{-i\phi}\rangle\langle\alpha e^{-i\phi}|
+|\alpha e^{i\phi}\rangle\langle\alpha e^{i\phi}|
\right),
$$

not a pure cat. A cavity cat is obtained only after a Ramsey rotation and atomic measurement or postselection:

$$
|C_\pm\rangle\propto
|\alpha e^{-i\phi}\rangle
\pm e^{i\vartheta}|\alpha e^{i\phi}\rangle.
$$

The branch separation is

$$
|\delta\alpha|^2=4|\alpha|^2\sin^2\phi,
$$

so photon loss gives the coherence magnitude

$$
|C(t)|
=\exp\!\left[-2|\alpha|^2\sin^2\phi\left(1-e^{-\gamma t}\right)\right]
$$

and the short-time rate

$$
\Gamma_{\mathrm{dec}}=2\gamma|\alpha|^2\sin^2\phi.
$$

**Assessment:** **[FATAL]** The current code defines `catRho` directly and never evolves the atom-cavity joint state, applies Ramsey pulses, samples atomic outcomes, or computes the stated two-atom correlation $\eta$. Thus it demonstrates decay of an inserted cat, not “making a cat with an atom” or reading it with two atoms.

**Required repair.** Build the joint state, implement the Ramsey analysis pulse and first-atom postselection, let the cavity decay, send the second probe, and compute the conditional probabilities entering $\eta$. Keep trace distance only as a secondary theoretical proxy, not as a replacement for the claimed observable.

### 5. Dispersive cavity readout

The equation is

$$
\dot\rho=-i[\epsilon(a+a^\dagger)+\chi a^\dagger a\,\sigma_z,\rho]
+\gamma\mathcal D[a]\rho.
$$

**Already explored.** The two conditional pointer trajectories are derived correctly, the joint qubit-cavity master equation shows coherence loss, the Fock cutoff is checked, and conditional Wigner functions expose the connection between pointer distinguishability and interference. This is physically rich and mostly good.

**Most important results and consequences.** Since $[\sigma_z,H]=0$, the ideal dispersive model is quantum nondemolition with respect to $\sigma_z$: it dephases but does not flip the qubit. For $s=\pm1$,

$$
\alpha_s(t)=\frac{-i\epsilon}{\gamma/2+is\chi}
\left(1-e^{-(\gamma/2+is\chi)t}\right).
$$

The steady pointer separation is

$$
|\delta\alpha_{\mathrm{ss}}|
=\frac{2|\epsilon\chi|}{(\gamma/2)^2+\chi^2}.
$$

For an ideal monitored output, the asymptotic measurement-induced dephasing rate is

$$
\Gamma_\phi^{\mathrm{meas}}
=\frac\gamma2|\delta\alpha_{\mathrm{ss}}|^2.
$$

During ring-up, the reduced qubit coherence contains both the overlap of the two intracavity coherent states and the overlap of the emitted fields. The accumulated output information scales as

$$
\mathrm{SNR}^2(t)\propto
\gamma\int_0^t|\alpha_+(s)-\alpha_-(s)|^2\,ds,
$$

with convention-dependent factors set by homodyne versus heterodyne detection and efficiency.

**Issue.** The IQ clouds in the catalog sample a one-time Gaussian distribution centered at the instantaneous intracavity amplitude. The formula

$$
P_{\mathrm{err}}=\frac12\operatorname{erfc}
\left(\frac{|\delta\alpha(t)|}{2}\right)
$$

is correct for that terminal coherent-mode discrimination in the chosen normalization. It is **not** the error curve for an integrated output record. Once the cavity reaches steady state, the plotted error saturates, whereas a continued ideal output integration keeps accumulating SNR.

**Assessment:** **[MISLEADING]** Relabel the present curve as terminal intracavity-mode discrimination, or simulate

$$
dJ(t)=\sqrt{\eta\gamma}\,q_s(t)\,dt+dW(t)
$$

and classify the integrated, optimally filtered output. The latter is the actual circuit-QED readout problem. See [Gambetta et al., *Measurement-induced dephasing and number splitting in circuit QED*](https://arxiv.org/abs/cond-mat/0602322).

### 6. High-temperature quantum Brownian motion

The Caldeira-Leggett equation used is

$$
\dot\rho=-i[H,\rho]-i\gamma[x,\{p,\rho\}]
-2\gamma Mk_BT[x,[x,\rho]].
$$

**Already explored.** The catalog correctly separates drift from diffusion, evolves the covariance, identifies thermal equilibration, warns that the high-temperature equation is not generally completely positive, compares it with a Lindblad thermal oscillator, and shows separation-dependent cat decoherence. This is a strong section.

**Most important results and consequences.** In this convention the momentum damping rate is $2\gamma$:

$$
\dot{\langle x\rangle}=\frac{\langle p\rangle}{M},
\qquad
\dot{\langle p\rangle}=-M\Omega^2\langle x\rangle-2\gamma\langle p\rangle.
$$

The underdamped poles are

$$
\lambda_\pm=-\gamma\pm i\sqrt{\Omega^2-\gamma^2}.
$$

The stationary covariance is the classical thermal covariance

$$
\Sigma_{xx}^{\mathrm{ss}}=\frac{k_BT}{M\Omega^2},
\qquad
\Sigma_{pp}^{\mathrm{ss}}=Mk_BT,
\qquad
\Sigma_{xp}^{\mathrm{ss}}=0.
$$

In the position representation,

$$
\partial_t\rho(x,x')\supset
-2\gamma Mk_BT(x-x')^2\rho(x,x'),
$$

so spatial coherences decay at

$$
\Gamma_{\mathrm{dec}}(\Delta x)=2\gamma Mk_BT(\Delta x)^2.
$$

For the completely positive thermal Lindblad model, a coherent-state separation initially decoheres at

$$
\Gamma_{\mathrm{dec}}simeq
\frac\gamma2(2n_T+1)|\delta\alpha|^2.
$$

**What to improve.** State explicitly that the ground covariance is $\mathbb 1/2$ only in the adopted dimensionless units $M\Omega=1$. Keep the Caldeira-Leggett and thermal-Lindblad coefficients visibly separate; they are related physical limits, not the same generator.

### 7. Photon counting

For the leak $c$, the jump unravelling has click probability

$$
\Pr(dN=1)=\lambda(t)dt,
\qquad
\lambda(t)=\langle c^\dagger c\rangle_c,
$$

and no-click evolution generated by

$$
H_{\mathrm{eff}}=H-\frac{i}{2}c^\dagger c.
$$

**Already explored.** Click records, state jumps, non-Hermitian no-click evolution, waiting-time statistics, antibunching, and ensemble recovery of the master equation are all covered. The reset-and-rebuild physical narrative is very good.

**Most important results and consequences.** For the resonantly driven atom,

$$
R_{\mathrm{ss}}
=\gamma\rho_{ee}^{\mathrm{ss}}
=\frac{\gamma\Omega^2}{\gamma^2+2\Omega^2}.
$$

After an emission resets the atom to $|g\rangle$, the survival probability and waiting-time density are

$$
S(\tau)=\left\|e^{-iH_{\mathrm{eff}}\tau}|g\rangle\right\|^2,
\qquad
w(\tau)=-\frac{dS}{d\tau}
=\left\|c e^{-iH_{\mathrm{eff}}\tau}|g\rangle\right\|^2.
$$

Thus $w(0)=0$ and $g^{(2)}(0)=0$: a second photon cannot be emitted immediately after the first. Drive-induced rebuilding produces damped Rabi structure in the waiting-time law and in $g^{(2)}(\tau)$.

**What to add.** Show $\int_0^\infty w(\tau)d\tau=1$, the steady flux, and $g^{(2)}(0)=0$ explicitly. The Bernoulli step requires $\lambda dt\ll1$ and is first-order; add a step-size check.

### 8. Homodyne monitoring

The diffusive SME is

$$
d\rho_c=\mathcal L\rho_c,dt
+\mathcal H[c]\rho_c,dW,
\qquad
dJ=\langle c+c^\dagger\rangle_cdt+dW.
$$

**Already explored.** The observed/inferred ledger, singular current, integrated record, state estimate, innovation, autocorrelation, and ensemble recovery form one of the catalog's best pedagogical sequences.

**Most important results and consequences.** The demonstrated atom starts in a pure state and all of its emission is monitored ideally. Its conditioned state therefore does not “purify smoothly”; it **remains pure**. A mixed initial state would purify. For a local-oscillator phase $\phi$,

$$
dJ_\phi=
\sqrt\eta\left\langle e^{-i\phi}c+e^{i\phi}c^\dagger\right\rangle_cdt+dW,
$$

so the chosen phase determines which output-field quadrature is directly observed and how the measurement backaction correlates with the atomic motion. Inefficiency inserts an unobserved fraction of the channel and permits conditional mixing.

Innovation whiteness is a necessary consistency check, but it is not by itself proof that a filter is correct or optimal: a wrong model can pass limited autocorrelation tests, and simulated data generated by the same filter are not an independent validation.

**Assessment:** **[IMPRECISE]** Replace “purifies smoothly” with “diffuses continuously while remaining pure in the ideal example; a mixed initial state would purify.” Add efficiency and local-oscillator phase as the two principal experimental controls.

### 9. Heterodyne monitoring

The two real monitored channels are

$$
c_I=\frac{c}{\sqrt2},
\qquad
c_Q=\frac{ic}{\sqrt2},
$$

and indeed

$$
\mathcal D[c_I]+\mathcal D[c_Q]=\mathcal D[c].
$$

**Already explored.** The two records, two innovations, separate whiteness checks, cross-correlation check, pure conditioned histories, and common unconditioned master equation are handled well.

**Most important results and consequences.** Heterodyne detection divides the output between two conjugate field quadratures. Each carries half the signal power of an optimally phased homodyne channel—the familiar $3\,$dB simultaneous-quadrature penalty. With the catalog's convention,

$$
dZ=dW_I+i,dW_Q,
\qquad
dZ^2=0,
\qquad
dZ\,dZ^*=2dt.
$$

The alternative convention $dZ=(dW_I+i,dW_Q)/\sqrt2$ gives $dZ,dZ^*=dt$; either is valid if kept consistent.

**Issue.** Homodyne and heterodyne monitor field quadratures, not atomic $\sigma_x$ and $\sigma_y$ as two simultaneously projected observables. A homodyne filter still estimates the full density matrix from its one record; heterodyne supplies two noisier field-quadrature records. Tightening this language would prevent an incorrect uncertainty-principle interpretation.

### 10. Continuous-measurement quantum Zeno crossover

The equation is

$$
\dot\rho=-i\left[\frac\Omega2\sigma_x,\rho\right]
+k\mathcal D[\sigma_z]\rho.
$$

**Already explored.** Weak- and strong-measurement trajectories, the deterministic Liouvillian eigenvalues, the exceptional point, the turnover of the slow rate, and the $1/k$ Zeno law are correctly connected. This section is in good shape.

**Most important results and consequences.** The ensemble variables obey

$$
\ddot z+2k\dot z+\Omega^2z=0,
$$

with poles

$$
\lambda_\pm=-k\pm\sqrt{k^2-\Omega^2}.
$$

Thus the exceptional point is $k=\Omega$. Deep in the Zeno regime,

$$
\Gamma_Z=k-\sqrt{k^2-\Omega^2}
\simeq\frac{\Omega^2}{2k}.
$$

This is the decay rate of the **ensemble mean**. A symmetric conditional telegraph process with per-direction switching rate $r$ has mean decay $2r$, hence

$$
r\simeq\frac{\Omega^2}{4k}.
$$

**What to add.** Distinguish the ensemble slow rate from the conditional dwell-time switching rate. “Measure more often” should be phrased as “increase the continuous measurement strength” in this model. The trajectory crossover can also be described in stages rather than as an instantaneous onset of telegraph behavior; see [Snizhko, Kumar, and Romito](https://arxiv.org/abs/2003.10476).

### 11. Charge qubit monitored by a QPC

The equation is the same Zeno generator with

$$
\dot\rho=-i\left[\frac\Omega2\sigma_x,\rho\right]
+\frac\kappa4\mathcal D[\sigma_z]\rho,
\qquad k=\frac\kappa4.
$$

**Already explored.** The weak Rabi and strong telegraph regimes, the quantum-regression construction, normalization of spectral weight, relation to a detector-current noise floor, and the distinction between measurement strength and detector transparency are valuable.

**Most important results and consequences.** For the maximally mixed steady state,

$$
C_{zz}(0)=1,
\qquad
\dot C_{zz}(0)=0,
\qquad
\ddot C_{zz}+2k\dot C_{zz}+\Omega^2C_{zz}=0.
$$

The two-sided symmetrized spectrum is therefore

$$
S_{zz}(\omega)
=\frac{4k\Omega^2}
{(\Omega^2-\omega^2)^2+4k^2\omega^2},
\qquad
\int_{-\infty}^{\infty}\frac{d\omega}{2\pi}S_{zz}(\omega)=1.
$$

Its nonzero maxima occur at

$$
\omega_{\mathrm{peak}}=\pm\sqrt{\Omega^2-2k^2},
\qquad k<\frac{\Omega}{\sqrt2}.
$$

They merge into a single central maximum at

$$
k=\frac{\Omega}{\sqrt2},
$$

whereas the Liouvillian poles cease oscillating and coalesce only at

$$
k=\Omega.
$$

**Assessment:** **[FATAL]** The catalog identifies the spectral maxima with the pole imaginary parts $\pm\sqrt{\Omega^2-k^2}$ and says the peaks merge at the exceptional point. Both claims are false. Damping shifts maxima away from pole frequencies, and the spectral merger precedes eigenvector coalescence. At the plotted ratio $k/\Omega=0.75$, the exact spectrum already has one central maximum even though the Liouvillian is still underdamped.

Deep in the Zeno regime the central width approaches $\Omega^2/(2k)$, so the catalog's asymptotic narrowing claim is correct. The observed current spectrum has the form

$$
S_I(\omega)=S_0+(\Delta I)^2\times
\text{(normalization)}\times S_{zz}(\omega).
$$

For an ideal weak continuous detector, the Korotkov-Averin signal-to-noise ratio is bounded by four; see [Korotkov and Averin](https://arxiv.org/abs/cond-mat/0002203).

**Required repair.** Replace the generic modal calculation here with the exact $S_{zz}$ above, label pole locations and spectral maxima separately, and show the spectral-merger point and exceptional point as two distinct markers.

### 12. Measurement-induced localization

With no Hamiltonian and ideal $\sigma_z$ monitoring,

$$
d\rho=k\mathcal D[\sigma_z]\rho\,dt
+\sqrt{k}\,\mathcal H[\sigma_z]\rho\,dW.
$$

**Already explored.** The trajectory fan, Born-rule frequencies, long-time record discrimination, and distinction between a flat record integral and the Bayesian state estimate are all well chosen. This case is strong.

**Most important results and consequences.** For a pure state,

$$
dz=2\sqrt{k}(1-z^2)dW.
$$

There is no drift: $z(t)$ is a bounded martingale. Its limiting probabilities therefore obey

$$
\Pr(z_\infty=+1)=\frac{1+z_0}{2},
\qquad
\Pr(z_\infty=-1)=\frac{1-z_0}{2},
$$

which is the Born rule. The log-likelihood coordinate $q=\operatorname{artanh}z$ obeys

$$
dq=4k\tanh q\,dt+2\sqrt{k}\,dW,
$$

making the runaway toward one of the two hypotheses explicit. Meanwhile the ensemble average obeys pure dephasing even though each ideal conditional trajectory remains pure.

**What to improve.** The poles $z=\pm1$ are reached asymptotically in the continuum SDE, not generally at a finite deterministic collapse time. Say “localizes arbitrarily close” rather than “is absorbed” unless an explicit threshold is defined.

### 13. Quantum Kalman filter

The conditional covariance obeys

$$
\dot\Sigma=A\Sigma+\Sigma A^T+D-\Sigma C^TC\Sigma,
$$

with the matrices given in the catalog.

**Already explored.** Deterministic covariance evolution, conditional squeezing, purity through $\det\Sigma=1/4$, ellipse geometry, noisy mean tracking, comparison with a hidden reference trajectory, and innovation whiteness are all covered well.

**Most important results and consequences.** Writing

$$
s=\sqrt{1+4k^2},
$$

the positive stabilizing Riccati fixed point is

$$
\Sigma_{xx}=\frac{1}{\sqrt{2(1+s)}},
\qquad
\Sigma_{xp}=\frac{s-1}{4k},
\qquad
\Sigma_{pp}=s\,\Sigma_{xx}.
$$

It satisfies

$$
\det\Sigma=\frac14.
$$

Thus position is conditionally squeezed below $1/2$, but momentum broadens and the ellipse tilts so the Heisenberg bound is exactly respected. The covariance is deterministic because the Gaussian filtering equations close; the conditional mean remains stochastic. Without conditioning, measurement backaction heats the oscillator rather than producing a stationary squeezed state.

**What to add.** Put the analytic fixed point next to the numerical `Solve` result. Replace “whiteness is the working definition of an optimal filter” with the weaker and more accurate statement that innovation whiteness is an important necessary diagnostic of model/filter consistency.

### 14. Markovian measurement feedback

For unit efficiency the feedback master equation used is

$$
\dot\rho=G\mathcal D[\sigma_z-i\sigma_y]\rho.
$$

**Already explored.** The discrete record-to-rotation loop, dark-state check, ensemble comparison with the feedback master equation, causal record/control/state view, and efficiency sweep are all useful and basically correct.

**Most important results and consequences.** Since

$$
\sigma_z-i\sigma_y=2|+x\rangle\langle-x|,
$$

the channel optically pumps the qubit toward $|+x\rangle$. Starting from the maximally mixed state,

$$
x(t)=1-e^{-4Gt},
\qquad y(t)=z(t)=0.
$$

With the catalog's inefficiency term,

$$
\dot\rho=G\mathcal D[\sigma_z-i\sigma_y]\rho
+G\frac{1-\eta}{\eta}\mathcal D[\sigma_y]\rho,
$$

the stationary polarization and purity are

$$
x_{\mathrm{ss}}=\frac{2\eta}{1+\eta},
\qquad
P_{\mathrm{ss}}=\frac12\left[1+left(\frac{2\eta}{1+\eta}\right)^2\right].
$$

**What to add.** These closed forms make the rate and efficiency ceiling immediately legible. State that the discrete “measure then rotate” construction approaches the Markovian feedback equation as $dt\to0$; delay, bandwidth, and finite-step ordering are physical assumptions, not mere implementation details.

### 15. Thermalization of an oscillator

The equation is

$$
\dot\rho=-i[\omega a^\dagger a,\rho]
+\gamma(n_T+1)\mathcal D[a]\rho
+\gamma n_T\mathcal D[a^\dagger]\rho.
$$

**Already explored.** Attraction from hot and cold initial states, the whole geometric Fock distribution, detailed balance, and the logarithmic population slope are all handled very well.

**Most important results and consequences.** The occupation satisfies the exact scalar law

$$
\frac{d}{dt}\langle n\rangle=-\gamma(\langle n\rangle-n_T),
$$

so

$$
\langle n(t)\rangle=n_T+igl(\langle n(0)\rangle-n_T\bigr)e^{-\gamma t}.
$$

The first field moment obeys

$$
\frac{d}{dt}\langle a\rangle
=-\left(i\omega+\frac\gamma2\right)\langle a\rangle.
$$

The unique stationary distribution is

$$
p_n=(1-r)r^n,
\qquad
r=\frac{n_T}{n_T+1}=e^{-\beta\omega}.
$$

**What to improve.** **[MISSING PAYOFF]** Add the exact transient occupation and coherence law. The displayed equation includes $\omega a^\dagger a$, but the numerical evolution uses a zero Hamiltonian; this is harmless for the diagonal initial states and populations shown, but it does not demonstrate phase rotation or coherence decay. In a finite Fock cutoff, compare with the **normalized truncated** geometric distribution and include an explicit cutoff-convergence check.

### 16. Rapid purification by adaptive measurement

The mean impurity equation is

$$
\mathbb E[dL]
=-4kL\left(\sin^2\theta+2L\cos^2\theta\right)dt,
\qquad
L=\frac12(1-|\mathbf a|^2).
$$

**Already explored.** Fixed and crosswise strategies, deterministic exponential purification, fixed-time distributions, first-passage distributions, and the reversal of which strategy is “faster” under different cost functions are all covered exceptionally well.

**Most important results and consequences.** Crosswise measurement, $\theta=\pi/2$, removes the stochastic term in $dL$ and gives

$$
L(t)=L(0)e^{-4kt}.
$$

At late times this doubles the average fixed-time purification exponent relative to an uncontrolled fixed-axis measurement. But a fixed-axis protocol has rare favorable records; its mean time to hit a prescribed high-purity threshold is asymptotically half that of the deterministic protocol. The catalog correctly treats these as different optimization objectives. These interpretations agree with [Jacobs](https://arxiv.org/abs/quant-ph/0301056) and [Wiseman and Ralph](https://arxiv.org/abs/quant-ph/0603062).

**What to improve.** Explain that the implementation adaptively rotates the **measurement axis** rather than applying a unitary feedback Hamiltonian. Physical realization assumes fast, accurate state estimation and control bandwidth. A step-refinement check should quantify how closely finite-$dt$ paths follow the exact deterministic law.

### 17. Measurement-based feedback cooling

The catalog uses a Gaussian conditional filter for the oscillator together with thermal damping and the feedback force

$$
f(t)=-G\langle p\rangle_c.
$$

**Already explored.** Conditional covariance, measurement heating, decomposition of unconditional energy into conditional covariance plus spread of conditional means, closed-loop trajectories, a gain optimum, and efficiency dependence are all developed carefully. This is one of the most complete cases.

**Most important results and consequences.** With the loop open, the unconditional occupation obeys

$$
\dot n=-\gamma(n-n_T)+\frac{k}{2},
$$

so measurement alone produces

$$
n_{\mathrm{open}}^{\mathrm{ss}}=n_T+\frac{k}{2\gamma}.
$$

The $k/2$ term is measurement-backaction heating. Feedback changes the drift of the conditional mean but not the conditional covariance—the quantum version of the separation principle. The cooled occupation has the structure

$$
n_{\mathrm{cool}}+rac12
=\frac12\operatorname{Tr}\left(V_c+V_{\bar x}\right),
$$

where $V_c$ is the Riccati fixed point and $V_{\bar x}$ solves the closed-loop Lyapunov equation. Too little gain leaves thermal and measurement energy; too much gain converts estimation noise into motion, so a finite optimum is expected.

**What to improve.** State explicitly that the code has set $\omega=1$ and uses ideal instantaneous, noiseless actuation with no delay or saturation. The residual energy is estimation/feedback noise within that ideal LQG model; a real actuator adds further noise and bandwidth limits. Add the simple open-loop formula above as a strong analytic check on the longer covariance calculation.

### 18. Linear versus nonlinear stochastic evolution

Under a reference probability measure $Q$, the unnormalized state satisfies

$$
d|\tilde\psi\rangle
=\left(-iH-\frac12c^\dagger c\right)|\tilde\psi\rangle dt
+c|\tilde\psi\rangle dW_Q.
$$

**Already explored.** The unnormalized ensemble reproducing the master equation, normalization of the mean weight, weight histograms, one path's direction and log weight, and effective-sample-size collapse are all important and well chosen.

**Most important results and consequences.** The weight

$$
w_t=\langle\tilde\psi_t|\tilde\psi_t\rangle
$$

is the Radon-Nikodym derivative

$$
w_t=\frac{dP}{dQ}\bigg|_{\mathcal F_t},
\qquad
\mathbb E_Q[w_t]=1.
$$

Consequently,

$$
\mathbb E_Q\!left[|\tilde\psi_t\rangle\langle\tilde\psi_t|\right]
=\rho(t).
$$

The phrase “unweighted average” is correct only in the Monte Carlo sense that no **additional** scalar is applied: each outer product already contains its likelihood in its norm.

**Issue.** Normalizing a path drawn under $Q$ gives the correct state direction associated with that canonical record, but an unweighted collection of those normalized paths is not distributed according to the physical measure $P$. The nonlinear physical trajectory requires both normalization and the Girsanov change of measure, equivalently replacing the reference noise by the physical innovation. Thus “the normalized state is the physical trajectory” needs the qualifier “after the associated measure change.”

The effective sample size

$$
N_{\mathrm{eff}}=\frac{(\sum_j w_j)^2}{\sum_jw_j^2}
$$

correctly diagnoses importance-sampling degeneracy.

**Assessment:** **[MISLEADING]** The code and ensemble identity are good; repair the pathwise-versus-distribution language and add a direct nonlinear-ensemble comparison under the physical measure.

### 19. Driven atom, averaged and watched

The general equation is

$$
\dot\rho=-i\left[\frac\Delta2\sigma_z+\frac\Omega2\sigma_x,\rho\right]
+\gamma(n_T+1)\mathcal D[\sigma_-]\rho
+\gamma n_T\mathcal D[\sigma_+]\rho
+\gamma k_d\mathcal D[\sigma_z]\rho.
$$

**Already explored.** Only the resonant, zero-temperature, zero-extra-dephasing corner is actually calculated. The catalog then compares one ideal pure homodyne trajectory with the mixed unconditional steady state. That comparison is correct but largely repeats Cases 2 and 8; it does not justify calling the general equation “fully worked.”

**Most important results and consequences.** Define

$$
\Gamma_1=\gamma(2n_T+1),
\qquad
\Gamma_2=\frac{\Gamma_1}{2}+2\gamma k_d,
\qquad
z_{\mathrm{eq}}=-\frac{1}{2n_T+1}.
$$

Then

$$
\begin{aligned}
\dot x&=-\Gamma_2x-\Delta y,\\
\dot y&=\Delta x-\Gamma_2y-\Omega z,\\
\dot z&=\Omega y-\Gamma_1(z-z_{\mathrm{eq}}).
\end{aligned}
$$

With

$$
D=\Gamma_1(\Gamma_2^2+\Delta^2)+\Omega^2\Gamma_2,
$$

the steady state is

$$
x_{\mathrm{ss}}=-\frac{\gamma\Omega\Delta}{D},
\qquad
y_{\mathrm{ss}}=\frac{\gamma\Omega\Gamma_2}{D},
\qquad
z_{\mathrm{ss}}=-\frac{\gamma(\Gamma_2^2+\Delta^2)}{D}.
$$

These formulas expose thermal inversion loss, detuning suppression, power broadening, and extra transverse linewidth. At $\Omega=0$, $z_{\mathrm{ss}}=z_{\mathrm{eq}}$; at strong drive the state saturates toward the center of the Bloch ball.

**Assessment:** **[MISLEADING]** The sentence “restore any of $\Delta$, $n_T$, or $k_d$ and the conditioned state becomes mixed” is false. Detuning is a Hamiltonian term and preserves purity. Thermal and dephasing channels cause conditional mixing only when their information is unmonitored or inefficiently monitored; they too can have pure-state unravellings when every relevant channel is monitored ideally.

**Required repair.** Either implement the general equation, plot the analytic steady state over $\Delta$, $n_T$, and $k_d$, and specify which channels are monitored, or relabel the case as a synthesis of the resonant ideal corner.

### 20. Mollow triplet and resonance-fluorescence squeezing

The resonant master equation is the same as Case 2, and the inelastic spectrum is

$$
S_{\mathrm{inel}}(\mu)
\propto\operatorname{Re}\int_0^\infty
e^{i\mu t}
\langle\delta\sigma_+(t)\delta\sigma_-(0)\rangle_{\mathrm{ss}}dt.
$$

**Already explored.** Weak/strong spectra, the triplet heat map, elastic weight, sideband finding, nonnegativity, and phase-sensitive homodyne squeezing are exactly the right physical targets. Qualitatively, these are the important phenomena.

**Most important results and consequences.** The relevant Bloch poles are

$$
-\frac\gamma2,
\qquad
-\frac{3\gamma}{4}\pm i\sqrt{\Omega^2-\frac{\gamma^2}{16}}.
$$

In the strong-drive limit the central line has half-width approximately $\gamma/2$, the sidebands have half-width approximately $3\gamma/4$, and their centers approach $\pm\Omega$. The steady coherent dipole and total excited population are

$$
|\langle\sigma_-\rangle_{\mathrm{ss}}|^2
=\frac{\Omega^2\gamma^2}{(\gamma^2+2\Omega^2)^2},
\qquad
\rho_{ee}^{\mathrm{ss}}
=\frac{\Omega^2}{\gamma^2+2\Omega^2}.
$$

Hence the coherent fraction of the emitted power is

$$
f_{\mathrm{el}}
=\frac{|\langle\sigma_-\rangle_{\mathrm{ss}}|^2}
{\rho_{ee}^{\mathrm{ss}}}
=\frac{\gamma^2}{\gamma^2+2\Omega^2},
$$

which tends to one for weak drive and zero for strong drive. The integrated inelastic fluctuation power is

$$
\langle\delta\sigma_+\delta\sigma_-\rangle_{\mathrm{ss}}
=\rho_{ee}^{\mathrm{ss}}-|\langle\sigma_-\rangle_{\mathrm{ss}}|^2.
$$

**Assessment:** **[FATAL]** for the present numerical spectra. All curves use the faulty `regressionSpectrum` helper. The qualitative Mollow-triplet claim is standard and correct, but the displayed line shapes, minima, and quantitative squeezing depth are not validated by the current code. Near exceptional points, the modal method is also numerically ill-conditioned even after changing `Transpose[big]` to `big`.

The phrase “spectrum from the record alone” is also inaccurate: the code computes a model-based quantum-regression spectrum, not a periodogram or correlation estimate from a simulated detector record.

**Required repair.** Use a trace-projected resolvent or direct time-domain regression; verify the spectrum against its equal-time sum rule and positivity; compare strong-drive widths with $\gamma/2$ and $3\gamma/4$; and obtain the same homodyne spectrum independently from long simulated records. Only then retain the squeezing-depth claim.

---

## Recommended revision order

1. Repair the regression/resolvent machinery, then regenerate Cases 11 and 20.
2. Rebuild Case 4 as an actual joint atom-cavity Ramsey protocol.
3. Correct the QPC distinction among pole frequency, spectral maximum, spectral merger, and exceptional point.
4. Fix Case 19's monitored-purity statement and either work the general model or narrow its title.
5. Relabel the dispersive IQ calculation or replace it with an integrated output-record simulation.
6. Add the compact analytic payoffs listed above to the otherwise sound cases.
7. Add step-size, ensemble-size, and cutoff convergence checks wherever an SDE, trajectory ensemble, or Fock truncation is used.

## Compact coverage ledger

| Case | What is already in good shape | Most important missing or incorrect item |
|---:|---|---|
| 1 | Coherence loss, Bloch contraction, steady manifold | $T_2$, purity law, physical treatment of degenerate steady states |
| 2 | Bloch equations, fixed point, modes, geometry | Saturation, flux, limits; explicitly realize the proposed ensemble |
| 3 | Exact loss channel, Wigner fringes, cutoff check | Energy law and final-vacuum clarification |
| 4 | Decay of a pre-inserted cat | Actual atom-cavity entanglement, Ramsey postselection, two-atom correlation |
| 5 | Pointer dynamics and joint dephasing | Integrated output SNR; current IQ-error interpretation is wrong |
| 6 | Covariance, diffusion, CP caveat, cat scaling | Explicit damping/decoherence rates and unit conventions |
| 7 | Click paths, waiting times, antibunching story | Flux, normalization, $g^{(2)}(0)$, timestep control |
| 8 | Record/state/innovation ledger | Pure start remains pure; phase/efficiency dependence |
| 9 | Two innovations and ensemble equivalence | Field quadratures are not direct simultaneous atomic $X,Y$ measurements |
| 10 | EP and Zeno turnover | Ensemble decay versus telegraph switching rate |
| 11 | Device interpretation and spectral normalization | Exact spectrum; merger occurs at $k=\Omega/\sqrt2$, not the EP |
| 12 | Born rule and record discrimination | Asymptotic rather than finite-time absorption |
| 13 | Riccati dynamics, tracking, conditional purity | Closed-form fixed point; whiteness is necessary, not sufficient |
| 14 | Dark-state pumping and efficiency sweep | Analytic rate and efficiency ceiling; finite-delay assumptions |
| 15 | Attractor, full thermal distribution, detailed balance | Exact transient $n(t)$ and $\langle a(t)\rangle$; Hamiltonian omitted in code |
| 16 | Both optimization objectives and distributions | Physical bandwidth and finite-step qualification |
| 17 | Full conditional/unconditional cooling bookkeeping | Open-loop analytic heating law; ideal-control assumptions |
| 18 | Linear ensemble, likelihoods, ESS | Reference versus physical probability measure |
| 19 | Ideal resonant trajectory versus average | General steady state; false detuning/purity statement |
| 20 | Correct target phenomena | Spectral implementation must be replaced and revalidated |

## Primary sources consulted for interpretation

- P. Rouchon and J. F. Ralph, [*Efficient quantum filtering for quantum feedback control*](https://arxiv.org/abs/1410.5345).
- J. Gambetta et al., [*Measurement-induced dephasing and number splitting in circuit QED*](https://arxiv.org/abs/cond-mat/0602322).
- K. Snizhko, P. Kumar, and A. Romito, [*Quantum Zeno effect appears in stages*](https://arxiv.org/abs/2003.10476).
- A. N. Korotkov and D. V. Averin, [*Continuous weak measurement of quantum coherent oscillations*](https://arxiv.org/abs/cond-mat/0002203).
- K. Jacobs, [*How to project qubits faster using quantum feedback*](https://arxiv.org/abs/quant-ph/0301056).
- H. M. Wiseman and J. F. Ralph, [*Reconsidering Rapid Qubit Purification by Feedback*](https://arxiv.org/abs/quant-ph/0603062).
