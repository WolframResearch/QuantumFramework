# Fresh Assessment of *Continuous Measurement in a Potential*

Source assessed: [`measurement-in-a-potential.md`](measurement-in-a-potential.md)

Source SHA-256: `01baf9ac52e6a448ed6057c8435de86baaa573006b3a92b6a6073537cba7c937`

Assessment date: August 4, 2026

This report supersedes the August 3 assessment. It evaluates the revised source from the beginning rather than assuming that an earlier objection still applies.

## Assessment brief

**Goal:** Determine whether the revised essay now gives a plain and accurate account of the physics, keeps recorded quantities separate from model estimates and simulator-only variables, and supports its formulas and simulations.

**Reader:** The author and serious students with basic quantum mechanics but no assumed background in stochastic equations, quantum filtering, or control theory.

**Success test:** At every point, the reader should know whether a quantity was recorded by electronics, obtained by calibration, inferred through a model, or created only by a simulation. Every continuum claim based on a finite calculation should have an invariant, limiting, or convergence check.

## Verdict

The revision has repaired the main physics failures. I found no remaining fatal error in the ideal-model equations. The reciprocal-squeezing claim is corrected, inefficient detection is treated correctly, the conditional covariance is distinguished from the wandering conditional mean, finite-grid failure at long times is demonstrated, and the agreement between the grid and the five-variable filter is now described as an internal same-model check.

The essay is nevertheless not ready in its present form. The remaining problems are:

1. Several sentences still call a model estimate “the real position” or imply that a Gaussian filter does not hold a quantum state.
2. The glossary says residual whiteness tests that the model is right. Whiteness is a necessary diagnostic, not proof of model validity.
3. The final dictionary says every row after the detector record is an inference, although its final row is explicitly simulator-only bookkeeping.
4. The essay is not concise or jargon-light. It is now 10,355 words with 90 Wolfram Language blocks, and the opening presents specialist names before the physical measurement chain is clear.
5. The displayed simulations are correct for the chosen parameters, but several reusable functions do not validate malformed inputs, extreme outcomes, or violations of their documented potential interface.
6. A few numerical and plotting claims are semantically mislabeled or lack visible acceptance criteria.

The current draft is technically strong and pedagogically overgrown. It needs a focused correction pass and substantial cutting, not another expansion.

## The plain physical account the essay should preserve

A detector records a noisy electrical signal. Calibration converts that signal into a position-sensitive time series. A mathematical model combines the calibrated record with the assumed dynamics and initial state to estimate the particle's conditional state. The estimate predicts future records and uncertainty; it is not a second detector reading.

A simulation creates both a synthetic record and a conditional state from the same model. Agreement between two algorithms using that record checks their internal consistency. It does not establish that the model describes a laboratory system. Experimental credibility requires calibration and tests using information that was not used to produce the estimate.

The revision now says most of this explicitly in line 28. That paragraph should become the organizing principle for the entire essay.

## What is recorded, calibrated, inferred, and simulated

| Quantity | Correct standing | Precise meaning |
|---|---|---|
| Photodetector voltage or digitized current | **Recorded** | The electronic time series produced by a particular detector |
| \(\Delta Y_k\) | **Calibrated from the record** | A scaled, integrated detector signal for one time bin |
| \(\bar x_k=\Delta Y_k/(2\sqrt{\lambda}\Delta t)\) | **Calibrated from the record** | A noisy position-sensitive number; it is not the particle's position during that bin |
| Innovation \(\Delta W_k=\Delta Y_k-2\sqrt{\lambda}\langle x\rangle_{k|k-1}\Delta t\) | **Computed residual** | The calibrated record minus the model's predicted conditional mean before that update |
| Conditional mean \(\langle x\rangle_c\) | **Model estimate** | The filter's estimate of the mean position conditioned on the record |
| Conditional covariance | **Model estimate** | The filter's stated residual uncertainty; deterministic for a Gaussian initial state in the ideal linear harmonic model, but still model dependent |
| Conditional state \(\psi_c\) or \(\rho_c\) | **Model estimate** | A compact probability assignment conditioned on the record and assumptions |
| Unconditional state \(\rho\) | **Ensemble prediction** | The state predicted when records are discarded or averaged over |
| Synthetic detector record | **Simulation output** | A model-generated counterpart of a calibrated laboratory record, not laboratory data |
| Grid state and `xGrid` | **Synthetic model estimate** | Model-generated conditional state and mean; their laboratory counterparts would also be inferred, not directly observed |
| Latent draw \(X\) in `drawExact` | **Simulation only** | One device for sampling the synthetic outcome distribution; no trajectory-wise laboratory counterpart |
| Random seed and underlying random-number stream | **Simulation only** | Reproducibility bookkeeping, not hidden laboratory variables |

### Remaining violations of this distinction

| Severity | Source | Current wording or implication | Correction |
|---|---|---|---|
| **MISLEADING** | Lines 56 to 58 | \(dY\) and \(\bar x\) are called “directly observed.” | The raw electronic current is directly recorded. \(dY\) and especially \(\bar x\) are calibrated and rescaled forms of it. |
| **MISLEADING** | Line 58 | The innovation is “recoverable,” and its whiteness tests that the model “is right.” | The innovation is computed from the record and the model estimate. Whiteness can reject a bad model but cannot prove that the model is uniquely correct. |
| **MISLEADING** | Lines 567 and 581 to 585 | The filter recovers the “real position,” and the position is “buried in noise.” | The record contains a signal whose expected value depends on the conditional mean. The filter estimates that mean under the model. |
| **MISLEADING** | Line 567 | The experiment “never holds the quantum state.” | In the Gaussian model, the two means and three covariances are the conditional quantum-state estimate. The computer avoids a grid wavefunction; it does not avoid representing a quantum state. |
| **MISLEADING** | Lines 642 and 652 | The filter is said to recover the mean from the record “alone” or from “nothing but” the readings. | The filter also requires a calibrated measurement scale, an initial state, a Hamiltonian, a measurement rate, and a noise model. The record selects one conditional estimate within that assumed model. |
| **IMPRECISE** | Lines 652 and 656 | The five-number calculation is described as what a real experiment runs and “the real thing.” | The cited experiment uses the same linear-Gaussian filtering principle, but a calibrated discrete-time state-space model with finite efficiency, external noise, feedback, and sometimes an augmented state for colored noise. |
| **INTERNAL CONTRADICTION** | Lines 863 to 871 | “Only the first row is data; every later row is inference,” followed by a row labeled “no counterpart.” | Say: rows two through four are model estimates; the last row is simulator-only bookkeeping. |
| **IMPRECISE** | Line 873 | “Pure conditional trajectories” are presented without immediately restoring the efficiency qualifier. | Pure state-vector trajectories apply to the unit-efficiency model. With missed information, the conditional trajectory is generally mixed. |

## 1. Jargon and concision

### Overall judgment

The revision improved precision by adding many explanations, derivations, and safeguards. It did not meet the requested standard of concise, low-jargon exposition.

The abstract alone asks a new reader to absorb “stochastic Schrödinger equation,” “SDE,” “Bayesian updating,” “Trotter-Suzuki,” “quantum Kalman filter,” “conditional position,” “quadratic world,” and “grid engine.” Most of those names are not needed to understand what is physically recorded.

The main text then adds “Wiener increment,” “Kraus operator,” “prior predictive,” “Lindblad,” “Strang,” “Riccati,” “innovation,” “a priori,” and “a posteriori.” Each term is legitimate. Their density is the problem.

### Recommended plain-language hierarchy

Use the plain description first, then the formal name once:

| Formal term | First explanation |
|---|---|
| continuous weak measurement | many rapid, individually noisy readings that also disturb the particle |
| Wiener increment | an ideal independent Gaussian noise increment with variance \(dt\) |
| innovation | the newest calibrated record minus the filter's prediction |
| conditional state | the state estimate after using the record observed so far |
| unconditional state | the prediction obtained when the record is ignored |
| Kraus operator | the position-dependent multiplier applied after one reading |
| Riccati equation | the deterministic update law for the uncertainty matrix |
| Strang splitting | half a potential step, one kinetic step, then another half potential step |

### Material that belongs in technical notes or an appendix

The following is correct but interrupts the central lesson:

- the full taxonomy of first-, second-, and higher-order operator splittings;
- the characteristic-polynomial proof that all covariance eigenvalues share a real part;
- the long finite-grid demonstration after the cutoff is reached;
- the full inefficient-detector symbolic mapping and asymptotic threshold analysis;
- the double-well double-commutator code;
- implementation details of `Sow` and `Reap`.

Keep the results in the main essay. Move the proofs and extended diagnostics to “Technical checks” sections or a companion notebook. A useful target is to cut the narrative by roughly one third and reduce the main sequence to the code cells that directly establish a physical claim.

## 2. Formula audit

### Core formulas that pass

Here \(\mu=2\hbar\lambda/(m\omega^2)\) is the dimensionless measurement strength used in the harmonic formulas, and \(\eta\) is the information efficiency.

The following were independently re-derived and recomputed:

1. The normalized unit-efficiency conditional state equation.
2. The record convention

   \[
   dY=2\sqrt{\lambda}\,\langle x\rangle_c\,dt+dW.
   \]

3. The expected log-likelihood-ratio, equivalently KL-divergence, rate \(2\lambda(\Delta x)^2\) for two fixed position hypotheses.
4. The Gaussian finite-step likelihood and amplitude update.
5. The record-averaged master equation

   \[
   \dot\rho=-\frac{i}{\hbar}[H,\rho]-\frac{\lambda}{2}[x,[x,\rho]].
   \]

6. The measurement contribution \(d\langle p^2\rangle/dt=\lambda\hbar^2\).
7. The free conditional steady width \((\hbar/4\lambda m)^{1/4}\).
8. The harmonic conditional mean and covariance equations.
9. The ideal harmonic steady covariance, purity identity, free limit, and squeezing ratio.
10. The correction to reciprocal squeezing:

    \[
    \frac{\Sigma_x}{\Sigma_x^{\mathrm{zp}}}
    \frac{\Sigma_p}{\Sigma_p^{\mathrm{zp}}}
    =1+\left(\frac{2C}{\hbar}\right)^2.
    \]

11. The fixed point's stability and common decay rate \(4\lambda\Sigma_x^{\mathrm{ss}}\).
12. The inefficient-detector covariance equations, determinant

    \[
    \Sigma_x\Sigma_p-C^2=\frac{\hbar^2}{4\eta},
    \]

    squeezing ratio, and threshold

    \[
    \mu^*=\frac{2\sqrt{1-\eta}}{\eta^{3/2}}.
    \]

13. The continuum heating law

    \[
    \frac{dE}{dt}=\frac{\lambda\hbar^2}{2m}
    \]

    for \(H=p^2/(2m)+V(x)\), with no damping or other environment.
14. The statement that the harmonic conditional covariance can settle while the conditioned center continues to diffuse.
15. The failure of five-moment closure for a generic nonlinear force.

No fatal algebraic or physical error was found in these results.

### Remaining formula and interpretation issues

#### IMPRECISE: the likelihood is written only up to proportionality

The omitted normalization is harmless when the update is immediately renormalized. It is not harmless when \(P(\bar x)\) is called a probability density. The normalized likelihood is

\[
P(\bar x\mid x)
=\sqrt{\frac{2\lambda\Delta t}{\pi}}\,
  e^{-2\lambda\Delta t(\bar x-x)^2}.
\]

State once that the constant is omitted only because it cancels from the conditional state update.

#### IMPRECISE: “distinguishability rate” needs a definition

The value \(2\lambda(\Delta x)^2\) is correct for the expected log-likelihood-ratio rate, equivalently the KL-divergence rate, between the two equal-variance record laws. “Distinguishability” by itself is not a unique metric; other statistical distances have other normalizations. Name the measure when the rate first appears.

#### IMPRECISE: “the middle term only keeps the norm”

The nonlinear drift is the normalization correction in the Itô equation, but it also contributes to the deterministic change of the normalized state. Calling normalization its origin is accurate; saying that normalization is its only visible effect is too narrow.

#### MISSING DERIVATION: why five moments fail in the double well

The numerical mismatch of a harmonic filter does not by itself prove that no five-variable exact filter exists. The proof is one line. For

\[
V(x)=\frac{(x^2-4)^2}{8},
\]

the mean momentum obeys

\[
\frac{d\langle p\rangle}{dt}
=-\frac{1}{2}\left(\langle x^3\rangle-4\langle x\rangle\right)+\text{measurement term}.
\]

The third moment is not determined by the five Gaussian moments once the state is non-Gaussian. Put this equation before the simulation. Then the double-well run illustrates a closure failure already established analytically.

#### NUMERICAL: unstable small-\(\mu\) squeezing formula

The defined function

\[
\frac{\sqrt{2(\sqrt{1+\mu^2}-1)}}{\mu}
\]

is algebraically correct but loses precision by subtracting nearly equal numbers. At machine \(\mu=10^{-8}\), the current function returns \(0\) instead of \(1\). Use

\[
\sqrt{\frac{2}{\sqrt{1+\mu^2}+1}}
\]

for numerical evaluation, and retain the original form only for symbolic presentation.

#### IMPRECISE: continuum density versus probability per grid point

Line 79 correctly states that the stored entry is \(\sqrt{\Delta x}\,\psi(x_j)\), so \(|\psi_j|^2\) is probability in one grid bin. Later:

- line 295 labels the stored amplitude as \(\psi(x)\);
- lines 310 to 317 call \(|\psi_j|^2\) a density and label the vertical axis \(|\psi|^2\).

For a continuum amplitude, plot \(\psi_j/\sqrt{\Delta x}\). For a continuum density, plot \(|\psi_j|^2/\Delta x\). The current plots have the correct shape but the wrong vertical meaning and scale.

#### IMPRECISE: “each factor is exact”

Line 246 should say that the dephasing subflow and the kinetic and potential subflows are exact separately. The combined Hamiltonian step is a second-order approximation when \(T\) and \(V\) do not commute. The full symmetric deterministic composition is still second order.

## 3. Simulation and Wolfram Language audit

### Fresh-kernel result

In my audit run, all 90 Wolfram Language blocks were extracted from the source hash named above, numbered, and evaluated in document order in a fresh Wolfram Language 15.0 kernel started without initialization files. That run produced no source-cell message, timeout, unevaluated load-bearing call, or `Indeterminate` in the displayed examples. The table below is the compact audit trail of the load-bearing outputs, not a block-by-block execution transcript.

The important outputs were:

| Claim tested | Fresh result |
|---|---:|
| amplitude update versus likelihood-times-prior | maximum difference \(1.73\times10^{-17}\) |
| free-run maximum edge and momentum-tail masses | \(1.58\times10^{-29}\), \(7.94\times10^{-31}\) |
| free steady-width bias at \(\Delta t=0.005,0.0025\) | \(1.7773\times10^{-3}\), \(8.9010\times10^{-4}\) |
| free unitary and \(p^2\) Hermiticity errors | \(5.55\times10^{-16}\), \(1.78\times10^{-15}\) |
| free momentum endpoint error | \(4.26\times10^{-13}\) |
| initial cat variance | \(36.99999945\) |
| final seeded cat mean | \(-7.05317743\) |
| maximum two-point-posterior versus grid-lobe difference | \(0.02980769\) |
| ideal position variance ratio at \(\mu=4\) | \(0.62481053\) |
| ideal momentum variance ratio at \(\mu=4\) | \(2.57615983\) |
| covariance fixed-point eigenvalues | \(-2.49924,\,-2.49924\pm3.20097i\) |
| harmonic grid covariance | \((0.31317677,0.39168865,1.28815428)\) |
| Riccati covariance | \((0.31240527,0.39038820,1.28807991)\) |
| grid uncertainty determinant | \(0.25000000\) |
| filter-grid maximum gap at \(\Delta t=0.002,0.001\) | \(0.00535910,\;0.00248550\) |
| inefficient variance ratio at \(\eta=0.34,\mu=4\) | \(1.28947844\) |
| threshold at \(\eta=0.34\) | \(8.19565348\) |
| trapped energies at \(t=0,0.5,1,1.5,2\) | \(0.5,1.00000004,1.50000027,2.00000071,2.50000119\) |
| long-run momentum-tail mass, start and end | \(1.41\times10^{-16},0.0724324\) |
| double-well variance range after the transient | \(0.151859\) to \(0.412712\) |
| double-well harmonic-filter error and mean range | \(0.519781\) over \(2.87997\) |
| double-commutator checks, harmonic and double well | \(-1.00000000,-1.00000000\) |
| four double-well energy slopes | \(0.999999,\;0.999995,\;1.000003,\;1.000001\) |

These results support the essay's present ideal-model claims.

### Simulation assumptions and how far they are checked

| Assumption | Standing in the essay | Audit judgment |
|---|---|---|
| One nonrelativistic degree of freedom with \(H=p^2/(2m)+V(x)\) | Defines the model | Legitimate scope choice, not something the simulation can verify experimentally |
| Markovian position monitoring with independent white Gaussian increments | Defines the record and backaction | Internally implemented consistently; a laboratory must test bandwidth, correlations, and the innovation spectrum |
| Perfect detection for every state-vector trajectory | Needed for purity | Correctly used in the trajectory code, but the concluding pure-trajectory language needs the \(\eta=1\) qualifier; the inefficient case is treated only at moment level |
| No damping, thermal bath, technical force noise, or feedback in the ideal examples | Needed for the stated heating law and steady conditional covariance | Mostly explicit; it sharply limits direct comparison with an experiment |
| Known initial state and exact parameters | Required by the synthetic filter comparison | Assumed, not inferred or stress-tested; the phrase “record alone” is therefore false |
| Periodic finite grid with finite momentum bandwidth | Numerical approximation | Edge and momentum-tail diagnostics are computed, spatial refinement is good in the harmonic control, and the long-run breakdown is correctly exposed |
| Vectorized, finite, real-valued potential on the grid | Code interface assumption | Works for every displayed potential but is documented rather than enforced |
| Split measurement and Hamiltonian substeps | Time-discretization assumption | Time-step refinement supports the claimed first-order trajectory bias; the Hamiltonian Strang substep itself is second order |
| Gaussian initial state and quadratic Hamiltonian for five-moment closure | Analytical restriction | Verified against the grid in the harmonic case; it is not valid for the cat or double well |
| One seeded path can illustrate nonlinear behavior | Presentation choice | Adequate for a qualitative example, not for ensemble or convergence claims |
| The latent \(X\) used by `drawExact` is only a sampling device | Algorithmic device | Correct; treating it as a hidden physical position would add an unsupported ontology |

The assumptions are therefore mostly coherent and increasingly explicit. What is missing is a compact declaration of them before the simulations, plus a visible statement beside each result saying whether it is an exact model identity, a finite-step numerical check, or a single synthetic example.

### Independent controls not shown in the essay

At fixed \(L=24\) and \(\Delta t=0.002\), the final harmonic variance was

| \(n\) | Final variance |
|---:|---:|
| 128 | \(0.313176774503953\) |
| 256 | \(0.313176774503950\) |
| 512 | \(0.313176774503950\) |

This confirms that the visible harmonic offset is not a spatial-grid error at these resolutions.

At fixed \(n=256,L=24\), time refinement gave

| \(\Delta t\) | Final variance |
|---:|---:|
| 0.004 | \(0.3139630817\) |
| 0.002 | \(0.3131767745\) |
| 0.001 | \(0.3127848449\) |
| analytic steady value | \(0.3124052669\) |

The result converges toward the analytic covariance as claimed. This table should appear in the essay if “independent of grid spacing” and “first-order splitting bias” remain in the prose.

The double-well qualitative conclusions also survived \(n=128,256,512\) at the displayed time step. Its single-trajectory numbers do not converge monotonically under a change of \(\Delta t\) when the simulations use the same integer seed but not a coupled physical noise path. For stochastic pathwise convergence, use coupled Wiener increments or compare ensemble statistics rather than unrelated seeded trajectories.

### Code defects and robustness gaps

| Severity | Code | Finding | Improvement |
|---|---|---|---|
| **CORRECTNESS OUTSIDE DEMO** | `grid` | The function silently requires positive box length and a positive even integer \(n\). `grid[5,10.]` fails. | Add patterns, validation, and messages, or implement odd grids. |
| **CORRECTNESS OUTSIDE DEMO** | `measUpdate` | A sufficiently extreme outcome underflows every likelihood weight. `Normalize` then silently returns the zero vector. | Subtract the maximum log weight before exponentiation and reject a zero norm. |
| **API LIMITATION** | `potentialStep` | The source explicitly requires a vectorized potential, so rejecting `v[x_?NumericQ] := x^4` is documented behavior. The function does not detect a contract violation and instead returns a nonnumeric expression. | Keep the fast vectorized contract, but validate it and issue a clear message; alternatively map scalar potentials over the grid. |
| **NUMERICAL** | `sqRatio` | Machine evaluation near zero suffers catastrophic cancellation. | Use the stable rationalized form. |
| **MISSING VALIDATION** | public functions | No checks enforce \(\lambda>0,\Delta t>0,s>0\), normalized nonzero state vectors, numeric finite potential values, or compatible dimensions. | Validate inputs at the public boundary. |
| **MISSING ACCEPTANCE TESTS** | numerical demonstrations | The prose interprets raw numbers, but the code does not state tolerances or return pass/fail results. | Use `VerificationTest` or explicit associations containing value, target, error, tolerance, and status. |
| **STATE RISK** | document globals | Many examples depend on global `hbar`, `mass`, grids, rates, and time steps. The new separation of `λR` and `λH` fixes the earlier worst collision, but out-of-order reevaluation remains risky. | Localize each complete example or pass a parameter association. |
| **API RISK** | `trajectory` | Untagged `Sow`/`Reap` makes the return shape depend on the caller's wrapper. | Return an association with explicit `"States"` and `"Record"` keys. |
| **IMPRECISE NAME** | `drawExact` | “Exact” refers only to sampling the discrete finite-step outcome law, not the continuum dynamics or an experimental truth. | Rename to `drawSyntheticOutcome`. |
| **GLOBAL COUPLING** | `Emeas`, `EmeasV` | Both are tied to global grids and precomputed operators. | Pass the grid and operators explicitly. |
| **PLOT TIMING** | record plots | The first reading is plotted at \(t=0\), although it belongs to the first interval ending at \(t=\Delta t\). | Plot readings at \(\Delta t\{1,\ldots,n_t\}\). |
| **INCOMPLETE DENSITY TEST** | density-matrix table | Positivity and trace are shown, but Hermiticity is not printed and taking only real parts can hide a complex residual. | Report trace, Hermiticity norm, minimum eigenvalue, and tolerances together. |

### What the simulations can and cannot establish

The present simulations support:

- correct implementation of the stated ideal equations on the displayed and independently tested inputs;
- agreement between grid and moment representations in the Gaussian harmonic case;
- the expected time-step bias;
- the distinction between conditioned and record-averaged evolution;
- failure of a harmonic Gaussian filter in the displayed nonlinear example;
- the finite grid's eventual failure to reproduce continuum heating.

They do not establish:

- that a real detector has white, instantaneous, unit-efficiency noise;
- that experimental parameters and initial conditions are known exactly;
- that the filtered conditional covariance is correct in a particular run;
- that innovation whiteness uniquely validates a model;
- that one seeded double-well trajectory represents ensemble behavior;
- that simulator-only latent draws correspond to hidden physical positions.

## 4. Comparison with the cited experiment

The primary paper, [*Real-time optimal quantum control of mechanical motion at room temperature*](https://arxiv.org/abs/2012.15188), supports these statements:

- the physical oscillator is described by position and momentum means plus a conditional covariance;
- a discrete-time Kalman estimator ran on an FPGA with \(32\,\mathrm{ns}\) sampling;
- the reported conditional position standard deviation was \(1.30\,z_{\mathrm{zpf}}\);
- the information efficiency was about \(0.34\);
- innovation spectra, calibrated noise models, and out-of-loop Raman sideband thermometry were used to test the state estimate;
- colored noise required an extended state-space model.

The revision correctly converts the \(1.30\) standard-deviation ratio to about \(1.7\) in variance and now explains independent validation. What remains too strong is the claim that the laboratory literally runs this essay's five-number filter. The shared object is the linear-Gaussian filtering structure. The laboratory implementation includes calibration, discrete hardware, feedback, finite efficiency, environmental forces, and an augmented colored-noise model.

## Part-by-part assessment

### Abstract and Part I

**Strong:** Line 28 gives the clearest account in the essay of data, inference, and synthetic records. The rate convention and discretization warning are now explicit.

**Remaining:** The abstract is still formalism-first and jargon-heavy. The glossary calls calibrated quantities directly observed and overstates the evidential power of innovation whiteness.

**Prescription:** Begin with line 28's distinction in plainer language, then introduce the record equation, and only then the conditional state equation.

### Part II

**Strong:** The grid convention, periodic boundary, Fourier convention, discrete-amplitude normalization, edge mass, momentum-tail mass, vectorized-potential contract, and split-order scope are now mostly correct and explicit.

**Remaining:** Plot labels later ignore the stated \(\sqrt{\Delta x}\) convention. Input contracts are documented but not enforced. “Any potential” remains broader than the accepted function interface.

**Prescription:** Make the public functions safe, correct the continuum-density plotting scale, and return structured diagnostics with tolerances.

### Part III

**Strong:** The free width bias is separated into transient and time-step effects. The momentum operator is verified. The cat inference is honestly introduced as a two-point approximation and compared with the grid lobe mass.

**Remaining:** “Collapse” and “measurement forces a definite choice” remain interpretation-heavy. The approximation differs from the grid lobe weight by as much as \(0.0298\), which should be printed rather than described only as “almost nothing.”

**Prescription:** Use “conditional weight concentrates in one lobe,” print the maximum approximation error, and label every record and state plot as synthetic.

### Part IV

**Strong:** This is the best section. The ideal covariance, nonreciprocal fixed-axis squeezing, stability, grid cross-check, finite-step filter convergence, inefficient-detector mapping, purity loss, and threshold are correct.

**Remaining:** The section still says “real position,” denies that the five variables represent a quantum state, and describes the real experimental filter too literally. The squeezing function is numerically unstable near zero.

**Prescription:** Replace “real position” with “conditional mean,” say “five variables represent the full Gaussian conditional state,” use the stable ratio, and describe the experiment as a richer implementation of the same filtering principle.

### Part V and the double well

**Strong:** The potential-independent continuum heating law, finite-grid cutoff failure, and conditioned-center heating check. The displayed double-well outputs reproduce and their qualitative behavior survives spatial refinement. The essay now clearly separates short-time continuum agreement from long-time finite-grid failure.

**Remaining:** The analytic reason for moment-hierarchy failure is asserted rather than shown. The nonlinear demonstration uses one main trajectory, and the time-step comparison of stochastic paths is not controlled by shared noise increments.

**Prescription:** Add the \(\langle x^3\rangle\) equation, label the double-well numbers as one synthetic record, and use ensembles for quantitative nonlinear claims.

### Part VI and conclusion

**Strong:** The final dictionary now labels the detector record, state estimate, conditional mean, conditional variance, and simulator-only latent variables.

**Remaining:** Its opening sentence contradicts its last row. It omits the calibrated-record layer between electronics and \(\bar x\), and it uses “a priori” and “a posteriori” after the essay's stated goal of plain language.

**Prescription:** Replace the table with the four-level evidence table near the beginning of this report, using “recorded,” “calibrated,” “model estimate,” and “simulation only.”

## Priority revision plan

### Priority 1: fix remaining conceptual language

1. Replace every “real position” or “position buried in noise” claim with a statement about the conditional mean signal and its model-based estimate.
2. Say explicitly that five Gaussian moments are a representation of the conditional quantum state.
3. Remove “record alone”: state the required calibration, model, parameters, and initial condition.
4. Change innovation whiteness from “the model is right” to “a necessary residual diagnostic.”
5. Fix the contradiction in Part VI.
6. Qualify pure trajectories by \(\eta=1\).

### Priority 2: correct numerical semantics and robustness

1. Use the stable squeezing formula.
2. Correct continuum-amplitude and density plot scales.
3. Add input validation and log-stable measurement weights.
4. Add explicit numeric tolerances and pass/fail outputs.
5. Print the harmonic spatial-refinement table and the cat approximation error.
6. Add the one-line moment-hierarchy proof before the double-well simulation.

### Priority 3: cut the essay

1. Move lengthy splitting taxonomy, fixed-point stability proof, long-cutoff example, and extended efficiency algebra to technical notes.
2. Keep one physical question, one minimal derivation, one simulation, and one verification table per main section.
3. Replace repeated claims that a plot “shows” something with a number and tolerance.
4. Remove specialist names from the abstract unless the abstract defines them in plain language.

## Publication checklist

- [ ] Raw current, calibrated record, model estimate, and simulator-only quantities are separated.
- [ ] \(\bar x\) is never called a direct position measurement without calibration.
- [ ] The innovation is described as a model-dependent residual and whiteness as necessary but insufficient.
- [ ] The five Gaussian variables are correctly called a quantum-state representation.
- [ ] “Real position” is removed from the filtering section.
- [ ] The filter is never said to use the record alone; its calibration, model, parameters, and initial condition are named.
- [ ] The unit-efficiency qualifier accompanies every pure-trajectory statement.
- [ ] Continuum amplitude and density plots include the proper \(\Delta x\) factors.
- [ ] The small-\(\mu\) formula is numerically stable.
- [ ] Public functions validate their domains and cannot silently return a zero state.
- [ ] Numeric checks state tolerances and pass/fail status.
- [ ] Harmonic spatial refinement is visible in the essay.
- [ ] The double-well closure failure is established analytically before being illustrated numerically.
- [ ] Quantitative stochastic claims use coupled paths or ensemble statistics.
- [ ] The cited experiment is presented as a richer state-space filter, not literally the essay's five-variable code.
- [ ] The main narrative is substantially shorter and jargon is defined only after the physical need appears.

## Final opinion

The revised essay is no longer confused about its central mathematics. Its ideal continuous-measurement model is coherent, and the simulations now test many of the right invariants and limits. The earlier false squeezing and efficiency claims are gone.

The remaining weakness is a mismatch between the essay's stated pedagogical goal and its form. The piece wants to teach what is measured, what is inferred, and what is only simulated, but it still sometimes calls an estimate “the real position,” and it now surrounds the lesson with more formal machinery than a new reader needs.

The next pass should be subtractive. Correct the few remaining operational statements, harden the code, move extended proofs out of the main path, and let the detector-to-estimate chain carry the narrative:

\[
\text{electronic record}
\longrightarrow
\text{calibrated record}
\longrightarrow
\text{model-based state estimate}
\longrightarrow
\text{independent statistical test}.
\]

The simulation belongs beside that chain:

\[
\text{assumed model}
\longrightarrow
\text{synthetic record and state}
\longrightarrow
\text{consistency, invariants, and convergence checks}.
\]

With those changes, the essay can become both technically reliable and genuinely clear.
