# Comparing Stochastic Integrators for Continuous Quantum Measurement

## Question 1

> Consider a stochastic differential equation, such as a state diffusion representing a continuous measurement of an observable. Suppose I have different solvers or integrators. How can I compare their performance, and how can I verify that they behave as expected?

For a stochastic state-diffusion solver, "performance" has at least three distinct meanings. All mathematical expressions below use standard TeX delimiters: `\( ... \)` inline and `\[ ... \]` for display equations.

1. Does it converge to the correct stochastic process?
2. Does it preserve the physical structure of the state?
3. How much computation does it require for a given accuracy?

Suppose the conditional quantum state satisfies an Itô stochastic master equation

\[
d\rho_t=\mathcal L(\rho_t)\,dt+
\sqrt{\eta}\,\mathcal H[c](\rho_t)\,dW_t,
\]

with measurement record

\[
dY_t=\sqrt{\eta}\,
\operatorname{Tr}\!\left[(c+c^\dagger)\rho_t\right]dt+dW_t.
\]

Here is a practical validation framework.

### 1. Establish a reference solution

For each realization, generate Wiener increments on a very fine grid and compute a reference trajectory using:

- a substantially smaller time step;
- preferably a higher-order method;
- or an exact solution when one is available.

All solvers must use the same noise path. If a coarse interval contains fine increments, construct

\[
\Delta W_{\mathrm{coarse}}
=
\sum_j \Delta W_{\mathrm{fine},j}.
\]

Generating independent noise for each step size makes pathwise error comparisons meaningless.

Also refine the reference solution once more to verify that its remaining error is negligible.

### 2. Measure strong convergence

Strong convergence asks whether individual trajectories are reproduced:

\[
E_{\mathrm{strong}}(h)
=
\mathbb E\!\left[d\!\left(\rho_T^{(h)},\rho_T^{\mathrm{ref}}\right)\right].
\]

Possible state distances include

\[
\|\rho-\sigma\|_F,
\qquad
\frac{1}{2}\|\rho-\sigma\|_1,
\qquad
1-F(\rho,\sigma).
\]

Repeat over many shared noise paths and fit

\[
E_{\mathrm{strong}}(h)\approx C h^p.
\]

A log-log plot of error against \(h\) should have slope \(p\). For ordinary Itô SDEs, Euler-Maruyama usually has strong order \(1/2\); Milstein can reach order \(1\) under suitable assumptions. A quantum-state projection or renormalization step can alter these orders, so measure them rather than assuming them.

### 3. Measure weak convergence

Often the goal is not to reproduce a particular trajectory but to reproduce expectation values and distributions. Test several physically relevant functions:

\[
E_{\mathrm{weak}}^\phi(h)
=
\left|
\mathbb E[\phi(\rho_T^{(h)})]
-
\mathbb E[\phi(\rho_T^{\mathrm{ref}})]
\right|.
\]

Useful choices are

\[
\phi(\rho)=\operatorname{Tr}(A\rho),
\qquad
\operatorname{Tr}(\rho^2),
\qquad
[\operatorname{Tr}(A\rho)]^2,
\]

as well as threshold-crossing probabilities, collapse times, and final measurement probabilities.

Report Monte Carlo confidence intervals. Otherwise, sampling noise can easily be mistaken for discretization error. Increase the number of trajectories until the statistical uncertainty is substantially smaller than the numerical bias being measured.

### 4. Check quantum-state constraints

At every step, monitor

\[
|\operatorname{Tr}\rho-1|,
\qquad
\|\rho-\rho^\dagger\|,
\qquad
\max(0,-\lambda_{\min}(\rho)).
\]

Also check that:

- purity remains in its allowed range;
- expectation values remain inside the observable's spectral bounds;
- conserved quantities remain conserved when the model requires it;
- pure states remain pure for an ideal, fully efficient measurement when the underlying stochastic evolution has that property.

Do not look only at averages. Record maximum violations and violation frequencies. A solver can have accurate ensemble averages while producing many unphysical individual states.

Blindly projecting negative eigenvalues to zero may hide instability. If projection is used, measure how often it occurs and how large the correction is.

### 5. Check measurement-specific identities

Continuous-measurement equations provide unusually strong consistency tests.

The unconditional average should satisfy the corresponding deterministic master equation:

\[
\bar\rho_t=\mathbb E[\rho_t],
\qquad
\frac{d\bar\rho_t}{dt}=\mathcal L(\bar\rho_t).
\]

Compare the trajectory average against a high-accuracy deterministic solution.

The reconstructed innovations

\[
\Delta W_n
=
\Delta Y_n
-
\sqrt{\eta}\,
\operatorname{Tr}\!\left[(c+c^\dagger)\rho_n\right]\Delta t
\]

should have approximately

\[
\mathbb E[\Delta W_n]=0,
\qquad
\operatorname{Var}(\Delta W_n)=\Delta t,
\]

with negligible temporal autocorrelation. Systematic drift or colored innovations often reveal an inconsistent record update, incorrect normalization, or an Itô/Stratonovich mistake.

For a quantum nondemolition measurement, also test that:

- trajectories collapse toward the appropriate eigenspaces;
- final outcome frequencies agree with the initial Born probabilities;
- the relevant conditional populations have the expected martingale behavior.

### 6. Use limiting and exactly solvable cases

A good test suite should contain:

- zero measurement strength, reducing the problem to a deterministic master equation;
- zero Hamiltonian or zero dissipation;
- an additive-noise SDE with a known solution;
- a commuting or QND measurement model;
- a linear unnormalized filtering equation, when available;
- short-time moment expansions that can be checked analytically.

These cases isolate different pieces of the implementation.

### 7. Test stability separately from convergence

Increase the integration time and test stiff parameter regimes. Record:

- numerical blow-ups;
- negative eigenvalues;
- excessive renormalization;
- sensitivity to step size;
- biased steady states;
- incorrect rare-event or collapse-time distributions.

A method can converge asymptotically as \(h\to0\) yet be unusable at practical step sizes.

### 8. Compare accuracy against cost

For each solver, measure wall-clock time, the number of drift and diffusion evaluations, memory use, and accepted or rejected steps. Then plot quantities such as

\[
\text{strong error versus runtime},
\qquad
\text{weak error versus runtime},
\qquad
\text{physical violations versus runtime}.
\]

This is more informative than error versus step size alone.

| Property | Suggested metric |
|---|---|
| Pathwise accuracy | Mean trace or Frobenius distance using shared noise |
| Statistical accuracy | Bias in observables and distributions |
| Convergence order | Slope of \(\log(\text{error})\) versus \(\log(\text{step size})\) |
| Trace preservation | Maximum and RMS trace error |
| Positivity | Minimum eigenvalue and violation rate |
| Ensemble consistency | Distance from deterministic master equation |
| Record consistency | Innovation mean, variance, and autocorrelation |
| Efficiency | Error versus wall-clock time |
| Robustness | Failure rate across long and stiff runs |

The strongest evidence is the combination of the expected strong and weak convergence slopes, shrinking physical violations, correct unconditional dynamics, white innovations, correct limiting-case behavior, and favorable accuracy per runtime. No single test is sufficient by itself.

---

## Question 2

> The solvers are designed in such a way that I cannot fit or provide the same generated Wiener process to all of them.

That changes the comparison strategy. If the solvers cannot consume the same realization of the Wiener process, then a direct trajectory-by-trajectory comparison is not valid. Two correct solvers driven by independent Wiener processes will generally produce very different trajectories.

### 1. With the current interfaces, compare laws rather than paths

Use independent ensembles and test weak convergence. For each solver, step size \(h\), and observable \(\phi\), estimate

\[
\widehat{\mu}_{h}
=
\frac{1}{N}\sum_{k=1}^{N}\phi\!\left(\rho_{T,k}^{(h)}\right).
\]

Compare these ensemble means with an analytic value, a deterministic unconditional master-equation solution, or a highly refined ensemble:

\[
B_h^\phi=
\left|\widehat{\mu}_{h}-\mu_{\mathrm{reference}}\right|.
\]

The Monte Carlo standard error is

\[
\operatorname{SE}(\widehat{\mu}_h)
=
\frac{s_h}{\sqrt N}.
\]

When comparing solvers \(A\) and \(B\), the uncertainty in their difference is approximately

\[
\operatorname{SE}(\widehat{\mu}_A-\widehat{\mu}_B)
=
\sqrt{\frac{s_A^2}{N_A}+\frac{s_B^2}{N_B}}.
\]

An observed solver difference is meaningful only when it is appreciably larger than this uncertainty.

Test several quantities, not just the mean state:

- \(\mathbb E[\operatorname{Tr}(A\rho_t)]\);
- variances and higher moments;
- purity distributions;
- final measurement probabilities;
- collapse or first-passage times;
- complete histograms or empirical cumulative distributions.

For scalar quantities, compare distributions using Wasserstein distance, the Kolmogorov-Smirnov statistic, or quantile differences. For quantum states, comparing a sufficiently rich collection of observable statistics is often more interpretable than comparing matrix entries directly.

The most important measurement-specific reference is

\[
\mathbb E[\rho_t]=\bar\rho_t,
\qquad
\dot{\bar\rho}_t=\mathcal L(\bar\rho_t).
\]

Solve this deterministic master equation accurately and compare each stochastic solver's ensemble mean against it.

### 2. Perform self-convergence for each solver

Even without shared noise between solvers, run each solver at

\[
h,\qquad h/2,\qquad h/4,\qquad h/8
\]

and examine weak errors in the same observables. You should see stabilization or a power law,

\[
B_h^\phi\sim C h^q.
\]

Because the ensembles are independent, Monte Carlo uncertainty may obscure this trend. Increase the number of trajectories as \(h\) decreases so that sampling uncertainty remains below the expected discretization error.

You cannot reliably compute a strong convergence order from independently generated paths. That limitation should be stated explicitly in any benchmark.

### 3. Apply path-independent correctness checks

These remain valid without coupled Wiener processes:

- trace preservation;
- Hermiticity;
- positivity and minimum eigenvalues;
- purity bounds;
- correct unconditional master-equation evolution;
- correct Born probabilities in QND examples;
- correct stationary distribution;
- correct measurement-record mean and variance;
- white reconstructed innovations;
- deterministic behavior when measurement strength is zero;
- robustness over long integrations.

These checks may reveal errors that an ensemble-average comparison misses.

### 4. Compare accuracy against cost

Because independent ensembles are required, compare the solvers using

\[
\text{weak error versus runtime},
\qquad
\text{distribution error versus runtime},
\qquad
\text{physical violations versus runtime}.
\]

For adaptive solvers, compare at common output times and several requested tolerances. Do not assume that equal nominal tolerances correspond to equal accuracy.

### If modifying the architecture is possible

The clean design is to separate the noise source from the integrator. The solver should request values from a common Brownian-path object rather than generate random numbers internally.

For adaptive solvers, the object can construct intermediate values with a Brownian bridge:

\[
W(t)\mid W(a),W(b)
\sim
\mathcal N\!\left(
W(a)+\frac{t-a}{b-a}[W(b)-W(a)],
\frac{(t-a)(b-t)}{b-a}
\right).
\]

Different integrators can then request different time points while remaining coupled to the same underlying random function \(W(t)\). Simply giving both solvers the same random seed is generally insufficient: if they request random numbers in different orders or take different numbers of steps, their paths immediately diverge.

For multidimensional, noncommuting noise, sharing only \(W(t)\) may be insufficient for higher-order methods; consistent iterated stochastic integrals or Lévy areas may be needed. For a single measurement channel, a shared Brownian path is normally enough.

With the current design, the defensible comparison is ensemble-based weak convergence plus physical consistency and cost. Strong pathwise comparison requires an architectural change that makes the noise driver external or replayable.

---

## Question 3

> Regarding the path-independent correctness checks, if the solvers generate the same stochastic average density matrix, should everything else also be the same?

No. Equality of the average density matrix is much weaker than equality of the stochastic dynamics.

Consider two qubit ensembles:

\[
\mathrm{A}:\qquad \rho=\frac{I}{2}\quad\text{on every trajectory},
\]

and

\[
\mathrm{B}:\qquad
\rho=
\begin{cases}
\lvert 0\rangle\!\langle 0\rvert,&p=\frac{1}{2},\\
\lvert 1\rangle\!\langle 1\rvert,&p=\frac{1}{2}.
\end{cases}
\]

Both have exactly the same average,

\[
\mathbb E[\rho]=\frac{I}{2}.
\]

But they are very different:

\[
\mathbb E[\operatorname{Tr}(\rho^2)]
=
\begin{cases}
\frac{1}{2},&\mathrm{A},\\
1,&\mathrm{B},
\end{cases}
\]

and, for \(z=\operatorname{Tr}(\sigma_z\rho)\),

\[
\operatorname{Var}(z)
=
\begin{cases}
0,&\mathrm{A},\\
1,&\mathrm{B}.
\end{cases}
\]

The average density matrix determines every ensemble-averaged quantity that is linear in the state:

\[
\mathbb E[\operatorname{Tr}(A\rho_t)]
=
\operatorname{Tr}\!\left(A\,\mathbb E[\rho_t]\right).
\]

It does not determine:

- purity and entropy distributions;
- variances of conditional expectation values;
- collapse probabilities and collapse times;
- measurement-record statistics;
- temporal correlations;
- state-record correlations;
- positivity of individual trajectories;
- nonlinear functions of the state.

This is central to continuous measurement: many distinct stochastic unravellings can produce the same unconditional master equation. Homodyne, heterodyne, and jump-type trajectories may share the same ensemble density operator while having very different individual paths and measurement records.

There is an additional subtlety: even if two solvers reproduce the full distribution of \(\rho_t\) separately at every time \(t\), they might still have different temporal dynamics. The joint distributions

\[
P(\rho_{t_1},\rho_{t_2},\ldots,\rho_{t_n})
\]

must agree to establish equality of the stochastic processes. Equal one-time marginals do not guarantee equal autocorrelations or first-passage-time distributions.

The appropriate hierarchy is:

1. **Same \(\mathbb E[\rho_t]\):** same means of linear observables.
2. **Same one-time state distributions:** same instantaneous statistics.
3. **Same multi-time joint distributions:** same state process in law.
4. **Same joint state-measurement-record distributions:** same continuous-measurement model.

If the only objective is to reproduce the unconditional master equation, checking \(\mathbb E[\rho_t]\) may be sufficient. If the solvers are intended to produce conditional states and measurement trajectories, it is only a necessary test, not a sufficient one.
