# SDE integration schemes for continuous quantum measurement: audit and recommendation

Scope: numerical schemes for the diffusive stochastic master equation (SME), both the
conditioned state and the homodyne/photocurrent readout, judged against the baseline
already in this folder: the positivity-preserving integrator of Rouchon and Ralph,
*Efficient Quantum Filtering for Quantum Feedback Control*, Phys. Rev. A **91**, 012118
(2015), arXiv:1410.5345, implemented QF-free as `\[ScriptCapitalT]` in
`manual-implementation-ito-qf-independent.wl` (see `AUDIT.md`).

Companion code written for this audit: `sde-schemes-prototype.wl` (the recommended
pure-state map `\[ScriptCapitalP]` and an Euler-Maruyama foil `\[ScriptCapitalE]`) and
`prototype-verify-benchmark.wls` (every number quoted in section 6 is produced by one run
of that script; all checks pass).

---

## 1. Verdict first

The baseline is a first-order, unconditionally completely-positive (CP) Kraus-map step:
`rho -> N[M rho M^dag + dt (dissipators)]`, a sum of `A rho A^dag` terms over a trace
normalization, so `rho` stays a valid density matrix for **any** step size `dt`. It is
dimension-general and already handles multiple channels, inefficient detection
`eta < 1`, unmonitored Lindblad operators `V`, and non-commuting operators. Its
ensemble mean reproduces the Lindblad master equation to `O(dt)`.

Nothing in the literature dominates it on all axes at once. What beats it, and where:

| Regime you are in | Use | Why it wins | Cost of switching |
|---|---|---|---|
| Perfect detection `eta=1`, every channel monitored, pure initial state | **Pure-state (state-vector) map** `\[ScriptCapitalP]` | Exactly the same conditioned trajectory, but a ket (O(d) memory) and matrix-vector steps (fewer/no O(d^3) matmuls). Measured ~1.4x at qubit scale rising to ~3.3x at `d=32` and still growing; memory O(d) not O(d^2). | None on accuracy or structure. Only valid while the conditional state is pure. |
| You need accurate **means / expectation values** (Lindblad reconstruction, feedback cost functions) | **Rouchon `\[ScriptCapitalT]` + Talay-Tubaro (Richardson) extrapolation** of the ensemble mean | Upgrades the first-order mean to weak order 2 at ~1.5x cost, keeping the underlying trajectories positive. No bespoke second-order scheme needed. | Two runs at `dt` and `dt/2`; the extrapolated *mean* is order 2, individual trajectories are unchanged. |
| You need the **second-order-accurate conditioned state itself** and still want guaranteed positivity | **W-map** (Wonglakhon-Wiseman-Chantasri 2024, arXiv:2408.14105) | Unconditionally CP by Kraus construction, average evolution correct to `O(dt^2)` (vs `O(dt)` for Rouchon-Ralph). | More work per step (a fourth-order bath expansion); per-trajectory (strong) order is *not* improved. |
| Quadratic `H`, linear measurement, Gaussian state (free particle, harmonic trap, linear optomechanics) | **Gaussian moment / Kalman-Bucy filter** | Exact closure in a handful of moments; cost is `O(1)` in the Hilbert dimension, not `O(d^2)` or `O(d^3)`. Dramatically cheaper than any Hilbert-space integrator. | Fails the moment the dynamics are nonlinear or the state is non-Gaussian. Already implemented for the free particle and harmonic trap in the companion essays. |
| The record is **photon counts**, not a continuous current | **Quantum-jump / MCWF** (Dalibard-Castin-Molmer) | It is the correct unraveling for a counting record; a diffusive integrator solves a different measurement. | Different physics, not a drop-in for the homodyne SME. |

Where **Rouchon-Ralph is already the right choice and nothing beats it**: a mixed
conditional state (`eta < 1`, unmonitored `V`, or a mixed seed) at small-to-moderate
`d`, where you want a guaranteed-physical state at any step size and per-trajectory
(strong) accuracy. No CP scheme improves the strong (per-trajectory) order over the
baseline; that remains an open problem. The only honest ways to "beat" it there are
(a) extrapolate for better *means* (Talay-Tubaro, still positive), or (b) accept the
W-map's heavier step for a second-order-accurate *conditioned state*.

Two non-starters, both confirmed by prior work in this course: generic higher-order SDE
solvers (Roessler SRK, Platen, higher Milstein) buy strong/weak order by giving up the
one thing the baseline is for, positivity, and go non-physical at moderate `dt`; and
WL-native `ItoProcess`/`RandomFunction`/stochastic `NDSolve` drift the norm, stiffen,
and time out on these systems.

---

## 2. The axes that matter

A scheme for this problem is judged on:

1. **Weak (average) order** `p_w`: `|E[f(rho_dt)] - E[f(rho_exact)]| = O(dt^{p_w})`. Governs how well the trajectory ensemble reproduces the Lindblad master equation and any expectation value or feedback cost.
2. **Strong (per-trajectory) order** `p_s`: `E[||rho_dt(T) - rho_exact(T)||] = O(dt^{p_s})`. Governs how faithfully a *single* trajectory tracks the true conditioned state for a given noise realization: the relevant quantity for a filter driven by a *real* measurement record.
3. **Unconditional structure preservation**: does the step keep `rho` positive, Hermitian, unit-trace for **any** `dt`, or only below a stability threshold?
4. **Cost and real-time / feedback readiness**: flops and memory per step; whether the step is cheap and branch-free enough to run inside a feedback loop.
5. **Generality**: many channels, `eta < 1`, unmonitored `V`, non-commuting operators, and scaling to higher Hilbert dimension `d`.
6. **Record fidelity**: is the simulated homodyne current `dy` the physically correct innovation, and is it reconstructed consistently with the state update?

The Rouchon-Ralph baseline scores: `p_w = 1`; `p_s` the usual multiplicative-noise
`1/2`-to-`1` (not its selling point); structure preservation **unconditional**; cost
`O(d^3)`/step (dominated in practice by WL dispatch at small `d`, see section 6);
generality **full**; record **faithful** (the same `dW` drives state and readout).

---

## 3. Higher-order structure-preserving CP schemes (the schemes on the baseline's own turf)

This is the family that tries to beat Rouchon-Ralph while keeping the Kraus/positivity
structure.

**Rouchon 2022 tutorial** (arXiv:2208.07416, *Annual Reviews in Control*). The
authoritative modern statement of the Kraus-map viewpoint: the discrete-time SME is a
CPTP map `rho -> sum_mu M_mu rho M_mu^dag` per step, and the continuous-time diffusive
and jump SMEs are recovered as its infinitesimal limit. The numerical schemes it
presents are the **first-order linear** ones, i.e. the Rouchon-Ralph family. This is the
correct conceptual home of the baseline and confirms its first-order status.

**Guevara-Wiseman 2020** (Phys. Rev. A **102**, 052217, arXiv:1909.12455). Motivated by
the observation that accumulated numerical error in a trajectory shows up as a deviation
from valid (CP) quantum evolution. They give a CP map that satisfies the "valid average
evolution" conditions to `O(dt^2)` but reproduces the Lindblad equation only to
`O(dt)`. A partial second-order result: better-behaved trajectories, same average order
as the baseline.

**Wonglakhon-Wiseman-Chantasri 2024 (the "W-map")** (Phys. Rev. A **110**, 062207,
arXiv:2408.14105). The genuine second-order CP scheme for the *conditioned* diffusive
SME. Built from a fourth-order expansion of the system-bath interaction unitary; the
resulting Kraus operator is CP **by construction, for any `dt`**, and reproduces the
average (Lindblad) evolution to `O(dt^2)`, one order better than both Rouchon-Ralph and
Guevara-Wiseman. Tested on dispersive (`z`-measurement) and dissipative (fluorescence
homodyne) qubits, it attains the smallest average trace distance to the exact
trajectories. The crucial caveat, stated plainly in the paper: there is **no improvement
in the per-trajectory (strong) scaling** with `dt` for *any* of these maps. The W-map
buys weak/average order 2 and keeps unconditional positivity; it does not make a single
trajectory converge faster, and it costs more per step.

**Cao-Lu 2021 / 2025** (arXiv:2103.01194, *J. Sci. Comput.*): "structure-preserving
numerical schemes for Lindblad equations", of **arbitrarily high order**, CPTP via a
Kraus representation with normalization. This is the result behind the folklore that
"arbitrary-order positivity-preserving schemes exist in finite dimension." It is
essential to state its scope precisely: it is for the **deterministic** Lindblad master
equation (the unconditional, ensemble-averaged density matrix), **not** the conditioned
stochastic equation. The related "Kraus is king" high-order CPTP low-rank method
(*J. Comput. Phys.* 2025) is likewise deterministic. So arbitrary order + unconditional
CP is a solved problem for the *average* dynamics, and an unsolved one for the
*conditioned trajectory*, where second order (the W-map) is the current ceiling.

**Amini-Mirrahimi-Rouchon** (arXiv:1407.7810, "Scheme 5") is the historical
positivity-preserving predecessor; it preserves positivity but is reported inaccurate and
slow in high dimension, superseded by the simpler Rouchon-Ralph linear scheme.

**Opinion.** For the conditioned SME, the reachable "better than Rouchon-Ralph while
staying unconditionally positive" is exactly one order of *weak* accuracy, delivered by
the W-map, at a heavier step and with no strong-order gain. If what you actually want is
an accurate *mean*, section 5 gets you weak order 2 far more cheaply. If you want the
accurate *conditioned state itself* with a positivity guarantee, the W-map is the tool,
and it is worth implementing only when the second-order state (not just its average)
matters. For everything else, the baseline is not improved on this axis.

---

## 4. The classical SDE schemes (the foils, and why they lose here)

The paper's own foil is Euler-Milstein. More generally: strong Taylor 1.5, Roessler
SRK (SRA/SRK), Platen, and higher Milstein raise the strong and/or weak order of a
generic SDE. They are the wrong tool for this SME for two structural reasons:

- **They do not preserve positivity.** The update is `rho + drift dt + diffusion dW`,
  not a Kraus sum, so nothing keeps the eigenvalues non-negative. The companion
  prototype makes this quantitative: an Euler-Maruyama step of the nonlinear SME on a
  qubit reaches a minimum eigenvalue of `-0.10` at `dt = 0.02` and `-0.52` at
  `dt = 0.5`, while the Rouchon step sits at machine zero (`~ -1e-16`) throughout. A
  filter that reports a "density matrix" with negative eigenvalues has left the state
  space; renormalizing the trace does not fix it.
- **Multiplicative, non-commuting noise caps the cheap order.** Strong order beyond
  `1/2` for several non-commuting `L_r` requires simulating iterated stochastic
  (Levy-area) integrals, which is expensive and rarely worth it. The Rouchon-Ralph
  Milstein-type term `1/2 sum_{r,s} sqrt(eta_r eta_s) L_r L_s (dy_r dy_s - delta_rs dt)`
  already captures the leading noise correction that makes the *mean* correct, without
  paying for Levy areas.

**Opinion.** Reaching for a higher-order generic SDE integrator trades away the property
that makes SME integration hard in the first place. Use them only if you have separately
guaranteed the state stays physical (small `dt` plus a projection), which usually costs
more than it saves.

---

## 5. Pure-state unravelings, and the recommended win

When the conditional state is **pure** (perfect detection `eta = 1`, every dissipation
channel monitored, pure initial state), propagating a `d x d` density matrix is
wasteful: the state is a `d`-vector. This is the single highest-value practical change
for large-`d` monitored systems, and it is what "efficient quantum filtering" in the
paper's title refers to.

The key point, verified symbolically and numerically in the prototype: the Rouchon step
restricted to a pure state **is** a state-vector map. With `rho = |psi><psi|`,

```
M rho M^dag / Tr[M rho M^dag]  =  |phi><phi| ,   phi = M psi / ||M psi|| ,
```

so `psi -> M psi / ||M psi||` produces exactly the same conditioned trajectory (a
rank-one `rho` at every step) and the same readout `dy = <psi|(L+L^dag)|psi> dt + dW`.
This is not a new discretization; it is the baseline scheme on the pure-state manifold
(Rouchon's papers give both forms). The gain is entirely computational:

- **Memory**: `O(d)` instead of `O(d^2)`.
- **Arithmetic**: the step never forms `M`; it applies the pieces to the ket,
  `M psi = Heff psi + (S psi) + 1/2 S(S psi) - Lcorr psi`, all matrix-vector products
  (`O(d^2)`), versus the density matrix's `M rho M^dag` (`O(d^3)`).

Beyond perfect detection there are two well-developed pure-state routes, neither a free
lunch:

- **Quantum state diffusion** (Gisin-Percival 1992) and the homodyne stochastic
  Schrodinger equation (Wiseman-Milburn) are the diffusive pure-state unravelings; the
  Rouchon pure-state map above is their positivity-trivial (norm-preserving) integrator.
  The numerical-analysis backbone is the norm-preserving **exponential scheme** of
  Mora (arXiv:math/0508490, weak order 1) and its exponential-integrator development
  (Phys. Scr., doi 10.1088/1402-4896/acd5b2): weak order 1, strong order 1/2, and (this is
  what row 2 of the verdict uses) **Talay-Tubaro extrapolation upgrades the mean to weak
  order 2 at ~1.5x cost**.
- **Mixed conditional states as a mixture of pure states** (stochastic interacting wave
  functions, *J. Comput. Phys.* 2018, arXiv:1709.08223): when `eta < 1` or `V != None`
  the conditional state is genuinely mixed, but it can still be written as an ensemble of
  interacting pure states each obeying an SSE. Cheaper in memory than the density matrix
  at large `d`, at the price of more trajectories. Situational: worth it when `d` is
  large and the mixedness is low rank.

**Opinion.** At `eta = 1` with full monitoring, do not propagate `rho`; propagate
`psi`. The recommendation is unconditional in that regime and the benefit grows with
`d`. When the state is mixed, stay on the density matrix at small `d`; consider the
stochastic-interacting-wavefunction route only when memory, not accuracy, is the binding
constraint at large `d`.

---

## 6. The prototype: measured results

`sde-schemes-prototype.wl` implements the pure-state map `\[ScriptCapitalP]` (matrix-vector
form, `O(d^2)`/step) and an Euler-Maruyama foil `\[ScriptCapitalE]`.
`prototype-verify-benchmark.wls` checks and benchmarks them against the baseline
`\[ScriptCapitalT]`; the following are from one run (all sections pass).

- **Equivalence (`eta = 1`)**: with a shared seed, `|psi><psi|` from `\[ScriptCapitalP]`
  matches the density-matrix trajectory from `\[ScriptCapitalT]` to `~1e-14`, and the
  Sow-ed readout to `~1e-16`, on a qubit, a spin-1 qutrit, and two-operator
  **non-commuting (`nL = 2`)** cases that exercise the double sum
  `1/2 sum_{r,s} dy_r dy_s L_r L_s` with `L_r L_s != L_s L_r`. The identity itself is
  proved symbolically (`2x2`, generic `M`) and spot-checked at `d = 5`. Same physics,
  same record; the numeric difference is floating-point associativity, not the scheme.
- **Positivity**: `\[ScriptCapitalP]` holds `||psi|| = 1` to `~5e-16` even at
  `dt = 5`, so the state is positive and unit-trace by construction at any step.
- **Lindblad mean**: the ensemble mean of `|psi><psi|` (K=800) matches the independent
  vectorized-Liouvillian reference to `0.0121` (below the Monte-Carlo floor
  `1/sqrt(K) = 0.035`) and reproduces the analytic amplitude-damping law
  `rho00(t) = e^{-gamma t}` to `0.005`. Weak order 1 of the mean, the defining property of
  the Rouchon-Ralph CP map for any monitored `L`, is established without Monte-Carlo noise
  by the multilevel coupled differences below and by the Mora / acd5b2 theorem. (A
  closed-form symbolic average of the *normalized* nonlinear one-step map is subtle: the
  naive `dt->eps^2`, `dW->eps g` moment expansion mis-handles the `1/Tr` normalization for
  some operators, so it is not claimed as a proof.)
- **Positivity trade (the foil)**: Euler-Maruyama reaches minimum eigenvalue `-0.10`,
  `-0.35`, `-0.52` at `dt = 0.02, 0.2, 0.5`; Rouchon stays at `~ -1e-16`.
- **Cost across dimension** (spin-`j`, `H = Jz`, `L = sqrt(0.3) J-`, 1000 steps, min over
  5 with `ClearSystemCache`), capped at `d = 32` where the timings are stable:

  | `d` | dense `\[ScriptCapitalT]` (s) | pure form-M (s) | pure matvec (s) | form-M speedup | matvec speedup |
  |---|---|---|---|---|---|
  | 2  | 0.0122 | 0.0086 | 0.0133 | 1.41x | 0.91x |
  | 4  | 0.0127 | 0.0092 | 0.0139 | 1.38x | 0.91x |
  | 8  | 0.0140 | 0.0099 | 0.0160 | 1.40x | 0.88x |
  | 16 | 0.0203 | 0.0127 | 0.0231 | 1.60x | 0.88x |
  | 32 | 0.081  | 0.0245 | 0.0483 | 3.32x | 1.68x |

  (Absolute timings are machine- and run-dependent; the speedup ratios and their growth
  with `d` are the invariant claim.) The density-matrix step is ~three `O(d^3)` matmuls
  (`s.s`, then `M rho`, then `. M^dag`). Forming `M` once and applying it to a ket
  ("form-M") removes two of them and wins at **every** `d`, from 1.4x at qubit scale to
  3.3x at `d = 32` and still climbing. The strictly `O(d^2)` matvec form removes the last
  matmul but overtakes form-M only past `d ~ 16-32`: at small `d`, WL evaluator dispatch,
  not arithmetic, dominates (the baseline `AUDIT.md` notes the same), so its extra
  per-step dispatches make it slower there. Beyond `d = 32` the dense trajectory is also
  `O(d^2)` in memory (a `d = 64` trajectory is ~65 MB vs the ket's ~1 MB), so its timing
  turns memory-bound and noisy while the pure-state maps stay light; the advantage only
  grows. Practically: use the form-M pure map (`\[ScriptCapitalF]`) for qubit/qudit-scale
  filtering; the matvec `\[ScriptCapitalP]` is the right asymptotic choice at large `d`.
  Both pure forms use `O(d)` memory.
- **Weak order and cheap extrapolation**: with common Brownian paths (a multilevel
  coupling, so the estimators are bias- not Monte-Carlo-limited), the base Rouchon mean
  is confirmed weak order 1 (coupled differences halve: ratios 2.04, 2.36), and
  Talay-Tubaro extrapolation removes the leading `O(dt)` error (the Richardson coupled
  difference is 2% of the base's). The residual `O(dt^2)` is below the Monte-Carlo floor
  at K=4000, consistent with the weak-order-2 theorem (acd5b2 Cor. 3.1) without
  over-claiming a numerically resolved slope.

---

## 7. Reconstructing the measurement record

The homodyne/photocurrent increment is `dy_r = sqrt(eta_r) <L_r + L_r^dag> dt + dW_r`,
with `dW_r` the innovation that also drives the state update. The baseline Sows exactly
this `dy` with the same `dW` that enters `M`, so the simulated record is already faithful:
state and current share one noise realization, which is the whole content of "the record
is the innovation." The `<L + L^dag>` is evaluated on the pre-step (Ito) conditional
state, the correct convention for the Ito record; a midpoint/Stratonovich average would
describe the same physics in a different calculus, not a more faithful record.

Two refinements exist but do not change the baseline's correctness:

- The **W-map** produces a slightly more accurate conditioned state at `O(dt^2)`, hence
  marginally better per-trajectory record statistics, at higher per-step cost.
- **Exact signal correlators** (Tilloy, Phys. Rev. A 2018) give analytic homodyne-current
  correlation functions without sampling trajectories at all, when correlators (not
  sample paths) are the goal.

For simulation and reconstruction of the record itself, the baseline is right and the
pure-state map inherits that correctness (it Sows the identical `dy`).

---

## 8. WL-native tooling and its limits

`ItoProcess` / `RandomFunction` and stochastic `NDSolve` are the wrong backbone for these
SMEs, confirmed by prior work in this course: they drift the state norm, stiffen on the
nonlinear measurement drift `(x - <x>)^2` (capping the usable step), and time out at
modest Hilbert dimension. A hand-rolled functional step (`FoldList` over a captured pure
`Function`, matrices/vectors) is both faster and structure-preserving. The one WL-native
success is the *linear-Gaussian* moment filter (section 1, row 4): there the state is a
handful of moments obeying a closed SDE with deterministic Riccati covariances, which
`ItoProcess` handles well and which is `O(1)` in the Hilbert dimension. Outside that
special structure, hand-rolled wins.

---

## References

- P. Rouchon, J. F. Ralph, *Efficient Quantum Filtering for Quantum Feedback Control*, Phys. Rev. A **91**, 012118 (2015); arXiv:1410.5345. (baseline)
- P. Rouchon, *A tutorial introduction to quantum stochastic master equations based on the qubit/photon system*, Annual Reviews in Control (2022); arXiv:2208.07416.
- I. Guevara, H. M. Wiseman, *Completely positive quantum trajectories with applications to quantum state smoothing*, Phys. Rev. A **102**, 052217 (2020); arXiv:1909.12455.
- N. Wonglakhon, H. M. Wiseman, A. Chantasri, *Completely positive trace-preserving maps for higher-order unraveling of Lindblad master equations*, Phys. Rev. A **110**, 062207 (2024); arXiv:2408.14105. (the W-map: unconditionally CP, weak order 2)
- Y. Cao, J. Lu, *Structure-preserving numerical schemes for Lindblad equations*, J. Sci. Comput. (2025); arXiv:2103.01194. (arbitrary-order CPTP, deterministic Lindblad only)
- H. Amini, M. Mirrahimi, P. Rouchon, *Models and feedback stabilization of open quantum systems*, arXiv:1407.7810. (Scheme 5, historical positivity-preserving predecessor)
- C. M. Mora, *Numerical solution of conservative finite-dimensional stochastic Schrodinger equations*, Ann. Appl. Probab. **15**, 2144 (2005); arXiv:math/0508490. (norm-preserving SSE, weak order 1)
- *On the rate of convergence of an exponential scheme for the non-linear stochastic Schrodinger equation with finite-dimensional state space*, Phys. Scr. (2023), doi 10.1088/1402-4896/acd5b2. (exponential norm-preserving scheme: weak order 1, strong order 1/2; Talay-Tubaro to weak order 2, Cor. 3.1)
- Stochastic interacting wave functions, *Numerical solution of stochastic master equations using stochastic interacting wave functions*, J. Comput. Phys. (2018); arXiv:1709.08223. (mixed-state unraveling as a mixture of pure SSEs)
- N. Gisin, I. C. Percival, *The quantum-state diffusion model applied to open systems*, J. Phys. A **25**, 5677 (1992). (QSD)
- J. Dalibard, Y. Castin, K. Molmer, *Wave-function approach to dissipative processes in quantum optics*, Phys. Rev. Lett. **68**, 580 (1992). (MCWF / quantum jump, for counting records)
- A. C. Doherty, K. Jacobs, *Feedback control of quantum systems using continuous state estimation*, Phys. Rev. A **60**, 2700 (1999). (Gaussian / Kalman-Bucy filter)
- H. M. Wiseman, G. J. Milburn, *Quantum Measurement and Control*, Cambridge University Press (2010).
- K. Jacobs, *Quantum Measurement Theory and its Applications*, Cambridge University Press (2014).
- A. Tilloy, *Exact signal correlators in continuous quantum measurements*, Phys. Rev. A **98**, 010104(R) (2018); arXiv:1712.05725. (analytic record correlators)
