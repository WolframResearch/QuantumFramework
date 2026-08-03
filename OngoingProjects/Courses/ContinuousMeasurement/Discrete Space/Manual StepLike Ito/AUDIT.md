# Audit: making `𝒯` (Rouchon SME integrator) QuantumFramework-independent, general, and fast

Scope: the `Manual evolution` section of `manual-implementation-ito-V2-community.nb`,
which implements the positivity-preserving stochastic-master-equation integrator of

> P. Rouchon and J. F. Ralph, *Efficient Quantum Filtering for Quantum Feedback Control*,
> Phys. Rev. A **91**, 012118 (2015), arXiv:1410.5345, Eqs. (3)-(4).

Deliverable: `manual-implementation-ito-qf-independent.wl` (new `\[ScriptCapitalT]` / `\[ScriptCapitalR]`, plain-matrix API).

## The scheme

For monitored operators `L_r` with efficiencies `eta_r` and unmonitored operators `V_j`,

```
rho(n+1) = num / Tr[num] ,
num      = M rho M^dag  +  dt Sum_j V_j rho V_j^dag  +  dt Sum_r (1-eta_r) L_r rho L_r^dag ,
M        = I - (i H + 1/2 Sum_j V_j^dag V_j + 1/2 Sum_r L_r^dag L_r) dt
             + Sum_r Sqrt[eta_r] L_r dy_r
             + 1/2 Sum_{r,s} Sqrt[eta_r eta_s] L_r L_s (dy_r dy_s - delta_rs dt) ,
dy_r     = Sqrt[eta_r] Tr[(L_r + L_r^dag) rho] dt + dW_r .
```

`num` is a sum of terms `A rho A^dag`, so the conditioned state stays positive
semidefinite for **any** step size. The readouts `dy_r` are `Sow`-ed each step.

Writing `S = Sum_r dy_r Sqrt[eta_r] L_r`, the double sum collapses to
`1/2 S^2 - 1/2 dt Sum_r eta_r L_r^2`, i.e. `M = Heff + S + 1/2 S^2 - Lcorr` with
`Heff = I - (iH + 1/2 Sum V^dag V + 1/2 Sum L^dag L) dt` and the Itô correction
`Lcorr = 1/2 dt Sum_r eta_r L_r^2`. That is exactly the new implementation.

## Correctness findings (all confirmed symbolically in a kernel)

The notebook computes `Heff = I - dt(iH + 1/2 Sum(L^2 + L^dag L) + 1/2 Sum(V^2 + V^dag V))`
and `M = Heff + Upsilon + 1/2 Upsilon^2`, `Upsilon = Sum_r Sqrt[eta_r] L_r dy_r`.
Comparing term by term against `M` above (single operator, general `eta`):

| # | Finding | Where it bites | Where it is invisible |
|---|---|---|---|
| 1 | **Dimension hardcoded**: `IdentityMatrix[2]`. | Any `d > 2`, *including the paper's own two-qubit (d=4) example*: `Heff` fails to form (`Thread::tdlen`) and the call crashes. | Every notebook example (all single-qubit, d=2). |
| 2 | **Itô correction missing the efficiency factor**: `M_notebook - M_paper = -1/2 (1-eta) dt L^2`. | `eta < 1` **and** `L^2 not in {scalar, 0}`. Biases the ensemble away from the (eta-independent) Lindblad average. | `eta = 1` (identical to the paper); or `L` a Pauli (`L^2 = I`, normalized away) or `L = J-` (`(J-)^2 = 0`). Every notebook example is one of these, so it never showed. |
| 3 | **Spurious `V^2` term** in `Heff`: `-1/2 dt V^2` that the paper's `M` has nowhere (V enters only as `V^dag V` and `V rho V^dag dt`). | `V != None` **and** `V^2 not in {scalar, 0}`. | `V = None` (every notebook example); or involutory / nilpotent `V`. |

The new implementation follows the paper exactly, so it is **identical to the
notebook on every case the notebook runs correctly** (eta=1, V=None) and additionally
correct for inefficient detection and unmonitored channels. Findings 2 and 3 are
latent precisely because a Pauli measurement operator squares to the identity and a
qubit `J-` squares to zero: the erroneous term is either normalized away or vanishes.
They reappear the moment `L` or `V` is a higher-spin ladder operator (e.g. spin-1 `J-`,
`(J-)^2 != 0`) or lives in `d > 2`.

## QuantumFramework independence

The notebook `𝒯` took QF objects and pulled matrices out with
`#["Computational"]["Matrix"]` / `["DensityMatrix"]`. The new `\[ScriptCapitalT]`
takes ordinary dense complex matrices directly and calls nothing from QF. Extract the
matrices once from whatever source you like and pass them in:

```wl
Get["manual-implementation-ito-qf-independent.wl"];
X = PauliMatrix[1]; Z = PauliMatrix[3];
traj = \[ScriptCapitalT][{{1.,0.},{0.,0.}}, 1.0 X, {Sqrt[2.0] Z}, 0.005, 20.];
{traj, {records}} = Reap @ \[ScriptCapitalT][rho0, H, Ls, dt, tf];   (* readouts via Reap *)
```

The call signature and the `Sow`-ed readouts are unchanged, so example code ports by
swapping QF objects for their matrices. The pattern guards (`_?MatrixQ`) make the
function fire only on matrix input, so a stray QF-object call stays unevaluated rather
than silently misbehaving.

## Performance

The hot loop is `FoldList` over `Floor[tf/dt]` steps; the notebook's two-point example
is 400,000 steps. Levers applied:

- **Everything step-invariant is precomputed once**: `Heff`, `Sqrt[eta_r] L_r`,
  `L_r + L_r^dag`, the Itô correction, and the dissipation channels. The notebook
  rebuilt `L_r + L_r^dag`, `Sqrt[eta] L`, and re-summed `L^dag L` inside every step.
- **The step is a captured pure `Function`**, not a pattern-matched downvalue, so no
  pattern matching happens 400,000 times.
- **The inefficiency / V dissipators are dropped entirely when absent** (`extra = 0&`),
  instead of the notebook's `dt Sum (1-eta_r) L rho L^dag` that multiplied full matrix
  products by zero every step when `eta = 1`.
- `S = dy . (Sqrt[eta] L)` and the numerator are single contractions on packed arrays.

Benchmark (pure WL both sides, `ClearSystemCache[]` before each, min over repeats):

| case | notebook algorithm | new `𝒯` | speedup |
|---|---|---|---|
| d=2, 1 op, 4,000 steps | 0.079 s | 0.045 s | 1.74x |
| d=2, 1 op, 40,000 steps | 0.971 s | 0.459 s | 2.11x |
| d=4, 2 ops, 5,000 steps, eta=0.85 | 0.130 s | 0.097 s | 1.34x |

The notebook's actual two-point run (d=2, 400,000 steps) with the new `𝒯`: **4.6 s**.

The win is ~1.7-2.1x when `eta = 1` (the inefficiency dissipator is dropped
entirely) and ~1.3x at `eta < 1` (both versions must then apply `L rho L^dag`).
The remaining cost is **WL evaluator overhead**, not arithmetic: 400,000 steps is
~1.6 million tiny (2x2) matrix products, each dominated by dispatch rather than
flops. The next lever, if a run ever needs it, is a `Compile`d hot loop specialized
to a fixed `(d, nL)` (roughly another order of magnitude), at the cost of the
generality and readability of the current functional version. It is deliberately
not taken here.

## Test results (QF-free harness, `tests.wls`)

- **Independent Lindblad reference** validated against the analytic amplitude-damping
  law `rho00(t)=(1/2)e^{-k t}`, `rho01(t)=(1/2)e^{-k t/2}` to ~1e-17.
- **Exact reproduction** of the notebook algorithm on all five examples: trajectory
  and readouts match to `0.` (bit-identical) on the `J-` examples and ~1e-14 on the
  `Z` examples (only `Total`-vs-`Sum` round-off).
- **Positivity + trace + Hermiticity**: `min eig >= -1e-15`, `|Tr - 1| ~ 2e-16`, and
  `|rho - rho^dag| ~ 1e-15` on every trajectory, still holding at an intentionally
  huge `dt = 0.5` (structure preservation).
- **Ensemble mean (K=800) -> Lindblad**: `max|diff| = 0.0121`, below the `1/Sqrt[K] ~ 0.035` Monte-Carlo scale.
- **Non-commuting monitored operators** (two operators with `L1 L2 != L2 L1`, d=2):
  reproduces the notebook algorithm to round-off and the ensemble mean matches Lindblad;
  the double-sum's `L_r L_s` ordering is also verified symbolically for two
  non-commuting operators with general efficiencies.
- **Public `\[ScriptCapitalR]`** round-trips the integrator's first sown readout (the
  standalone readout and the inlined one cannot drift).
- **Convergence**: the ensemble mean tracks Lindblad at both `dt = 0.04` and `0.01`
  (refinement does not degrade it).
- **Discriminators** (where the notebook's latent bugs would show; ensembles run on a
  common random stream so `|notebook - new|` isolates the bug free of Monte-Carlo noise):
  - d=4 two-qubit (paper's own example): notebook `IdentityMatrix[2]` cannot form
    `Heff`; new `𝒯` runs, stays physical, and its 300-trajectory mean matches the d=4
    Lindblad equation to 0.0075.
  - eta<1 qutrit (spin-1 `J-`, `eta=0.2`, `L^2 != 0`): new mean vs Lindblad `= 0.0014`,
    notebook (eta-unweighted) vs Lindblad `= 0.0355` (**25x** farther); MC-free
    systematic `|notebook - new| = 0.036`.
  - V!=None qutrit (spin-1 `J-` unmonitored): new mean vs Lindblad `= 0.0022`, notebook
    (spurious `V^2`) vs Lindblad `= 0.0188` (**8.6x** farther); systematic `= 0.019`.

## Notes / caveats

- The middle-optional signature `\[ScriptCapitalT][rho,H,L,eta,V,dt,tf]` is kept for
  drop-in parity with the notebook. Optionals between required args is a fragile
  pattern, so every slot is guarded: `rho`/`H` by `?MatrixQ`, `dt`/`tf` by `?NumericQ`,
  and `eta`/`V` by `(None | _List)`. A malformed call (a scalar where a list belongs,
  a non-number for `dt`) therefore stays unevaluated rather than mis-binding into a
  plausible-looking wrong trajectory. The 5-/6-/7-argument shapes and the fail-fast on
  each guarded slot are checked in `tests.wls` Section 9. A re-author using named
  options (`Efficiencies ->`, `Unmonitored ->`) would drop the fragile pattern entirely.
- Structure vs accuracy (they are different, and large `dt` separates them). The
  update is a sum of `A rho A^dag` terms over a trace normalization, so the state
  stays a valid density matrix (positive, Hermitian, trace 1) for **any** `dt`,
  tested up to `dt = 100`. Accuracy is a separate matter: the scheme is first-order
  (Rouchon-Ralph), so the ensemble mean carries an `O(dt)` bias. Small `dt` is
  accurate (bias below the Monte-Carlo floor at `dt ~ 0.02`); large `dt` is not (the
  bias grows to ~3x the MC floor by `dt = 1`). So a large step keeps `rho` physical
  but not correct in the mean: refine `dt` for accuracy, and do not read positivity
  at large `dt` as a sign the trajectory is right.
- Validation is at the top of the call (matrix guards), not in the loop, by design.
