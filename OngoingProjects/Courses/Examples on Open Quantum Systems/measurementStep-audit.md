# Audit: `measurementStep` and its three predicates

Target: `open-system-simulation-catalog.md`, section "The positivity-preserving measurement step": the `measurementStep` function and its three predicate helpers `realScalarQ`, `hermitianMatrixQ`, `densityOperatorQ`. This audit made no changes to the catalog. Evidence below is from `wolframscript` runs reproduced in this session.

Note on line numbers: the catalog was under concurrent edit while this audit ran (its content shifted 7-14 lines and the `steadyState` `NonUniqueSteadyState` branch was removed mid-audit), so all catalog line numbers here are approximate, as of the audit; anchor on the symbol names. As of the last check the predicates are at ~237-243, `measurementStep` at ~245-291, and the five call sites at ~297 (`fixedStep`), ~324 (`trajectory`), ~366 (`oneStepBias`), ~2390 (`steer`), ~2602 (`steerStep`). The `measurementStep` code itself matched verbatim through every check.

## 1. Verdict on the physics: faithful, bit-exact to Rouchon-Ralph

The implementation is a correct, term-by-term realization of the filter stated in the prose (arXiv 1410.5345). Every correspondence the essay claims holds:

| Formula (prose) | Code | Line (~) |
|---|---|---|
| $K = i\hat H + \tfrac12\sum_{\text{all}}\hat c^\dagger\hat c$ | `drift = I H + Fold[Plus, 0 H, ConjugateTranspose[#].# &/@ channels]/2`, `channels = Join[watched, unwatched]` | 255, 269 |
| $\tfrac{dt}{2}\sum_k\eta_k\hat c_k^2$ (watched only) | `corr = Fold[Plus, 0 H, MapThread[#1 (#2.#2)&, {realEffs, watched}]]/2`, applied as `corr stepSize` | 270, 285 |
| $s = \sum_k\sqrt{\eta_k}\,dJ_k\,\hat c_k$ | `sig = Fold[Plus, 0 H, MapThread[Sqrt[#2] Re[#1] #3 &, {dJ, realEffs, watched}]]` | 284 |
| $M = \mathbb 1 - K\,dt + s + \tfrac12 s^2 - \tfrac{dt}{2}\sum\eta\hat c^2$ | `M = id - drift stepSize + sig + sig.sig/2 - corr stepSize` | 285 |
| $\tilde\rho = M\rho M^\dagger + dt\sum_j\hat l_j\rho\hat l_j^\dagger + dt\sum_k(1-\eta_k)\hat c_k\rho\hat c_k^\dagger$ | `top = M.rho.ConjugateTranspose[M] + stepSize (…unwatched…) + stepSize (…(1-eff) watched…)` | 286-288 |
| $\rho' = \tilde\rho/\mathrm{Tr}\,\tilde\rho$ | `top/normalization`, `normalization = Re@Tr[top]` | 289-291 |

The `0 H` seed passed to every `Fold[Plus, …]` is an additive-zero matrix of the right `{d,d}` shape, so empty `unwatched` (or empty `watched`) contributes an exact zero matrix rather than the scalar `0`. That is the reason the three sums stay matrix-valued at the boundary cases the essay actually calls (`unwatched = {}` at every site).

Evidence:

- **Bit-exact against an independent reimplementation of the prose formula.** A separate `refStep` built straight from the equations, run on 20 random cases with two watched channels at $\eta<1$, one unmonitored leak, random Hermitian operators and a random density matrix: `max |essay - reference| = 0.` A cold-kernel reproduction strengthened this with general **non-Hermitian** channel operators (which would expose any `ConjugateTranspose`-vs-`Transpose` slip, invisible when $\hat c^\dagger=\hat c$) and still got `0.` The deviation is exactly zero when the reference builds the sums with the same structure as the essay; any residual comes from summing the terms in a different order or association, which at three or more machine-precision summands perturbs the result at $\sim10^{-16}$ (floating-point non-associativity), not a formula error. The physics fidelity is solid either way.
- **The two asserted facts reproduce.** Naive update: `Min@Eigenvalues = -0.140`, a negative probability, off the state space. Positivity step `measurementStep[0 id2,{Z},{1.},{},0.1][densityMatrix[plus],{0.5}]`: `{min eig, trace} = {-3.7e-17, 1.0}`, i.e. on the boundary, positive to numerical precision, trace one.
- **The Ito correction is right.** Averaging the normalized single step over the Gaussian record and comparing to $\rho_0 + \mathcal L\rho_0\,dt$ gives a Frobenius gap of $3.2\times10^{-5}$ at $dt=10^{-3}$, consistent with Monte-Carlo scatter plus $O(dt^2)$ bias. The $-\tfrac{dt}{2}\sum\eta\hat c^2$ term is what cancels the deterministic part of $\tfrac12 s^2$ so the ensemble mean reproduces the Lindblad generator instead of double-counting the measurement back-action.
- **The upper efficiency bound is the positivity hypothesis; the lower bound is a separate physicality requirement.** Forcing $\eta=1.5$ drives the numerator's `min eig = -0.038 < 0`: the $(1-\eta)\hat c\rho\hat c^\dagger$ term flips sign and the Kraus-form positivity argument collapses. A scan confirms the numerator stays PSD for every $\eta\le1$, touches zero at $\eta=1$, and goes negative only for $\eta>1$. The lower bound is **not** a positivity boundary: at $\eta=-0.5$ the numerator's `min eig = +0.127 > 0` (still PSD), because $M\rho M^\dagger$ is positive-semidefinite for *any* $M$ and $(1-\eta)>0$ when $\eta<0$. What $\eta<0$ breaks is the *reality* of the measurement kick $s=\sqrt{\eta}\,dJ\,\hat c$ (an imaginary $\sqrt{\eta}$ is not a physical homodyne back-action), which is a modeling requirement, not the positivity theorem. So positivity holds for all $\eta\le1$; the full range $0\le\eta_k\le1$ is the physical range of a detection efficiency, and keeping the guard on both bounds is right, but only the upper bound is what the section's positivity claim rests on.

No correctness defect. The heavy error handling is the only thing in question, and it is a style/pedagogy question, not a physics one.

## 2. Two-axis grade of every check

Two facts frame the grades. First, **no call site inspects the result for `Failure`**: all five callers (`fixedStep`, `trajectory`, `oneStepBias`, `steer`, `steerStep`) pass known-valid arguments and use the return directly. A fired `Failure` would not teach anything; it would propagate as an opaque object into a `FoldList` or a plot. The objection is not to `Failure` as such: the essay uses it legitimately as a return value elsewhere (for example `regressionSpectrum`'s `DefectiveEigenbasis`). It is that these thirteen `measurementStep` tags are read by no caller and named by no prose in the section, so a fired one informs no one. Second, **the per-step revalidation is measurably expensive and provably redundant**: it multiplies the per-step cost by a factor of order **2 to 4** (length- and machine-dependent; 3.8x at n=8000 and 3.5-3.9x at n=1600 on this machine, ~2.1x at n=40000 on another cold-kernel run) for a **bit-identical** trajectory (deviation 0). `densityOperatorQ`'s per-step `Eigenvalues` is the **dominant** slice: on a leave-one-out at n=8000 (min over 8, `ClearSystemCache[]`), removing that one check alone closes about 90% of the full-to-lean gap, the other four per-step checks together the remaining ~10%.

| # | Tag | Stage | Physically load-bearing? | Noise for this file? | Verdict |
|---|---|---|---|---|---|
| 1 | InvalidChannelContainers | setup | No: a type guard, no physics | Yes: the signature already says lists | **DROP** |
| 2 | InvalidHamiltonian | setup | Weak: `d` is read from `H` | Mostly: a non-square `H` makes `Dot` error clearly and locally | **DROP** → precondition comment |
| 3 | NonHermitianHamiltonian | setup | Weak: non-Hermitian `H` is wrong physics but does **not** break positivity ($M\rho M^\dagger$ is PSD for any `M`) | Mostly: no site passes a non-Hermitian `H` | **DEMOTE** → comment |
| 4 | ChannelEfficiencyLengthMismatch | setup | Weak: encodes the watched↔effs pairing | Mostly: a mismatch makes `MapThread` error clearly | **DROP** → comment |
| 5 | ChannelDimensionMismatch | setup | Weak: conformability | Mostly: a wrong `{d,d}` makes `Dot` error clearly | **DROP** → comment |
| 6 | **InvalidEfficiency** ($0\le\eta\le1$) | setup | **Yes: $\eta\le1$ is the exact positivity hypothesis** (section 1); $\eta\ge0$ keeps the $\sqrt\eta$ kick physical | No | **KEEP** as one `ConfirmAssert` |
| 7 | InvalidTimeStep ($dt>0$) | setup | Moderate: a step must advance forward | Partly: no site passes $dt\le0$ | **KEEP** as one `ConfirmAssert` (the lighter of the two; droppable to a comment) |
| 8 | StateDimensionMismatch | per-step | No: the fold only ever feeds back a valid `{d,d}` state | Yes | **DROP** |
| 9 | **InvalidDensityOperator** | per-step | No: the update **provably preserves** positivity/trace, so re-checking it every step contradicts the section's own thesis | Yes, the dominant slice (~90%) of the 2-4x cost | **DROP** (strongest) |
| 10 | InvalidRecordContainer | per-step | No | Yes | **DROP** |
| 11 | RecordLengthMismatch | per-step | Weak: `MapThread` pairing | Yes | **DROP** |
| 12 | InvalidRecord | per-step | No | Yes | **DROP** |
| 13 | **NonPositiveNormalization** | per-step | Yes: guards a real divide-by-nonpositive | No: nearly free, `Tr[top]` is computed anyway | **KEEP** as a one-line guard |

Net: of 13 branches, **3 survive** (efficiency bound, positive $dt$, normalization guard), all others drop or demote to a comment. The three predicates `realScalarQ`, `hermitianMatrixQ`, `densityOperatorQ` have **no other consumer in the file** (verified by grep: their only uses are the checks above), so all three are deleted with the checks they served. The surviving efficiency bound needs no predicate: `VectorQ[effs, 0 <= # <= 1 &]` is self-contained.

## 3. Recommended leaner `measurementStep`

The preconditions live in one comment. Two `ConfirmAssert`s guard the two conditions a reader could plausibly violate while experimenting and that fail *silently* if violated (a non-positive "state", or time running backward). Dimension and Hermiticity mismatches are left to surface as ordinary WL errors, which in a live notebook point at the exact operation that failed. The stepper carries no per-step validation; its only runtime check is the division guard.

```wl
(* measurementStep: one positivity-preserving (Rouchon-Ralph, arXiv:1410.5345) step.
   Preconditions, enforced once at setup, not per step: H Hermitian and d*d; every
   channel operator d*d; 0 <= eta_k <= 1 (eta_k <= 1 is the bound under which the
   Kraus-form update stays positive; eta_k >= 0 keeps the sqrt(eta) kick physical);
   dt > 0. Inside the returned stepper the update provably preserves
   Hermiticity, positivity, and trace, so the state and the record are trusted. *)
measurementStep[H_, watched_, effs_, unwatched_, dt_] := Enclose[
  Module[{d = Length[H], id = IdentityMatrix[Length[H]],
    channels = Join[watched, unwatched], drift, corr},
   ConfirmAssert[VectorQ[effs, 0 <= # <= 1 &], "efficiencies must lie in [0,1]"];
   ConfirmAssert[TrueQ[dt > 0], "time step must be positive"];
   drift = I H + Fold[Plus, 0 H, ConjugateTranspose[#] . # & /@ channels]/2;
   corr = Fold[Plus, 0 H, MapThread[#1 (#2 . #2) &, {effs, watched}]]/2;
   Function[{rho, dJ}, Module[{sig, M, top, nrm},
     sig = Fold[Plus, 0 H, MapThread[Sqrt[#2] #1 #3 &, {dJ, effs, watched}]];
     M = id - drift dt + sig + sig . sig/2 - corr dt;
     top = M . rho . ConjugateTranspose[M] +
       dt Fold[Plus, 0 H, # . rho . ConjugateTranspose[#] & /@ unwatched] +
       dt Fold[Plus, 0 H, MapThread[(1 - #2) #1 . rho . ConjugateTranspose[#1] &, {watched, effs}]];
     nrm = Re@Tr[top];
     If[TrueQ[nrm > 0], top/nrm,
       Failure["NonPositiveNormalization", <|"Normalization" -> nrm|>]]]]]];
```

The three predicates `realScalarQ`, `hermitianMatrixQ`, `densityOperatorQ` (~237-243) are **deleted outright**. They exist only to feed the dropped checks and appear nowhere else in the catalog.

Rationale, per surviving element:

- **`ConfirmAssert[… 0 <= # <= 1 …]`** keeps the physical range of a detection efficiency; its upper end $\eta\le1$ is the precondition the section's positivity claim rests on, and an $\eta>1$ otherwise yields a silently non-positive matrix that still looks like a state. Guarding both ends is right (an efficiency outside $[0,1]$ is meaningless), even though the lower end is a physicality bound rather than a positivity one.
- **`ConfirmAssert[TrueQ[dt > 0]]`** rejects a non-advancing or backward step. This is the lighter of the two; a reader wanting the absolute minimum can move it into the comment and keep only the efficiency guard.
- **`If[TrueQ[nrm > 0], …, Failure[…]]`** guards the sole fragile operation, the division. With valid inputs `Tr[top]` is a sum of traces of PSD matrices, hence $\ge 0$, and $>0$ except in degenerate cases, so this essentially never fires; it costs nothing because `nrm` is needed for the division regardless.

The `Re[...]` coercions on `effs`, `dJ`, and `dt` in the original are dropped: every call passes real efficiencies, real homodyne records, and a real step, so the coercions were defensive and change no result (the bit-exact check below passes without them).

Dropped tags, and why the essay stays correct without each:

- **InvalidChannelContainers, InvalidHamiltonian, ChannelEfficiencyLengthMismatch, ChannelDimensionMismatch** — every call passes lists, a square numeric `H`, and matching, conformable `{d,d}` operators. A violation surfaces immediately as a `Dot`/`MapThread` error at the offending line, which is clearer in a live notebook than a wrapped `Failure`.
- **NonHermitianHamiltonian** — non-Hermitian `H` is wrong physics but not a positivity break; no site passes one. Named in the precondition comment.
- **StateDimensionMismatch, InvalidDensityOperator, InvalidRecordContainer, RecordLengthMismatch, InvalidRecord** — all per-step. The stepper is only ever driven by a `FoldList` (or a direct valid call) whose carried state is the previous output, which the update keeps Hermitian, unit-trace, and PSD by construction. Re-deriving that every step via `Eigenvalues` is the dominant slice (~90%) of the 2-4x tax for zero information; it also contradicts the section's thesis that the step *is* the thing that keeps $\rho$ a state.

## 4. What each call site sees

With valid arguments (all five), the setup returns the stepper exactly as before, so nothing downstream changes. The behavior on *bad* input shifts in a way worth stating:

- Bad **efficiency** or **non-positive dt**: `measurementStep[…]` now returns `Failure["ConfirmationFailed", …]` from the `Enclose` (verified: `FailureQ` true for $\eta=1.5$ and for $dt=-0.1$), where the original returned `Failure["InvalidEfficiency"|"InvalidTimeStep", …]`. Same shape (a `Failure`), now carrying the precondition text.
- Bad **dimensions / non-Hermitian / mismatched lengths**: instead of a structured `Failure`, the caller now sees an ordinary WL error message at the failing `Dot`/`MapThread`. For a notebook that is a more locatable failure, not a less informative one.
- Bad **per-step state or record**: no longer intercepted. This is the intended change: those inputs cannot arise from the essay's own folds, and intercepting them cost the 2-4x with no payoff.

Confirmations for this rewrite: bit-identical to the original over 50 random valid inputs (`max |orig - lean| = 0`); FACT 2 reproduces (`{-3.7e-17, 1.0}`); the `trajectory` `FoldList` shape yields a Hermitian, trace-one, PSD final state; the `steer` measure-then-rotate shape yields trace one and a PSD state.

## Summary

The measurement step is a faithful, bit-exact implementation of the Rouchon-Ralph filter, and its positivity guarantee rests on one supplied precondition, $\eta_k\le1$ (the upper end of the physical efficiency range $0\le\eta_k\le1$ the guard should keep), while the rest of the wall should shed. Of thirteen `Failure` branches and three predicates, three checks are load-bearing (the efficiency bound, a positive step, and the division guard) and everything else is defensive plumbing that no caller reads and no prose teaches; the per-step density-operator revalidation is not merely noise but a near-quadrupling of the trajectory cost to re-prove an invariant the update already guarantees. The leaner function above is the same physics in a third of the lines, still refusing the efficiency above one that would silently break positivity.
