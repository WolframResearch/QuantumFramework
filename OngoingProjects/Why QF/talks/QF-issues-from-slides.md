# QF issues found while revising the slide decks

Temporary list from the 10-slide pilot and the sweep over the remaining 129 slides of Modeling-Slides.md and Continuous-Measurement-Slides.md (2026-11-10, QF checkout 48be806e, WL 15.0.1).

Every entry marked "re-run" was reproduced again in a fresh kernel from its reproduction file, with the reported result. Reproduction files live in `QF-sweep/batchNN/checks/issue-*.wl` (run from the batch folder with `wolframscript -file checks/run.wls <file> < /dev/null`); my own checks are in `QF-sweep/mine/checks/`. Each entry was also searched for in QF-Defect-Ledger.md and BUGS.md; none is listed there unless stated.

## Likely bugs

1. **`"PauliDecompose"` is undefined on an operator built from a state.** `QuantumOperator[qs]["PauliDecompose"]` gives `QuantumOperator::undefprop` for a two-qubit density-matrix state; `QuantumOperator[qs["DensityMatrix"], {1, 2}]["PauliDecompose"]` works. Re-run (`mine/checks/pilot-issues.wl`). Source: pilot, CHSH slide.

2. **Substitution into a QF state silently does nothing.** On WL 15.0.1, `ReplaceAll` no longer enters a `SparseArray` (`SparseArray[{a, 0, b}] /. a -> 1` returns it unchanged), so `qs /. a -> 3/10` returns the same state with no message, for plain, evolved, tensor-product and partial-trace states alike. `qs[<|a -> 3/10|>]` also returns it unchanged when `a` was not declared with `"Parameters"`; only `QuantumState[qs, "Parameters" -> {a}][3/10]` substitutes. The audit gives the broken route as advice (`audit/mistakes.md`, "Property cache has no auto-invalidation": `qs2 = qs /. v -> v2`). Re-run (`mine/checks/own.wl`, `batch05/checks/issue-replaceall-*.wl`, `batch04/checks/issue-replaceall-state.wl`). Source: batches 4 and 5.

3. **`QuantumEvolve` returns a raw unevaluated `NDSolveValue` after its own error.** `QuantumEvolve[w QuantumOperator["X"], None, {t, 0, 1}]` (a free symbol, numeric times) prints `QuantumEvolve::error` and then returns `NDSolveValue[...]` with QF private symbols in it instead of a `Failure`; property calls on it then fail again. Re-run (`batch04/checks/issue-evolve-symbolic-numeric-range.wl`). Source: Gaussian-pulse slide.

4. **`QuantumEvolve[None, jumps -> rates, state, t]` stays unevaluated with no message**, while `"Liouvillian"[None, ...]` accepts `None`. Re-run (`batch01/checks/issue-evolve-none.wl`). Source: sudden-death slide.

5. **`QuantumLinearSolve` prints `StandardDeviation::shlen` twice when b has one nonzero entry.** The answer is right (9/8, 3/8 for the slide's system). It also prints its progress lines in a headless run. Re-run (`batch06/checks/issue-qls-stdev.wl`). Source: HHL slide.

6. **`QuantumNaturalGradientDescent` returns parameters as `x + 0.*I`.** The numeric metric is complex-typed, so every step after the first is `Complex`. Re-run (`batch06/checks/issue-qngd-complex.wl`). Source: VQE slide.

7. **`qo["Eigenvalues"]` and `qo["Eigensystem"]` print `Eigensystem::arh` on a machine-number operator of matrix size 256.** QF passes its stored `SparseArray` to a full eigen-solve; the values are right and a dense operator of the same size is silent. This is the machine-number version of QF-B46 (`Eigensystem::arhm` on exact sparse operators). Re-run (`batch11/checks/issue-eigenvalues-arh.wl`). Source: Young-diagram slide.

8. **`QuantumPartialTrace` of a pure state is far slower than the reduced matrix needs.** Keeping 9 qubits: 0.31 s at n = 10, 0.94 s at n = 11, 5.7 s at n = 12, against 0.011 s for a reshape and product with the same matrix (difference 1e-18). The cause is already in `audit/mistakes.md` (the full density matrix is built first); the timings are new. Batch 7 also saw the n = 14 call run about 100 s past a `TimeConstrained` limit of 120 s (not re-run). Re-run up to n = 12 (`mine/checks/own.wl`). Source: Ising slide.

9. **An unknown option to `BosonicNormalOrder` gives `OptionValue::nodef` and the call continues** instead of being rejected. Re-run (`batch09/checks/issue-fieldvariables.wl`). Source: beam-splitter slide.

## Open questions (values right; behavior may be intended)

10. **Symbolic rates come out as `Sqrt[γ] Conjugate[Sqrt[γ]]`** in `QuantumEvolve` and in named channels (amplitude damping gives Bloch z = `-(Sqrt[1 - γ] Conjugate[Sqrt[1 - γ]]) + Sqrt[γ] Conjugate[Sqrt[γ]]`). On the sudden-death slide this makes the concurrence a 48583-leaf expression, off by 1.6e-8 at machine precision though exact at 30 digits. QF does not assume parameters real; see the ParameterAssumptions work. Re-run (`batch01/checks/issue-evolve-conjugate.wl`, `batch04/checks/issue-symbolic-bloch-conjugate.wl`).

11. **Symbolic `"VonNeumannEntropy"` at an endpoint.** For `{{p, 0}, {0, 1 - p}}` it returns `((-1 + p) Log[1 - p] - p Log[p])/Log[2]` bits, which gives `Indeterminate` with `Infinity::indet` at p = 0; the numeric state at p = 0 gives 0. This is how WL evaluates p Log p at 0, so it is a limit-versus-substitution question, not a wrong formula. Re-run. Source: pilot, capacity slide.

12. **Symbolic `"ProbabilitiesList"` wraps each entry in `Abs` and divides by the total**: `{Abs[p]/(Abs[1 - p] + Abs[p]), ...}` for `{{p, 0}, {0, 1 - p}}`; `Simplify[..., 0 < p < 1]` recovers `{p, 1 - p}`. Probably a guard for numeric input; it blocks closed forms. Re-run (`batch01/checks/issue-probabilities-abs.wl`, `batch04/checks/issue-symbolic-probabilities.wl`).

13. **A Lindblad solution fails QF's own physicality test.** On 32 Fock levels the evolved state has an eigenvalue of -1.6e-8 from integration error, so `"PhysicalQ"` is False and `QuantumDistance` against it prints `::notphysical`. The ledger records strict machine precision as an accepted decision for the entanglement functions; this shows it also reaches `QuantumDistance` on QF's own output. Re-run (`batch02/checks/issue-evolve-psd.wl`, `issue-distance-notphysical.wl`).

14. **A measured circuit returns a `QuantumState` or a `QuantumMeasurement` depending on whether every record controls a later gate.** Re-run (`batch06/checks/issue-record-head.wl`). Source: teleportation slide.

15. **`"SymbolicMeasure"` stores a determined outcome as a nested `Boole[Mod[...]]` chain** that grows with every earlier measurement; the truth table is right (XOR). Re-run (`batch06/checks/issue-symbolic-outcome.wl`).

16. **Probabilities of a circuit that ends in measurements cost about 50x** the same circuit without them, because each measurement adds a pointer qubit (N = 21 Shor: 9.4 s against 0.18 s). This is the growth studied in `OngoingProjects/Measurement Dilation/`. Re-run (`batch07/checks/issue-qpe-measurement-cost.wl`; values agree to 6e-17).

17. **`QuantumInterferometer["Fourier"[m]]` is machine-precision**: exact zero probabilities come out near 1e-32, while `"BS"` stays exact. Re-run (`batch08/checks/issue-fourier-numeric.wl`).

18. **Exact zero probabilities of the 5-mode Fourier network do not simplify to 0**: `PossibleZeroQ` and `FullSimplify` give 0, but none of the 100 dark patterns is `=== 0`; m = 6 gives exact 0. Re-run (`batch09/checks/issue-fourier5-zero.wl`).

19. **Exact Clements settings keep unreduced `Root` objects**, and three phases keep `Arg[Root[...]]` after `FullSimplify`. Not re-run. Source: batch 9.

20. **The symbolic-photon-number probability has no range guard**: substituting an impossible pattern leaves `Abs[WignerD[...]]^2` unevaluated, or gives `Indeterminate` at numeric angles, where the numeric-pattern route returns `QuantumInterferometer::occ`. Re-run (`batch08/checks/issue-symbolic-range.wl`).

21. **No qudit T gate.** `QuantumOperator["T"[3]]` fails with `QuantumOperator::invalidArgs`, and `"P"[θ, 3]` puts the phase on the last level only, so the qutrit T = diag(1, ζ, ζ⁻¹) has to be typed as a matrix, though `"X"[d]`, `"Z"[d]` and `"H"[d]` exist. Re-run (`batch10/checks/issue-qudit-T.wl`).

22. **Bad input is accepted or slow to fail.** `QuantumState` takes strings as amplitudes with no message (re-run, `batch05/checks/issue-string-amplitudes.wl`). `QuantumState[{"Register", {3, 3}, 4}]` runs about 3 minutes into `$IterationLimit::itlim` instead of failing at once (not re-run).

## Intended changes and behavior (not bugs)

23. **`FieldVariables` and `BosonicNormalOrder` changed in 9a5fa392** (Bruno, 2026-11-05): `FieldVariables` returns `{{a†}, {a}}`, the `"Scalars"` option and `BosonicRelations` are gone. The commit message says QF's tests and docs still use the old forms, and `audit/function-index.md` describes the old return. Two slides break (see below).

24. **`"Realignment"` prints `QuantumEntanglementMonotone::realignment` on every call.** The message explains that realignment is a separability test, not a monotone, so it is deliberate.

## Errors in the slides themselves

All five are fixed in Modeling-Slides.md (2026-11-11; none of these slides is in Continuous-Measurement-Slides.md). The edited cells run without messages in `QF-sweep/mine/fixcheck/` (the qutrit cell's only messages are the headless front-end ones of its `Manipulate`).

- **Qutrit Wigner function under dephasing:** `"PhaseSpace"` is indexed [p, q], so after `Transpose` the axis labeled q shows p. Fix: swap the two `FrameLabel` entries. Values and death times are unaffected.
- **NV center in a tilted field:** the slide reads `"Matrix"` of `"JX"`, `"JY"`, `"JZ"`, which QF stores diagonal in each operator's own eigenbasis, so the Hamiltonian is not the textbook one (the spectrum and displayed values happen to be right). Fix: `"MatrixRepresentation"`.
- **Balanced homodyne detection, and the beam splitter in the ladder-operator algebra:** broken by item 23 (empty curves, a False identity, about 20 messages). Fix in `batch09/proposals.md`.
- **Young-diagram slide:** loads TensorNetworks by context only, so it fails on a machine without TensorNetworks installed.
- **Transmon spectrum slide:** the exact charge dispersion reaches only 0.87 of the Koch asymptote at E_J/E_C = 80, so "merges with the exact levels" holds by eye on the log scale only.

## Audit entries that are now wrong

- `audit/mistakes.md`, "Property cache has no auto-invalidation": recommends `qs /. v -> v2`, which does nothing on WL 15.0.1 (item 2).
- `audit/function-index.md`, `FieldVariables`: describes the pre-9a5fa392 return form (item 23).
