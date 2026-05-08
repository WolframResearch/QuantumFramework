# Formula_Test — strict verification of stabilizer-formalism formulas against QF

A self-contained audit of the `Wolfram`QuantumFramework`` stabilizer subsystem
against the formulas extracted from 28 papers in
`OngoingProjects/Stabilizer/tex/`.

## Files

| File | Purpose |
|---|---|
| [`stabilizer-formulas.md`](stabilizer-formulas.md) | The 25-section formula reference, distilled from the 28 stabilizer papers |
| [`stabilizer-formulas-test.wls`](stabilizer-formulas-test.wls) | 132 strict `VerificationTest`s mapped 1:1 to the reference; runs on a fresh kernel via `wolframscript` |
| [`findings-report.md`](findings-report.md) | Three confirmed real-bug findings (F1, F2, F3) with fresh-kernel reproductions and **concrete kernel patches** (Patch P1 / P2 / P3) |

## Running

From a fresh kernel:

```bash
wolframscript -f OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas-test.wls
```

## Status (2026-05-07, post-fix)

**132 / 132 passing.** All three findings (F1 / F2 / F3) have landed as
kernel patches (P1 / P2 / P3). The .wls runs green, and the same battery
also runs as part of the regular `Tests/RunTests.wls` runner via
[`Tests/Stabilizer/Formulas.wlt`](../../../Tests/Stabilizer/Formulas.wlt)
(top-level `VerificationTest` mirror generated from the .wls).

**Failing-bug history.** The TestIDs listed below carry an `F<n>` infix
as a permanent breadcrumb pointing at the finding that surfaced them.
After the patches, every one of these now passes; the names stay so a
future regression is easy to triage:

| TestID | Finding | Patch |
|---|---|---|
| `S3-F1-5QubitCode-Got97-16termSum-uptoPhase` | F1 (named codes returned `\|+_L⟩` not `\|0_L⟩`) | P1 → `NamedCodes.m` |
| `S3-F1-SteaneCode-LogicalZero-8amps-STRICT` | F1 | P1 |
| `S5-F2-PostMeas-State-IsEigenstate-STRICT` | F2 (Pauli-string measurement didn't collapse) | P2 → `PauliMeasure.m` |
| `S7-F1-Steane-CSS-Codeword-8terms-STRICT` | F1 | P1 |
| `S16-F3-CliffordChannel-StatePrep-Yields-Zeros-STRICT` | F3 (`cc[ps]` state-prep dispatch missing) | P3 → `CliffordChannel.m` |
| `S24-F1-SteaneCode-LogicalZero-8Amps-STRICT` | F1 | P1 |

The patches themselves are documented in
[`findings-report.md` § Patches](findings-report.md#patches-apply-to-kernel--out-of-scope-for-the-formula_test-folder).

## Quick map

| Test section | Formula reference §                                 | Status |
|---|---|---|
| S0  | §0  Smoke tests / Pauli identities                       | 10/10 |
| S1  | §1  xz, GF(4) encoding, symplectic IP                    | 5/5  |
| S2  | §2  Heisenberg gate table                                | 11/11 |
| S3  | §3  Standard test states + named codes                   | 15/15 |
| S4  | §4  Density-matrix from stabilizer                       | 3/3  |
| S5  | §5  Measurement (single-qubit + Pauli string)            | 8/8  |
| S6  | §6  Inner products                                       | 7/7  |
| S7  | §7  Sum-of-stabilizers projection                        | 3/3  |
| S8  | §8  CSS structure                                        | 2/2  |
| S9  | §9  MacWilliams identity                                 | 1/1  |
| S10 | §10 Quantum bounds                                       | 4/4  |
| S11 | §11 Magic states / Clifford hierarchy                    | 5/5  |
| S12 | §12 Symplectic ↔ Clifford                                | 3/3  |
| S13 | §13 BSM / teleportation / fusion                         | 3/3  |
| S14 | §14 Local complementation                                | 3/3  |
| S15 | §15 Canonical forms                                      | 2/2  |
| S16 | §16 CliffordChannel (Yashin25)                           | 4/4  |
| S17 | §17 Qudit sanity                                         | 1/1  |
| S18 | §18 Pauli tracking                                       | 3/3  |
| S19 | §19 Universal one-line identities                        | 10/10 |
| S20 | §20 Round-trips                                          | 6/6  |
| S21 | §21 Stim-style                                           | 3/3  |
| S22 | §22 Cross-package fixtures                               | 1/1  |
| S23 | §23 Caveats                                              | 3/3  |
| S24 | §24 Concrete numeric fixtures                            | 3/3  |
| S25 | §25 Stabilizer Olympics 13-item conformance              | 13/13 |

## Findings (see [`findings-report.md`](findings-report.md))

All three findings are now resolved on this branch by Patches P1 / P2 / P3
in [`findings-report.md`](findings-report.md). The descriptions below
record what the bug was, for cross-reference and future regression triage.

- **F1** (resolved by P1) — `PauliStabilizer["5QubitCode" | "SteaneCode" | "9QubitCode"]`
  returned `|+_L⟩`, not `|0_L⟩` as documented (the n-th stabilizer was
  logical X̄ instead of logical Z̄). Fix: replace `XXXXX` / `XXXXXXX` / `XXXXXXXXX`
  with the corresponding Z̄ string in
  [`Kernel/Stabilizer/NamedCodes.m`](../../../QuantumFramework/Kernel/Stabilizer/NamedCodes.m).
- **F2** (resolved by P2) — `ps["M", pauli_string]` did **not** correctly
  collapse the state to a +1 eigenstate of the measured Pauli when the
  operator anti-commuted with the stabilizer; the single-qubit
  `ps["M", q_integer]` path was already correct. Fix: write the pauli
  vector into the tableau on the random-outcome branch in
  [`Kernel/Stabilizer/PauliMeasure.m`](../../../QuantumFramework/Kernel/Stabilizer/PauliMeasure.m).
- **F3** (resolved by P3) — `CliffordChannel[ps_state][ps_state']`
  (state-prep channel applied to a same-dim PauliStabilizer) emitted
  `CliffordChannel::stateevol` and fell back to dense materialization,
  contradicting the documented `nA == 0` dispatch case. Fix: add the
  `nA == 0` arm to `cliffordChannelToPauliStabilizer` dispatch in
  [`Kernel/Stabilizer/CliffordChannel.m`](../../../QuantumFramework/Kernel/Stabilizer/CliffordChannel.m).
