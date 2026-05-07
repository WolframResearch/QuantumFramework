# Formula_Test — strict verification of stabilizer-formalism formulas against QF

A self-contained audit of the `Wolfram`QuantumFramework`` stabilizer subsystem
against the formulas extracted from 28 papers in
`OngoingProjects/Stabilizer/tex/`.

## Files

| File | Purpose |
|---|---|
| [`stabilizer-formulas.md`](stabilizer-formulas.md) | The 25-section formula reference, distilled from the 28 stabilizer papers |
| [`stabilizer-formulas-test.wls`](stabilizer-formulas-test.wls) | 132 strict `VerificationTest`s mapped 1:1 to the reference; runs on a fresh kernel via `wolframscript` |
| [`findings-report.md`](findings-report.md) | Three confirmed findings (F1, F2, F3) with reproductions and suggested fixes |

## Running

From a fresh kernel:

```bash
wolframscript -f OngoingProjects/Stabilizer/Formula_Test/stabilizer-formulas-test.wls
```

Current status: **126 pass / 6 fail**. The 6 failures collapse into 3 distinct
findings that are real package issues, not test bugs — see the report.

## Quick map

| Test section | Formula reference §                                 | Status |
|---|---|---|
| S0  | §0  Smoke tests / Pauli identities                       | 10/10 |
| S1  | §1  xz, GF(4) encoding, symplectic IP                    | 5/5  |
| S2  | §2  Heisenberg gate table                                | 11/11 |
| S3  | §3  Standard test states + named codes                   | 13/15 (F1) |
| S4  | §4  Density-matrix from stabilizer                       | 3/3  |
| S5  | §5  Measurement (single-qubit + Pauli string)            | 7/8 (F2) |
| S6  | §6  Inner products                                       | 7/7  |
| S7  | §7  Sum-of-stabilizers projection                        | 2/3 (F1 downstream) |
| S8  | §8  CSS structure                                        | 2/2  |
| S9  | §9  MacWilliams identity                                 | 1/1  |
| S10 | §10 Quantum bounds                                       | 4/4  |
| S11 | §11 Magic states / Clifford hierarchy                    | 5/5  |
| S12 | §12 Symplectic ↔ Clifford                                | 3/3  |
| S13 | §13 BSM / teleportation / fusion                         | 3/3  |
| S14 | §14 Local complementation                                | 3/3  |
| S15 | §15 Canonical forms                                      | 2/2  |
| S16 | §16 CliffordChannel (Yashin25)                           | 3/4 (F3) |
| S17 | §17 Qudit sanity                                         | 1/1  |
| S18 | §18 Pauli tracking                                       | 3/3  |
| S19 | §19 Universal one-line identities                        | 10/10 |
| S20 | §20 Round-trips                                          | 6/6  |
| S21 | §21 Stim-style                                           | 3/3  |
| S22 | §22 Cross-package fixtures                               | 1/1  |
| S23 | §23 Caveats                                              | 3/3  |
| S24 | §24 Concrete numeric fixtures                            | 2/3 (F1 downstream) |
| S25 | §25 Stabilizer Olympics 13-item conformance              | 13/13 |

## Findings (see [`findings-report.md`](findings-report.md))

- **F1** — `PauliStabilizer["5QubitCode" | "SteaneCode" | "9QubitCode"]` returns
  `|+_L⟩`, not `|0_L⟩` as documented (the n-th stabilizer is logical X̄
  instead of logical Z̄).
- **F2** — `ps["M", pauli_string]` does **not** correctly collapse the
  state to a +1 eigenstate of the measured Pauli when the operator
  anti-commutes with the stabilizer. The single-qubit `ps["M", q_integer]`
  path works correctly.
- **F3** — `CliffordChannel[ps_state][ps_state']` (state-prep channel applied
  to a same-dim PauliStabilizer) emits `CliffordChannel::stateevol` and
  falls back to dense materialization, contradicting the documented
  `nA == 0` dispatch case.
