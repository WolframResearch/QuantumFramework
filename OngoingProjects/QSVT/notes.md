# QSVT — current state

Brief report cross-referencing `OngoingProjects/QSVT/QSVT.nb` with the audit findings under `audit/`. Anchored at `cbbe9368`.

## What `QSVT.nb` shows

Pedagogical / exploratory notebook with four sections:

1. **Code** — source-level copy of all QSVT helpers (re-defines `BlockEncode`, `PCPhase`, `QSVT`, `QSP2QSVTangles`, `PauliDecompose`); identical to [QuantumFramework/Kernel/QSVT.m](QuantumFramework/Kernel/QSVT.m).
2. **Block encoding** — embeds a normalized $A$ into the unitary $\bigl[\begin{smallmatrix}A & \sqrt{I-AA^\dagger}\\ \sqrt{I-A^\dagger A} & -A^\dagger\end{smallmatrix}\bigr]$, then does ancilla post-selection (`"0"^†`). Then redoes the same encoding via **LCU**: `{αk,Uk}=PauliDecompose[A]`, build `PREP` from $\sqrt{α_k}$, `SEL=Multiplexer@@Uk`, sandwich as `{PREP, SEL, PREP^T}` to recover $A/\sum_k|α_k|$.
3. **QSVT** — two demos:
   - 3rd-degree Legendre polynomial via 4 pre-optimized angles `{-0.204, -0.912, 0.912, 0.204}`.
   - Inverse $s/x$ with $κ=4,\ s=0.10145775$, 44 phi angles consumed via `"QSP" -> True` (exercises `QSP2QSVTangles`).
4. **Quantum LinearSolve** — runs PennyLane externally (`StartExternalSession[..., "EnvironmentName" -> pennylane, ...]`) to optimize 51 phi angles via `qml.AdagradOptimizer`, then plugs them into a hand-built "real part" LCU `realCirc = {"H", C-QSVT[A^T, phi], C-QSVT[A, phi]^†, "H"}` to invert $A$. Cross-checked against `Normalize@LinearSolve[A,b]`. Mid-notebook caveat from the author: *"this even-odd circuit … is just unnecessary over-complicating the example … And I couldn't make it work for smaller matrix"* — `sumEvenOddCirc` cell is left as known-broken.

## Backing kernel

[QuantumFramework/Kernel/QSVT.m](QuantumFramework/Kernel/QSVT.m) — 49 LOC, 2024-04-12, **zero public exports**. All five symbols (`BlockEncode`, `PCPhase`, `QSP2QSVTangles`, `PauliDecompose`, `QSVT`) are `PackageScope` (`QSVT.m:3-7`). Nothing in `PacletInfo.wl` Symbols list. The notebook redefines them at top-level so it can run without going through the Wolfram\`QuantumFramework\`PackageScope\` ` context.

## Audit findings that touch this code path

- **`PauliDecompose` is defined twice** ([audit/CHANGELOG.md:70](audit/CHANGELOG.md), [audit/mistakes.md:476-482](audit/mistakes.md)).
  - `QuantumOperator/Properties.m:626-649` — operator-property path, returns Association `"XYZ" → coefficient`, accessed via `qo["PauliDecompose"]`.
  - `QSVT.m:33-41` — QSVT-internal helper, returns the `{coefficients, paulis}` **pair**.
  - The notebook's LCU section uses the pair form (`{αk, Uk} = PauliDecompose[A]`, line ~95 of converted .md) — i.e., the QSVT-internal one. The duplicate is benign here because both files load and the QSVT.m definition wins for matrix args.
- **Activity status: 🟢 stable / less-active** ([audit/activity-status.md:70-71](audit/activity-status.md)) — recent edits are rare; "PackageScope-only utility for QSVT circuits."
- **`problem-patterns.md` entry is a stub** ([audit/problem-patterns.md:98-99](audit/problem-patterns.md)) — only the file path is listed; no recipe / variants / doc-example pointer filled in. This notebook is a candidate source for filling it.
- **No companion doc page.** `audit/CHANGELOG.md` Phase-3 entry confirms 25 reference pages + 10 tutorials read; QSVT is not among them. The `OngoingProjects/QSVT/QSVT.nb` is currently the only end-user-facing material.

## Open items

1. **No public surface.** `QSVT[A, angles]` is reachable only via `PackageScope` or the redefinition trick the notebook uses. Either expose it (Symbols + Usage) or move it under an existing public umbrella (e.g., `QuantumCircuitOperator[{"QSVT", A, angles}]`).
2. **`sumEvenOddCirc` is broken** for the small-matrix case (per author's own mid-notebook note); the working path is `realCirc`. Either fix or delete `sumEvenOddCirc`.
3. **PennyLane dependency for angle generation.** The inverse-$1/x$ angles in the notebook are produced by external PennyLane optimization. There is no native WL angle-finder. The variational `QuantumLinearSolve` in [Kernel/QuantumOptimization.m](QuantumFramework/Kernel/QuantumOptimization.m) ([audit/activity-status.md:79-80](audit/activity-status.md)) is a *separate* HHL-style path (NMinimize on Fidelity), not a QSVT angle generator — so the two coexist without one obviating the other.
4. **Pauli-decomposition collision.** Either rename the QSVT.m helper (e.g., `QSVTPauliDecompose`) or have the QSVT path consume `qo["PauliDecompose"]` and reshape, eliminating the duplicate.
5. **Fill `problem-patterns.md` "Quantum signal processing / QSVT"** with the three recipes the notebook exercises: block-encoding-by-square-root, LCU-via-PauliDecompose, and QSVT polynomial approximation (give Legendre and inverse-$1/x$ as canonical examples).
