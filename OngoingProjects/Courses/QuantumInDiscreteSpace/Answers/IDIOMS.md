# Minimal-idiom catalog (WL and QF)

The answer key's standing quality bar is not "correct" but "correct **and** the minimal faithful
construction." A verified sledgehammer still passes every correctness check, so this catalog names the
lightest idiom for each recurring task and lists the over-engineering anti-patterns to hunt for. It is
the rubric for the authoring idiom-search gate and for the `/wl-verify` minimal-idiom review role.

**How to use.** (1) Authoring: before finalizing a cell, enumerate >=2 constructions, read the primitive's
own doc page for a lighter form, pick the minimal faithful one; if you keep a heavier tool, the prose (or
a note here) must say why the lighter one is unfaithful. (2) Review: for each cell, propose a lighter
built-in equivalent and prove agreement in a fresh kernel. This catalog is a *starting rubric*, not
gospel: every proposed replacement is verified in-kernel before it is applied.

## Preferred minimal idiom by task

| Task | Minimal WL | Minimal QF | Avoid |
|------|-----------|-----------|-------|
| Pauli vector $\vec\sigma$ | `PauliMatrix[{1,2,3}]` (Listable) | -- | `Table[PauliMatrix[k],{k,3}]`, `PauliMatrix/@Range[3]` |
| Unit axis $\hat n(\theta,\phi)$ | `FromSphericalCoordinates[{1,\[Theta],\[Phi]}]` | -- | hand-typed `{Sin[\[Theta]]Cos[\[Phi]],...}` |
| Generic pure qubit | `{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}` | `QuantumState[{...}]` / named `QuantumState["Plus"]` | real `{Cos t, Sin t}` (pins $\langle Y\rangle$, phase) |
| Generic mixed qubit | `1/2 (IdentityMatrix[2] + {rx,ry,rz} . PauliMatrix[{1,2,3}])` | `QuantumState[rho]` | matrix-entry `{{a,c},{c^*,1-a}}` |
| Normalize | `Normalize[v]` | `state["Normalized"]` | `v/Sqrt[v.v]`, `v/Sqrt[d]` |
| Uniform superposition | `Normalize@ConstantArray[1,d]` | `QuantumState["UniformSuperposition", d]` | `ConstantArray[1,d]/Sqrt[d]` |
| Propagator $U(t)=e^{-iHt}$ | `MatrixExp[-I H t]` | `QuantumOperator[MatrixExp[-I H t]]` | building $\cos/\sin$ blocks by hand |
| **Apply propagator to a state** | **`MatrixExp[-I H t, \[Psi]0]`** (two-arg $e^m v$) | `QuantumEvolve[QuantumOperator[H], state, t]` | `DSolveValue[...]` (unless the ODE itself is the point); `MatrixExp[-I H t].\[Psi]0` (unless the matrix is shown) |
| Solve time-dependent ODE | `NDSolveValue`/`DSolveValue` with `vec[t] \[Element] Vectors[2]` | `QuantumEvolve` | component form `{a'[t],b'[t]}` (heavier than the `Vectors[2]` domain) |
| Pure-state expectation | `Conjugate[\[Psi]] . A . \[Psi]` | `QuantumMeasurementOperator[A][state]["Mean"]` | -- |
| Mixed-state expectation | `Tr[A . rho]` | `QuantumMeasurementOperator[A][state]["Mean"]` | `Conjugate[\[Psi]].A.\[Psi]` on a MIXED state (wrong) |
| Variance | `Tr[A.A.rho] - Tr[A.rho]^2` (or pure analog) | compute `<A^2>-<A>^2` | (no QF `"Variance"`) |
| Eigit-decomposition | `Eigensystem` / `Eigenvalues` / `Eigenvectors` | `["Eigenvalues"]`,`["Eigenvectors"]`,`["Eigensystem"]`,`["Diagonalize"]` | manual characteristic-poly |
| Spectral reconstruction | `EigenvalueDecomposition` | `["Diagonalize"]` | hand-built $\sum a\lvert a\rangle\langle a\rvert$ |
| Inner product / overlap | `Conjugate[a] . b` | `(bra["Dagger"] @ ket)["Scalar"]` | -- |
| Composite system | `KroneckerProduct[...]` | `QuantumTensorProduct[...]` / string `QuantumOperator["XX"]` | hand-indexed tensor |
| Commutator / anticommutator | `a.b-b.a` or `Commutator[a,b,Dot]` / `a.b+b.a` | `Commutator[qoA,qoB]` (reduces) / `a@b+b@a` | (QF `Anticommutator` does NOT reduce) |
| Born probabilities | `Abs[\[Psi]]^2` / `Tr[P.rho]` | `state["ProbabilitiesList"]` | dilating just to read a probability |
| Trotter step | -- | `QuantumCircuitOperator["Trotterization"[{...},...]]` | hand-alternated `MatrixPower` (unless the point) |
| Pauli expansion | `Tr[PauliMatrix[j].A]/2` | `qo["PauliDecompose"]` | -- |

## Anti-patterns to flag in review

1. `DSolveValue`/`NDSolveValue` used to get $e^{-iHt}\lvert\psi_0\rangle$ where `MatrixExp[-I H t, \[Psi]0]`
   is the one-call answer (the ODE machinery earns its place only when *solving the ODE* is the taught skill).
2. `Table[f[k],{k,...}]` / `Map[f,list]` where `f` is **Listable** (`f[range]`, `Sin[data]`).
3. Manual `/Sqrt[d]` or `v/Sqrt[v.v]` instead of `Normalize`.
4. Hardcoded matrices/vectors where a named state/gate or a built-in generator exists
   (`QuantumState["Plus"]`, `PauliMatrix`, `FourierMatrix`, `FromSphericalCoordinates`).
5. `\langle\psi\rvert A\lvert\psi\rangle` sandwich on a **mixed** state (must be `Tr[A rho]`).
6. Procedural `Do`/`For`/`While`/`AppendTo` where `Map`/`Table`/`Fold`/`Nest`/`Sow`-`Reap` is native.
7. `MatrixExp[m].v` where the two-arg `MatrixExp[m,v]` is meant (unless the matrix $m$ is itself displayed).
8. Extracting `["Matrix"]//Normal` as the QF *deliverable* when the rich object (that acts, composes,
   reports its spectrum) is the point (show the object; a matrix is one view).
9. Re-deriving a closed form and presenting it as the definition (compute from the definition instead).

## Audit log

2026-07-04, full-manual idiom audit (Parts 1-10). Method: grep for the mechanical anti-patterns across
all ten files, then a semantic grep (QuantumChannel / measurement-dilation / recomputation / object-flattening),
each candidate verified in a fresh kernel before applying. Result: the manual is largely minimal already;
6 genuine fixes, all of the "apply-a-propagator-by-dot" or "recompute/relist" class.

- **7.1b (Part 7)**: `DSolveValue`-only for the evolved state -> now leads with `MatrixExp[-I H t, \[Psi]0]`
  (minimal two-arg), keeps `DSolveValue` as the explicit ODE-integration contrast; three-way agreement with
  `QuantumEvolve` verified `=== True`. Cause: first-working-construction lock-in (shipped the handed
  `DSolveValue` literal without an idiom search).
- **7.4 (Part 7)**: `MatrixExp[-I (\[Omega]/2) PauliMatrix[3] t] . {1,1}/Sqrt[2]` -> `MatrixExp[..., {1,1}/Sqrt[2]]`
  (two-arg). Downstream Bloch vector unchanged (`{Cos[t\[Omega]],Sin[t\[Omega]],0}`).
- **7.7 (Part 7)**: `MatrixExp[-I h t] . ({1,1}/Sqrt[2])` -> `MatrixExp[-I h t, {1,1}/Sqrt[2]]` (two-arg).
- **2.2 (Part 2)**: `\[Sigma] = Table[PauliMatrix[j], {j, 3}]` -> `PauliMatrix[{1,2,3}]` (Listable; matches the
  Part-3+ convention). Verified `Table[PauliMatrix[j],{j,3}] === PauliMatrix[{1,2,3}]`.
- **8.x spin-coherent (Part 8)**: `MatrixExp[-I \[Beta] jy] . UnitVector[n,1]` -> `MatrixExp[-I \[Beta] jy, UnitVector[n,1]]`.
- **8.x total-spin (Part 8)**: `QuantumOperator[Sum[stot[[k]].stot[[k]],{k,3}]]["Eigenvectors"]` recomputed
  the already-bound `s2` -> reuse `QuantumOperator[s2]["Eigenvectors"]`.

Deliberately LEFT (not over-engineering): `Table[Conjugate[\[Psi]].PauliMatrix[j].\[Psi],{j,3}]` (the Bloch
vector *definition*, clearer than a rank-3 broadcast); `MatrixExp[A].MatrixExp[B]` in Trotter/Wigner-D
(matrix x matrix, two-arg N/A); `NDSolve` in 7.12 Landau-Zener (genuinely time-dependent H); state literals
`{1,1}/Sqrt[2]` etc.; `["Matrix"]//Normal` inside WL===QF comparison cells; `QuantumMeasurementOperator[A][state]["Mean"]`
for expectations (the established idiom). No `QuantumChannel` appears anywhere, so no channel-over-operator inflation.

VERIFIED (fresh-context wl-verifier, OPEN ISSUES 0): all six edits reproduced in a clean kernel, including
the exactness sanity `MatrixExp[m,v] === MatrixExp[m].v` (numeric + symbolic), the 7.1b three-way agreement,
and driver `errored/messaged={}` on Parts 2/7/8 with no em/en dashes. (Edit 5 is question 8.3 line 109; the
8.7 coherent-state chains at Part-08 lines 203/210 were deliberately left as `.UnitVector` and are out of scope.)
