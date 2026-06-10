# Symbolic two-qubit KAK: verified plan

Goal: extend `qo["KAK"]` to **structured / parametric** symbolic 2-qubit gates
(`RXX[θ]`, `CPhase[φ]`, `Can[a,b,c]`, products of commuting generators, etc.),
returning a closed-form `QuantumCircuitOperator`. This is the differentiator over
Qiskit, which cannot decompose a symbolic gate at all. Generic dense symbolic 4x4
is explicitly out of scope (its eigendecomposition is a quartic `Root` blowup).

All claims below are verified with `wolframscript` (see "Verification log").

## Why the current numeric path cannot be reused as-is

`TwoQubitKAK` is hard-wired to numerics in three places:
1. `mat = N @ qo["MatrixRepresentation"]` forces machine numbers.
2. `Eigensystem[Re[m2] + Sqrt[2.] Im[m2]]` uses a **machine** separator `Sqrt[2.]`
   and `Re`/`Im` that do not split on symbolic complex expressions.
3. branch logic `If[Re[Det[...]] < 0, ...]`, `Arg[...]` have no truth value on
   free symbols.

The fix is a parallel symbolic branch, not a patch of the numeric one.

## The symbolic algorithm (same magic-basis KAK, exact arithmetic)

For a 4x4 unitary `U` and assumptions `assum` (the parameters are real), inside
`Assuming[assum, ...]` with `sf = FullSimplify[#, assum] &`:

1. `su = U / Det[U]^(1/4)` (exact; the 1/4 root branch is absorbed by the global
   phase at the end, exactly as in the numeric path).
2. `uM = sf @ ComplexExpand[mbd . su . mb]`.
3. `m2 = sf @ ComplexExpand[Transpose[uM] . uM]` (complex-symmetric, unitary).
4. `reM = ComplexExpand[Re[m2]]`, `imM = ComplexExpand[Im[m2]]` — **`ComplexExpand`
   is the load-bearing builtin**: it is what makes `Re`/`Im` split once the
   parameters are declared real.
5. `P = Transpose[Orthogonalize[Eigensystem[reM + Sqrt[2] imM][[2]]]]` — **exact
   `Sqrt[2]`** (not `Sqrt[2.]`), and **`Orthogonalize`** to force a real-orthonormal
   eigenbasis even on degenerate eigenspaces.
6. `theta = sf[Arg[Diagonal[Transpose[P] . m2 . P]] / 2]`.
7. `K1 = sf @ ComplexExpand[uM . P . DiagonalMatrix[Exp[-I theta]]]`; SO(4) sign fix
   guarded by `sf[Re[Det[...]]] < 0` (decidable once `assum` is in scope).
8. `coords = sf[Transpose[$kakBellSigns] . theta / 4]`.
9. `corr = sf @ ComplexExpand[mb . DiagonalMatrix[Exp[I (theta - signs . coords)]] . mbd]`.
10. `left = kf[mb . K1 . mbd . corr]`, `right = kf[mb . Transpose[P] . mbd]`.

The angles end up as `Arg`/`ArcTan` expressions that `FullSimplify` reduces under
`assum` (e.g. `RXX[θ]` -> `c = {-θ/2, 0, 0}`).

## Builtins used (from the V15.0 matrix-decomposition audit), and why

| Step | Builtin | Why it is the right choice |
|---|---|---|
| Re/Im split | **`ComplexExpand`** | the only way to evaluate `Re`/`Im` of a symbolic complex matrix once parameters are real; without it `Re[m2]` stays inert and step 5 fails. |
| simultaneous diagonalization | **`Eigensystem`** of `reM + Sqrt[2] imM` | exact symbolic eigensystem of a real-symmetric matrix; the commuting-pair trick keeps it one eigensolve instead of two. |
| real-orthonormal eigenbasis | **`Orthogonalize`** | a symbolic `Eigensystem` returns un-normalized, possibly non-orthogonal eigenvectors on degenerate eigenspaces; `Orthogonalize` (Gram-Schmidt) repairs both. |
| Kronecker factoring | **`SingularValueDecomposition`** (reshuffle) | rank-1 realignment split; needs a zero-singular-value guard (see edge cases). |
| simplification | **`FullSimplify[..., assum]`** | collapses the `Arg`/`ArcTan` angle expressions to clean closed forms. |

`Det`, `Arg`, `MatrixExp`, `DiagonalMatrix`, `Transpose` are all exact-friendly.

## Edge cases to harden (found by verification)

1. **`kf` on a degenerate local (`Sqrt[0]`).** `CPhase[φ]` produced correct coords
   but a `1/Sqrt[0]` in the reshuffle SVD, because one local factor's reshuffle has
   a zero leading singular value pattern. Fix: in the symbolic `kf`, special-case a
   diagonal / rank-deficient reshuffle (read the factors off a nonzero pivot
   instead of dividing by `Sqrt[s[[1]]]`), or `Simplify` the SVD under `assum`.
2. **Branch decidability.** `Re[Det[P]] < 0` and `Arg[...]` need `assum` in scope;
   require the caller to pass real-parameter assumptions (default `Element[vars,
   Reals]`, vars from `qo["Parameters"]`).
3. **Generic dense symbolic.** Guard with a `LeafCount` / `TimeConstrained` bail:
   if `m2` or `coords` exceed a leaf budget, or the eigensystem times out, fall
   back to the numeric path on `N[...]` substitutions or message the user. (Generic
   symbolic is a `Root` blowup and is out of scope.)

## Surface and dispatch

- Detect symbolic at the top of `qo["KAK"]`: if `! MatrixQ[mat, NumericQ]`, route to
  the symbolic worker `TwoQubitKAKSymbolic[mat, assum]` instead of the current
  `kaksymbolic` bail. Assumptions come from `Element[qo["Parameters"], Reals]`
  (fall back to `Element[Variables[mat], Reals]`).
- Emission reuses the existing gadget path; the single-qubit `"U"` gates carry
  symbolic angles via `UnitaryAnglesWithPhase` (already symbolic-capable — it is
  what `ZYZ` uses on symbolic 1-qubit gates). The global phase is recovered
  symbolically (`Arg` of the projected trace, simplified under `assum`) rather than
  numerically.
- `TimeConstrained` + `LeafCount` guard for the generic-symbolic bail.

## Staging

1. `TwoQubitKAKSymbolic` worker (steps 1-10) with the `kf` zero-SV guard; verify
   exact reconstruction (`FullSimplify[recon - U] == 0`) for `RXX[θ]`, `RYY`, `RZZ`,
   `CPhase[φ]`, `Can[a,b,c]`, `RXX·RZZ`.
2. Symbolic emission (gadgets with symbolic angles; symbolic global phase).
3. Dispatch: symbolic 2-qubit -> symbolic worker; generic-symbolic / timeout ->
   numeric-substitution fallback or message.
4. Tests: add a symbolic section to `Tests/TwoQubitKAK.wlt`
   (`FullSimplify[circuitMatrix - U] == 0` for the structured battery).

## Verification log (wolframscript)

- Coords closed-form: `RXX[θ]` -> `{-θ/2, 0, 0}` (LeafCount 22); `CPhase[φ]` ->
  `{-φ/4, 0, 0}` (LeafCount 35); `RXX·RZZ` -> closed form (LeafCount 184).
- **Full reconstruction `RXX[θ]`: `FullSimplify[recon - U] == 0` (EXACT)** with
  `coords = {-θ/2, 0, 0}`, using exact `Sqrt[2]` + `ComplexExpand` + `Orthogonalize`
  + `Assuming[0 < θ < π]`.
- `CPhase[φ]`: correct coords, `kf` hit `Sqrt[0]` -> the documented edge-case guard.
