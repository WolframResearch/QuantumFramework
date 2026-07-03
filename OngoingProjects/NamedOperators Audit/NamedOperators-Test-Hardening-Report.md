# Named-Operator Test Hardening: Physics/Math Coverage and the Lowercase-Alias Fix

**Date:** 2026-07-02. **Scope:** all 108 entries of `$QuantumOperatorNames`
(`QuantumFramework/Kernel/QuantumOperator/NamedOperators.m:10-35`) and every call
signature reachable from them. **Deliverable:** six new `Tests/*.wlt` files
(292 tests), one source fix, and this report.

## Big picture

Every named operator in the framework is a physical object with a definite
matrix, algebra, and set of defining identities. The old test suite only checked
that each operator had the right **shape** (its `"Dimensions"`) and that it
survived a serialize/reconstruct round trip. Both checks are blind to a wrong
matrix of the right shape: a buggy $J_x$ that silently returns the diagonal
$J_z$ is still $2\times 2$ and still round-trips, which is exactly how that bug
lived undetected for a long time. This work replaces the shape-only net with one
that tests **mathematical content**: closed forms, the su(2) and Weyl-Heisenberg
algebras, unitarity and Hermiticity, generator-and-angle conventions for the
parametric gates, the snake and spider identities of the categorical operators,
and the trace-preserving structure of the channels. A component that returns the
wrong operator now fails a physics assertion, not just a shape check.

## What was audited

The operators were enumerated two ways and reconciled: the registry list
`$QuantumOperatorNames` (verified live at 108 names, no drift), and the actual
pattern definitions in `NamedOperators.m` read in full so that every call
**shape** per name was covered (bare name, `name[d]`, `name[params]`,
`name[params, order]`, rule-forms, and the doubled-picture arguments). Coverage
runs past the default qubit case: the generalized Paulis, rotations, `SUM`,
`Fourier`, `HeisenbergWeyl`, `Curry`/`Uncurry`, `Cup`/`Cap`, `Copy`, `Trace`,
`Marginal`, `WignerD`, and the ladder operators are all exercised at
$d = 3, 5$ and, where relevant, on multiple qudits.

The six files, organized by operator family:

| File | Tests | Covers |
|---|---:|---|
| `NamedOperators-Pauli.wlt` | 56 | X, Y, Z, PauliX/Y/Z, Shift, ShiftPhase, H, Hadamard, NOT, RootNOT, S, T, V, SX, "0", "1", lowercase aliases, Pauli strings |
| `NamedOperators-Rotations.wlt` | 57 | RX/RY/RZ, XRotation/YRotation/ZRotation, R, U/U1/U2/U3, Phase, P, PhaseShift, GlobalPhase, Diagonal, FlipSign, WignerD |
| `NamedOperators-TwoQubit.wlt` | 46 | C/Controlled, C0/Controlled0, CX/CY/CZ/CH/CS/CT/CPHASE/CNOT/C0NOT, SWAP, RootSWAP, CSWAP/Fredkin, C0SWAP, Toffoli, Deutsch, Multiplexer, Braid |
| `NamedOperators-Qudit.wlt` | 65 | qudit Pauli conventions, Weyl relation, HeisenbergWeyl, SUM, Fourier/InverseFourier, Hadamard at $d = 2,3,5$ |
| `NamedOperators-Structural.wlt` | 44 | Identity, Zero, Permutation, Curry/Uncurry, Spider family, Cup/Cap/Trace, Copy/Measure/Encode/Decohere/Marginal/Discard/Reset, Left/Right/Double |
| `NamedOperators-Channels.wlt` | 24 | ladder ($a$, $a^\dagger$), RandomUnitary/RandomHermitian, Switch, Channel, Measurement, Liouvillian, Hamiltonian |

The pre-existing `Tests/AngularMomentum.wlt` (94 tests) already covers the
$J_x, J_y, J_z, J_{\pm}$ family and is left in place as the pattern the new files
follow.

## The bug found and fixed: dead case-insensitive dispatch

The framework carries a case-folding layer so that a lowercase spelling resolves
to the registered operator (`$upperCasesOperatorNames` at `NamedOperators.m:933`).
It never fired. The guard demanded both that the name contain a lowercase letter
(`ToUpperCase[name] =!= name`) and that the name itself be a key of a map whose
keys are all uppercase (`KeyExistsQ[$upperCasesOperatorNames, name]`). Those two
conditions are mutually exclusive, so `QuantumOperator["cnot"]`,
`QuantumOperator["x"]`, `QuantumOperator["Rx"[0.3]]` all returned a
`QuantumOperator::invalidName` failure instead of the gate. This is the same
class of defect as the $J_x$ bug: an intended behavior silently broken, invisible
to a shape test.

The fix keys the lookup on `ToUpperCase[name]` and adds a `! MemberQ` guard so
that a mixed-case registered name (like `"Fourier"` or `"Identity"`) stays on its
own rule instead of looping back through the alias rule:

```wolfram
QuantumOperator[name_String, opts___] /;
    ToUpperCase[name] =!= name && ! MemberQ[$QuantumOperatorNames, name] &&
    KeyExistsQ[$upperCasesOperatorNames, ToUpperCase[name]] :=
    QuantumOperator[$upperCasesOperatorNames[ToUpperCase[name]], opts]
```

The `MemberQ` guard is load-bearing: without it, `"Fourier"` (which does contain
lowercase letters) would map to `$upperCasesOperatorNames["FOURIER"]` = `"Fourier"`
and recurse forever. Lowercase aliases now resolve, unknown names still fail with
`invalidName`, and no registered name changes behavior. The `alias-*` tests in
`NamedOperators-Pauli.wlt` guard the fix and fail if it is reverted.

## Convention quirks locked with characterization tests

Where the "correct" answer is a convention choice rather than a universal truth,
the current behavior is pinned by a test with a present-tense comment, so the
framework can choose to keep or change the convention deliberately rather than
drift by accident.

- **Qudit rotations are non-unitary for $d > 2$.** `RX`/`RY`/`RZ` are built as
  $\exp(-i\theta X_d/2)$ from the generalized Pauli generators, which for
  $d > 2$ are unitary but not Hermitian, so the exponential is not unitary and no
  message is raised. A genuine qudit rotation needs a Hermitian generator (the
  spin-$j$ operators `"JX"/"JY"/"JZ"` are the natural choice). The tests lock
  unitarity at $d = 2$ and characterize the non-unitarity at $d = 3$; if the
  generator is ever made Hermitian, the three `NonUnitary-*-d3` tests flip and
  announce the change. This is a real defect, surfaced but not silently patched,
  because fixing it is a convention decision (which Hermitian generator) rather
  than an unambiguous correction.

- **The clock is the inverse clock.** $Z_d = \mathrm{diag}(1, \omega^{-1},
  \dots, \omega^{-(d-1)})$ with $\omega = e^{2\pi i/d}$, the complex conjugate of
  the textbook clock $\mathrm{diag}(1, \omega, \dots)$. It agrees at $d = 2$.

- **The Weyl relation is $X_d Z_d = \omega\, Z_d X_d$.** Verified live at
  $d = 3, 5$. (The earlier named-operator audit stated the opposite sign,
  $XZ = \omega^{-1} ZX$; re-derivation in the kernel shows $XZ = \omega ZX$, and
  the test locks the sign that actually holds.)

- **`HeisenbergWeyl` uses the standard clock.** `HeisenbergWeyl[d,0,1]` is
  $\mathrm{diag}(1, \omega, \dots)$, so it is the conjugate of, and disagrees
  with, the framework's own `"Z"[d]`. Both facts are asserted.

- **`Fourier` diagonalizes the shift under $F^\dagger$.** $F^\dagger X_d F = Z_d$
  (QF's inverse clock); $F X_d F^\dagger$ gives the standard clock instead.

- **`Deutsch[theta]` is the half-angle gate.** Its doubly-controlled block is
  $i\,R_x(\theta)$, so the Toffoli is reached at $\theta = \pi$, not the
  literature's $\theta = \pi/2$.

- **Bare `"SX"` is the two-qubit $S \otimes X$, not $\sqrt{X}$.** It is the only
  registered name whose characters all lie in the Pauli-string alphabet, so the
  chain rule parses it as a tensor product; the single-qubit $\sqrt{X}$ is
  `"SX"[]` or `"V"`.

- **`"+"`/`"-"`/`Up`/`Down` are oscillator ladder operators.** `"+"[d]` is the
  truncated annihilation $a$ (superdiagonal $\sqrt{n}$), `"-"[d]` is $a^\dagger$,
  and $a^\dagger a$ is the number operator $\mathrm{diag}(0,1,\dots,d-1)$. At
  $d = 2$, `"+"` is $\sigma_+ = |0\rangle\langle 1|$.

- **The `Switch` order.** `Switch[A,B]` is
  $|0\rangle\langle 0| \otimes AB + |1\rangle\langle 1| \otimes BA$ (control on
  qubit 1).

- **`"J+"/"J-"` are reversed versus the textbook labels** (already locked in
  `AngularMomentum.wlt`): in the highest-weight-first basis, QF's `"J+"` equals
  $J_x - i J_y$ and lowers $J_z$.

## A signature that is degenerate, documented not tested

The integer-control shorthand `Controlled[op, n]` (encoding the control pattern
in the bits of $n$) does not produce a valid operator: the target defaults to
qubit 1 while the controls are assigned qubits $1, 2, \dots$, so target and
control collide even when the target wire is given explicitly. The clean
multi-control form `Controlled["NOT", {1,2}]` with target order `{3}` works and is
tested (it reproduces the Toffoli). The degenerate integer form is left as a
finding rather than pinned, since locking a broken shape adds no value; it is a
candidate for a follow-up fix.

## Results

All 292 new tests pass. The full suite is `2001/2002`; the single failure is
`CrossPackage_QuantumClifford.wlt / QC-PackageSource-OnDisk`, a check for whether
the external QuantumClifford Julia package is installed on this machine. It is
pre-existing and environmental, unrelated to this work. The kernel suites that
exercise named-operator construction are all green, confirming the source fix
introduced no regression: `QuantumOperator.wlt 86/86`,
`QuantumCircuitOperator.wlt 97/97`, `AngularMomentum.wlt 94/94`,
`Qiskit.wlt 132/132`.

```
NamedOperators-Pauli.wlt        56/56
NamedOperators-Rotations.wlt    57/57
NamedOperators-TwoQubit.wlt     46/46
NamedOperators-Qudit.wlt        65/65
NamedOperators-Structural.wlt   44/44
NamedOperators-Channels.wlt     24/24
```

## How the tests are written

Each file compares the operator's computational-basis matrix
(`op["Computational"]["Matrix"]`) against a closed form or a reference built from
first principles, using a numeric comparator (`Chop` after `N`, tolerance
$10^{-10}$) for anything with radicals or exponentials, and a symbolic comparator
(`FullSimplify`) for the parametric gates so that the whole one-parameter family
is checked at a symbolic angle rather than one sampled value. Parametrized checks
run inside `Do` loops over dimension or spin with distinct `TestID`s, which
`TestReport` captures. The files live at repo-root `Tests/`, are auto-discovered
by `Tests/RunTests.wls`, and load the local `QuantumFramework/` source.
