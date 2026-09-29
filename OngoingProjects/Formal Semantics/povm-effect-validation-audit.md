# Should the effect-list `QuantumMeasurementOperator` validate that each effect is Hermitian and positive semidefinite?

Investigation and verdict. No kernel change was made; this is a report. Verified live against the worktree kernel at HEAD `923c7035` (paclet `2.1.1`, Wolfram Language 15), which is ahead of the observation anchor `a1799b61` cited in the mismatch record. The two effect-list branches under study are unchanged from that anchor: `Kernel/QuantumMeasurementOperator/QuantumMeasurementOperator.m:111` (`MatrixPower[#, 1/2] & /@ tensor`, effect-matrix list) and `:137` (`Sqrt /@ ops`, operator list). Context read first: item 3 of `implementation-mismatches.md`, and sections 1, 5, and 11 of `QF-formal-semantics.md`.

## Verdict

The genuine defect is **not** "the effect-list constructor is missing a physics check." It is that **the constructor performs a derived, partial computation, the matrix square root `M_m = sqrt(E_m)`, and on inputs where that computation is undefined it neither stores the input faithfully nor fails cleanly**: it returns an object whose head is still `QuantumMeasurementOperator` but which carries an unevaluated `MatrixPower`, leaks `System` messages, and (in two of the three cases) applies to finite but physically wrong numbers.

A Hermitian-and-positive-semidefinite guard on numeric effects is the right repair, and it does **not** contradict QF's shapes-not-physics stance. That stance is not "never inspect the numbers"; as QF actually practices it, it is "accept any shape-valid input, store it faithfully so every property stays computable, and expose the physical predicates (`"PhysicalQ"`, `"UnitaryQ"`, `"TracePreservingQ"`) as tests the user reads afterward." The effect-list form cannot honor that contract for a non-positive-semidefinite effect, because it does not store the input at all: it stores `sqrt(E_m)`, which for such an input either does not exist or no longer represents the requested effect. Guarding it is a domain check on a partial construction, the same kind of check QF already added as `QuantumChannel::emptyKraus` and the single-operator Kraus reroute, not an enforcement of physical validity on an otherwise well-formed object.

Precisely, a clean version checks, for each supplied effect that is a fully numeric matrix, that it is Hermitian and that its least eigenvalue is at least `-tol`, and otherwise emits a `QuantumMeasurementOperator::` message and returns a bare `Failure`. It leaves symbolic and parametric effects unchecked, and it does not check the completeness relation `sum E_m = I`. The details, and why the weaker candidates fail, are below.

## 1. The observation, reproduced

Load on its own line (`Needs` after `PacletDirectoryLoad`, so no QF symbol is created in `Global`` before the context loads):

```wl
PacletDirectoryLoad["<worktree>/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"]
```

Three inputs, each a list read as POVM effects. The head stays `QuantumMeasurementOperator` in every case, so an invalid object is returned rather than a clean failure.

**Single nilpotent effect** (`{{0,1},{0,0}}` has no square root):

```wl
QuantumMeasurementOperator[{{{0, 1}, {0, 0}}}]
```

Leaks `TensorRank::rect` three times and `General::stop` at construction. The stored operator carries `MatrixPower[{{0,1},{0,0}}, 1/2]` unevaluated inside a `HoldForm[DiagonalMatrix[...]]`; the object is so malformed that reading its `"Operators"` cascades further, `DiagonalMatrix::nosmat`, `AssociationThread::idim`, and `Values::invrl`.

**Nilpotent effect with the identity** (the two-element form, routed through the same rank-3 branch):

```wl
QuantumMeasurementOperator[{{{0, 1}, {0, 0}}, IdentityMatrix[2]}]
```

Constructs **silently**. Its stored `"Operators"` and `"POVMElements"` carry the same unevaluated `MatrixPower[{{0,1},{0,0}}, 1/2]`. Applying it to a state returns a `QuantumMeasurement` whose `"ProbabilitiesList"` is `Abs[MatrixPower[...]]^2 / (1 + Abs[MatrixPower[...]]^2)` and the like, symbolic expressions rather than numbers, and the application leaks `MatrixPower::nosol` three times and `General::stop`.

**Hermitian but negative effect** (`{{-1,0},{0,0}}` is Hermitian, eigenvalue `-1`):

```wl
QuantumMeasurementOperator[{{{-1, 0}, {0, 0}}, {{0, 0}, {0, 1}}}]
```

Constructs **silently**. Its square root `sqrt(diag(-1,0)) = diag(I,0)` exists but is complex; `M_1^dagger M_1 = diag(1,0)`, a positive projector that is not the negative effect the user supplied. `"POVMElements"` collapses to `ComplexInfinity` and `Indeterminate` (the rescale step divides by a zero it should not have reached, raising `Power::infy` and `Infinity::indet`), yet applying the object to `|0><0|` returns the finite, meaningless list `{1, 0}`: the complex root silently substituted a different, positive POVM, so an input that is not a POVM at all passes without any error.

The `TensorRank::rect`, `MatrixPower::nosol`, `Power::infy`, and `Infinity::indet` outputs are all leaked `System` messages, not QF diagnostics.

## 2. Core-design audit: where QF validates physics and where it does not

For each constructor: what it accepts silently, and what property lets a user test validity afterward. This is the backbone of the argument.

| Constructor | Physics checked at construction | Input stored faithfully? | Working validity readout afterward |
|---|---|---|---|
| `QuantumState[matrix]` | none: normalization, unit trace, positivity, Hermiticity all unchecked | yes, verbatim | `"PhysicalQ"` / `"Type"` (`PositiveSemidefiniteMatrixQ[N @ dm]`, `QuantumState/Properties.m:624`), `"NormalizedQ"`, `"PureStateQ"` |
| `QuantumOperator[matrix]` | none: unitarity unchecked | yes, verbatim | `"UnitaryQ"`, `"HermitianQ"` |
| `QuantumMeasurementOperator[obs]` (single matrix, rank 2, `:125`) | none: Hermiticity unchecked; stores the operator itself via `QuantumOperator[tensor]`, no square root | yes, verbatim | `"HermitianQ"`, plus the eigenvalues it reads out |
| `QuantumChannel[{M1,...}]` | none: trace preservation and complete positivity unchecked, but `{}` is rejected (`QuantumChannel::emptyKraus`, `QuantumChannel.m:27`) and a one-element list is rerouted to the single-operator constructor (`:40`) | yes, the Kraus operators are stored | `"TracePreservingQ"`; complete positivity needs no test, since a Kraus/Stinespring dilation is completely positive by construction |
| `QuantumMeasurementOperator[{E1,...}]` (rank-3, `:111`) and `[{qo1,...}]` (`:137`) | none | **no**: stores the derived `sqrt(E_m)`, not the effects | **none that works**: `"POVMQ"` is `qmo["Type"] === "POVM"` (`QuantumMeasurementOperator/Properties.m:197`), a structural tag of which branch built the object, and returns `True` on every effect-list object, valid or not |

Live confirmations behind the table:

- `QuantumState[{{2,0},{0,-1}}]` (Hermitian, not positive, eigenvalues `{2,-1}`) constructs with no message; `"DensityMatrix"` equals the input verbatim, `"Eigenvalues"` returns `{2, -1}`, and `"PhysicalQ"` is available to flag it. The state constructor stores raw input and stays queryable.
- `QuantumOperator[{{2,0},{0,3}}]` (not unitary) constructs; `"UnitaryQ"` returns `False`.
- `QuantumMeasurementOperator[{{0,1},{2,0}}]` (a single non-Hermitian matrix, read as an observable) constructs; its outcomes are the operator's eigenvalues (`{Sqrt[2], -Sqrt[2]}` here), read off the stored operator. No square root is taken; the object is faithful and queryable, exactly like a non-unitary operator.
- `QuantumChannel[{}]` returns `Failure["EmptyKraus", ...]` with the `QuantumChannel::emptyKraus` message; `QuantumChannel[{X}]` builds the channel `rho -> X rho X^dagger` with `"TracePreservingQ"` `True`.
- `QuantumMeasurementOperator[...]["POVMQ"]` returns `True` for a valid projective set, a valid trine POVM, the Hermitian-negative input, the non-Hermitian input, and the nilpotent input alike.

The effect-list form is the only row that stores a derived quantity instead of the input, and the only row whose validity property does not actually test validity.

## 3. The load-bearing distinction: faithful storage versus a partial derived construction

Section 11 of the model says QF "checks shapes, not physics ... it does not check normalization, positivity, unitarity, or trace preservation, which are assumptions of the model, readable afterward but not enforced." Section 1 says "a head names a representation, not a physical object ... whether an object is a valid state, a unitary, or a trace-preserving channel is a separate question the constructor does not ask."

Read against the code, that principle rests on two conditions that hold for states, operators, observables, and channels, and fail for the effect list:

1. **The input is stored faithfully, so every property stays computable.** A non-normalized, non-positive `QuantumState` is a completely well-defined object: its density matrix is exactly what was supplied, its eigenvalues and traces are readable. The physical defect is an attribute of a valid object, not a broken object.
2. **A genuine predicate reads the physical status afterward.** `"PhysicalQ"` runs `PositiveSemidefiniteMatrixQ` on the state; `"UnitaryQ"` tests unitarity; `"TracePreservingQ"` tests `sum M_m^dagger M_m = I`. The user, not the constructor, decides whether to treat the object as physical, and has a working test to do so.

The effect-list `QuantumMeasurementOperator` breaks both:

- It does **not** store the effects. Section 5 states the design: it "takes the `E_m` to be POVM effects and forms each measurement operator as the matrix square root `M_m = sqrt(E_m)`, the canonical choice. The stored measurement operators are these `sqrt(E_m)`." The square root is a nonlinear, derived computation, and it is **partial**: it is undefined on some shape-valid inputs. For the nilpotent effect it has no solution, so `MatrixPower[..., 1/2]` stays unevaluated and the object is half-built. For the Hermitian-negative and non-Hermitian effects it produces a complex or non-Hermitian operator whose recomputed effect `M_m^dagger M_m` is not the effect that was supplied. There is no faithfully-stored object to fall back on.
- Its validity predicate does not work. `"POVMQ"` is a type tag (`qmo["Type"] === "POVM"`), true for any object the effect-list branch produced. It returns `True` on the nilpotent, negative, and non-Hermitian inputs. So even a user who knows to check has no working test on the object.

This is why the case is different in kind from a non-normalized state. The genuine defect is a **partial derived construction that returns an invalid object and leaks a `System` message, with no faithful storage and no working readout to recover from**, rather than a missing physics check on an otherwise valid object. QF already recognizes this category and fails cleanly in it: `QuantumChannel::emptyKraus` guards the operator sum against an empty Kraus set ("a channel's operator-sum runs over a nonempty Kraus set, so an empty list is not a channel"), and the single-operator Kraus reroute guards against a dilation that "would ... leave the output wire on the non-positive environment label and the channel invalid." The non-positive effect is the same situation: a shape-valid input on which the intended construction is undefined or produces an invalid object.

A useful contrast within QF's own recent work sharpens the boundary. The measures audit added `numericStateNotPSDQ` (`QuantumDistance.m:19`) and per-head `::notphysical` messages, and deliberately **did not** guard the `QuantumState` constructor, because a non-positive `QuantumState` is a legitimate intermediate: the partial transpose inside the negativity path is exactly such a state, stored faithfully and used correctly. There, non-positive input is meaningful, the object is valid, and the right move is to warn at the measure, not fail at the constructor. For an effect, non-positive input has no legitimate reading (there is no "partial transpose of an effect" that a downstream computation wants), and the construction is undefined on it, so the right move is the opposite: fail at the constructor, as `emptyKraus` does. The two decisions are consistent, not in tension: warn where the object is valid and the input is a legitimate intermediate; fail where the construction is undefined and no legitimate input produces it.

## 4. Boundary probes

What the effect-list constructor does across the space of inputs, all verified live. "Numbers" means `"ProbabilitiesList"` is a numeric vector; "clean" means no leaked message.

| Input | Hermitian? | PSD? | Root exists? | Result |
|---|---|---|---|---|
| `{P0, P1}` (projective) | yes | yes | yes | clean; `"POVMElements"` `{P0,P1}`; probabilities `{1,0}` |
| `{2 P0, 2 P1}` (PSD, not normalized) | yes | yes | yes | clean; `"POVMElements"` rescaled to `{P0,P1}` (unit mean diagonal); probabilities `{1,0}` |
| trine `(2/3)|psi_k><psi_k|`, sums to `I` | yes | yes | yes | clean; `"POVMElements"` the trine; probabilities `{2/3,1/6,1/6}` |
| `{P0, P0}` (PSD, does not sum to `I`) | yes | yes | yes | clean; accepted; probabilities renormalized to `{1/2,1/2}` |
| `{{-1,0},{0,0}}, {{0,0},{0,1}}` (Hermitian, negative) | yes | no | yes (complex) | silent build; `"POVMElements"` `ComplexInfinity`/`Indeterminate`; probabilities finite `{1,0}` but for a different, substituted POVM |
| `{{1,1},{0,2}}` (non-Hermitian) | no | (WL says yes, see below) | yes (real) | silent build; finite probabilities; `M_m^dagger M_m` is not the supplied matrix, so not a POVM |
| `{{0,1},{0,0}}` (nilpotent) | no | no | no | leaks; unevaluated `MatrixPower`; object invalid |

Three readings from this map:

- **Every valid positive-semidefinite input works, including non-projective and non-normalized ones.** The model's "rescaled so their sum has unit mean diagonal" step (section 5) is real: `{2 P0, 2 P1}` comes back as `{P0, P1}`. A set that does not sum to the identity is accepted and its probabilities are renormalized. These behaviors are shapes-not-physics working as intended and must be preserved: the completeness relation `sum E_m = I` should **not** be enforced.
- **"The square root exists" is strictly weaker than "the effects form a POVM."** The non-Hermitian input `{{1,1},{0,2}}` has a real square root and sails through to finite numbers, but its recomputed effect is not the supplied matrix, so it is not a measurement. A guard that only rejects non-square-rootable input would let this through.
- **The only failures are non-positive-semidefinite inputs**, and they fail in two distinct ways: loudly (nilpotent, leaked messages, unusable object) and silently (Hermitian-negative and non-Hermitian, finite wrong numbers). The silent failures are the more dangerous, and they are exactly the ones a "fix the leak" repair would miss.

### The discriminating predicate

`HermitianMatrixQ[E] && PositiveSemidefiniteMatrixQ[E]` separates the valid inputs from the invalid ones, and `PositiveSemidefiniteMatrixQ` **alone does not**:

| Effect | `HermitianMatrixQ` | `PositiveSemidefiniteMatrixQ` | conjunction |
|---|---|---|---|
| `{{0,1},{0,0}}` nilpotent | False | False | reject |
| `{{-1,0},{0,0}}` Hermitian-negative | True | False | reject |
| `{{1,1},{0,2}}` non-Hermitian | False | **True** | reject |
| `{{1,0},{0,0}}` `P0` | True | True | accept |
| `{{2/3,0},{0,0}}`, trine effect | True | True | accept |

`PositiveSemidefiniteMatrixQ[{{1,1},{0,2}}]` returns `True` because Wolfram Language tests the symmetric (Hermitian) part `(m + m^dagger)/2`, which is positive here, even though the matrix itself is not Hermitian. So the check must be the conjunction: Hermitian first, then positive.

### The symbolic caveat

A parametric effect list constructs correctly today and must keep doing so:

```wl
QuantumMeasurementOperator[{{{p, 0}, {0, 1 - p}}, {{1 - p, 0}, {0, p}}}]
```

builds a valid parametric measurement whose `"POVMElements"` are the supplied matrices and whose probabilities are the expected functions of `p`. But `HermitianMatrixQ[{{p,0},{0,1-p}}]` is `False` and `PositiveSemidefiniteMatrixQ` is `False` for a free `p` (the kernel cannot prove either without knowing `p` is real and in `[0,1]`). A hard predicate guard would wrongly reject this. So the guard must fire only on fully numeric effects (`MatrixQ[E, NumericQ]`) and pass symbolic or parametric ones through, exactly as `numericStateNotPSDQ` does for states.

## 5. Candidate resolutions

### (A) Hermitian-and-positive-semidefinite guard at the constructor, numeric only. Recommended.

For each supplied effect that is a fully numeric matrix, confirm it is Hermitian and its least eigenvalue is at least `-tol` (the `-1.*^-8` tolerance `numericStateNotPSDQ` already uses); leave symbolic and parametric effects unchecked; do not check `sum E_m = I`. On the first failing numeric effect, emit a `QuantumMeasurementOperator::` message and return a bare `Failure`, matching `QuantumChannel::emptyKraus`. This covers both entry points, `:111` (effect matrices) and `:137` (operators).

- Rejects the nilpotent (not Hermitian), the Hermitian-negative (fails the eigenvalue test), and the non-Hermitian (not Hermitian) inputs, the exact three failure cases.
- Admits every valid input in the boundary map, including non-normalized and non-projective POVMs and parametric effect lists.
- Sits in the category QF already fails cleanly in: a domain check on a partial derived construction, not a physics check on a valid object, so it is consistent with sections 1 and 11 rather than in tension with them. The precondition it enforces, `E_m` Hermitian and positive, is precisely the condition under which the model's stated `M_m = sqrt(E_m)` exists as a measurement operator.

### (B) No physics guard; fix only the leak and expose a validity property. Rejected.

The leak is the loudest symptom, not the defect. The Hermitian-negative and non-Hermitian inputs leak nothing at construction and return finite numbers; suppressing `TensorRank::rect` would leave them returning wrong probabilities silently. The "expose a validity property" half is already present and already broken: `"POVMQ"` returns `True` on every invalid input, because it tests the object's type tag, not the effects. So this candidate leaves silent physical corruption in place with no working detector, and is strictly worse than the states case it is modeled on, where the object really is valid and `"PhysicalQ"` really does test it. (Repairing `"POVMQ"` to test the effects would help a diagnostic-minded user, but it does not stop the constructor from returning a broken object, and is better added alongside (A) than in place of it.)

### (C) Guard only that the square root exists. Rejected.

"Square-rootable" is strictly weaker than "valid POVM," proven by the non-Hermitian input `{{1,1},{0,2}}`, which is square-rootable and passes to finite numbers while not being a measurement, and by the Hermitian-negative input, whose complex root exists and yields a substituted positive POVM. This candidate would catch only the loud nilpotent case and let both silent cases through. It fixes the crash and leaves the corruption.

## 6. Literature context

The recommendation is the standard definition of the object, not an added constraint. In the POVM formalism a measurement is a set of effects `E_m` that are positive semidefinite operators (hence Hermitian) summing to the identity, and the canonical, minimal-disturbance measurement operators are `M_m = sqrt(E_m)`; the square root is well-defined precisely because each `E_m` is positive. This is the textbook presentation (Nielsen and Chuang, section 2.2.6). QF's model adopts exactly this in section 5. So requiring `E_m` to be Hermitian and positive before taking the root does not narrow the accepted class of measurements below the standard one; it declines inputs that are not measurements in the first place.

## 7. What is now true

QF's design does not call for a positivity test on effects in the sense section 11 forbids, and it does call for one in the sense sections 1 and 11 permit. The shapes-not-physics rule is a promise to store shape-valid input faithfully and expose the physical status afterward; the effect-list `QuantumMeasurementOperator` cannot keep that promise for a non-positive effect, because it stores `sqrt(E_m)` rather than the effect, and that root either does not exist or misrepresents the input. The genuine defect is a partial derived construction that returns an invalid object with a leaked `System` message and no working `"POVMQ"` readout, which is the category `QuantumChannel::emptyKraus` and the single-operator Kraus reroute already handle by failing cleanly. The clean version is the Hermitian-and-positive-semidefinite guard the mismatch record already names, made precise in three ways this investigation establishes: the predicate is Hermitian **and** positive (positive alone admits the non-Hermitian `{{1,1},{0,2}}`); it fires only on numeric effects, so parametric POVMs still build; and it checks each effect, not the completeness relation, so non-normalized and non-projective POVMs still build. It rejects exactly the nilpotent, Hermitian-negative, and non-Hermitian inputs, and admits everything that is a POVM.
