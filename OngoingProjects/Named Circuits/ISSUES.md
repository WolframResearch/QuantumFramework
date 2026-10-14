# Named Circuits: issues log

Running dump of signature/behavior problems found while exercising the named
`QuantumCircuitOperator` constructions against QF 2.1.3. Each entry: what was
run, what happened, what was expected.

Triage 2026-10-13 (working tree at `ad590d7e`, paclet 2.1.3): each entry carries
a **Verdict**. Confirmed bugs: Association varSpec, BooleanOracleR MapAt,
Deutsch/DeutschJozsa constant, Decohere recursion, Spider recursion, XSpider/ZSpider
color, registry aliases. Not bugs: JX/JY (wrong property read) and GrayOracle
(works as typed). Channel-on-ket pending API intent. No fixes applied.

## Association varSpec rejected across BooleanOracle, PhaseOracle, Grover

**Verdict: BUG (severity 6/8).** Confirmed live: all three return `Failure` with
`Join::incpt` / `Lookup::invrl`; list form builds. Loud failure of a documented
input form.

The association form of the variable spec (`<|var -> wire|>`) fails wherever it
is accepted positionally, with `Join::incpt` / `Lookup::invrl` from
`OrderBooleanVariables` trying `Join[<|...|>, {}]`. The list form works
everywhere. This contradicts the signatures: `PhaseOracle` and `BooleanOracle`
both document an `_Association` varSpec, and `Grover` forwards to them.

Repro (all three return `Failure`):

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && \[FormalX][2], <|\[FormalX][1] -> 1, \[FormalX][2] -> 2|>]}]
QuantumCircuitOperator[{"PhaseOracle"[\[FormalX][1] && \[FormalX][2], <|\[FormalX][1] -> 2, \[FormalX][2] -> 3|>]}]
QuantumCircuitOperator["Grover"[\[FormalX][1] && \[FormalX][2], <|\[FormalX][1] -> 2, \[FormalX][2] -> 3|>]]
```

Works (list form):

```wl
QuantumCircuitOperator[{"BooleanOracle"[\[FormalX][1] && \[FormalX][2], {1, 2}]}]
```

Root cause is in `OrderBooleanVariables` (Kernel/QuantumCircuitOperator/NamedCircuits.m):
the association branch builds `Join[rules, Thread[...]]` where one side is an
association and the other a list. Expected: accept the association and order the
variables by it. (Corrects an earlier note here that claimed `PhaseOracle`
accepts the association, it does not.)

## BooleanOracleR negative-formula branch mangles the rotation gate

**Verdict: BUG (severity 2/8).** Confirmed live: `MapAt::partw` fires yet the call
still returns a `QuantumCircuitOperator`, so it silently builds the wrong oracle
(the negation gate is dropped). Wrong result, only a secondary message as a tell.

When the formula's negation has fewer ESOP terms than the formula (so the
oracle builds from the negation), `BooleanOracleR` emits
`MapAt::partw: Part {2} of RZ[Pi] does not exist` and silently fails to negate
the rotation angle. The rotation gate has the form `"RZ"[angle]` (angle at part
1), but the code negates part 2.

Repro (emits the message):

```wl
QuantumCircuitOperator[{"BooleanOracleR"[\[FormalX][1] || \[FormalX][2] || \[FormalX][3]]}]
QuantumCircuitOperator[{"BooleanOracleR"[Implies[\[FormalX][1], \[FormalX][2]], Automatic, Automatic, "RY"[Pi/2]]}]
```

Clean (formula whose negation is not shorter, e.g. a single product term):

```wl
QuantumCircuitOperator[{"BooleanOracleR"[\[FormalX][1] && \[FormalX][2]]}]
```

Root cause: `MapAt[Minus, rotationGate, {2}]` in the `isNegative` branch
(Kernel/QuantumCircuitOperator/NamedCircuits.m, `BooleanOracleR` rule). Should be
`{1}` (the angle), or `MapAt[Minus, rotationGate, 1]`.

## Deutsch / DeutschJozsa fail on a constant given as True/False

**Verdict: BUG (severity 7/8).** Confirmed live: `Deutsch[True]` and
`DeutschJozsa[True, 3]` return `Failure` via `Range::range` + `Thread::tdlen`;
integer form builds. Loud failure of the constant case.

A constant function passed as `True` or `False` (no Boolean variables) breaks
with `Range::range` and `Thread::tdlen` internal errors instead of building the
constant oracle. The integer-index form works for the constant cases.

```wl
QuantumCircuitOperator[{"Deutsch"[True]}]        (* => errors, FAIL *)
QuantumCircuitOperator[{"DeutschJozsa"[True, 3]}] (* => errors, FAIL *)
QuantumCircuitOperator[{"Deutsch"[1]}]            (* => builds, constant f=0 *)
```

Expected: a variable-free constant should build the constant oracle (identity on
the inputs, with the target set by the constant), or fail with a clear message.

## Decohere[p] blows the recursion limit

**Verdict: BUG (severity 4/8).** Confirmed live: returns
`TerminatedEvaluation[RecursionLimit]` and aborted the enclosing script even under
`CheckAbort`. Loud and uncatchable; takes down the surrounding evaluation.

`QuantumCircuitOperator[{"Decohere"[0.1]}]` (and the same operator generally with
a numeric argument) exceeds `$RecursionLimit` (depth 1024) and returns
`TerminatedEvaluation[RecursionLimit]`, aborting the enclosing evaluation; it is
not catchable by `TimeConstrained`. `Decohere[]` with no argument builds but
returns a non-normalized object (trace 1/2 density matrix). Needs a signature /
base-case fix.

## Channel: flat two-argument form silently returns a wrong result

**Verdict: BUG (severity 3/8), scoped to the flat call form.** (Corrects an
earlier claim here that the named element itself was broken, and the triage note
above that it does not decohere: both used the flat two-argument form. The
documented nested form works.)

The documented nested form decoheres correctly in a circuit:
`Channel["Depolarizing"[3/10]]` on `|+⟩` gives `{{1/2, 7/20}, {7/20, 1/2}}`,
purity 149/200, matching `QuantumChannel["Depolarizing"[3/10]]`. The flat
two-argument form,
`QuantumCircuitOperator[{"Channel"["Depolarizing", 0.3]}][QuantumState["+"]]`,
instead returns a purity-1 state with a density matrix full of
`Missing[NotAvailable, DiscardExtraQudits]` (and a stray `Set::nosym`). Minor:
the flat (wrong-arity) call should error rather than silently produce a
malformed result.

## JX and JY return the Jz matrix (angular momentum broken)

**Verdict: NOT A BUG (wrong property read).** The raw `"Matrix"` is each J-operator
in its own eigenbasis, diagonal by design (QF audit 2.2). Read on
`"MatrixRepresentation"` instead: there `JX =!= JZ` and `[JX, JY] === I JZ` holds
(both confirmed live). The physics is correct; the entry tested the stored basis
form. This entry should be dropped, not fixed.

`QuantumOperator["JX"[j]]`, `QuantumOperator["JY"[j]]`, `QuantumOperator["JZ"[j]]`
all return the SAME diagonal matrix (the Jz matrix, diag(-j..j)). Verified at
j=1: all three equal diag(-1,0,1), and `[JX,JY] = i JZ` fails. Jx and Jy must be
the off-diagonal (tridiagonal) components. JZ is correct; `J+`/`J-` return proper
ladder operators (off-diagonal). Only the x and y Cartesian components are wrong.
This is a correctness bug in the angular-momentum operators, not a documentation
choice, so JX/JY are excluded from the qudits page.

```wl
Normal@QuantumOperator["JX"[1]]["Matrix"] === Normal@QuantumOperator["JZ"[1]]["Matrix"]  (* True; should be False *)
```

## Spider[phase] blows the recursion limit

**Verdict: BUG (severity 5/8).** Confirmed live: returns
`TerminatedEvaluation[RecursionLimit]`; `Spider[]`, `ZSpider[phase]`,
`XSpider[phase]` build. Same loud, uncatchable abort as `Decohere[p]`.

`QuantumCircuitOperator[{"Spider"[Pi/2]}]` exceeds `$RecursionLimit` and returns
`TerminatedEvaluation[RecursionLimit]`, aborting the evaluation (same failure mode
as `Decohere[p]`). `Spider[]` with no argument builds. `ZSpider[phase]` and
`XSpider[phase]` build fine, so the crash is specific to the generic `Spider`
with a phase argument.

## XSpider and ZSpider coincide at one leg each

**Verdict: BUG (severity 1/8), worse than written.** This is not a single-leg doc
limitation. The rule form builds the multi-leg spiders
(`QuantumOperator["XSpider" -> {1} -> {1, 2}]`), and there X and Z still return the
*identical* matrix. `XSpider[π]` is `diag(1, -1) = Z`, not `X`. The X/Z color is
silently collapsed at every arity: a wrong operator returned with no error, the
most dangerous failure mode here.

At the default 1-input/1-output arity, `XSpider[a]` and `ZSpider[a]` return the
identical matrix diag(1, e^{i a}) (checked at a = π and a = π/2). The X/Z color
distinction only appears with multiple legs, but the multi-leg spider signatures
(`ZSpider[phase, inLegs, outLegs]`, `ZSpider[phase, {dims}, {dims}]`, etc.) all
fail with `invalidArgs`, so there is no working way found to build a genuine
multi-leg colored spider from the named form. Documenting the color distinction
is therefore not possible from the single-leg form alone.

## GrayOracle needs a pure function, not a symbol

**Verdict: NOT A BUG (works as typed).** The signature is `_Function | _List`, so a
bare `Sin` is correctly rejected with `invalidArgs`; a pure function or angle list
builds (both confirmed live). At most an error-message usability nit, not a defect.

`GrayOracle`'s first argument is typed `_Function | _List`, so a bare symbol
like `Sin` does not match and the call fails with `invalidArgs`. A pure function
(`#^2 &`) or an angle list works. Minor, but the error message does not hint at
the pure-function requirement.

```wl
QuantumCircuitOperator[{"GrayOracle"[Sin]}]       (* => Failure, invalidArgs *)
QuantumCircuitOperator[{"GrayOracle"[(#^2 &)]}]   (* => builds *)
```

## Registry aliases not bare-callable

**Verdict: BUG / inconsistency (severity 8/8).** Confirmed live: the synonyms are
absent from `$QuantumCircuitOperatorNames`, so a bare string fails `invalidName`
while the `"Name"[...]` call-form builds. Clear error and a workaround exists, so
lowest severity.

The amplitude/operator synonyms (`"GroverAmplification"`, `"GroverAmplification0"`,
`"GroverPhaseAmplification"`, `"GroverPhaseAmplification0"`, `"GroverOperator"`,
`"GroverOperator0"`, `"GroverPhaseOperator"`, `"GroverPhaseOperator0"`,
`"Multiplexor"`) have working dispatch rules but are absent from
`$QuantumCircuitOperatorNames`, so a bare string fails the name guard while the
`"Name"[...]` form builds.

Repro:

```wl
QuantumCircuitOperator["GroverOperator"]          (* => Failure: invalid name *)
QuantumCircuitOperator["GroverOperator"[QuantumOperator["X", {4}]]]  (* => builds *)
```

Decision needed: either add these to the registry (so they are first-class
synonyms) or drop the rules (so there is one name per construction). Right now
the behavior is split and undocumented.
