# Rebuild plan for steps 4 to 6 (2026-11-09)

Source of truth: `design-note.md`. Gap list: `audit-steps-1-6.md`; the IDs below (A4, G2, ...) are its row IDs. Every MISSING, PARTLY and DIFFERENT row appears in exactly one phase, or under "Decisions for Mads" when the design leaves it open. Nothing here is new design: where a line goes beyond the design note, it says so.

## Status and execution to the end (2026-11-09, agreed with Mads)

| Phase | State |
|---|---|
| 0 case set | done; 56 cases after 6 were added for phases 4 to 7 |
| 1 wrong answers, undecided answers | committed d82b6720 |
| 2 one Refine, zero and sign tests, messages | committed f84cca0b |
| 3 assumed signs | committed 4a97d2bb |
| speed against main (decision 10) | committed a6529453, 1895d6e1, 304c43fa |
| merge of main | committed 5fe4f845 |
| 4 unitarity from structure, carried norm | built 2026-11-10, not committed (see below) |
| 5 facts reach monotones, rates, phase-space forms, distances | committed 1895d6e1, 6206c351 |
| 6 rotations unsimplified when built | committed 7efc14a0; other unasked Simplify not changed (see below) |
| 7 eigensystems on the facts | committed 7efc14a0 |
| fixes and speed after phase 7 | committed 307130c5 |
| end | full suite and timing against main, then the report |

How it runs: no questions to Mads. Each phase: tests first, build, the symbol-domain tests only, a local commit (no push, no PR). At the end: merge main, the full suite and the phase 0 cases once; a failure is traced to its phase by testing the per-phase commits. The adversarial verification loop is not run (it needs Mads's yes); the report says so. The result is one private report page: what each phase changed with an example, the 50 cases main against the branch (answer, size, time), the suite count, every default chosen and every departure from this plan, and what is left.

Defaults chosen for the open points:
- Phase 4: unitarity is read from structure at the time it is asked (named gate with real angles, dagger, products, tensor products, controlled gates of unitaries) before any entry check; the carried norm is kept in a table keyed by the state, filled when a unitary is applied, read first by "Norm", "Weights" and probabilities, and never stored in the basis (a basis that differs only by a norm would no longer compare equal to its twin).
- Phase 5: every conjugation site listed goes through `conjugateOnDomains`; every `ComplexExpand` site passes the record's complex symbols; Lindblad rates get the role [0, Infinity); the entanglement monotones read the record (decision 2).
- Phase 6: each unasked `Simplify` or `FullSimplify` becomes one helper, `Refine` on the facts then a time-limited `Simplify` kept only when smaller; build-time Simplify of rotations removed (decision 1).
- Phase 7: `Eigensystem` with QF's zero test under the facts; triangular and block matrices split first; a count of eigenvectors that differs from a multiplicity gives Indeterminate; no conditions at special values (decision 6).

### Phases 4 to 7 and the speed work, as built (2026-11-10)

Speed (decision 10). Measured against main with alternating runs (fastest of 7 per case, main's spread as the allowance). What the branch had added to every call, and what replaced it:
- the yes/no message scope opened on every property call: now only on the five yes/no properties;
- every basis carried an empty "SymbolDomains" entry and every merge read records through property calls: the entry is kept only when it holds a fact, and records are read from the data;
- three new basis rules were tried on every basis: one rule now, which exits at once on an empty record; QF's own records skip the checks a user's declaration needs; a record's and a symbol's validity are checked once;
- the basis validity check used assertions: now a plain chain of tests, faster than main;
- state properties tested the state before the property name, so a state with facts paid on every property: the name is matched first;
- "UnitaryQ" and "HermitianQ" built the matrix twice: once now;
- PhysicalQ of rational real matrices: no Hermitian test needed, one exact Reduce for the signs;
- probabilities: weights factor by factor (|e^y|^2 = e^(2 Re y) on the facts), conjugation read off the structure without Refine, the norm's 0 and 1 coefficients by an exact zero test, the norm reused from the weights;
- Refine runs only for what the structure leaves (conjugates of real expressions, square roots of squares and Abs of a known sign are settled first).

Phase 4, first measured and left out; built after a stricter timing (2026-11-10). The earlier claim that no case was slower compared the branch with main's noisiest run; against main's quietest run two cases were slower, the 3-qubit "ProbabilitiesList" (about 20%) and "HermitianQ" with undeclared symbols (about 40%). As built:
- facts proven about one object are kept with that object, keyed on the object itself (an expression store), never in its basis: a True or False yes/no answer with its messages, and "UnitaryQ" for a gate unitary by construction;
- a named qubit gate built in a circuit with arguments real on its facts is marked unitary as built (a fixed list of names, each checked against its matrix by a test); a label is never trusted, since a user can give any matrix the label "H" and "RX" on a qutrit carries the qubit label but is not unitary; the dagger, the inverse and moving to other wires keep the mark;
- a circuit whose operators are all marked keeps the norm of the state vector it acts on (computational bases only); the output's norm and probabilities read it from the input, so nothing is expanded;
- a yes/no check notes the entries it could not cancel, and the "true when real" message is read from them (Conjugate[s] -> s), so the check no longer runs twice;
- an entry that is a polynomial in symbols with no fact and their conjugates is decided by expanding it, before any point is tried;
- the symbols of an expression are listed by first collecting its atoms (the same list, about 5 times faster), which was most of the cost of "Weights";
- the norm's size limit (LeafCount <= 4 n + 8, set from timings) is replaced by a rule: a sum of numbers times squares is kept, anything else is simplified on the facts. $weightSimplifyLeafCount stays as a named setting.
Not built: marks for gates built directly with QuantumOperator, for qudit gates, and for "Fourier", "InverseFourier", "Fredkin", "Toffoli" and "Deutsch", which a circuit builds through other paths; these are checked from their matrices as before.

Phase 5, built: monotones refined on the state's facts (negativity of Werner[p] holds no Conjugate[p]); jump rates in [0, Infinity); the Wigner, Husimi and s-ordered functions in x and p treat as real only x, p and the facts' real symbols (the Wigner function of |0> + z|1> at z = 1 + I was complex, a wrong answer on main too); distances are refined on the joined facts, and the Bloch vector reads a state already in the computational basis as it is (Bloch vectors and trace distances of real-angle states hold no Conjugate). "Computational" itself still drops declared parameters, as on main: keeping them made a parameter bound to a complex value after an interferometer give a wrong state (QuantumInterferometer test "QI-loss-state-parameters-survive-and-bind-after-propagation"). Not built: the scattered conjugation copies in stabilizer, optimization and second-quantization code, and QuantumEvolution's discontinuity search, since no case shows a wrong or messy answer there.

Phase 6, built: no rotation is simplified when built (decision 1); the RZPhases and RotationWires tests changed to the stored, unsimplified form. Not built: replacing the other unasked Simplify calls, since none shows in the case set as slower than main or wrong.

Phase 7, built: eigensystems solved on the matrix refined on its facts ({{1, Sqrt[p^2] - p}, {0, 1}} with p >= 0 gave a zero eigenvector, on main too). Not built: the multiplicity check giving Indeterminate, since a defective matrix legitimately has fewer eigenvectors and QF's own matrices gave none wrong in the cases.

Tests changed by decisions, not by regressions: DeclaredReal-NegativeParameter-Sqrt (decision 9: Sqrt[theta] of a bare parameter now assumes theta >= 0), DeclaredReal-FormalSymbol-StaysComplex (weight written |r|^2 Cos^2 instead of |r Cos|^2, the same value), RZPhases and RotationWires (decision 1).

## What the work must achieve

1. **Right answers on symbolic objects.** A gate's dagger, a state's norm and probabilities, and the yes/no properties ("UnitaryQ", "HermitianQ", "PhysicalQ", "TracePreservingQ") are correct for the values the facts allow, in short form.
2. **Honest answers.** Each yes/no property gives True, False, or `Indeterminate` with a message saying why; a result that rests on a fact the record lacks says which fact.
3. **No cost to numeric work, bounded cost to symbolic work.** Numeric objects run as on main. A symbolic operation is never much slower than on main and never hangs: every costly step has a time and size limit.
4. **Nothing surprising.** User symbols are never changed, private symbols never escape, and the facts a user states win over the defaults.

## Questions the design must answer

Weakest premise first: the questions that the most later work rests on, and that are least checked, come first. Each phase below names the questions it answers; a phase that answers none is suspect. Answers are recorded here, in QF's own terms, as they are settled.

1. **Which cases decide success, and by what measure?** A fixed set of QF's own symbolic cases (from `bench.wls`, `../SymbolicForQuantum/bench/cases.wl` and `cases_hard.wl`, and the tests), each with main's output and time, and a target: correct, no longer than a stated size, and no slower than a stated multiple of main. Every later question is judged on this set. *Not yet answered.*
2. **What must QF add on top of `Refine` and `Simplify` run under one fixed fact expression?** `Refine` already applies sign facts factor by factor. Phase 2's cleanup function is worth building only for what that leaves undone or too slow on the case set. *Not yet answered; decides the scope of phase 2.*
3. **Which facts does QF carry, and how do they pass through operations?** Facts about symbols only, or also about operators ("unitary for every allowed value", "Hermitian"), passed through products, tensor products, controlled gates, direct sums and exponentials of Hermitian generators. *Partly answered: phase 4 takes unitarity from structure first.*
4. **What does a symbol with no recorded fact mean?** Today it is complex and unconstrained, so many questions stay undecided. A solver instead answers for generic values and reports conditions only on request. Which default, and when a condition is reported, decides how often QF answers at all. *Not yet answered.*
5. **What counts as proof?** An exact form proves; a value at an allowed point disproves "zero everywhere" but never proves "nonzero everywhere"; a timeout proves nothing. *Answered: phase 2 zero tests, and decision 5 (`Indeterminate` with a message, naming where the answer holds when the answer depends on the values).*
6. **How is cost bounded?** Per-stage caps inside one budget per operation, size estimated before a costly step, no cache entry from a timeout, nothing done twice in one operation. *Answered in outline (phase 2 revisions); the limits are set on the case set of question 1.*
7. **Which quantities are carried instead of recomputed, and when do they go stale?** The norm through unitary operators; unitarity from structure; a scalar factor through `Abs`. *Answered in outline (phase 4).*
8. **How is a substituted value checked against the facts?** Exact, inexact and symbolic values; only a violation gives a message (step 7). *Answered in outline in `design-note.md`, not built.*
9. **How do eigensystems behave under the facts?** Structure first, multiplicity checked, generic values by default (decision 6). *Answered in outline (phase 7).*

## Decisions on the answers (Mads, 2026-11-09)

- **Question 2, accepted: `Refine` early.** `Refine` runs on amplitudes and matrix entries as they are conjugated, not on finished results: at 8 qubits it takes 28 ms on the amplitudes against 600 ms on the final probabilities (7.4 million leaves), and `Simplify` on top 10 s (`phase0/scale.wls`). `Simplify` runs only where a short form is the point and the result is below a size limit. The collected complex-amplitude norm stays. The private real copies and the per-entry sign settling of step 5 go.
- **Question 4, construction steps, accepted: generic values with one message.** A construction step treats an expression that is not identically zero as nonzero, says once that a special-value condition was dropped, and an option keeps the conditions instead.
- **Speed target and operator facts, accepted.** Numeric work no slower than main beyond noise; a symbolic case at most about 25% slower than main, and only when its result is shorter or more correct. Only two operator facts, unitary and Hermitian, derived from structure when an operator is built; any operation without a rule drops them; each rule is tested, including the cases that must drop the fact (a sum of unitaries, multiplying by I).
- **Question 4, undeclared symbols, accepted: a complex number, with three additions.** Entries of QF arrays are already treated as numbers by QF and by `ConjugateTranspose` (`decisions/scalar.wls`); real by default would make |α|^2 into α^2 for an amplitude, silently. (1) When a yes/no answer is `Indeterminate` only because symbols lack facts, its message names them ("declare t, Δ, Ω real to decide this"). (2) One option marks every symbol of an object real (`"SymbolDomains" -> Reals`), expanded into the symbols present when the object is built, so no pattern is ever stored. (3) An entry that is not a number at all (a noncommutative product, a matrix symbol, a QF object) gets no facts, and a yes/no property on it is `Indeterminate` with a message naming the entry; today `QuantumOperator[{{a ** b, 0}, {0, 1}}]` gives the dagger entry `Conjugate[a ** b]` silently. Roles are never guessed from how an object is built. Goes in phase 1 with the record changes.

## Answers from phase 0 and the second study (2026-11-09, proposed, awaiting Mads)

Measured with `phase0/` (cases, `run.wls`, `plain.wls`; results `main.tsv`, `branch.tsv`, `plain.tsv`), main at `a3dc4e2c`, branch at `bd72b1e2`.

1. **Question 1.** The case set is the 50 cases of `phase0/cases.wl`. Each value result is checked against a plain linear-algebra computation of the same quantity at an allowed point: all pass on main and on the branch. Target per case: correct; no `Conjugate`, `Abs`, `Re`, `Im`, `Arg` or `Sign` that the facts remove; no larger than the smaller of the branch's result and plain `Refine` then `Simplify`; no slower than main beyond noise, numeric cases unchanged.
2. **Question 2.** One `Refine` under a fixed fact expression, applied to main's raw output, gives exactly the branch's result (same size, no hard heads) on every conjugation, density-matrix and channel case, in 0.1 to 90 ms. `Refine` then `Simplify` gives the branch's norm 1 and its probabilities on every real-angle case, and shorter results than the branch on reduced states, purity, Bloch vectors and negativity, in 10 to 360 ms. `Simplify` alone is slower and once hits its 20 s limit. Neither does what the branch's own norm does with complex amplitudes (3 qubits: 45 leaves against 1,234), and neither decides a yes/no property. So QF adds only: the fact expression and one `Refine` per result; `Simplify` after it where a short form is the point (norm, probabilities, reduced states); the norm and probabilities collected by complex amplitudes; and the yes/no tests (zero test, sign provers). The private real copies and the per-entry sign settling are replaced by `Refine` for displayed results. Phase 2 shrinks to this.
3. **Question 3.** QF carries facts about symbols and, computed once when an operator is built or composed, about operators: unitary for every allowed value, Hermitian. Rules: `Exp` of an anti-Hermitian generator is unitary (`-I θ H` with θ real and H Hermitian is anti-Hermitian); products, tensor products, controlled gates and direct sums of unitaries are unitary; `c U` is unitary only for |c| = 1 and `c H` Hermitian only for real c; sums keep Hermitian, never unitary; `U` next to `U^dagger` cancels. A matrix symbol cannot be declared unitary (`design-questions.md`, answer 8), so these rules are QF's own.
4. **Question 4.** A symbol with no recorded fact stays complex; a symbol given an inequality is real (a complex mark together with an inequality is a contradiction, reported once when the record is built). Two rules: a yes/no property must hold for every allowed value (otherwise `Indeterminate`, decision 5); a construction step (dividing by a norm, normalizing, separating eigenvalues) treats an expression that is not identically zero as nonzero, and drops conditions that hold only at special values (equalities, x != c on a continuous symbol), keeping inequalities and every condition on an integer symbol. One option with levels (drop special-value conditions, keep them all) controls phases 3 and 7 alike.

Further changes from the study, to fold into the phases when each is built:
- `Conjugate` passes through `Sqrt`, `Log`, non-integer powers and inverse trig only when the argument is proven off the branch cut (`Sqrt`: argument >= 0; `Log`: > 0; `ArcSin`, `ArcCos`: in [-1, 1]), never merely real.
- Trig zero tests use c = Cos, s = Sin with c^2 + s^2 = 1 (after reducing angles to integer combinations of base angles), not the half-angle tangent, which misses θ = π.
- Sample points are accepted only where the expression is defined (denominators, `Log`, `Sqrt` and power arguments, `Abs` and piecewise boundaries checked exactly); a bounded integer symbol with at most 100 values is enumerated, which proves; closed ends of a range are tested as separate points.
- The sign entry answers positive, negative, varies (with two witness points of opposite sign, which also gives `PhysicalQ` its False), or unknown with a reason; samples choose which direction to prove; a one-symbol rational expression is decided exactly from end values and critical points.
- A stage result that brings in `Indeterminate` or an infinity is refused; terms are not cancelled across a denominator not proven nonzero.
- Messages from nested QF calls are collected and issued once by the outermost operation, naming the subexpression in the user's symbols.
- Internal `Reduce` and `Solve` calls pass their options explicitly, take only the facts linked to the question's symbols, and turn "cannot decide" messages into unknown.

## How each phase closes

1. Tests first: each item gets a test that fails on the current branch, written and run before the code changes.
2. The code.
3. A checklist in this file: each audit ID of the phase marked built, with file:line, or not built, with the reason.
4. Tests that check which stage answered, wherever an item is a stage (E9).
5. Full suite; timing against main on `bench.wls` plus the phase's own cases. Each cleanup stage can be switched off, and a stage stays only if the timing shows it helps (I1). The tests of each phase are shown to fail on main as well as on the current branch (C5).
6. Commit only on Mads's yes. `/wl-verify` and `/wl-quality` only when he asks.

## Phase 0: the case set and the experiment behind phase 2

Answers questions 1 and 2. No kernel change.

1. Gather the cases of question 1 into one file, each run on main (`a3dc4e2c` or later) and on this branch: output, `LeafCount`, time (minimum of 5 runs, `ClearSystemCache[]`, property cache off), and whether the output is correct for the allowed values.
2. Set the target per case: correct, size no larger than a stated bound, time no more than a stated multiple of main.
3. On the same cases, run plain `Refine` and `Simplify` under one fixed fact expression per object, with `TimeConstraint` set, and record what they leave undone or do too slowly. That list is the scope of phase 2.

## Phase 1: the two wrong answers

Answers question 5.

| ID | Fix | Test that fails today |
|---|---|---|
| E6, G6a | "PhysicalQ" of a symbolic density matrix with declared "Parameters" runs the same positivity check as one without; the old `! UnknownQ` answer stays only for numeric and non-matrix states | diag(p, 2 - p) with p declared on [0, 3] gives False; with p on [0, 2] gives True |
| G6b, G5 | "Type" asks purity on the matrix after stage 1 (private real copies); an undecided purity is not read as "Mixed" | pure rho written with Conjugate[Sin[t]], t real, gives "Pure"; an undecided case gives neither "Pure" nor "Mixed" |
| G6c | summary boxes show Undefined, not `$Aborted`, when the symbolic check times out | slow symbolic "UnitaryQ" in the box shows Undefined |
| B3, K | basis "Options" lists "SymbolDomains" only when not empty (decision 3) | `QuantumBasis[]["Options"]` has no "SymbolDomains"; ApplyFixedCost's three option tests back to main's form |

### Phase 1 checklist (built 2026-11-09, not committed)

Tests: 11 new or changed in `Tests/SymbolDomains.wlt`, each failing on the code before phase 1; `Tests/ApplyFixedCost.wlt` back to main's option expectations. Kernel paths relative to `QuantumFramework/Kernel/`.

| Item | Status | Where |
|---|---|---|
| E6, G6a: "PhysicalQ" with declared parameters checked against the facts | built | `QuantumState/Properties.m:186` (`symbolicDensityStateQ`, no arity condition); test PhysicalQ-DeclaredParameters |
| G6b: "Type" asks purity after the real symbols are taken real; undecided purity is neither Pure nor Mixed | built | `QuantumState/Properties.m:196, 775`; test Type-PureWithConjugate |
| G6c: summary boxes show `Indeterminate`, with no message | built | the five `Formatting.m` files (QuantumOperator, QuantumChannel, QuantumMeasurementOperator, QuantumHamiltonianOperator) |
| Decision 5: `Indeterminate` with a message per reason (time, tests, fails, region) | built | `QuantumBasis/QuantumBasis.m:656` (`zeroVerdict`), 717 (`whereNonnegative`), 773 (`yesNo`), 796 (`reportVerdict`); test PhysicalQ-Split |
| Decision 5: the region named when a time-limited `Reduce` finds it | built | `QuantumBasis/QuantumBasis.m:717` |
| Decision 5: an `Indeterminate` is never cached | built | `Utilities.m:1690` |
| Decision 5: internal uses go through `TrueQ` and say nothing | built | `QuantumState/Properties.m:328, 775, 1190` |
| Decision 7: a user's plain Reals or Complexes replaces a role's range; an interval meets it | built | `QuantumBasis/QuantumBasis.m:170, 193`; test ExplicitDomainReplacesRole |
| Question 4: complex mark beside a range or integer domain is a contradiction | built | `QuantumBasis/QuantumBasis.m:68, 187`; test ComplexWithRange |
| Question 4: hint naming symbols to declare | built, for False answers | `QuantumBasis/QuantumBasis.m:787`; test Undeclared-Hint |
| Question 4: `"SymbolDomains" -> Reals` for every symbol of an object, met with role ranges | built | `QuantumBasis/QuantumBasis.m:200, 456`, `QuantumState/QuantumState.m:124`; test AllReal |
| Question 4: an entry that is not a number leaves yes/no answers undecided | built for yes/no answers; conjugation unchanged until phase 2 | `QuantumBasis/QuantumBasis.m:773`; test NotANumber |
| B3, K: basis "Options" omit an empty "SymbolDomains" | built | `QuantumBasis/Properties.m:66`; test OptionsHideEmpty |

Two refinements found while building, confirmed by Mads (2026-11-09):
- A yes/no property that is an equality (unitary, Hermitian, unital, trace preserving) answers False when it fails at an allowed point, since an expression that is not identically zero vanishes only at special values; `Indeterminate` with a region applies to the inequality properties ("PhysicalQ"). So `QuantumOperator[MatrixExp[-I t (Δ Z + Ω X)]]["UnitaryQ"]` is False, with the message that it is True when t, Δ, Ω are real.
- Messages are issued by the property asked; collecting them once per outermost operation stays with phase 2.

## Phase 2: the cleanup function (step 5 as designed)

Answers questions 2, 5 and 6. Its scope waits on question 2.

One function, called once per user operation on the whole result. Every existing caller (conjugation, "HermitianQ", "UnitaryQ", "UnitalQ", "PhysicalQ", weights, norm) goes through it, so one "ProbabilitiesList" call conjugates once, not twice (G1).

| ID | Item |
|---|---|
| G1 | the single function and its callers |
| G2 | split a matrix into blocks of rows and columns that never meet through a nonzero entry; clean each block; one cache keyed on expression and facts, shared by every entry, dropped when the call ends |
| F2, G3 | head scan before any work: no Conjugate, Re, Im, Abs, Arg, Sign, Sqrt or Log of a recorded symbol means no swap, no block, no Refine; the scan decides which stages run |
| G3, F4, E8 | every stage has a size limit and a time limit, keeps its result only if no larger; one time budget per user operation, nested calls get what is left; every Simplify in these paths limited, returning its input on timeout, never `$Failed` (the entanglement callers included) |
| G3.1, E0 | stage 1, private real copies, used for every computation in these paths, replacing `ComplexExpand` plus `Refine` in `simplifyOnReals`, `refineOnParameters` and the sign proof |
| G3.2, E4, F3, E2 | stage 2: the `Sqrt` rule (`Sqrt[p^2]` becomes `Abs[p]`, `Abs[p]` becomes p only for p known nonnegative); each sign question asked for the claim and its negation; three cheap provers (intervals, corner evaluation for expressions linear in each symbol, a rule per head); `Refine` on what is left inside one `$Assumptions` block per operation, holding only the facts of the symbols present, on terms with Conjugate, Abs, Sqrt or Log |
| F9 | the private copy of a positive symbol carries `Positive[q] = True` |
| G3.3 | stage 3: common denominator, numerator and denominator cleaned separately |
| G3.4, F10 | stage 4: `FullSimplify` on pieces below a size limit, with a time limit and a short `TransformationFunctions` list (trig of real angles, square roots of probabilities, rational functions); KAK's open-ended `FullSimplify` uses it |
| G4 | on timeout, the best form found so far for each piece |
| E6, F6 | zero-test sample points: distinct values inside the record, angles over a full period, the failure computed exactly, not at 30 digits against 10^-10; each cached answer records which test produced it |
| F11, G5 | `PossibleZeroQ` only where the exact test does not apply, its answer marked likely, never replayed as proven |
| E10 | test that no output of any tested property holds a private or removed symbol, not only that no names are left |
| E9 | tests that name the stage: intervals versus Refine in the sign work, intervals versus Simplify in the nonnegative proof, each cleanup stage |

### Phase 2a checklist (built 2026-11-09, not committed)

| Item | Status | Where |
|---|---|---|
| Question 2: one `Refine` under the basis's fixed fact expression replaces the private real copies and the per-entry sign settling | built | `QuantumBasis/QuantumBasis.m` (`refineOnFacts`, `conjugateOnDomains`, `conjugateTransposeOnDomains`) |
| `Refine` on amplitudes and entries as they are formed, with its own `TimeConstraint`, inside a block replacing the user's `$Assumptions` | built | same |
| Conjugate distributed over sums, products and integer powers before `Refine` (the private copies did this silently; the collected norm needs it) | built | `conjugatedEntries`, `conjugateDistributed` |
| Each matrix refined before it is compared with its conjugate (Hermitian, unitary, semidefinite checks; purity) | built | `QuantumOperator/Properties.m` (`hermitianVerdict`, `unitaryVerdict`), `positiveSemidefiniteVerdict`, `symbolicPureQ` |
| Same results as phase 1 | checked | every phase 0 case passes with the same size; times within noise (`phase0/phase2a.tsv`); full suite 4815 pass |

Not in 2a: `refineOnParameters` (the weights of amplitudes holding complex symbols) still uses `simplifyOnReals`; tests pin its output forms, so it moves with 2b.

### Phase 2b checklist (built 2026-11-09, not committed)

| Item | Status | Where |
|---|---|---|
| Two questions kept apart: zero for every allowed value (`zeroVerdict`), nonzero for every allowed value (`nonzeroVerdict`, only a proof says yes, a zero at an allowed point says no) | built | `QuantumBasis/QuantumBasis.m` |
| Zero test: proof on the entry as given first; then, for entries left, one allowed point refutes; then the full test: `Refine` on the facts, a product zero when a factor is, a positive power when its base is, the exponential form over one denominator, and trig of linear angles in c = Cos, s = Sin of base angles reduced with s^2 = 1 - c^2 | built | `zeroTest`, `provenZeroTest`, `trigZeroQ` |
| `PossibleZeroQ` not used: its "likely" never changed an answer, and on exact values it raised precision messages | dropped | |
| Sample points seeded by the question (expression and ranges); five per symbol over equal parts of its range ([-3, 3] when unbounded, so an angle covers a period), each finite end with the others inside, both corners; a point is kept only if the facts hold there and every denominator, `Log`, `Abs`, `Sign` and non-integer power argument is nonzero there | built | `questionPoints`, `allowedPoints`, `endPoints`, `valueAt` |
| Values at points computed with the coordinates at 40 digits and the error tracked (no 10^-10 threshold; an exact point made the kernel stall on quantities exactly 0 or 1); a zero reading counts as a zero only when it reads zero again at 80 digits | built | `valueAt`, `numberSign`, `pointSign` |
| Bounded integer symbols with at most 100 values in total are taken at every value, which proves | built | `questionPoints` |
| Sign entry: positive, nonnegative, negative, nonpositive, zero, varies with two allowed points of opposite sign, not real with a point, or undecided; interval arithmetic first (also on the factored form of a rational expression), then the samples choose the direction, then `Reduce` for one-symbol algebraic expressions, then `Simplify` | built | `signTest`, `intervalSign`, `provenSign`, `signProof` |
| Nonnegativity: negative everywhere gives False; negative, varying or not real somewhere sends the question to `Reduce` for the region | built | `nonnegativeVerdict`, `whereNonnegative` |
| `Reduce` and `Simplify` read only the facts given: `$Assumptions` replaced, every `Reduce` option passed, only the facts of the question's symbols, `Reduce` only on algebraic input (so no "cannot solve" message) | built | `reduceOnFacts`, `simplifyOnFacts`, `linkedFacts`, `algebraicQ` |
| A cleanup result that brings in `Indeterminate` or an infinity is refused | built | `keptStage`, used by `refineOnFacts` and the weights |
| Weights off `simplifyOnReals`: `Refine`, then \|x\|^2 = x Conjugate[x] (and Re, Im likewise) with the conjugate on the facts, simplified under them | built | `QuantumState/Properties.m` (`valueOnFacts`) |
| Tests | 6 new, 1 changed (Log[Exp[x]] - x is now proven zero; Gamma[x + 1] - x Gamma[x] stands for the undecided case) | `Tests/SymbolDomains.wlt`: ZeroOnFacts, NonzeroEverywhere, SignAnswer, IgnoresUserAssumptions, WeightOfSumOnFacts, StageRefusesIndeterminate |
| Bug found: a user's global `$Assumptions` changed `PhysicalQ` (diag(p, 2 - p), p in [0, 3], gave True under p <= 2) and the wrong True was cached | fixed | IgnoresUserAssumptions |
| Same results as 2a on the phase 0 cases, no case over 25% slower | checked | `phase0/phase2b.tsv` |

Not in 2b: the rule per head for signs (`Log`, odd powers, `Sin`/`Cos` on half-periods) and the `Conjugate` branch-cut rules, since interval arithmetic, `Reduce` and `Simplify` decide the whole test set; "no cancelling across a denominator not proven nonzero", which concerns cleanup results, not the zero proof. An empty record (an integer symbol between 1/3 and 2/3) was already caught when the record is built.


### Phase 2c checklist (built 2026-11-09, not committed)

| Item | Status | Where |
|---|---|---|
| The outermost property call on a QF object owns the messages of every check run for it: inside it a message is recorded, and on return each distinct message is issued once, under the head of the object the user asked | built | `verdictScope`, `verdictMessage` in `QuantumBasis/QuantumBasis.m`; the dispatchers of QuantumState, QuantumOperator, QuantumChannel, QuantumMeasurementOperator, QuantumCircuitOperator, QuantumMeasurement, QuantumHamiltonianOperator |
| Measured first: no operation repeated a message (internal checks were already silent); the visible change is that a circuit's or measurement operator's message names that object | checked | `phase2b/msgs.wls` |
| Tests | 3 new | MessageUnderCircuit, MessageUnderMeasurementOperator, MessagesOncePerOperation |
| A check answered and cached gives no message on later calls, as before | kept | |

### Phase 2 closed

Full suite 4824 pass, 0 fail (`phase2c-suite.txt`); every phase 0 case gives the same result as 2a, none more than 25% slower (`phase0/phase2c.tsv`). 2a, 2b and 2c await one commit.

## Phase 3: conditional answers

Answers question 4 (when a condition is reported).

| ID | Item |
|---|---|
| E7, F5 | when a step needs a sign the record lacks, return `ConditionalExpression[result, fact]` with the smallest such fact, attached to the whole object, and warn. Test: `{Sqrt[1 - q], Sqrt[q]}` with q declared but unbounded gives the density matrix under 0 <= q <= 1 and a message |


### Phase 3 checklist (built 2026-11-09, not committed)

Decisions (Mads, 2026-11-09): 8. An assumed sign lives on the object, as a fact; values stay plain, never wrapped in `ConditionalExpression`. 9. A symbol with any stated domain (the user's, a role's, `"SymbolDomains" -> Reals`) is never given an assumed sign, so decision 7 holds; only a symbol named in `"Parameters"` with no domain gets one.

| Item | Status | Where |
|---|---|---|
| When a state or operator is built, a square root, other non-integer power or log of a real expression in one bare parameter, whose sign the facts do not give, gets one: argument >= 0 (> 0 under a log or negative power) when the facts allow, else the opposite; the range this leaves joins the record; one message names the signs and the condition | built | `assumedSignRecord`, `assumeSign`, `conditionInterval` in `QuantumBasis/QuantumBasis.m`; constructor hook in `QuantumState/QuantumState.m` |
| The condition travels with the facts (tensor products, operators, every merge) and later steps read it as a fact | built | the record |
| A half-line (q <= 1, g >= 0) is accepted as a domain, since one sign leaves one | built | `symbolDomainQ`, `domainFact` |
| Example: `{Sqrt[1 - q], Sqrt[q]}` with q a parameter gives q in [0, 1], density matrix with 1 - q and q, norm 1, probabilities {1 - q, q} | checked | AssumedSign |
| Tests | 3 new, 1 changed (Sqrt[x] with x explicitly Reals keeps Conjugate) | AssumedSign, NothingAssumed, AssumedSignTravels, ConjugateKeptWithoutFacts |

Differs from the plan as written: the sign is chosen nonnegative first when the facts allow, not from a sample point, so the choice is predictable; a sign on two symbols at once is not assumed (no range of one symbol states it); `Abs` and `Log` met later in a computation, not in the amplitudes, get no assumption; the option with levels waits for phase 7, where its second level is needed.

## Phase 4: the carried norm (design 3 to 5, order of work 3 and 4)

Answers questions 3 and 7.

| ID | Item |
|---|---|
| A3 | named rotations answer "UnitaryQ" from their name when their angles are real; others by the symbolic check, asked once per operator even for parametric operators (stored on the object, since the property cache skips them) |
| B1 | the norm lives in the state's data; every construction path keeps it or drops it |
| A4, C3 | applying an operator unitary on the domain passes the input norm to the output, through both the Schrodinger and the tensor-network paths |
| A5, C4 | probabilities read the carried norm; "Weight" reuses "Norm" instead of recomputing; the recomputation is retired wherever a carried value exists. Test: a probe shows no norm recomputation on a unitary circuit |

## Phase 5: remaining step 4 items

Answers question 3 (where the facts reach).

| ID | Item |
|---|---|
| E1 | Lindblad rates recorded nonnegative in `QuantumEvolve` and the Lindblad named operators |
| D2, E3 | `projector` and the scattered conjugation copies (QuantumDistance, InnerProduct, QuantumEntanglement, StabilizerFrame, QuantumOptimization, SecondQuantization/Utils, Stabilizer/Constructors) go through the cleanup function's conjugation |
| E11, F13 | the six `ComplexExpand` calls that take every symbol real pass the record's complex symbols (PhaseSpaceRepresentations 171, 216, 321; QuantumMeasurement 51; QuantumOptimization 410; QuantumEvolution 395); KAK's call too |
| A2 | after phases 2 and 5 the norm, weights and sign work use one helper, not five routes |
| C2 | the entanglement monotones (concurrence, negativity and the rest) read the object's record; `QuantumEntangledQ`, KAK and the stabilizer phase check stay as they are (decision 2). Test: negativity of "Werner"[p] holds no Conjugate, Re or Ramp of p |

## Phase 6: step 6 redone

Answers question 6.

| ID | Item |
|---|---|
| H2 | every unasked Simplify or FullSimplify becomes the cleanup function without its last stage |
| H1 | all of inventory group 3, not 4 of about 20: eigensystem wrapper, "Bipartition", measurement probabilities, distribution, entropy and amplitudes, measurement-operator eigenbases, entanglement monotones, stabilizer zero tests, Multiway, phase-space and Wigner transforms, "Werner", CHSH circuits, bosonic algebra |
| H6 | the Schrodinger route's final `Simplify` becomes the cleanup function, with its time limit |
| H5 | exponentials of a sum of noncommuting terms (Rabi) go through the cleanup function |
| H3 | RZ phases and qudit RX, RY, RZ lose their build-time Simplify (decision 1) |

## Phase 7: eigensystems

Answers question 9.

| ID | Item |
|---|---|
| F7 | triangular and block-diagonal cases first; for a Hermitian matrix the eigenvector count must equal each multiplicity, else undecided; coincidences at special parameter values returned as conditions |

## Revisions after study (2026-11-09)

These refine the phases above; where they disagree, this section holds. Kernel behavior they rely on was checked with `decisions/claims.wls`.

**Phase 1**
- The undecided answer cannot be `Undefined`: it passes through `And` and `Or` (`Undefined && True` is `Undefined`) and reads as a false condition, not as unknown. It becomes `Indeterminate` with a message (decision 5).
- An explicit `Reals` or `Complexes` from the user drops the role's bounds (decision 7).

**Phase 2, the facts and Refine**
- One fact expression per object, built once from the record in a fixed order, written only as `Element[s, Reals | Integers]` and plain inequalities, and reused unchanged by every operation on that object; inside its block, `Refine` and `Simplify` are called without a second argument. This replaces "only the facts of the symbols present", which would rebuild the prepared facts on every call.
- A symbol the user marks complex sends no domain fact at all: a real fact sent beside it wins (`Element[x, Complexes] && Element[x, Reals]` makes `Im[x] == 0` True).
- Time: `Refine` gets its own `TimeConstraint` (it keeps every part it proved and leaves the rest), `Simplify` and `FullSimplify` get the pair `{per step, total}`, never a `TimeConstrained` wrapper. Each costly stage has a fixed cap below the operation's remaining budget, so one stage cannot starve the cheap ones after it.
- Signs: `Refine` applies sign facts factor by factor (`Abs[p z]` gives `p Abs[z]` for p >= 0); QF does not duplicate those rules. QF's cheap provers (intervals, corners, a rule per head) decide only whether `Refine` needs to run, through one entry that returns positive, negative or unknown, in strict and weak forms. The rule per head adds: `Abs` is nonnegative, positive when its argument is proven nonzero; `Log` has the sign of its argument minus 1; an odd power has its base's sign; `Sin` and `Cos` are nonnegative on the matching half-period. A power is real for a real base and an integer exponent only if the exponent is positive or the base nonzero, and for a nonnegative base with a positive exponent. `Conjugate` passes freely through trig, hyperbolic, `ArcTan`, `ArcCot` and `Exp`, but through `Sqrt`, `Log`, non-integer powers, `ArcSin` and `ArcCos` only for a proven real argument. `Conjugate[s] s` becomes `Abs[s]^2` for a complex s before any sign question.
- The head scan adds non-integer `Power`, `Ramp`, `Min`, `Max`, `Norm`, `UnitStep`, `Unitize` and `RealAbs`, restricted to arguments holding a recorded symbol.
- The stages run as a loop, not one pass: a stage's result is accepted only if strictly smaller by a cost that weighs numbers by their digits (not `LeafCount`), or if it removes `Conjugate`, `Abs`, `Re`, `Im`, `Arg` or `Sign` heads with bounded growth; after an acceptance the loop restarts at the cheapest stage, with an iteration cap. The accepted form and the working form passed to the next stage are kept apart; the working form may grow up to a fixed factor (2 for combining over one denominator, 10 for expanding).
- Size is estimated before a costly step (terms after expansion, denominators after combining, integer multiples inside trig) and the step is skipped above the limit; above an overall size only the per-entry cleanup runs. Size limits are options, not constants.
- Common denominator: a form without polynomial gcds first, then `Together` with `Trig -> False`; trig rewriting belongs to the zero test only.
- Block splitting moves out of the cleanup function (it does nothing for entry-by-entry cleanup) into the whole-matrix checks and phase 7, using one permutation for rows and columns, read from literal zeros after stage 1, never from a zero test.
- Cache: kept across calls, keyed on expression, fact expression, enabled stages and step limit; an answer that came from a timeout is never stored; each stored answer records which test produced it.
- A flag set for the duration of a QF operation stops nested QF calls from opening a second block or splitting again.
- Contradictions: before any `Refine`, an integer symbol's bounds are tightened (a strict bound moves to the next integer) before they are compared, and `Refine[1, facts]` is checked once when the record is built, with its message caught by `Check` and reported. Under contradictory facts every claim comes back proven, so this check is required.

**Phase 2, zero tests**
- Two questions, kept apart: "is this identically zero" (for `U^dagger U - 1`, `H - H^dagger`), where a nonzero value at an allowed point disproves; and "is this nonzero for every allowed value" (for denominators, a norm QF divides by, distinct eigenvalues), which only an exact proof answers.
- Before the exponential form: a product is zero when one factor is; a positive rational power is zero when its base is; one-angle trig at small integer multiples goes through the half-angle tangent substitution, a smaller real rational function.
- `PossibleZeroQ` ignores the record (it answers False for `Sqrt[(1 - p)^2] - (1 - p)`, which is zero for every p <= 1): only its True is used, as likely; its False is never read as nonzero.
- Sample points: drawn inside `BlockRandom`, seeded by the question itself (expression and facts), so the same question always gets the same points whatever the call order; exact interior rationals for disproof, algebraically independent constants for a pass that counts as likely.

**Phase 3**
- One list of assumed conditions and one sample point per operation: each missing sign is chosen so it holds at that point and is appended, and the object carries their conjunction, so two steps cannot assume x > y and y > x.
- Each condition is reduced against the record first (known-sign factors drop out, so `Sqrt[q (1 - q)]` with q >= 0 known needs only q <= 1), conditions implied by earlier ones are dropped, and generic conditions of the form x != 0 are not attached.
- Later cleanup of the object uses its condition as a fact. Undecided answers say whether time ran out or the question could not be decided.

**Phase 4**
- Unitarity from structure first: an exponential `Exp[-I θ H]` with H Hermitian and θ real is unitary; products, tensor products, controlled gates and direct sums of unitaries are unitary; U next to U^dagger cancels. The entry check of `U^dagger U - 1` runs only on gates whose structure says nothing. A named gate never reaches the zero test.
- A scalar factor c multiplies a carried norm by |c|, through the sign-aware `Abs` cleanup.

**Phase 5**
- `PowerExpand` takes an `Assumptions` option: the KAK and stabilizer-phase calls pass the record's facts through it.
- `ComplexExpand` runs only when `Re`, `Im`, `Conjugate` or `Arg` is present, with `TargetFunctions -> {Re, Im}`, indexed complex symbols listed as written (`θ[1]`, or `θ[_]`), and its result accepted only if smaller. Stage 2 runs before it, since it reads no inequality facts. A private copy carries `Positive[q] = True` only for the duration of that one call.

**Phase 7**
- `Eigensystem` is called with `ZeroTest ->` QF's exact test inside the fact block. Triangular and block structure are read from literal zeros. Fewer or more eigenvectors than an eigenvalue's multiplicity, or a zero vector in the output, give undecided. Conditions at special parameter values are not attached by default (decision 6).

## Decisions (Mads, 2026-11-09)

Each was measured or checked before deciding; scripts in `decisions/`.

1. **H3, rotations: no build-time simplification for any rotation**, qubit or qudit, as the design says. Measured (`decisions/rot.wls`): for qubit RZ with a compound angle, Simplify changes nothing (about 2 ms either way); for d = 3 and 4 it makes the build 1.3 to 8 times slower (qutrit RY 133 ms against 16), shrinks the matrix by at most a quarter, leaves apply and dagger costs unchanged, and makes "UnitaryQ" 2 to 3 times slower; at d = 5 it pays only for RX (matrix 10 times smaller, "UnitaryQ" 2.3 s against 6.2 s) and costs over 1 s more for RY. The rare large case is left to the cleanup function when a property needs it. Goes in phase 6; the RZPhases and RotationWires tests change to the unsimplified form.
2. **C2, which functions read the record: the entanglement monotones only.** Checked (`decisions/facts.wls`): concurrence and negativity of "Werner"[p] and of a bit-flipped Bell state come back with Conjugate[p], Re and Ramp, because they do not know p is real; reading the record cleans them. `QuantumEntangledQ` does not need it: it already counts only the values that make a density matrix, so [0, 1] adds nothing. KAK and the stabilizer phase check take every symbol real already and use nothing more. Goes in phase 5.
3. **K, work the design did not ask for.** Keep the role ranges (Depolarizing up to 4/3, ResetError, affine slots, one-component Bloch vectors, PhaseShift integers) and the new yes/no answers (symbolic "UnitalQ", Pure/Mixed "Type" after the phase 1 fix). Change: basis "Options" lists "SymbolDomains" only when it is not empty, so objects without facts look as on main and main's three ApplyFixedCost expectations return to their original form; goes in phase 1. Keep substitution dropping a symbol's facts until step 7 replaces it with a check.
4. **Order: as listed.** Phases 1 and 2 first, since every later phase calls the cleanup function; phase 4 (carried norm) before step 7.

5. **The undecided answer: `Indeterminate`, with a message**, as `QuantumEntangledQ` does (`Indeterminate && False` is False, `Indeterminate || True` is True; `Undefined` absorbs both). Three cases, each with its own message: time ran out; the answer differs across the allowed values, with the message naming where it holds (for example p <= 2 for diag(p, 2 - p), p in [0, 3]) when a time-limited `Reduce` finds it; or the question is outside what the tests decide. An `Indeterminate` that came from a timeout is never stored in the property cache. Every internal `If` or `Which` on these checks goes through `TrueQ` (the "Physical" property's `If[qs["PhysicalQ"], ...]` stays unevaluated today), and the summary boxes show the same value instead of `$Aborted`. Goes in phase 1.
6. **Eigensystem conditions at special values: off by default**, available through an option; QF answers for generic parameter values, as `Eigensystem` does. Goes in phase 7.
7. **An explicit `Reals` or `Complexes` from the user drops a role's bounds**, so a channel can be studied outside its physical region on purpose; a user interval still meets the role's interval. Goes in phase 1 with the record changes.

The numeric side (tolerance-aware checks on machine-number objects, and one nearest-state repair: "Physical" clips and rescales while the entanglement code uses the simplex projection) is a separate task, `nearest-state-fix`, landing on main.
