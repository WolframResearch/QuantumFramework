# Structured arrays in QF: task briefs

QF stores every operator as a `SparseArray`, so when an operator's structure already fixes an answer (the exponential of a diagonal is its entries' exponentials, the inverse of a unitary is its conjugate transpose, the spectrum of a permutation comes from its cycles), QF still computes that answer from scratch. The survey `../structured-arrays-in-qf.md` and the table `../plan-structured-arrays-benchmark.md` measured ten such operations in September 2026. Each brief here turns one remaining operation into a task someone can pick up without any other context.

## Where the survey stands

| Survey row | Operation | State |
|---|---|---|
| 1 | Exponential of a diagonal operator | Landed 2026-10-05 (`191f1ac3`), see `../exp-diagonal-shortcut-report.md` |
| 9 | Exponential of a dephasing Liouvillian | Landed with row 1 and `f6299b43` (the superoperator stays a `SparseArray`) |
| 2 | f of a diagonal operator (`Cos[qo]`, `Sqrt[qo]`) | Already entry by entry in `f[qo]`; the September note that `Cos[qo]` is wrong no longer holds |
| 6 | Change of basis into the Fourier basis | Brief 1, `01-unitary-change-of-basis.md` |
| 5 | Spectrum of a diagonal operator | Brief 2, `02-diagonal-spectrum.md` |
| (open question 6 of the report) | Building a Liouvillian | Brief 3, `03-liouvillian-build.md`, investigate first |
| 3, 4 | Spectrum of a permutation and of the QFT | Brief 4, `04-permutation-and-qft-spectra.md`, investigate first |
| 7 | QFT applied to a state | Brief 5, `05-qft-as-fft.md`, investigate first |
| 8 | Controlled evolution | Not briefed: only a constant factor, from about eleven qubits |
| 10 | H on every qubit | Not briefed: WL has no structured form that applies it fast |

The survey and the table describe QF as it was before 2026-10-05; their rows 1, 2 and 9 are out of date. Every number in the briefs was measured on main at `f6299b43` (2026-10-05) by `baseline.wls`, `supplement.wls` and `supplement2.wls` in this folder, whose raw output is in `baseline-out/`.

Briefs 1 and 2 are ready to implement. Briefs 3 to 5 start with an investigation whose result decides the change.

## How to do any of these tasks

**Load QF from the checkout**, not from an installed paclet:

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
```

**Compare before and after on two copies.** Export the current main twice and edit one copy:

```bash
mkdir -p before after && git archive HEAD QuantumFramework Tests | tar -x -C before && git archive HEAD QuantumFramework Tests | tar -x -C after
```

Run every script on both copies, one kernel at a time, and keep the raw output. A script loads a copy with `PacletDirectoryLoad[FileNameJoin[{dir, "QuantumFramework"}]]`.

**Run the tests** with the current runner only:

```bash
wolframscript -file Tests/RunTests.wls
```

```bash
wolframscript -file Tests/RunTests.wls DiagonalExp
```

The second form runs only the test files whose name contains the pattern. Read the `TOTAL` and `STATUS` lines; `STATUS GREEN` means every test passed and nothing printed a message. Never run a `Tests/RunTests.wls` taken from an older commit: the runners from `16df487b` to `ab5ba368` uninstall the installed QF paclet. Run one full suite at a time.

**Write the tests first, and prove they test the change.** A new test file goes in `Tests/`, named after the area (`QuantumOperatorDiagonalExp.wlt` is the model to copy). Every test file shares one `Global`` context, so give helper symbols a prefix no other file uses. Run the new file on the unchanged copy: the tests that cover the change must fail there and pass on the changed copy, and the regression tests must pass on both.

**House rules for the code** (the QF maintainers enforce these):

- No `Quiet`. A test that expects a message names it in `VerificationTest`'s third argument.
- No `Print`. `wolframscript -file` prints nothing on its own, so scripts write with `WriteString[$Output, ...]`.
- `With` and `Block` over `Module`; `Enclose`/`Confirm` for failure in new code; pattern-matched definitions over `If` cascades.
- A `SparseArray` stays sparse, and a machine result stays a packed array.
- Comments describe what the code does, never what it replaced.

**Comparing numbers.** `===`, `==` and `Hash` treat machine numbers that differ in the last bits as equal, and a test that demands agreement to the last bit breaks when WL changes an algorithm. Compare machine results against a stated bound (`Max[Abs[a - b]] <= 2 $MachineEpsilon`). To show that two results are identical to the last bit, compare `Hash[BinarySerialize[Developer`FromPackedArray[Normal[x]]]]`.

**Timing.** Call `ClearSystemCache[]` before each run, take the fastest of three when one run is short, and put every slow case inside `TimeConstrained`. Other kernels on the machine move timings by tens of percent: to judge a difference that small, time both versions alternately in separate kernels, or both functions in one kernel.

**Landing.** Commit straight to `main`, one commit per independent part, each with its own tests, so each can be reverted alone. The message says what is now true: `Area: one sentence`, then a short body. No `Co-Authored-By` line. Run the full suite in the checkout before pushing, and push only when Mads says so.

## Traps met on this work

- **Arrays paclet and older QF.** With the `Wolfram/Arrays` paclet installed on 2026-10-05 (1.4.0), QF older than `bd321ccb` fails in `ArrayContract` whenever an operator acts on a mixed state. Base every comparison on current main.
- **The `"Diagonal"` constructor's label.** `QuantumOperator["Diagonal"[list], order]` labels the operator with the whole list, and arithmetic on that label can cost more than the operation being timed. Time with `"Label" -> None`.
- **Shadowed symbols.** `wolframscript -print all` and `wolframscript -code "Needs[...]; ..."` parse the whole input before `Needs` runs, so QF symbols land in `Global`` and stay inert. Use `-file` with one top-level expression after another.
- **A killed kernel** makes `wolframscript` print "license error". Check memory use before suspecting the license.
- **Internal symbols** such as `matrixFunction` and `eigensystem` live in ``Wolfram`QuantumFramework`PackageScope` ``; tests reach them because the runner puts that context on the path, scripts must spell it out.
