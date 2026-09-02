# `PauliStabilizer["M", q]` disagrees with `["Expectation", ...]` on a deterministic measurement

Found on 2026-08-31 while validating the QEC syndrome-extraction circuit
(`OngoingProjects/QEC/QECCore/Circuit.wl`) against the engine.

## Summary

`ps["M", q]` returns an outcome that depends on **which generating set the state
was constructed from**, not on the state. Two `PauliStabilizer` objects that are
the same physical state (identical state vector, identical expectation values on
every Pauli) give different measurement outcomes.

`["Expectation", ...]` is correct in every case checked; `["M", q]` is the one
that is wrong.

## Minimal reproduction

```wl
Needs["Wolfram`QuantumFramework`"];

a = PauliStabilizer[{"ZZ", "-IZ"}];        (* |11> *)
b = PauliStabilizer[2]["X", 1]["X", 2];    (* |11> *)

Normal[a["State"]["StateVector"]]          (* {0, 0, 0, 1} *)
Normal[b["State"]["StateVector"]]          (* {0, 0, 0, 1} *)

Table[op -> {a["Expectation", op], b["Expectation", op]},
      {op, {"ZI", "IZ", "ZZ", "XX"}}]
(* {"ZI" -> {-1, -1}, "IZ" -> {-1, -1}, "ZZ" -> {1, 1}, "XX" -> {0, 0}} *)

Table[First[Keys[a["M", k]]], {k, 2}]      (* {0, 1}  <-- wrong on qubit 1 *)
Table[First[Keys[b["M", k]]], {k, 2}]      (* {1, 1}  <-- right *)
```

Both measurements are deterministic (the returned Association has a single key),
so there is no sampling involved: one of the two answers is simply false.

## Pattern

On two qubits, `PauliStabilizer[{s1, s2}]["M", q]` returns the sign bit of the
`q`-th **input generator** rather than the eigenvalue of `Z_q`:

| input | `["M", {1,2}]` | `(1 - Expectation[Z_k])/2` |
|---|---|---|
| `{"ZZ", "IZ"}`   | `{0, 0}` | `{0, 0}` (agree) |
| `{"ZZ", "-IZ"}`  | `{0, 1}` | `{1, 1}` |
| `{"-ZZ", "IZ"}`  | `{1, 1}` | `{1, 0}` |
| `{"-ZZ", "-IZ"}` | `{1, 0}` | `{0, 1}` |

The two agree exactly when the input generators are already single-qubit `Z`
operators, which is the case in essentially every small hand-written test. That
matches the note in the `qf` skill that this class of Aaronson-Gottesman
bookkeeping bug is "invisible to single-measurement tests".

Three qubits, same shape:

```wl
PauliStabilizer[{"ZZI", "-IZZ", "-IIZ"}]
(* ["M", {1,2,3}]                 -> {0, 1, 0} *)
(* (1 - Expectation[Z_k])/2       -> {0, 0, 1} *)
```

## Where it bites

Any workflow that builds a state from its stabilizer generators (rather than by
applying gates to `PauliStabilizer[n]`) and then measures. Syndrome extraction is
exactly that workflow: the encoded state is naturally specified by the code's
generators.

## Workaround in use

`OngoingProjects/QEC/` reads deterministic check outcomes with

```wl
(1 - ps["Expectation", pauliString]) / 2
```

which is the route `codePhysicalCorrectionCycle` already used and which the QEC
test suite exercises. No use of `["M", q]` on a generator-constructed state.

## Suggested next step

Compare `Kernel/Stabilizer/Measurement.m`'s deterministic branch against the
canonicalisation done by the `PauliStabilizer[{__String}]` constructor: the
constructor appears to store the given rows verbatim (`["Stabilizers"]` returns
them unchanged) while the measurement's deterministic path assumes a form in
which row `q` is the one carrying `Z_q`.
