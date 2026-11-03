# BosonicNormalOrder on a list of expressions

Verdict: unsupported input that should raise a message instead. A list of two or more expressions is not a documented call, but the argument pattern accepts it and the default method returns a wrong value with no message. The wrong value is not a mis-ordered list: it is the internal Groebner reduction record of the second list element only.

## Minimal reproduction

```wl
{a, ad} = FieldVariables[{1}];
BosonicNormalOrder[{a ** ad, ad ** a}, {a, ad}]
(* {{0 &}, SuperDagger[\[FormalA]1] ** \[FormalA]1}, no message *)
Map[BosonicNormalOrder[#, {a, ad}] &, {a ** ad, ad ** a}]
(* {1 + SuperDagger[\[FormalA]1] ** \[FormalA]1, SuperDagger[\[FormalA]1] ** \[FormalA]1}, correct *)
```

## Source line

`QuantumFramework/Kernel/SecondQuantization/BosonicAlgebra.m` at commit fbebae87:

- line 74: `BosonicNormalOrder[expr_, vars_List, opts : OptionsPattern[]]` takes `expr_`, so a List gets through.
- line 82: the default method (`Automatic | "GrobnerBasis"`) passes `expr` straight to `grobnerNormalOrder`.
- line 48 (line 47 when `"Scalars"` is non-empty): `NonCommutativePolynomialReduce[expr, rels, NonCommutativeAlgebra[...]][[2]]`.

For one polynomial, `NonCommutativePolynomialReduce` returns `{{f1, ...}, r}` and `[[2]]` is the remainder `r`, the normal-ordered form. For a list of polynomials it returns one such pair per element (the probe calls it directly and gets `{{{0 &}, ...}, {{-#1 &}, ...}}`), so `[[2]]` picks the whole pair for the second element. The `{0 &}` in the reported output is that element's quotient list (`ad ** a` is already normal ordered, so its quotient is zero), and the second entry is its remainder. Nothing in the chain is a sum-matching pattern or a pure function mapped over the list.

`BosonicAntinormalOrder` (lines 118, 126) goes through the same `grobnerNormalOrder` call and fails the same way. `bosonicReduce` (line 149) has the same `[[2]]` but is internal.

## Is list input supported?

No sign that it is:

- `Attributes[BosonicNormalOrder]` is `{}` (not `Listable`).
- The kernel usage message names a single `expr` in all three forms.
- No page under `QuantumFramework/Documentation` mentions `BosonicNormalOrder`.
- No test under `Tests` passes a list as the first argument (searched `Tests/SecondQuantization.wlt` and `Tests/SecondQuantization/*.wlt`).
- Only two definitions exist: `expr_, {}` and `expr_, vars_List`.

## What comes back (local checkout and installed 2.1.3 are identical)

| Call | Result | Message |
|---|---|---|
| `a ** ad` (single) | `1 + ad ** a` | none |
| `{a ** ad}` | unevaluated `Part[...]` | `Part::partw` |
| `{a ** ad, ad ** a}` | `{{0 &}, ad ** a}` (second element's record) | none |
| `{ad ** a, a ** ad}` | `{{#1 &}, 1 + ad ** a}` | none |
| `{a ** ad, a ** a ** ad, ad ** a ** a ** ad}` | `{{#1 ** a + a ** #1 &}, 2 a + ad ** a^2}` (second element only; third is dropped) | none |
| `{a ** ad + 2 a, 3 + a ** a ** ad ** ad}` | `{{...}, 5 + ad^2 ** a^2 + 4 ad ** a}` (second element only) | none |
| `{a ** ad ** b, b ** bd ** a}`, vars `{a, ad, b, bd}` | `{{6 quotient functions}, a + bd ** a ** b}` (second element only) | none |
| `{a ** ad, ad ** a}`, `Method -> "Blasiak"` | fails (value not recorded, `Check` caught the message) | `ParseBlasiakMonomial::badfactor` |
| `{a ** ad, ad ** a}` with no vars | unevaluated | none |
| `a ** ad` with no vars | unevaluated | none |
| `BosonicAntinormalOrder[{ad ** a, a ** ad}, {a, ad}]` | `{{0 &}, a ** ad}` (second element's record) | none |

Here `a, ad, b, bd` stand for `\[FormalA]1, SuperDagger[\[FormalA]1], \[FormalA]2, SuperDagger[\[FormalA]2]`. Every remainder shown is the correct normal (or antinormal) form of the second element; the first and any later elements are lost. Lists of length 1 and the Blasiak method do emit a message; the default method on lists of length 2 or more never does.

The missing-vars case is not part of this defect: there is no one-argument definition, so the call stays unevaluated whether or not the argument is a list.

## Installed QF 2.1.3

The run without the load line resolves to the installed paclet repository copy (version 2.1.3) and returns exactly the same results for every probe. The local checkout also reports version 2.1.3 in its PacletInfo, so the version string alone does not distinguish them; the paclet location does.

## Files

- `load.wl`: the two `Needs` lines.
- `probe.wl`: every probe above; writes the results as InputForm to stdout.
- `run-installed.wls`: loads `load.wl` and `probe.wl` from this folder, so `Needs` picks the installed paclet. To probe the local checkout, use a private runner whose first line is `PacletDirectoryLoad` on the checkout's `QuantumFramework/` folder, then `Get` the same two files.

Run with `wolframscript -file run-installed.wls` from this folder.
