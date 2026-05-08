# wl/ style guide

Conventions for everything under `wl/`: paclet kernel sources (`*.wl`),
test specs (`*.wlt`), runner scripts (`*.wls`), and the LibraryLink C
bridge (`CSource/*.c`).

## Rules the user has explicitly called out

These are non-negotiable. Strip on sight.

### No `Print`

Never call `Print` in a `.wl`, `.wlt`, or `.wls` file. If textual output
is needed, define a local helper `print` that wraps either
`WriteString[$Output, ...]` (notebook / kernel context) or
`WriteString["stdout", ...]` (`wolframscript -c` / `-f` context) and
call that.

```wolfram
print[args___] := WriteString["stdout", StringJoin @@ Map[ToString, {args}], "\n"]
```

`Print` also trips an IDE lint warning ("Suspicious use of session
symbol Print") and is noisy in batch test runs.

### No em dashes (`—`, U+2014)

Don't write em dashes in source files, docs, comments, or commit
messages. Use a plain hyphen (`-`), a colon, or a sentence break.

### No Unicode box-drawing characters

Don't use `─` (U+2500) or related box-drawing chars (`┌ ┐ │ ┴ ─`) as
visual section banners in comments. Use plain ASCII:

```c
// === lifecycle ===     // good
// --- lifecycle ---     // good
// ─── lifecycle ───     // BAD
```

Same applies in `.wl` (`(* === lifecycle === *)`) and `.md`
(`## Lifecycle`).

### No decorative Unicode arrows in source

Use `->` (or `>` in shell prompt strings) instead of `→`, `←`, etc.
Same spirit as the rules above: ASCII-only in source files unless the
character carries meaning (mathematical typesetting in a comment is
fine when relevant).

### Dark mode + Standard colors

Always design WL output for both light and dark Wolfram themes from
the start.

- Use **Standard color names** (`StandardBlue`, `StandardRed`,
  `StandardGreen`, `StandardOrange`, `StandardYellow`,
  `StandardPurple`, `StandardGray`) instead of `RGBColor[...]` or
  `Darker[Blue, 0.4]`. They read correctly on light and dark
  backgrounds.
- Use `LightDarkSwitched[lightValue, darkValue]` when a value
  (color, opacity, thickness) genuinely needs to differ between
  modes. Avoid hard-coded `White`/`Black` for `Background`,
  `FaceForm`, etc.; either let the theme drive it, or wrap in
  `LightDarkSwitched`.
- Don't unconditionally pin `Background -> White` on a Graph or
  Graphics. If a static PNG export needs an explicit background,
  set it on the Export call, not on the Graph.

## Indentation and spacing

- 4-space indentation. Never 2.
- Spaces around infix operators and pattern tests: `t_ ? NumericQ`,
  not `t_?NumericQ`.
- Prefer structural indentation over column alignment.
- One semantic unit per line in long forms (`Which`, `Switch`,
  `Association`, `Table`, option lists).
- For multi-line `If`, put a space after the opening bracket so the
  test argument lines up with the branches:

  ```wolfram
  If[ Length[dirs] === 0,
      debugPrint["no examples to run"];
      Exit[1]
  ]
  ```

  Single-line `If[cond, then, else]` does not need the leading space.

### Never split a binary operator's operands across lines

Operators (`+`, `-`, `*`, `/`, `.`, `&&`, `||`, etc.) must have both
operands on the same line.  The Wolfram IDE flags split operands as
`DifferentLine` and the resulting form reads worse than a longer
single line.

```wolfram
(* BAD: trips DifferentLine lint *)
TSet[wTen, wTen - lrHat * mTen
                  / (Sqrt[vTen] * invSqrtB2cor + eps)];

(* GOOD: keep the whole arithmetic chain on one line, even if it's
   long; or factor an intermediate into a named binding first. *)
TSet[wTen, wTen - lrHat * mTen / (Sqrt[vTen] * invSqrtB2cor + eps)];

(* GOOD: factor when the line gets unreadably wide *)
denom = Sqrt[vTen] * invSqrtB2cor + eps;
TSet[wTen, wTen - lrHat * mTen / denom];
```

The lint rule is `wolfram lint(DifferentLine)`; if you see it,
either join the line or factor.

### Optional arguments

For functions that take Wolfram-style options, use
`opts : OptionsPattern[]`, not `opts___ ? OptionQ`.  The former
integrates with `Options[fn]` declarations and `OptionValue[fn, key]`,
which is the canonical Wolfram pattern.

```wolfram
Options[myFn] = {GraphLayout -> "LayeredDigraphEmbedding"};

myFn[args___, opts : OptionsPattern[]] :=
    With[{layout = OptionValue[GraphLayout]},
        ...
    ]
```

## Definitions

Prefer `Block` for local workspaces unless `Module`'s unique-symbol
guarantee is actually required. Don't add a trailing `;` to
`SetDelayed` definitions.

```wolfram
f[x_] := x + 1

g[args__] := Block[{
    x = ...,
    y
},
    body
]
```

The `Block`/`Module` variable list opens on the same line as the
head (`Block[{`), variable bindings are indented one level deeper,
and `}`, body, and the closing `]` all return to the column where
the line that opens `Block[` starts. This keeps the variable list
visually distinct from the body and matches IDE auto-formatting.

## Control flow

### No `For` loops

Don't write `For[i = 1, i <= n, i++, ...]`. It reads like C and
forces explicit counters and mutation. Use:

- `Do[body, {i, n}]` for side effects with a counter.
- `Table[expr, {i, n}]` to build a list.
- `Nest`, `NestList`, `NestWhile`, `Fold`, `FoldList` when there's
  iterative state to thread through a fixed transition.

```wolfram
(* BAD *)
For[i = 1, i <= n, i++,
    AppendTo[acc, f[i]]
]

(* GOOD *)
acc = Table[f[i], {i, n}]
```

### No `Head[expr] === Foo` -- use `MatchQ`

`Head[x] === Foo` is a structural test that doesn't compose with
patterns. `MatchQ[x, _Foo]` (or `MatchQ[x, Foo[args...]]` for a
shape check) is the canonical form, and it composes with patterns
when you need them.

```wolfram
(* BAD *)
If[ Head[expr] === Inactive[Equal] && expr[[1]] === expr[[2]],
    ...
]

(* GOOD *)
If[ MatchQ[expr, Inactive[Equal][x_, x_]],
    ...
]
```

## Mutation

### No `AppendTo` (or other list-growing mutation)

Don't grow a list with `AppendTo`/`PrependTo` inside a loop.
`AppendTo[xs, y]` is `xs = Append[xs, y]`, which copies on every
step. Build the result with `Table`, `Map`, or `Fold` so the
final list is allocated once.

```wolfram
(* BAD *)
out = {};
Do[ AppendTo[out, f[i]], {i, n} ];
out

(* GOOD *)
out = Table[f[i], {i, n}]

(* BAD: state threaded by mutating a local *)
acc = init;
Do[ acc = step[acc, x], {x, xs} ];
acc

(* GOOD: state threaded by Fold *)
acc = Fold[step, init, xs]
```

The same goes for incrementally building an Association: use
`Association[Table[k -> v, ...]]`, `AssociationMap`, or
`Fold[Append, <||>, kvPairs]` rather than `assoc[k] = v` inside
a `Do`.

### Bracket alignment

A closing `]` (or `}`, `|>`, `)`) goes on its own line, indented
to the same column as the opening head. Never end a multi-line
form with `...]` on the last expression's line.

```wolfram
(* GOOD *)
Module[{x, y, z},
    body1;
    body2
]

(* BAD *)
Module[{x, y, z},
    body1;
    body2]
```

This makes block boundaries scan-readable and matches what the IDE
auto-formatter expects.

## Naming

- Public symbols: `CamelCase` (e.g. `TLam`, `THeap`, `TWnf`).
- Internal helpers: `lowerCamelCase` (e.g. `debugPrint`, `loadFn`).
- Don't prefix internal helpers with `i...`.

## Tests

Use `VerificationTest` from the standard testing framework. Test specs
live in `*.wlt` files; the runner is a `.wls` script invoked by
`make wl-test`.

```wolfram
VerificationTest[
    expression,
    expectedOutput,
    TestID -> "human-readable id"
]
```

Don't put `Print` inside test bodies. The runner reports outcomes via
`TestReport`.

## File structure

When practical, organize a `.wl` file in this order:

1. Short file comment if needed.
2. `BeginPackage` declarations and public `::usage` strings.
3. `Begin["`Private`"]`.
4. Small general helpers (e.g. `debugPrint`).
5. Domain-specific helpers.
6. Main entry-point definitions near the end.
7. `End[]; EndPackage[];`

## Comments

- Comment non-obvious behavior, quirks, or external format
  constraints (LibraryLink calling conventions, paclet layout
  expectations, etc.).
- Don't narrate obvious code.
- Prefer one short section comment over many tiny inline comments.

## CSource (LibraryLink bridge)

The `CSource/*.c` bridge files follow the C `src/` style guide
([../AGENTS.md](../AGENTS.md)) but with one extra rule: **scalar in,
scalar out**. Every exported `EXTERN_C DLLEXPORT int <name>(...)`
function takes only `Integer` arguments and returns a single `Integer`
via `MArgument_setInteger(res, ...)`.

Higher-level constructors that bind names (`TLam[x, body]` with
HoldAll, `TDup[body, {dp0, dp1} |-> ...]` returning a continuation
result) are synthesized on the WL side from these scalar primitives.
This keeps the C surface tiny and the bridge testable from C alone.
