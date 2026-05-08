# Per-Function Doc Page Runbook (any Wolfram paclet)

Hand this file back to Claude with:

> *"Follow `<path-to-this-file>` for function `F` in the paclet at `<PACLET>`."*

Replace:
- `F` with any public symbol of the paclet (e.g. `QuantumChannel`, `MyFunction`).
- `<PACLET>` with the absolute path of the paclet root — the directory containing `PacletInfo.wl`.

Optional extras the user can supply at invocation time:
- **Primary context override** (default: read from `PacletInfo.wl`).
- **Sub-contexts** to also scan (default: read from `PacletInfo.wl` Kernel extension).
- **Reports directory** (default: `<PACLET>/scratch/docs-audit/` — created if missing).
- **Project-specific deep skill** name (e.g. `qf` for `QuantumFramework`). Use if present.
- **Exclusion list** — symbols that are public but should not get doc pages.
- **Project-specific kernel→file lookup table** — for paclets where symbol names don't match their declaring filenames; otherwise grep.

The runbook handles both cases:
- **Branch A** — `F` already has `<PACLET>/Documentation/English/ReferencePages/Symbols/F.nb`. Audit it against the kernel and fill gaps.
- **Branch B** — no such file exists. Triage `F`, then (if confirmed public) write `::usage` and create the page.

> **Read carefully, do not guess.** Inventory first, file:line anchors on every claim, verbatim excerpts. No summaries of summaries. If a kernel symbol's behavior depends on a file you have not yet read, stop and read it.

## 0. Session preamble (always perform first)

### 0.1 Collect or derive the paclet configuration

Resolve these variables before doing anything else. Do **not** assume defaults; verify by reading the paclet manifest.

| Variable | How to resolve |
|---|---|
| `$PACLET` | User input — absolute path to the paclet root. Confirm `$PACLET/PacletInfo.wl` exists. |
| `$KER` | Read `"Kernel"` extension's `"Root"` field in `PacletInfo.wl`; default `$PACLET/Kernel`. |
| `$REF` | `$PACLET/Documentation/<Lang>/ReferencePages/Symbols/`; `<Lang>` is the `"Language"` field of the `Documentation` extension (typically `English`). Confirm the directory exists. |
| `$GUIDES` | `$PACLET/Documentation/<Lang>/Guides/` — list `.nb` files; if exactly one, that's the canonical guide. |
| `$TUTORIALS` | `$PACLET/Documentation/<Lang>/Tutorials/`. May be empty. |
| `$DOCS` (reports dir) | User-provided or default `$PACLET/scratch/docs-audit/`. **Create if missing.** Holds the verification harness, per-page reports, and the issues report. |
| `<PrimaryContext>` | Read `"PrimaryContext"` field in `PacletInfo.wl` (e.g. `Wolfram\`QuantumFramework\``). |
| `<Subcontexts>` | Read the `"Context"` list in the Kernel extension; filter to those that actually export symbols (skip `Init`-style internal contexts). |
| `<PacletName>` | Read `"Name"` field (e.g. `Wolfram/QuantumFramework`). Used in `paclet:<PacletName>/ref/F` links. |

If `PacletInfo.wl` is unreadable or absent, stop and ask the user to confirm the path.

Run a quick sanity check:
```
wolframscript -code 'PacletDirectoryLoad["<PACLET>"]; Needs["<PrimaryContext>"]; Length[Names["<PrimaryContext>*"]]'
```
should print a positive integer. If zero, the paclet didn't load — investigate before continuing.

### 0.2 Read the user's persistent memory once at session start

Path: `~/.claude/projects/<project-slug>/memory/MEMORY.md` and any referenced files. Honor whatever feedback rules are encoded there. Common rules to expect (project-independent):
- Use `wolframscript` for kernel evaluation, not the WolframLanguageEvaluator MCP. Long runs go via `run_in_background`.
- Never use `Quiet` to suppress messages. Capture them or surface them.
- For any constructor/accessor pair the FIRST round-trip test must be exact equality `A[C[x]] === x`.
- Don't add `Co-Authored-By` trailers to commits.

### 0.3 Skill stack — invoke as the work demands

These are not optional dressing; they encode conventions you would otherwise reinvent badly.

| Skill | When to use |
|---|---|
| `<DeepSkill>` (if user provides one) | **Always active for the paclet's domain.** E.g. `qf` for `QuantumFramework` encodes paclet-specific pitfalls. Trumps broader skills. Skip this row if no project-specific skill exists. |
| `documentation-writing` | **Primary skill for writing or substantively editing any reference page.** Knows WFR/paclet Reference Page conventions, Guide page structure, and Tech Note tutorial format. Invoke at the start of Branch A §A.5 (filling gaps) and Branch B §B.6 (creating a new page). |
| `nb-writer-v2` | Preferred for low-level cell construction inside the page that `documentation-writing` produces. Kernel-driven serialization handles BoxData, RowBox, CellGroupData, FormBox, Greek/special characters, and inline-math escapes correctly. |
| `nb-writer` | Fallback if `nb-writer-v2` is unavailable or overkill for a one-cell edit. |
| `nb-reader` | For reading existing `.nb` files (existing reference page on Branch A, sibling templates on Branch B). Faster and more accurate than raw `Read` on serialized notebooks. |
| `nb-inline-math` | When adding inline math (e.g. `\(\rho\)`) inside a Text cell or ExampleText cell. Prevents bracket-nesting and Greek-character pitfalls in BoxData. |
| `wl-native-style` | When writing example code. Prefer functional forms (`Map`, `Fold`, `Table`, `Association`) over `Do`/`For`/`While`/`AppendTo`. |
| `wl-broadcasting` | When example code multiplies arrays of mismatched rank. WL is row-wise, not NumPy-style. |
| `learning-by-computing` | Optional — only if a Reference Page's Notes section grows long enough to warrant explanation-then-computation flow. Most pages do not need this. |
| Domain-area skills (e.g. `wl-quantum-computing`, `wl-symbolic-computation`, `wl-machine-learning`, `wl-visualization`, `wl-data-and-knowledge`) | Background reference when a specific example draws from those areas. Do not let them override `<DeepSkill>` for paclet-specific behavior. |

If you find yourself hand-rolling `.nb` cell expressions instead of using `documentation-writing` + `nb-writer-v2`, stop and switch to the skills.

### 0.4 Verification harness

Path: `$DOCS/verify-doc-page.wls`. **If missing, build it before any work that needs verification — see §6.** If present, smoke-test it on one existing page and confirm it produces a markdown report.

### 0.5 Issues report file

Path: `$DOCS/issues-report.md`. If missing, create it with a header skeleton (see §7).

## 1. Determine the branch

### 1.1 Check that `F` is a public paclet export

```
wolframscript -code 'PacletDirectoryLoad["<PACLET>"];
  Needs["<PrimaryContext>"];
  ctxs = {"<PrimaryContext>", <each subcontext>};
  Select[Flatten[Names[#<>"F"]&/@ctxs], !StringEndsQ[#, "Private`F"]&]'
```
If no context contains `F` as a public symbol, stop and tell the user `F` is not a public export.

### 1.2 Check the user-provided exclusion list (if any)

If `F` is excluded, stop and tell the user this symbol is intentionally undocumented.

### 1.3 Check existence of `$REF/F.nb`
- Exists → **Branch A** (§3, "Update existing").
- Does not exist → **Branch B** (§4, "Create new").

## 2. Glossary

- **Form**: a distinct constructor/call signature recognized by `F`. E.g. `F[name]`, `F[<|...|>]`, `F[a, b]` are three forms.
- **Option**: a key in `Options[F]`. Listed via `Options[F]` after the paclet is loaded; declared in kernel as `Options[F] = {...}`.
- **Property**: a string `F[obj]["propname"]` understands. Lives in dispatch tables, often in `Properties.m` or near the symbol's main downvalues.
- **Named instance**: a named-string form like `F["SomeName"]`. Often lives in `Named*.m` companion files.
- **Sibling symbol**: a symbol declared in the same kernel file. Often shares behavior.

## 3. Branch A — Update existing reference page

### A.1 Read every kernel file that touches F

If the user supplied a project-specific lookup table mapping symbols to declaring files, use it first. Otherwise grep:
```
grep -ln 'F::usage\|^F\[\|F\[.*\] :=\|F\[.*\] =' $KER/*.m $KER/*/*.m
```

Read **every** match. Declarations may be split across the family directory plus shared files (commonly: `Init.m`, `Usage.m`, `Utilities.m`, `Formatting.m`, plus any cross-cutting modules the paclet uses).

For each file you read, write down:
- The line ranges defining `F` downvalues.
- Every `OptionValue["..."]` callsite (reveals undeclared options).
- Every property dispatch line (e.g. `obj["X"] := ...`, `Switch[prop, "X", ...]`).
- Every named-instance handler.

### A.2 Inventory the four coverage axes

For `F` produce four lists with file:line anchors:

A.2.a **Forms** — every distinct downvalue signature. Cross-check against `F::usage` text; if usage lags behind code, the code wins.

A.2.b **Options** — run `wolframscript -code 'Needs["<PrimaryContext>"]; Options[F]'` for the canonical declared set. Then grep for `OptionValue` callsites in F's files to catch any options consumed but not declared (rare bug pattern).

A.2.c **Properties** — list every string `F` accepts in property position. Grep patterns to use:
- `"\"[A-Z][A-Za-z]*\""` near `prop_String`, `MemberQ[$Properties`, `Switch[prop`, `Replace[prop`, or `obj[("name1"|"name2"|...)] :=`.
- Don't miss properties added in `Formatting.m` or visualization modules.

A.2.d **Named instances** — for symbols whose `Named*.m` companion file exists, list every named class. Group by what they share.

### A.3 Read the existing F.nb

**Use the `nb-reader` skill** to read `$REF/F.nb`. (Raw `Read` on a serialized .nb is brittle.) Inventory:
- Existing input forms shown in the **Usage** cell (count `ModInfo` rows).
- Existing examples in **Primary Examples** (count Input cells).
- Existing **Scope** subsections (per uncovered form).
- Existing **Options** subsections (one ExampleSubsection per option).
- Existing **Properties & Relations** examples.
- Existing **Possible Issues** notes.
- Whether **Generalizations & Extensions**, **Applications**, **Interactive Examples**, **Neat Examples** sections exist (often empty placeholders).

### A.4 Diff: produce a gap list

For each axis (Forms, Options, Properties, Named instances), list items in §A.2 but not §A.3. Anchor each gap with the kernel file:line where it's defined.

If the gap list is empty, jump to §A.6 (verification only).

### A.5 Edit F.nb to fill gaps

**Invoke the `documentation-writing` skill before any edits.** It knows the paclet Reference Page conventions and will guide section ordering, link styles, and prose voice. Use the **NB schema (§8)** as the structural cross-check. For the actual cell-byte construction inside the skill's plan, use `nb-writer-v2`. For inline math inside ExampleText/Notes cells, use `nb-inline-math`. Specific rules:

A.5.a **For each missing form**: add an `ExampleText` describing the form, then an `Input` cell exercising it, then expected `Output`. Place under **Scope** (`ExampleSection` with `CellTags -> {"Scope", "ExtendedExamples"}`).

A.5.b **For each missing option**: add an `ExampleSubsection` titled with the option name under the existing **Options** section. Body: short `ExampleText` explaining the option, an `Input` showing default vs explicit value, and the resulting `Output`.

A.5.c **For each missing property**: add either (a) a brief `ExampleText` + `Input` + `Output` triple under **Properties & Relations**, or (b) an entry in the **Notes** 2ColumnTableMod table at the top of the page if it's better summarized than demonstrated. Default: at least one demonstration in Examples for every property.

A.5.d **For each missing named instance class**: add at least one `Input` example per group under **Scope**.

A.5.e **Style discipline**:
- Cell style names exactly: `"PrimaryExamplesSection"`, `"ExampleSection"`, `"ExampleSubsection"`, `"ExampleText"`, `"Input"`, `"Output"`, `"ExampleDelimiter"`. Do not invent new style names.
- Every cell needs `CellID` (use a fresh integer; do not collide with existing IDs in the file) and `ExpressionUUID` (use `CreateUUID[]`).
- Inter-example delimiter cell: `Cell[BoxData[InterpretationBox["\t", $Line = 0; Null]], "ExampleDelimiter"]`.
- Links to other paclet symbols use `ButtonBox[..., BaseStyle->"Link", ButtonData->"paclet:<PacletName>/ref/<Other>"]`.
- StyleDefinitions at file end: whatever the paclet's existing pages use (commonly `FrontEnd\`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]` for paclets that follow the Wolfram convention). Read a sibling page to confirm; do not change it.

A.5.f Skill cascade: `documentation-writing` → `nb-writer-v2` → `nb-writer`. If `documentation-writing` defers cell-level bytes to you, use `nb-writer-v2`. If that cannot match the existing schema exactly, fall back to `nb-writer`, or hand-construct cells by reading 2–3 sibling pages with `nb-reader` and copying their cell expressions verbatim, then editing the BoxData payload.

A.5.g Save F.nb. Do **not** delete or reorder existing cells unless they are clearly broken. If you must delete, anchor the removal in the issues report.

### A.6 Verify all examples in F.nb

Run the verification harness (§6):
```
wolframscript -file "$DOCS/verify-doc-page.wls" "$REF/F.nb"
```

The harness produces `$DOCS/reports/F.report.md` with per-cell results. Iterate until clean:
- A cell that errors → fix the example (or, if the kernel itself is broken, file a kernel-bug note in the issues report and add a Possible Issues warning).
- A cell that emits an unexpected message → tighten the example to avoid the message, or use `VerificationTest` semantics (with the message expected) inside Notes — never `Quiet`.
- A cell whose re-run output diverges from the stored Output → re-run interactively to confirm; if the new output is correct, paste it; if not, debug.

### A.7 Append to issues report

Open `$DOCS/issues-report.md` and append a section (skeleton in §7). For every gap filled, every gap *not* filled, and every verification finding: cite kernel file:line and the offending Input verbatim.

Per the user's no-Quiet rule: never silently swallow issues. Surface them at the end of your turn.

### A.8 Commit

One commit per function (typical). Message format:
```
Doc page: audit and update F

- Add examples for options: <list>
- Add examples for properties: <list>
- Add examples for forms: <list>
- All examples verified via wolframscript
```

Do **not** push unless the user asks. Do **not** add `Co-Authored-By` trailers.

## 4. Branch B — Create new reference page

### B.1 Sub-triage: is F really a public, documentable symbol?

Read the kernel definition of `F` and decide:

- **Promote (continue with B.2)**: `F` is genuinely user-facing API; stable signatures and a clear use case.
- **Demote**: `F` looks like an internal helper accidentally exported. Stop; recommend moving `F` into `<PrimaryContext>PackageScope\`` or another private context. **Do not move it yourself unless explicitly asked.**
- **Defer**: `F` is experimental/unstable (e.g. lives in an `Experimental\`` sub-context). Stop; recommend either an "Experimental" tag in the page or skipping for this round.

Common red flags for "should be Demote":
- Used only inside other kernel files (no external signatures).
- Pattern-matches on internal-only argument shapes.
- Has no `::usage` and no test coverage.

If unsure, ask the user before continuing.

### B.2 Confirm or write F::usage

If `F::usage` exists, read it and confirm the listed signatures match the kernel downvalues you'll inventory in B.4. Update wording only if egregiously wrong.

If `F::usage` is missing: write one matching the existing house style. Open 2–3 existing usage strings (typically in a centralized `Kernel/Usage.m` file or directly above the symbol's downvalues) and copy the BoxData/RowBox structure. A typical entry:
```
Cell[BoxData[GridBox[{
  {"", Cell[TextData[{
    Cell[BoxData[RowBox[{TemplateBox[List[Cell[TextData["F"]],
        "paclet:<PacletName>/ref/F",
        "Wolfram Package Symbol"], "PackageLink",
        Rule[BaseStyle, "InlineFormula"]],
      "[", "args", "]"}]],
      "InlineFormula", Rule[FontFamily, "Source Sans Pro"]],
    " \[LineSeparator]<description>"}]]},
  ...
}]]]
```
List one row per form. Place the new `F::usage = "..."` either inside `F`'s declaring kernel file (above the downvalues) or in the centralized `Usage.m` if the family centralizes there.

### B.3 Read every kernel file that touches F

Same protocol as A.1.

### B.4 Inventory all forms / options / properties / named instances

Same protocol as A.2.

### B.5 Pick template pages

**Use the `nb-reader` skill** to read 2–3 existing reference pages with similar shape (similar number of forms, has options, has properties dispatch) in full. Pick by shape:
- **Constructor + properties + named instances** → a sibling like the paclet's main object type.
- **Constructor + small property set** → a focused utility symbol.
- **Predicate (returns Boolean)** → any `*Q` symbol.
- **Transform/operator with no named instances** → any binary or transform symbol.
- **Cluster wrapping built-ins** → any symbol that delegates to system functions.

Use the picked templates as scaffolding for cell types and section ordering.

### B.6 Generate F.nb

**This is the canonical use case for the `documentation-writing` skill — invoke it first.** It will plan the page (sections, prose, examples) following paclet Reference Page conventions. Then construct the cell bytes using `nb-writer-v2` (preferred; kernel-driven serialization handles BoxData, Greek/special characters, FormBox, and inline-math escapes correctly) or `nb-writer` if `nb-writer-v2` is unavailable. Use `nb-inline-math` for any inline math in ExampleText/Notes cells. Cell sequence (mandatory order; see §8 for full schema):

1. **Categorization** (Closed CellGroup): EntityType="Symbol", Paclet=`<PacletName>`, Context=`<PrimaryContext>` (or sub-context if applicable), URI=`<PacletName>/ref/F`.
2. **Keywords** (Closed): space-separated descriptive keywords.
3. **Syntax Templates** (Closed): four placeholder Template cells (leave empty if not applicable).
4. **ObjectName**: header "F" + Usage cell (one row per form) + Notes cells (with 2ColumnTableMod tables for option/property summaries).
5. **Tech Notes**: link to relevant tutorial(s) under `paclet:<PacletName>/tutorial/<Name>`. If none, leave the section with placeholder "XXXX".
6. **Related Demonstrations**: usually placeholder "XXXX".
7. **Related Links**: usually placeholder "XXXX".
8. **See Also**: TextData with ButtonBox links to sibling functions, separated by `\[EmptyVerySmallSquare]`.
9. **More About**: link to `paclet:<PacletName>/guide/<MainGuide>` (filename of the canonical Guide page in `$GUIDES`).
10. **Examples Initialization**: `Needs["<PrimaryContext>"]` (or sub-context if needed).
11. **Primary Examples**: at least 3 examples covering the most-used forms. Format: `ExampleText` → `Input` → `Output` per item, with `ExampleDelimiter` between them.
12. **Extended Examples** (CellGroupData tagged "ExtendedExamples"): subsections in this order:
    - **Scope** — one ExampleSubsection per uncommon constructor form.
    - **Generalizations & Extensions** — usually placeholder.
    - **Options** — one ExampleSubsection per option.
    - **Applications** — at least one realistic-use example if possible.
    - **Properties & Relations** — one example per property (or grouped logically).
    - **Possible Issues** — known footguns.
    - **Interactive Examples** — placeholder unless `F` integrates with `Manipulate`/etc.
    - **Neat Examples** — at least one visual or otherwise-impressive demonstration.

Required global cell options (lift from a sibling page):
- `WindowToolbars -> "MultipurposeBar"`, `Magnification -> 1.5 Inherited`, `CellContext -> CellGroup`.
- `TaggingRules -> {... "Paclet" -> <PacletName>, ... "NotebookTemplate" -> "Reference.nb", ...}`.
- `StyleDefinitions` whatever the sibling pages use.

### B.7 Verify all examples

Same as A.6.

### B.8 Update auxiliary files

- **`$PACLET/AutoCompletionData/specialArgFunctions.tr`** (if the paclet has one): if `F` accepts a small fixed set of named-string arguments (e.g. property names, named instances), add an entry. Format: see existing entries.
- **`$GUIDES/<MainGuide>.nb`**: link `F` from the appropriate guide section. Open the guide in nb-reader, find the closest sibling group, add a `ButtonBox`/link entry mirroring the existing style.

> **Do not modify `PacletInfo.wl`.** The Symbols list there is owned by the paclet maintainer; doc work does not touch it.

### B.9 Append to issues report

Open `$DOCS/issues-report.md` and append a section with:
- Heading `## F (Branch B — new page)`.
- Triage outcome and justification.
- Whether `::usage` was newly written.
- Verification results (per-cell).
- Auxiliary files updated.

### B.10 Commit

One commit per function. Message format:
```
Doc page: add reference page for F

- Wrote F::usage (was missing)  [omit line if usage already existed]
- Added Documentation/<Lang>/ReferencePages/Symbols/F.nb
- All examples verified via wolframscript
```

## 5. Project-specific kernel-file lookup table (optional, user-supplied)

If the user supplies a table mapping symbols to their declaring kernel files, use it before grepping. Symbols often live in unexpectedly-named files, and a project-specific table saves time and prevents missed cross-cuts.

Format: a markdown table with two columns — `Symbol | Required kernel files`. Every row should be exhaustive (include every file that touches the symbol, not just the primary declaration). See §11 (worked example) for the format applied to QuantumFramework.

If no table is supplied, always grep first (as in §A.1).

## 6. Verification harness specification

If `$DOCS/verify-doc-page.wls` does not exist, build it with these properties:

**Inputs**: one or more `.nb` paths (positional), optional `--out <dir>` flag (default: `$DOCS/reports/`).

**Behavior** (per page):
1. Spawn a fresh kernel via `wolframscript -script` (or by `LaunchKernels[]` in a parent session — fresh kernel discipline is essential to avoid cross-page state leakage).
2. Inside the fresh kernel: `PacletDirectoryLoad["<PACLET>"]; Needs["<PrimaryContext>"]`.
3. Get the .nb as a WL expression: `nb = Get[path]`.
4. Evaluate the **ExampleInitialization** cell first.
5. Walk `nb` for `Cell[BoxData[boxes_], "Input", ___]`. Convert each via `ToExpression[ToString[boxes, StandardForm], StandardForm, Hold]`.
6. For each held expression, evaluate inside `CheckAll[Block[{$MessageList = {}}, ReleaseHold[expr]], $MessageList]` (or use `Internal\`InheritedBlock` + `$Messages` to capture). Record:
   - The result (or `$Failed`).
   - All emitted messages.
   - The timing.
   - The cell index (for traceability).
7. Optional cross-check: find the next sibling `Cell[..., "Output"]` and compare via `MatchQ` or string comparison after `Chop` for numeric noise. Mark mismatches as `MISMATCH` (do not gate; surface in report).

**Output** (per page): `<out>/<basename>.report.md` with:
```
# Verification report for F.nb

- Total Input cells: N
- Errors: K
- Messages emitted: M
- Output mismatches: P
- Total runtime: T

## Cell <idx>
Input:
```
<verbatim>
```
Result: <stringified result, truncated>
Messages: <list>
```

**Summary** (across all pages): write `<out>/_summary.csv` with columns: page, cells, errors, messages, mismatches, seconds.

**Constraints**:
- Do **not** wrap in `Quiet`. Capture messages via `$MessageList` or `MessageHandler`.
- Implementation runs only via `wolframscript`, not the MCP evaluator.
- Parameterize the harness on `<PACLET>` and `<PrimaryContext>` so it's reusable across paclets.

## 7. Issues report format

File: `$DOCS/issues-report.md`.

Skeleton at start of project:
```
# <PacletName> — Doc Page Audit Issues Report

Branch: <branch-name>
Started: <date>
Convention: each function gets a section. Each issue cites kernel file:line and a verbatim excerpt of the offending Input cell.

## <Function name> (Branch A or B — date)
### Coverage gaps filled
- Form: <signature>  → kernel anchor: <file:line>  → page section: <Scope/Options/etc.>
- Option: <name> → kernel anchor → page section
- Property: <name> → kernel anchor → page section
- Named instance: <name> → kernel anchor → page section

### Verification findings
- Cell <idx>: <verbatim Input> → <error/message verbatim>
  Kernel anchor: <file:line>
  Resolution: <fixed | filed kernel bug | TODO with reason>

### Open issues
- <description>  (TODO marker: KERNEL-BUG, DOC-TODO, etc.)
```

Append per session. Never edit prior sessions' findings — add follow-ups as new entries with timestamps.

## 8. Reference NB schema (canonical paclet ref-page structure)

This schema mirrors what Wolfram-style paclets emit (`FunctionPageStylesExt.nb` based). Confirm by reading 3+ sibling pages from the target paclet first; if the paclet uses a custom stylesheet, lift from those siblings instead.

**Top-level Notebook options**:
- `WindowToolbars -> "MultipurposeBar"`
- `WindowSize -> {989, 963}` (or sibling-page value)
- `Magnification -> 1.5 Inherited`
- `CellContext -> CellGroup`
- `TaggingRules -> {... "Paclet" -> <PacletName>, ... "NotebookTemplate" -> "Reference.nb", "GeneratedNotebookOptions" -> ..., "Author" -> ..., "CreationDate" -> ...}`
- `StyleDefinitions -> <whatever the paclet uses>` (commonly `FrontEnd\`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]`)
- `FrontEndVersion -> "<current version>"`

**Cell sequence** (in order):

1. **History cell** with HistoryData (New in / Modified in / Obsolete in / Excised in).
2. **AuthorDate cell** (if applicable).
3. **Categorization** (Closed CellGroup): one row each for Entity Type ("Symbol"), Paclet Name (`<PacletName>`), Context (`<PrimaryContext>`), URI (`<PacletName>/ref/F`).
4. **Keywords** (Closed CellGroup).
5. **Syntax Templates** (Closed CellGroup): 4 empty placeholder Template cells.
6. **ObjectName** (CellGroup): the symbol name as header. Followed by:
   - **Usage** cell — `Cell[BoxData[GridBox[{...}]], "Usage"]`. One GridBox row per form: a `ModInfo` cell + a TextData with PackageLink to the symbol + parens with arg names + line separator + description.
   - **Notes** cells — explanatory paragraphs and `2ColumnTableMod` GridBoxes for properties/options summaries (each row: `ModInfo` + `"name"` + `TableText` description).
7. **Tech Notes** (CellGroup): list of tutorial links via `ButtonBox` to `paclet:<PacletName>/tutorial/<Name>`. If none, "XXXX".
8. **Related Demonstrations** (CellGroup): usually "XXXX".
9. **Related Links** (CellGroup): usually "XXXX".
10. **See Also** (CellGroup): `TextData[{ButtonBox[...], " \[EmptyVerySmallSquare] ", ButtonBox[...], ...}]`.
11. **More About** (CellGroup): `ButtonBox` to `paclet:<PacletName>/guide/<MainGuide>`.
12. **Examples Initialization** (CellGroup): TemplateBox info icon + `Cell[BoxData[RowBox[{"Needs", "[", "\"<PrimaryContext>\"", "]"}]], "ExampleInitialization", InitializationCell -> True]`.
13. **PrimaryExamplesSection** (CellGroup): header + `InterpretationBox[...]` "More Examples" link + multiple `ExampleText` → `Input` → `Output` triples separated by `ExampleDelimiter`.
14. **ExtendedExamples** (CellGroup, tagged `"ExtendedExamples"`): the eight subsections — **Scope**, **Generalizations & Extensions**, **Options**, **Applications**, **Properties & Relations**, **Possible Issues**, **Interactive Examples**, **Neat Examples**. Each is a `Cell[..., "ExampleSection"]` wrapping a `CellGroupData[{...ExampleSubsection cells...}]`. Subsections that are empty contain the placeholder text "XXXX".

**Cell-level options always present**:
- `CellID -> <fresh integer>`
- `ExpressionUUID -> "<CreateUUID[]>"`

## 9. Excluded symbols

User-supplied list (if any) of symbols that are public exports but should not be documented on this branch. If `F` is on the list, stop and tell the user.

The list typically includes:
- Internal helpers that should eventually be moved to PackageScope but for now remain public.
- Experimental symbols whose API is too unstable to document.
- Symbols already documented under a different track (e.g. their own ongoing-project subdirectory).

## 10. Common pitfalls

10.1 **Symbols may not live in the file you'd guess.** Always grep first if no project-specific lookup table is supplied. Examples (paclet-dependent): a `*Q` predicate often lives in the file that defines the type it tests, not in its own file; a transform symbol may share a file with its companion transforms.

10.2 **Memoize wrappers** in `Init.m`-style files can produce stale results during verification. If the paclet uses `Memoize` (or any caching layer), clear it between verification runs by killing the kernel (per-page fresh kernel handles this automatically — §6).

10.3 **Symbols that are pure delegators**. If `F` simply calls another function, its page can be templated heavily, but verify the kernel actually delegates as you assume (e.g. argument order, units like radians vs degrees, default options).

10.4 **Stray pages for non-paclet symbols.** If a `.nb` exists in `$REF` but the symbol is a system built-in (or comes from a different paclet), treat it as audit-only: read the page, decide whether it adds paclet-specific value, and recommend keep / remove / convert-to-tutorial. Do not generate a new ref page for a built-in unless the user explicitly directs.

10.5 **Never `Quiet`** — never. If a verification run produces messages, capture them. If a kernel has a real bug producing an unwanted message, file it in the issues report and add a Possible Issues note explaining it.

10.6 **Round-trip discipline** — for any constructor/accessor pair `C`/`A`, the FIRST property example on the page should be `A[C[x]] === x`. Phase-oblivious / set-wise comparisons come after, opt-in via a named helper.

10.7 **Never use `Co-Authored-By: Claude` trailers.** Author commits with the user's identity only.

## 11. Worked example — invoking the runbook for QuantumFramework

This shows how a user would specialize the runbook for a specific paclet.

**Invocation:**

> *"Follow `<path>/per-function-runbook.md` for function `QuantumChannel` in the paclet at `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework`. Project-specific deep skill: `qf`. Reports directory: `OngoingProjects/Improving doc pages/`. Excluded symbols and project lookup table: see §11 of the runbook."*

**Resolved configuration** (auto-derived from `PacletInfo.wl`):
- `$PACLET` = `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework`
- `$KER` = `$PACLET/Kernel`
- `$REF` = `$PACLET/Documentation/English/ReferencePages/Symbols/`
- `$GUIDES` = `$PACLET/Documentation/English/Guides/`
- `<PrimaryContext>` = `Wolfram\`QuantumFramework\``
- `<Subcontexts>` = `Wolfram\`QuantumFramework\`Gates\``, `Wolfram\`QuantumFramework\`Experimental\``, `Wolfram\`QuantumFramework\`QuantumOptimization\``, etc.
- `<PacletName>` = `Wolfram/QuantumFramework`
- `<MainGuide>` = `WolframQuantumComputationFramework`
- `<DeepSkill>` = `qf`
- `$DOCS` = `$PACLET/../OngoingProjects/Improving doc pages/`

**Project-specific exclusion list** (do not document on this branch):
- `Kernel/QuantumEvolution.m`: `HamiltonianTransitionRate`, `LindbladTransitionRates`
- `Kernel/ZX.m`: `ZXTensorNetwork`, `ZXTensorNetworkQuantumCircuit`, `ZXExpression`
- `Kernel/Multiway.m`: `QuantumCircuitMultiwayCausalGraph`, `QuantumCircuitPathGraph`, `QuantumCircuitTokenEventGraph`
- `Kernel/QuantumCircuitOperator/TensorNetwork.m`: `TensorNetworkQuantumCircuit`
- `Kernel/Experimental.m`: `DecomposedQuantumStateProbabilities`, `QuantumBeamSearch`, `QuantumDiagramProcess`, `QuantumMPO`, `QuantumMPSApply`, `FeynmanBacktracking`
- `PauliStabilizer` (deferred to the Stabilizer track at `OngoingProjects/Stabilizer/`)

**Project-specific kernel-file lookup table:**

| Symbol | Required kernel files |
|---|---|
| QuantumState | `QuantumState/QuantumState.m`, `QuantumState/Properties.m`, `QuantumState/Formatting.m`, `QuantumState/NamedStates.m`; plus `Init.m`, `Labels.m`, `Visualization.m`, `Utilities.m` |
| QuantumOperator | `QuantumOperator/QuantumOperator.m`, `Properties.m`, `Formatting.m`, `NamedOperators.m`; plus `Gates.m`, `Labels.m`, `Visualization.m` |
| QuantumBasis | `QuantumBasis/QuantumBasis.m`, `Properties.m`, `Formatting.m` |
| QuantumChannel | `QuantumChannel/QuantumChannel.m`, `Properties.m`, `Formatting.m`, `NamedChannels.m` |
| QuantumCircuitOperator | `QuantumCircuitOperator/QuantumCircuitOperator.m`, `Properties.m`, `Formatting.m`, `NamedCircuits.m`, `Draw.m`, `TensorNetwork.m`, `Qiskit.m`, `QuEST.m`, `Classiq.m` |
| QuantumMeasurement | `QuantumMeasurement.m`; plus `QuantumMeasurementOperator/*` for sibling behavior |
| QuantumMeasurementOperator | `QuantumMeasurementOperator/QuantumMeasurementOperator.m`, `Properties.m`, `Formatting.m`, `NamedMeasurementOperators.m` |
| QuantumMeasurementSimulation | `QuantumStateEstimation.m`; plus `QuantumMeasurementOperator/*` |
| QuantumStateEstimate / QuantumStateEstimation / QuantumStateSampler | `QuantumStateEstimation.m` |
| QuditBasis | `QuditBasis/QuditBasis.m`, `Properties.m`, `NamedBases.m`, `QuditName.m` |
| QuditName | `QuditBasis/QuditName.m` |
| QuantumDistance / QuantumSimilarity | `QuantumDistance.m` (shared file) |
| QuantumEntangledQ / QuantumEntanglementMonotone | `QuantumEntanglement.m` |
| QuantumEvolve | `QuantumEvolution.m` |
| QuantumPartialTrace | `QuantumPartialTrace.m`; plus per-family `Properties.m` for delegation |
| QuantumPartialTranspose | `QuantumPartialTranspose.m` |
| QuantumTensorProduct | `QuantumTensorProduct.m`; plus per-family `Properties.m` for delegation |
| QuantumWignerTransform / QuantumWeylTransform | `QuantumWignerTransform.m`; plus `QuantumPhaseSpaceTransform.m` (shared dispatch) |
| QuantumWignerMICTransform | `MIC.m` |
| QuantumPhaseSpaceTransform / QuantumPositiveTransform | `QuantumPhaseSpaceTransform.m`; plus `QuantumWignerTransform.m` |
| QuantumCircuitMultiwayGraph | `Multiway.m` |
| QuantumShortcut | `Labels.m` |
| QuantumMPS | `Experimental.m` |
| ImportQASMCircuit / QiskitCircuit / ImportQPY | `QuantumCircuitOperator/Qiskit.m` |
| 27 Gate symbols + QuantumGateQ | `Gates.m` |
| QuantumOptimization symbols (`ParametrizedLayer`, `EntanglementLayer`, `GenerateParameters`, `FubiniStudyMetricTensor`, `FubiniStudyMetricTensorLayers`, `SPSRGradientValues`, `ASPSRGradientValues`, `ClassiqSetup`, `GradientDescent`, `QuantumLinearSolve`, `QuantumLockingMechanism`, `QuantumUnlockingMechanism`, `QuantumNaturalGradientDescent`) | `QuantumOptimization.m` (single file; symbol-by-symbol dispatch — read whole file) |

**Project-specific pitfalls** (in addition to the generic §10):

- The 27 gate symbols (`XGate`, `YGate`, ..., `RZZGate`, `U1Gate`–`UGate`) all delegate to `QuantumOperator` with named matrices — pages can be templated.
- `EinsteinSummation.nb` exists in `$REF` but is a built-in WL symbol, not a paclet export — see §10.4.
- `Wolfram\`QuantumFramework\`Init\`` uses `Memoize[FromOperatorShorthand]` — relevant when verifying memoization-sensitive examples.
