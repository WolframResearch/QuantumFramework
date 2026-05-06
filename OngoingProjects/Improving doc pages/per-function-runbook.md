
# QuantumFramework — Per-Function Doc Page Runbook

Hand this file back to Claude with: *"Follow `OngoingProjects/Improving doc pages/per-function-runbook.md` for function `F`."* Replace `F` with a paclet symbol (e.g. `QuantumChannel`, `QuantumStateEstimation`, `XGate`).

The runbook handles both cases:
- **Branch A** — `F` already has `Documentation/English/ReferencePages/Symbols/F.nb`. Audit it against the kernel and fill gaps.
- **Branch B** — no such file exists. Triage `F`, then (if confirmed public) write `::usage` and create the page.

> **Read carefully, do not guess.** Per `feedback_reading_protocol.md` (the user's auto-memory): inventory first, file:line anchors on every claim, verbatim excerpts. No summaries of summaries. If a kernel symbol's behavior depends on a file you have not yet read, stop and read it.

## 0. Session preamble (always perform first)

0.1 Confirm working tree:
- Branch should be `cleaning-and-editing-doc-pages` (created from `origin/main`). Run `git status` and `git rev-parse --abbrev-ref HEAD` to confirm. If on a different branch, ask the user before switching.
- Paclet root: `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework` (referred to below as `$PACLET`).
- Reference pages dir: `$PACLET/Documentation/English/ReferencePages/Symbols/` (referred to as `$REF`).
- Kernel dir: `$PACLET/Kernel/` (referred to as `$KER`).

0.2 Read the user's auto-memory once at session start: `/Users/mohammadb/.claude/projects/-Users-mohammadb-Documents-GitHub-QuantumFramework/memory/MEMORY.md` and any referenced files. Important rules to honor:
- Use `wolframscript` for kernel evaluation (not the WolframLanguageEvaluator MCP). Long runs go via `run_in_background`.
- Never `Quiet` to suppress messages. Capture them or surface them.
- For any constructor/accessor pair the FIRST round-trip test must be exact equality `A[C[x]] === x`; only then add phase-oblivious checks.

0.3 Confirm `wolframscript` works: `wolframscript -code 'PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"]; Needs["Wolfram`QuantumFramework`"]; Length[Names["Wolfram`QuantumFramework`*"]]'` should print a number ≥ 41.

0.4 **Skill stack — invoke these as the work demands. They are not optional dressing; they encode conventions you would otherwise reinvent badly.**

| Skill | When to use |
|---|---|
| `qf` | **Always active for this work.** Encodes the QuantumFramework paclet's pitfalls (paclet 1.6.5 anchor SHA `cbbe9368`, source-audit reference). Trumps the broader `wl-quantum-computing` skill. |
| `documentation-writing` | **Primary skill for writing or substantively editing any reference page.** It knows WFR/paclet Reference Page conventions, Guide page structure, and Tech Note tutorial format. Invoke at the start of Branch A §A.5 (filling gaps) and Branch B §B.6 (creating a new page). |
| `nb-writer-v2` | Preferred for low-level cell construction inside the page that `documentation-writing` produces or instructs. Kernel-driven serialization handles BoxData, RowBox, CellGroupData, FormBox, Greek/special characters, and inline-math escapes correctly. Use when `documentation-writing` defers to you for the actual notebook bytes. |
| `nb-writer` | Fallback if `nb-writer-v2` is unavailable or overkill for a one-cell edit. |
| `nb-reader` | For reading existing `.nb` files (the existing reference page on Branch A, sibling templates on Branch B). Faster and more accurate than raw `Read` on serialized notebooks. |
| `nb-inline-math` | When adding inline math (e.g. `\(\rho\)`) inside a Text cell or ExampleText cell. Prevents the bracket-nesting and Greek-character pitfalls in BoxData. |
| `wl-native-style` | When writing example code. Prefer functional forms (`Map`, `Fold`, `Table`, `Association`) over `Do`/`For`/`While`/`AppendTo`. Per the user's `CLAUDE.md`. |
| `wl-broadcasting` | When example code multiplies arrays of mismatched rank. WL is row-wise, not NumPy-style. |
| `learning-by-computing` | Optional — only if a Reference Page's Notes section grows long enough to warrant explanation-then-computation flow. Most pages do not need this. |
| `wl-quantum-computing`, `wl-symbolic-computation`, `wl-machine-learning`, `wl-visualization`, `wl-data-and-knowledge` | Use as background reference when a specific example draws from those areas. Do not let them override `qf` for QF-specific behavior. |

If you find yourself hand-rolling `.nb` cell expressions instead of using `documentation-writing` + `nb-writer-v2`, stop and switch to the skills.

0.5 Confirm the verification harness exists. Path: `$PACLET/../OngoingProjects/Improving doc pages/verify-doc-page.wls`. **If missing, build it before any work that needs verification — see §6.** If present, smoke-test it on one existing page (e.g. `QuantumDistance.nb`) and confirm it produces a markdown report.

0.6 Confirm the issues report file exists. Path: `$PACLET/../OngoingProjects/Improving doc pages/issues-report.md`. If missing, create it with a header skeleton (see §7).

## 1. Determine the branch

1.1 Check that `F` is a public paclet export:
```
wolframscript -code 'PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
  Needs["Wolfram`QuantumFramework`"];
  ctxs = {"Wolfram`QuantumFramework`","Wolfram`QuantumFramework`Gates`",
          "Wolfram`QuantumFramework`Experimental`",
          "Wolfram`QuantumFramework`QuantumOptimization`"};
  Select[Flatten[Names[#<>"F"]&/@ctxs], #=!="Wolfram`QuantumFramework`F"&]'
```
(Replace `F` with the actual symbol.) If no context contains `F`, stop and tell the user `F` is not a public export.

1.2 Check the **Excluded list (§9)**. If `F` is excluded, stop and tell the user this symbol is intentionally undocumented.

1.3 Check existence of `$REF/F.nb`:
- Exists → **Branch A** (§3, "Update existing").
- Does not exist → **Branch B** (§4, "Create new").

## 2. Glossary (terms used below)

- **Form**: a distinct constructor/call signature recognized by `F`. E.g. `QuantumState[name]`, `QuantumState[<|...|>]`, `QuantumState[v, basis]` are three forms.
- **Option**: a key in `Options[F]`. Listed via `Options[F]` after the paclet is loaded; declared in kernel as `Options[F] = {...}`.
- **Property**: a string `F[obj]["propname"]` understands. Lives in dispatch tables, usually in `Kernel/<Family>/Properties.m`.
- **Named instance**: a named-string form like `QuantumState["Bell"]`, `QuditBasis["Wigner"]`, `QuantumChannel["Depolarizing"]`. Lives in `Kernel/<Family>/Named*.m`.
- **Sibling symbol**: a symbol declared in the same kernel file. Often shares behavior. See §5 lookup table.

## 3. Branch A — Update existing reference page

### A.1 Read every kernel file that touches F

Use the **Kernel lookup table (§5)**. If `F` is not in the table, run:
```
grep -ln 'F::usage\|^F\[' $KER/*.m $KER/*/*.m
```
to locate every file that mentions `F`. Read all of them — declarations may be split across the family directory plus shared files (`Init.m`, `Usage.m`, `Utilities.m`, `Labels.m`, `Visualization.m`, `Graphics.m`, `Formatting.m`, `Multiway.m`, `MIC.m`).

For each file you read, write down (mentally or in a scratchpad):
- The line ranges defining `F` downvalues.
- Every `OptionValue["..."]` callsite (reveals undeclared options).
- Every property dispatch line (e.g. `obj["Amplitudes"] := ...`, `Switch[prop, "Eigenvalues", ...]`).
- Every named-instance handler (e.g. `QuantumState["Bell", ...] := ...`).

### A.2 Inventory the four coverage axes

For `F` produce four lists with file:line anchors:

A.2.a **Forms** — every distinct downvalue signature. Cross-check against `F::usage` text; if usage lags behind code, the code wins.

A.2.b **Options** — run `wolframscript -code 'Needs["Wolfram`QuantumFramework`"]; Options[F]'` for the canonical declared set. Then grep for `OptionValue` callsites in F's files to catch any options consumed but not declared (rare bug pattern).

A.2.c **Properties** — list every string `F` accepts in property position. Grep patterns to use:
- `"\"[A-Z][A-Za-z]*\""` near `prop_String`, `MemberQ[$Properties`, `Switch[prop`, `Replace[prop`, or `obj[("name1"|"name2"|...)] :=`.
- For families with a `Properties.m` file, that's the primary source. Don't miss properties added in `Formatting.m` or `Visualization.m`.

A.2.d **Named instances** — for symbols whose `Named*.m` file exists, list every named class. Group by what they share (e.g. all 1-qubit gates, all 2-qubit Bell-style entanglers).

### A.3 Read the existing F.nb

**Use the `nb-reader` skill** to read `$REF/F.nb`. (Raw `Read` on a serialized .nb is brittle; `nb-reader` converts to Markdown and is the right tool.) The file is large; read it in slices if needed. Inventory:
- Existing input forms shown in the **Usage** cell (count `ModInfo` rows).
- Existing examples in **Primary Examples** (count Input cells).
- Existing **Scope** subsections (per uncovered form).
- Existing **Options** subsections (one ExampleSubsection per option).
- Existing **Properties & Relations** examples.
- Existing **Possible Issues** notes.
- Whether **Generalizations & Extensions**, **Applications**, **Interactive Examples**, **Neat Examples** sections exist (often empty placeholders).

### A.4 Diff: produce a gap list

For each axis (Forms, Options, Properties, Named instances), list items that are in §A.2 but not §A.3. Anchor each gap with the kernel file:line where it's defined.

If the gap list is empty, jump to §A.6 (verification only).

### A.5 Edit F.nb to fill gaps

**Invoke the `documentation-writing` skill before any edits.** It knows the paclet Reference Page conventions and will guide section ordering, link styles, and prose voice. Use the **NB schema (§8)** as the structural cross-check. For the actual cell-byte construction inside the skill's plan, use `nb-writer-v2` (it handles BoxData, Greek characters, and inline math correctly). For inline math inside ExampleText/Notes cells, use `nb-inline-math`. Specific rules:

A.5.a **For each missing form**: add an `ExampleText` describing the form, then an `Input` cell exercising it, then expected `Output`. Place under **Scope** (`ExampleSection` with `CellTags -> {"Scope", "ExtendedExamples"}`).

A.5.b **For each missing option**: add an `ExampleSubsection` titled with the option name under the existing **Options** section. Body: short `ExampleText` explaining the option, an `Input` showing default vs explicit value, and the resulting `Output`.

A.5.c **For each missing property**: add either (a) a brief `ExampleText` + `Input` + `Output` triple under **Properties & Relations**, or (b) an entry in the **Notes** 2ColumnTableMod table at the top of the page if it's better summarized than demonstrated. Default: at least one demonstration in Examples for every property.

A.5.d **For each missing named instance class**: add at least one `Input` example per group under **Scope**.

A.5.e **Style discipline**:
- Cell style names exactly: `"PrimaryExamplesSection"`, `"ExampleSection"`, `"ExampleSubsection"`, `"ExampleText"`, `"Input"`, `"Output"`, `"ExampleDelimiter"`. Do not invent new style names.
- Every cell needs `CellID` (use a fresh integer; do not collide with existing IDs in the file) and `ExpressionUUID` (use `CreateUUID[]`).
- Inter-example delimiter cell: `Cell[BoxData[InterpretationBox["\t", $Line = 0; Null]], "ExampleDelimiter"]`.
- Links to other paclet symbols use `ButtonBox[..., BaseStyle->"Link", ButtonData->"paclet:Wolfram/QuantumFramework/ref/<Other>"]`.
- StyleDefinitions at file end: `FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]`. Do not change this.

A.5.f Skill cascade: `documentation-writing` → `nb-writer-v2` → `nb-writer`. If `documentation-writing` produces a plan but defers cell-level bytes to you, use `nb-writer-v2`. If `nb-writer-v2` cannot match the existing schema exactly, fall back to `nb-writer`, or hand-construct cells by reading 2–3 sibling pages with `nb-reader` and copying their cell expressions verbatim, then editing the BoxData payload.

A.5.g Save F.nb. Do **not** delete or reorder existing cells unless they are clearly broken (e.g. a Notes table row pointing to a property that no longer exists). If you must delete, anchor the removal in the issues report.

### A.6 Verify all examples in F.nb

Run the verification harness (§6):
```
wolframscript -file "$PACLET/../OngoingProjects/Improving doc pages/verify-doc-page.wls" "$REF/F.nb"
```

The harness produces `OngoingProjects/Improving doc pages/reports/F.report.md` with per-cell results. Iterate until clean:
- A cell that errors → fix the example (or, if the kernel itself is broken, file a separate kernel-bug note in the issues report and skip that example with an `ExampleText` warning).
- A cell that emits an unexpected message → either tighten the example to avoid the message, or use `VerificationTest` semantics (with the message expected) inside Notes — never `Quiet`.
- A cell whose re-run output diverges from the stored Output → re-run interactively to confirm; if the new output is correct, paste it; if not, debug.

### A.7 Append to issues report

Open `OngoingProjects/Improving doc pages/issues-report.md` and append a section with:
- Heading `## F (Branch A — updated)`.
- For each cell that errored or emitted messages during verification: cell index in the .nb, the offending Input string verbatim, the message/error verbatim, file:line anchors back to the kernel definition that controls the behavior.
- For each gap filled: a one-liner saying which axis and which item, with kernel file:line.
- For each gap *not* filled (with reason): a one-liner with reason and a TODO marker.

Per `feedback_no_quiet.md`: never silently swallow issues. If something is broken in the kernel and you skip the example, surface it explicitly to the user at the end of your turn.

### A.8 Commit

One commit per function (typical). Message format:
```
Doc page: audit and update F

- Add examples for options: <list>
- Add examples for properties: <list>
- Add examples for forms: <list>
- All examples verified via wolframscript
```

Do **not** push unless the user asks. Do **not** add `Co-Authored-By` trailers (per the user's `CLAUDE.md`).

## 4. Branch B — Create new reference page

### B.1 Sub-triage: is F really a public, documentable symbol?

Read the kernel definition of `F` and decide:

- **Promote (continue with B.2)**: `F` is genuinely user-facing API; it has stable signatures and a clear use case.
- **Demote**: `F` looks like an internal helper accidentally exported. Stop; report to user with a recommendation to move `F` into `Wolfram`QuantumFramework`PackageScope`` or another private context. **Do not move it yourself in this branch unless explicitly asked.**
- **Defer**: `F` is experimental/unstable (e.g. lives in `Wolfram`QuantumFramework`Experimental``). Stop; report to user with a recommendation to add a "Experimental" tag in the page or skip for this round.

Common red flags for "should be Demote":
- Used only inside other kernel files (no examples, no externally facing signatures).
- Pattern-matches on internal-only argument shapes.
- Has no `::usage` and no test coverage.

If you're not sure, ask the user before continuing.

### B.2 Confirm or write F::usage

If `F::usage` exists, read it and confirm the listed signatures match the kernel downvalues you'll inventory in B.4. Update wording only if egregiously wrong.

If `F::usage` is missing: write one matching the existing house style. Open 2–3 existing usage strings (e.g. from `$KER/Usage.m` or directly above the symbol's downvalues in its own file) and copy the BoxData/RowBox structure. Each usage entry has the form:
```
Cell[BoxData[GridBox[{
  {"", Cell[TextData[{
    Cell[BoxData[RowBox[{TemplateBox[List[Cell[TextData["F"]],
        "paclet:Wolfram/QuantumFramework/ref/F",
        "Wolfram Package Symbol"], "PackageLink",
        Rule[BaseStyle, "InlineFormula"]],
      "[", "args", "]"}]],
      "InlineFormula", Rule[FontFamily, "Source Sans Pro"]],
    " \[LineSeparator]<description>"}]]},
  ...
}]]]
```
List one row per form. Place the new `F::usage = "..."` either inside `F`'s declaring kernel file (above the downvalues) or in `Kernel/Usage.m` if the family's other usage strings are centralized there.

### B.3 Read every kernel file that touches F

Same protocol as A.1.

### B.4 Inventory all forms / options / properties / named instances

Same protocol as A.2.

### B.5 Pick template pages

**Use the `nb-reader` skill** to read 2–3 existing reference pages with similar shape (similar number of forms, has options, has properties dispatch) in full. Examples by shape:
- **Constructor + properties + named instances** → `QuantumState.nb`, `QuantumOperator.nb`.
- **Constructor + small property set** → `QuantumDistance.nb`.
- **Predicate (returns Boolean)** → `QuantumEntangledQ.nb`.
- **Transform/operator with no named instances** → `QuantumPartialTrace.nb`, `QuantumTensorProduct.nb`.
- **Cluster wrapping Wolfram built-ins** → `QuantumWignerTransform.nb`.

Use the picked templates as scaffolding for cell types and section ordering.

### B.6 Generate F.nb

**This is the canonical use case for the `documentation-writing` skill — invoke it first.** It will plan the page (sections, prose, examples) following paclet Reference Page conventions. Then construct the cell bytes using `nb-writer-v2` (preferred; kernel-driven serialization handles BoxData, Greek/special characters, FormBox, and inline-math escapes correctly) or `nb-writer` if `nb-writer-v2` is unavailable. Use `nb-inline-math` for any inline math in ExampleText/Notes cells. Cell sequence (mandatory order; see §8 for full schema):

1. **Categorization** (Closed CellGroup): EntityType="Symbol", Paclet="Wolfram/QuantumFramework", Context="Wolfram`QuantumFramework`" (or sub-context if applicable), URI="Wolfram/QuantumFramework/ref/F".
2. **Keywords** (Closed): space-separated descriptive keywords.
3. **Syntax Templates** (Closed): four placeholder Template cells (leave empty if not applicable).
4. **ObjectName**: header "F" + Usage cell (one row per form) + Notes cells (with 2ColumnTableMod tables for option/property summaries).
5. **Tech Notes**: link to relevant tutorial(s) under `paclet:Wolfram/QuantumFramework/tutorial/<Name>`. If none, leave the section with placeholder "XXXX".
6. **Related Demonstrations**: usually placeholder "XXXX".
7. **Related Links**: usually placeholder "XXXX".
8. **See Also**: TextData with ButtonBox links to sibling functions, separated by `\[EmptyVerySmallSquare]`.
9. **More About**: link to `paclet:Wolfram/QuantumFramework/guide/WolframQuantumComputationFramework`.
10. **Examples Initialization**: `Needs["Wolfram`QuantumFramework`"]` (or sub-context if needed).
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
- `TaggingRules -> {"Author" -> "Wolfram Research", "CreationDate" -> "...", "Paclet" -> "Wolfram/QuantumFramework", ...}`.
- `StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]`.

### B.7 Verify all examples

Same as A.6.

### B.8 Update auxiliary files

- **AutoCompletionData/specialArgFunctions.tr**: if `F` accepts a small fixed set of named-string arguments (e.g. property names, named instances), add an entry. Format: see existing entries.
- **Documentation/English/Guides/WolframQuantumComputationFramework.nb**: link `F` from the appropriate guide section. Open the guide in nb-reader, find the closest sibling group, add a `ButtonBox`/link entry mirroring the existing style.

> **Do not modify `PacletInfo.wl`.** The Symbols list there is owned by the paclet maintainer; doc work does not touch it.

### B.9 Append to issues report

Open `OngoingProjects/Improving doc pages/issues-report.md` and append a section with:
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
- Added Documentation/English/ReferencePages/Symbols/F.nb
- All examples verified via wolframscript
```

## 5. Kernel file lookup table

Use this to know which files to read for `F`. **If F is not listed here, grep first; do not guess.**

| Symbol | Required kernel files |
|---|---|
| QuantumState | `QuantumState/QuantumState.m`, `QuantumState/Properties.m`, `QuantumState/Formatting.m`, `QuantumState/NamedStates.m`; plus `Init.m`, `Labels.m`, `Visualization.m`, `Utilities.m` |
| QuantumOperator | `QuantumOperator/QuantumOperator.m`, `QuantumOperator/Properties.m`, `QuantumOperator/Formatting.m`, `QuantumOperator/NamedOperators.m`; plus `Gates.m`, `Labels.m`, `Visualization.m` |
| QuantumBasis | `QuantumBasis/QuantumBasis.m`, `QuantumBasis/Properties.m`, `QuantumBasis/Formatting.m` |
| QuantumChannel | `QuantumChannel/QuantumChannel.m`, `QuantumChannel/Properties.m`, `QuantumChannel/Formatting.m`, `QuantumChannel/NamedChannels.m` |
| QuantumCircuitOperator | `QuantumCircuitOperator/QuantumCircuitOperator.m`, `Properties.m`, `Formatting.m`, `NamedCircuits.m`, `Draw.m`, `TensorNetwork.m`, `Qiskit.m`, `QuEST.m`, `Classiq.m` |
| QuantumMeasurement | `QuantumMeasurement.m`; plus `QuantumMeasurementOperator/*` for sibling behavior |
| QuantumMeasurementOperator | `QuantumMeasurementOperator/QuantumMeasurementOperator.m`, `Properties.m`, `Formatting.m`, `NamedMeasurementOperators.m` |
| QuantumMeasurementSimulation | `QuantumStateEstimation.m` (declared here); plus `QuantumMeasurementOperator/*` |
| QuantumStateEstimate | `QuantumStateEstimation.m` (sibling of `QuantumStateEstimation`, `QuantumStateSampler`) |
| QuantumStateEstimation | `QuantumStateEstimation.m` |
| QuantumStateSampler | `QuantumStateEstimation.m` |
| QuditBasis | `QuditBasis/QuditBasis.m`, `Properties.m`, `NamedBases.m`, `QuditName.m` |
| QuditName | `QuditBasis/QuditName.m` |
| QuantumDistance | `QuantumDistance.m` (also defines `QuantumSimilarity`) |
| QuantumSimilarity | `QuantumDistance.m` (shared file) |
| QuantumEntangledQ | `QuantumEntanglement.m` |
| QuantumEntanglementMonotone | `QuantumEntanglement.m` |
| QuantumEvolve | `QuantumEvolution.m` |
| QuantumPartialTrace | `QuantumPartialTrace.m`; plus per-family `Properties.m` for delegation |
| QuantumPartialTranspose | `QuantumPartialTranspose.m` |
| QuantumTensorProduct | `QuantumTensorProduct.m`; plus per-family `Properties.m` for delegation |
| QuantumWignerTransform | `QuantumWignerTransform.m`; plus `QuantumPhaseSpaceTransform.m` (shared dispatch); also defines `QuantumWeylTransform` |
| QuantumWeylTransform | `QuantumWignerTransform.m` |
| QuantumWignerMICTransform | `MIC.m` |
| QuantumPhaseSpaceTransform | `QuantumPhaseSpaceTransform.m` (also defines `QuantumPositiveTransform`); plus `QuantumWignerTransform.m` |
| QuantumPositiveTransform | `QuantumPhaseSpaceTransform.m` |
| QuantumCircuitMultiwayGraph | `Multiway.m` |
| QuantumShortcut | `Labels.m` |
| QuantumMPS | `Experimental.m` |
| ImportQASMCircuit | `QuantumCircuitOperator/Qiskit.m` |
| QiskitCircuit | `QuantumCircuitOperator/Qiskit.m` |
| ImportQPY | `QuantumCircuitOperator/Qiskit.m` |
| Gates: XGate, YGate, ZGate, HGate, SGate, TGate, PhaseGate, SXGate, SXdgGate, IdentityGate, CXGate, CYGate, CZGate, SWAPGate, CCXGate, CSWAPGate, RXGate, RYGate, RZGate, RXXGate, RYYGate, RZZGate, RZXGate, U1Gate, U2Gate, U3Gate, UGate | `Gates.m` (lines 37–66); each delegates to `QuantumOperator` |
| QuantumGateQ | `Gates.m:69` |
| QuantumOptimization symbols (`ParametrizedLayer`, `EntanglementLayer`, `GenerateParameters`, `FubiniStudyMetricTensor`, `FubiniStudyMetricTensorLayers`, `SPSRGradientValues`, `ASPSRGradientValues`, `ClassiqSetup`, `GradientDescent`, `QuantumLinearSolve`, `QuantumLockingMechanism`, `QuantumUnlockingMechanism`, `QuantumNaturalGradientDescent`) | `QuantumOptimization.m` (single file; symbol-by-symbol dispatch — read whole file) |
| PauliStabilizer | `PauliStabilizer.m`; **but defer to the Stabilizer track** — see `OngoingProjects/Stabilizer/` for active doc work; do not duplicate |

## 6. Verification harness specification

If `OngoingProjects/Improving doc pages/verify-doc-page.wls` does not exist, build it with these properties:

**Inputs**: one or more `.nb` paths (positional), optional `--out <dir>` flag (default: `OngoingProjects/Improving doc pages/reports/`).

**Behavior** (per page):
1. Spawn a fresh kernel via `wolframscript -script` (or by `LaunchKernels[]` in a parent session — fresh kernel discipline is essential to avoid cross-page state leakage).
2. Inside the fresh kernel: `PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"]; Needs["Wolfram`QuantumFramework`"]`.
3. Get the .nb as a WL expression: `nb = Get[path]`.
4. Evaluate the **ExampleInitialization** cell first.
5. Walk `nb` for `Cell[BoxData[boxes_], "Input", ___]`. Convert each via `ToExpression[ToString[boxes, StandardForm], StandardForm, Hold]`.
6. For each held expression, evaluate inside `CheckAll[Block[{$MessageList = {}}, ReleaseHold[expr]], $MessageList]` (or use `Internal`InheritedBlock` + `$Messages` to capture). Record:
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
- Per `feedback_no_quiet.md`: do **not** wrap in `Quiet`. Capture messages via `$MessageList` or `MessageHandler`.
- Per `feedback_use_wolframscript.md`: implementation runs only via `wolframscript`, not the MCP evaluator.
- Reference: prior style for Stabilizer audit harnesses lives in `OngoingProjects/Stabilizer/`. Mirror their idioms where useful.

## 7. Issues report format

File: `OngoingProjects/Improving doc pages/issues-report.md`.

Skeleton at start of project:
```
# QuantumFramework — Doc Page Audit Issues Report

Branch: cleaning-and-editing-doc-pages
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

## 8. Reference NB schema (canonical structure)

Every QuantumFramework reference page (lifted from `QuantumState.nb`, `QuantumOperator.nb`, `QuantumDistance.nb`, `QuantumPartialTrace.nb`, `QuantumTensorProduct.nb`):

**Top-level Notebook options**:
- `WindowToolbars -> "MultipurposeBar"`
- `WindowSize -> {989, 963}` (or similar)
- `Magnification -> 1.5 Inherited`
- `CellContext -> CellGroup`
- `TaggingRules -> {... "Paclet" -> "Wolfram/QuantumFramework", ... "NotebookTemplate" -> "Reference.nb", "GeneratedNotebookOptions" -> ..., "Author" -> ..., "CreationDate" -> ...}`
- `StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]`
- `FrontEndVersion -> "<current version>"`

**Cell sequence** (in order):

1. **History cell** with HistoryData (New in / Modified in / Obsolete in / Excised in).
2. **AuthorDate cell** (if applicable).
3. **Categorization** (Closed CellGroup): one row each for Entity Type ("Symbol"), Paclet Name ("Wolfram/QuantumFramework"), Context, URI ("Wolfram/QuantumFramework/ref/F").
4. **Keywords** (Closed CellGroup).
5. **Syntax Templates** (Closed CellGroup): 4 empty placeholder Template cells.
6. **ObjectName** (CellGroup): the symbol name as header. Followed by:
   - **Usage** cell — `Cell[BoxData[GridBox[{...}]], "Usage"]`. One GridBox row per form: a `ModInfo` cell + a TextData with PackageLink to the symbol + parens with arg names + line separator + description.
   - **Notes** cells — explanatory paragraphs and `2ColumnTableMod` GridBoxes for properties/options summaries (each row: `ModInfo` + `"name"` + `TableText` description).
7. **Tech Notes** (CellGroup): list of tutorial links via `ButtonBox` to `paclet:Wolfram/QuantumFramework/tutorial/<Name>`. If none, "XXXX".
8. **Related Demonstrations** (CellGroup): usually "XXXX".
9. **Related Links** (CellGroup): usually "XXXX".
10. **See Also** (CellGroup): `TextData[{ButtonBox[...], " \[EmptyVerySmallSquare] ", ButtonBox[...], ...}]`.
11. **More About** (CellGroup): `ButtonBox` to `paclet:Wolfram/QuantumFramework/guide/WolframQuantumComputationFramework`.
12. **Examples Initialization** (CellGroup): TemplateBox info icon + `Cell[BoxData[RowBox[{"Needs", "[", "\"Wolfram`QuantumFramework`\"", "]"}]], "ExampleInitialization", InitializationCell -> True]`.
13. **PrimaryExamplesSection** (CellGroup): header + `InterpretationBox[...]` "More Examples" link + multiple `ExampleText` → `Input` → `Output` triples separated by `ExampleDelimiter`.
14. **ExtendedExamples** (CellGroup, tagged `"ExtendedExamples"`): the eight subsections — **Scope**, **Generalizations & Extensions**, **Options**, **Applications**, **Properties & Relations**, **Possible Issues**, **Interactive Examples**, **Neat Examples**. Each is a `Cell[..., "ExampleSection"]` wrapping a `CellGroupData[{...ExampleSubsection cells...}]`. Subsections that are empty contain the placeholder text "XXXX".

**Cell-level options always present**:
- `CellID -> <fresh integer>`
- `ExpressionUUID -> "<CreateUUID[]>"`

## 9. Excluded symbols (do NOT document on this branch)

Per the user's directive, these public exports are intentionally undocumented in this branch. If `F` is in this list, stop and tell the user:
- `Kernel/QuantumEvolution.m`: `HamiltonianTransitionRate`, `LindbladTransitionRates`
- `Kernel/ZX.m`: `ZXTensorNetwork`, `ZXTensorNetworkQuantumCircuit`, `ZXExpression`
- `Kernel/Multiway.m`: `QuantumCircuitMultiwayCausalGraph`, `QuantumCircuitPathGraph`, `QuantumCircuitTokenEventGraph`
- `Kernel/QuantumCircuitOperator/TensorNetwork.m`: `TensorNetworkQuantumCircuit`
- `Kernel/Experimental.m`: `DecomposedQuantumStateProbabilities`, `QuantumBeamSearch`, `QuantumDiagramProcess`, `QuantumMPO`, `QuantumMPSApply`, `FeynmanBacktracking`

## 10. Common pitfalls

10.1 Symbols not in the file you'd guess:
- `QuantumEntangledQ`, `QuantumEntanglementMonotone` are in `QuantumEntanglement.m`, not `QuantumEntangledQ.m`.
- `QuantumEvolve` is in `QuantumEvolution.m`.
- `QuantumStateEstimate`, `QuantumStateEstimation`, `QuantumStateSampler`, `QuantumMeasurementSimulation` all live in `QuantumStateEstimation.m`.
- `QuantumWeylTransform` is in `QuantumWignerTransform.m`.
- `QuantumPositiveTransform` is in `QuantumPhaseSpaceTransform.m`.
- `QuantumSimilarity` is in `QuantumDistance.m`.
- `QuantumWignerMICTransform` is in `MIC.m`.
- `QuantumShortcut` is in `Labels.m`.
- `QuantumMPS` is in `Experimental.m`.

10.2 `Wolfram`QuantumFramework`Init`` uses `Memoize[FromOperatorShorthand]` — anything that flows through `FromOperatorShorthand` may have memoized results; clear it with `Wolfram`QuantumFramework`Init`Private`*` if needed during verification.

10.3 The 27 gates are pure delegators to `QuantumOperator`. Their pages can be templated heavily — they all have the same "build a `QuantumOperator` with named matrix; show as a circuit; show its action on common states" structure. But verify each kernel definition actually delegates to the named operator you think it does (e.g. `RXGate[θ, ord]` expects the angle in radians; confirm in the kernel).

10.4 `EinsteinSummation.nb` exists in `$REF` but `EinsteinSummation` is a built-in WL symbol, not a paclet export. If the user asks you to handle `EinsteinSummation`, treat it as an audit-only case: read the page, decide whether it adds paclet-specific value (e.g. demonstrates use within QuantumFramework workflows), and recommend keep / remove / convert-to-tutorial. Do not generate a new ref page for a built-in symbol unless the user explicitly directs.

10.5 Never `Quiet` — never. If a verification run produces messages, capture them. If a kernel has a real bug that produces an unwanted message in the example, file it in the issues report and add a Possible Issues note explaining it.

10.6 Round-trip discipline — for any constructor/accessor pair `C`/`A`, the FIRST property example on the page should be `A[C[x]] === x`. Phase-oblivious / set-wise comparisons come after, opt-in via a named helper. Origin: Phase 5c Y miss in Stabilizer.

10.7 Never use the `Co-Authored-By: Claude` trailer. Author commits with the user's identity only.

