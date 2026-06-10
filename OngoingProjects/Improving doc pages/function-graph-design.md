# Function-relation graph — design for the paclet-doc-rewrite-template

> **Status: SUPERSEDED — preserved as design rationale.**
>
> The design described here was implemented and then promoted into a standalone paclet, **`Wolfram/PacletFunctionGraph`**, which lives in this repository at [`PacletFunctionGraph/`](PacletFunctionGraph/). The standalone scripts and shipped `Resources/` files this design described (`build-function-graph.wl`, `graph-helpers.wl`, `$PACLET/Resources/FunctionGraph.wl`, `$PACLET/Resources/GraphHelpers.wl`) **have been deleted**; their behavior is now provided by:
>
> | Original artifact | Replaced by |
> |---|---|
> | `build-function-graph.wl` (CLI script) | `BuildFunctionGraph[pacletPath, opts]` in the paclet |
> | `graph-helpers.wl` `SymbolNeighborhoodGraph` | `SymbolNeighborhoodGraph[graphData, symName]` in the paclet |
> | `graph-helpers.wl` `PacletOverviewGraph` | `PacletOverviewGraph[graphData]` in the paclet |
> | Shipped `Resources/FunctionGraph.wl` cache | Built on-demand by `BuildFunctionGraph` (a few seconds; no on-disk cache to invalidate) |
> | Schema with `Returns` / `Properties` / `KernelAnchors` | Trimmed to `Kind` / `Accepts` / `OperatesOn` / `ConsumedBy` — the four fields actually used by the renderers |
> | Manual `"PredicateHeads"` table for paclet-specific predicates | Auto-detected by walking each predicate's own DownValues |
>
> **Current entry point:** `Needs["Wolfram`PacletFunctionGraph`"]; graphData = BuildFunctionGraph["/path/to/Paclet"]`. The paclet's [README.md](PacletFunctionGraph/README.md) and the six shipped ref pages (under `PacletFunctionGraph/Documentation/English/`) are the live references.
>
> The text below is preserved as the original design rationale and may use legacy paths/names. Treat it as historical context, not current API.

## Context

The doc-rewrite template at [`Doc writing skill Claude/refpage-template.md`](Doc writing skill Claude/refpage-template.md) defines a **Properties & Relations** subsection under "More Examples". Before this feature, the runbook at [`per-function-runbook.md`](per-function-runbook.md) (§A.5.c, §B.6) asked the writer to fill that subsection by reading the kernel from scratch every time a doc page is rewritten — and the result was uneven. QF doc pages reflected this: `QuantumState.nb` had 15 cells in this section, `QuantumOperator.nb` 3, `QuantumChannel.nb` and `QuantumDistance.nb` were essentially empty stubs.

The feature introduces a **function-interaction graph** that captures, for each public symbol:

- What inputs it accepts, classified by **the actual WL `Head`** the user types — both built-in heads (`List`, `Association`, `String`, `Integer`, `SparseArray`) and paclet heads (`QuantumState`, `QuantumBasis`, `QuantumOperator`). Example: `QuantumState[{1,2}]` — `Head[{1,2}]` is `List`, so the "Accepts" entry for `QuantumState` records `List` (with the descriptor "vector amplitudes"), not the prose "vector".
- What it returns, again classified by `Head` of the result — sometimes a paclet head (`QuantumOperator` returns `QuantumOperator`; `QuantumOperator[qo][QuantumState[…]]` returns `QuantumState`), sometimes a built-in head (`QuantumDistance[…]` returns a `Real`; `QuantumEntangledQ[…]` returns a `Symbol` — `True | False`; `qs["DensityMatrix"]` returns a `List`).
- Which other symbols accept this one as input (reverse edges — e.g. `QuantumDistance`, `QuantumEntangledQ`, `QuantumChannel`, `QuantumMeasurement` all consume `QuantumState`).
- Property dispatch — `qs["DensityMatrix"]` returns a `List`; `qs["NormalizedState"]` returns a `QuantumState`; `qs["VonNeumannEntropy"]` returns a `Real`.

**Boundary — API-level only, not internal representation.** The graph captures the user-facing composition surface (the heads the user *types* and *sees*), not internal storage. We do NOT walk `FullForm` to record that `QuantumBasis[]` internally stores `Association[Rule["Input", QuditBasis[…]], …]`. The graph edge `QuantumOperator[QuantumState[…], …] → QuantumState[…]` is exactly the granularity we want; the underlying `Rule[…SparseArray[…]…]` storage is out of scope.

**Scope clarification.** The deliverable updates the **paclet-agnostic** template/runbook infrastructure at `OngoingProjects/Improving doc pages/` — `per-function-runbook.md` is already paclet-agnostic (see §11 worked example for QF, §10 generic pitfalls), and the new pieces must follow the same pattern. QuantumFramework is the **test bed**: we exercise the infrastructure on QF, refine it, and commit only when it works end-to-end on a real paclet.

The architectural choice — *generate on demand each time, or build once and cache* — has a clear answer: **build once per paclet, cache, invalidate on two triggers.** The graph is the same artifact for every ref page in a paclet (and useful again for the Guide page); recomputing it per page is waste, and hand-grepping is non-deterministic enough to produce inconsistent ref-page content across sessions.

**Rebuild policy** (two triggers, distinct behaviors):

1. **Kernel-SHA drift** — when SHA256 over the inspected kernel files changes, the cache is "stale-by-SHA". Behavior: **warn, do not auto-rebuild.** Rationale: kernel edits during a doc session would silently trigger a multi-minute rebuild, which is the wrong cadence — a writer mid-task does not want the cache spontaneously regenerating mid-session. The writer must pass `--rebuild` explicitly.
2. **Calendar age** — when `$Meta.GeneratedAt` is older than the configured max age (default **30 days**, configurable via `--max-age <days>`), the cache is "stale-by-age". Behavior: **auto-rebuild on next `--check`.** Rationale: a time-based trigger keeps the graph fresh even when nobody touches the kernel for a stretch ("at least monthly" is the design floor). The script logs the auto-rebuild loudly so the writer knows it happened. Opt out with `--no-auto-rebuild` if running in a script that should fail fast on stale caches.

The two triggers are independent: SHA drift alone never auto-rebuilds; age alone always auto-rebuilds (unless opted out); both at once is treated the same as age alone (the rebuild happens, the SHA drift is incidentally resolved).

Pre-existing assets in the runbook that the new pieces should compose with:

- `per-function-runbook.md:32-43` (§0.1) — already resolves `$PACLET`, `$KER`, `$REF`, `$GUIDES`, `<PrimaryContext>`, `<PacletName>`, `$DOCS` for any paclet. The graph lives at `$DOCS/function-graph.wl` and inherits this discovery.
- `per-function-runbook.md:332-339` (§5) — already specifies how a project supplies a kernel-file lookup table; the build script reads this table when present and falls back to grep otherwise (per `per-function-runbook.md:336`).
- `per-function-runbook.md:506-535` (§11) — QF's lookup table, the seed for the QF test run.
- `audit/function-index.md` (737 lines, anchor SHA `cbbe9368`) — QF-only partial prose index. Header says "Examples and pitfalls pending Phase 2 reads." Treat as a sanity check, not a source of truth.

## Recommended approach

Three deliverables, all paclet-agnostic; one test run on QF.

### 1. Paclet-agnostic generator: `build-function-graph.wl`

Path: `OngoingProjects/Improving doc pages/build-function-graph.wl`. Invoked via `wolframscript -file build-function-graph.wl --paclet <PACLET>` (per `feedback_use_wolframscript.md`).

**Parameters (all derive from `PacletInfo.wl` per `per-function-runbook.md:32-43`):**

- `--paclet <PATH>` — paclet root (required).
- `--lookup-table <PATH>` — markdown path to the project-specific symbol→files table (optional; default: grep mode per `per-function-runbook.md:118-120` / §A.1).
- `--exclude <SYMBOL>...` — symbols to skip (optional; mirrors `per-function-runbook.md:498-504`).
- `--rebuild` — force regeneration regardless of cache state.
- `--max-age <days>` — calendar-age rebuild threshold (default `30`).
- `--no-auto-rebuild` — disable the age-based auto-rebuild (default behavior is to auto-rebuild on `--check` when age exceeds `--max-age`).
- `--symbol <NAME>` — regenerate one entry only.

**Output:** `$DOCS/function-graph.wl` where `$DOCS` follows the runbook default (`$PACLET/scratch/docs-audit/`) or user override. The file is an Association keyed by symbol. **All type slots store actual WL `Head` symbols, not strings** — so a downstream consumer can do `MatchQ[expr, _ @@ entry["Returns","Head"]]` style checks. Schema:

```wl
<|
  "$Meta" -> <|
    "Paclet" -> <PacletName>,                     (* e.g. "Wolfram/QuantumFramework" *)
    "PacletVersion" -> <version>,
    "PrimaryContext" -> <PrimaryContext>,
    "KernelSHA" -> <sha256 over the union of kernel files inspected>,
    "GeneratedAt" -> DateObject[…],
    "MaxAgeDays" -> 30,
    "Anchor" -> <git short SHA at generation time>,
    "GenerationIssues" -> {…}                    (* edges the generator could not confidently classify *)
  |>,

  "QuantumState" -> <|
    "Kind" -> "Constructor",                       (* Constructor | Operator | Predicate | Transform | Distance | NamedInstance | Utility *)

    "Accepts" -> {                                 (* One entry per distinct downvalue form, classified by Head of first non-OptionsPattern arg *)
      <|"Head" -> List,          "Description" -> "vector amplitudes"|>,
      <|"Head" -> List,          "Description" -> "density matrix (matrix-shaped List)"|>,
      <|"Head" -> Association,   "Description" -> "named-amplitude association"|>,
      <|"Head" -> SparseArray,   "Description" -> "sparse vector or matrix"|>,
      <|"Head" -> String,        "Description" -> "named state, e.g. \"Bell\""|>,
      <|"Head" -> QuantumState,  "Description" -> "basis-change copy"|>,
      <|"Head" -> QuantumBasis,  "Description" -> "basis only (canonical state)"|>
    },

    "Returns" -> <|"Head" -> QuantumState, "Description" -> "constructed state"|>,

    "OperatesOn" -> None,                          (* QuantumState is not operator-like; for QuantumOperator this map would be non-empty *)

    "ConsumedBy" -> {QuantumOperator, QuantumChannel, QuantumMeasurement, QuantumDistance, QuantumEntangledQ, …},  (* Reverse edges — symbol Heads, populated in a final pass *)

    "Properties" -> <|                              (* Property dispatch: prop string -> Head of returned object *)
      "DensityMatrix"     -> <|"Head" -> List,         "Description" -> "n×n density matrix"|>,
      "NormalizedState"   -> <|"Head" -> QuantumState, "Description" -> "renormalized copy"|>,
      "VonNeumannEntropy" -> <|"Head" -> Real,         "Description" -> "entropy in nats"|>,
      "PureStateQ"        -> <|"Head" -> Symbol,       "Description" -> "True | False"|>,
      …
    |>,

    "KernelAnchors" -> <|
      "Constructor" -> "Kernel/QuantumState/QuantumState.m:28-40",
      "Properties"  -> "Kernel/QuantumState/Properties.m"
    |>
  |>,

  "QuantumOperator" -> <|
    "Kind" -> "Operator",
    "Accepts" -> {…},
    "Returns" -> <|"Head" -> QuantumOperator, "Description" -> "constructed operator"|>,
    "OperatesOn" -> <|                             (* Non-empty: this is what makes QuantumOperator operator-like *)
      QuantumState    -> <|"Head" -> QuantumState,    "Description" -> "applied to a state"|>,
      QuantumOperator -> <|"Head" -> QuantumOperator, "Description" -> "composed with another operator"|>,
      QuantumChannel  -> <|"Head" -> QuantumChannel,  "Description" -> "composed with a channel"|>
    |>,
    "ConsumedBy" -> {QuantumChannel, QuantumMeasurementOperator, QuantumCircuitOperator, QuantumEvolve, …},
    "Properties" -> <|…|>,
    "KernelAnchors" -> <|…|>
  |>,

  "QuantumDistance" -> <|
    "Kind" -> "Distance",
    "Accepts" -> {<|"Head" -> QuantumState, "Description" -> "first state"|>, <|"Head" -> QuantumState, "Description" -> "second state"|>, <|"Head" -> String, "Description" -> "metric name (e.g. \"Fidelity\")"|>},
    "Returns" -> <|"Head" -> Real, "Description" -> "distance value"|>,   (* Built-in head — not a paclet object *)
    "OperatesOn" -> None,
    "ConsumedBy" -> {},
    "Properties" -> <||>
  |>,
  …
|>
```

Notes on the schema:

- `Accepts` is a **list of entries**, one per downvalue form. Same `Head` may appear multiple times with different descriptors (`List` can mean "vector" OR "matrix") — that's intentional; it preserves the distinct calling patterns the writer must demonstrate.
- `Returns` and properties record `Head` of the **user-observed result** — not the internal representation. `qs["DensityMatrix"]` reports `Head -> List`, even though internally that List may be backed by a `SparseArray`.
- Boolean returns are recorded as `Head -> Symbol` (since `Head[True]` is `Symbol`). The descriptor distinguishes (`"True | False"`).
- `ConsumedBy` is a flat list of symbol heads (also actual symbols, not strings) — populated in a final reverse pass.

**Algorithm:**

1. Parse `PacletInfo.wl` to derive `$KER`, `<PrimaryContext>`, `<PacletName>`. Compute the symbol roster from the `"Symbols"` field minus `--exclude`.
2. If `--lookup-table` was supplied, parse the markdown table; otherwise grep `$KER` for each symbol's declarations per `per-function-runbook.md:118-120`.
3. **Cache freshness check** — compute SHA256 over the union of inspected kernel files; read existing cache's `$Meta.GeneratedAt` and `$Meta.KernelSHA`. Decision matrix:
   - Cache missing → regenerate.
   - `--rebuild` passed → regenerate.
   - Age > `--max-age` and `--no-auto-rebuild` not passed → log "STALE-BY-AGE, auto-rebuilding" and regenerate.
   - SHA differs → log "STALE-BY-SHA, run `--rebuild` to refresh" and **exit without writing** (per the explicit-only-on-SHA-drift choice).
   - Otherwise → log "cache fresh (age `<n>` days, SHA match)" and exit.
4. Load the paclet (`PacletDirectoryLoad`, `Needs[<PrimaryContext>]`). For each symbol:
   - **Forms** — programmatic from `DownValues[F]`. For each downvalue LHS, identify the first non-`OptionsPattern` argument and record `Head[<that arg, evaluated as a sample if necessary>]` — the actual WL head the user types. Build one `<|"Head" -> …, "Description" -> …|>` entry per distinct form. **Do not** walk into the FullForm of any constructed object — we record the head of the *argument*, not the head of any internal field.
   - **Properties** — grep the symbol's `Properties.m` (or property-dispatch lines per `per-function-runbook.md:139-141`) for `Switch[prop_, ...]`, `obj_[prop_String]`, `obj_[("a"|"b"|...)]` patterns. For each property, evaluate one sample call (`qs["DensityMatrix"]` etc. against a canonical fixture state) and record `Head[result]`. Sample fixtures: a tiny 2-dim canonical state per type family, defined once at the top of the script.
   - **OperatesOn** — find `F[args___][x_Head] := …` patterns in downvalues. These are the operator-application edges. Record `<input Head> -> <output Head>` where the output head is derived by evaluating the application against a canonical fixture (e.g. `someOp[someState]` → record `Head[result]`).
5. **Reverse-pass:** for each symbol `S`, derive `ConsumedBy` by scanning every other symbol's `Accepts` and `OperatesOn` entries for any `"Head" -> S`.
6. Write the Association via `Put`. Never use `Quiet` (per `feedback_no_quiet.md`); surface every kernel message in the run log and record ambiguities in `$Meta.GenerationIssues`.

**Out of scope (explicit):** the script never inspects `FullForm[expr]` of any paclet object to look at its internal Association/Rule/SparseArray storage. The graph is the public composition surface only. A reader looking at the graph should never see `Rule[…SparseArray[…]…]` — only paclet/built-in heads.

The script also exports two paclet-agnostic helper functions to a side file (`graph-helpers.wl`) so doc-page cells can render without re-invoking the generator:

- `SymbolNeighborhoodGraph[graph_, sym_]` — directed `Graph[]` showing one symbol's neighborhood (Accepts in, OperatesOn-targets and ConsumedBy out). Used by per-page Neat-Examples cells.
- `PacletOverviewGraph[graph_]` — directed `Graph[]` of all symbols, grouped by `"Kind"`. Used by the Guide page.

Edges colored by type (`Accepts` blue, `OperatesOn` green, `ConsumedBy` orange, property-dispatch dashed) so a reader can decode the picture without a legend.

### 2. Paclet-agnostic runbook + template updates

**`OngoingProjects/Improving doc pages/per-function-runbook.md`** — three localized edits, all worded paclet-agnostically:

- **Add §0.6** between §0.5 and §1 (after `:85`): "Ensure `$DOCS/function-graph.wl` is fresh. Run `wolframscript -file build-function-graph.wl --paclet $PACLET --check`. Behavior: cache-missing or older than `--max-age` (default 30 days) → auto-rebuilds and prints `STALE-BY-AGE`; SHA-drift only → prints `STALE-BY-SHA, run --rebuild` and exits without writing; otherwise prints `cache fresh`. This graph is the source for Properties & Relations content (§A.5.c / §B.6) and the Neat-Examples Graph cell (§A.5.h / §B.6.13)."

- **Rewrite §A.5.c** (`:170-171`) and the corresponding §B.6 Properties & Relations bullet (`:293`) so they require sourcing edges from `$DOCS/function-graph.wl` and demonstrate three canonical patterns when applicable:
  - **Round-trip:** `Accessor[Constructor[x]] === x` (per `feedback_user_facing_roundtrip_first.md` and the existing `per-function-runbook.md:473-474` §10.6).
  - **Operator application:** the `OperatesOn` edges.
  - **Consumed-by spotlight:** at least one downstream symbol from the `ConsumedBy` set.

- **Add §A.5.h** (new sub-step after §A.5.g at `:183`) and **§B.6 step 14** (new step after the §8 Neat Examples bullet at `:296`): "Append a Neat-Examples cell that calls `SymbolNeighborhoodGraph[Get[\"$DOCS/function-graph.wl\"], \"F\"]`. The cell evaluates at notebook-open time and shows F's immediate neighborhood."

**`OngoingProjects/Improving doc pages/Doc writing skill Claude/refpage-template.md`** — two localized edits:

- Extend the "Properties and Relations (N)" section (`:180-196`) to list the three canonical edge patterns from §A.5.c.
- Extend the "Neat Examples (N)" section (`:216-218`) with the `SymbolNeighborhoodGraph` one-liner pattern.

**Optional addition (no commitment yet) — `OngoingProjects/Improving doc pages/per-function-runbook.md` Guide-page note**: a short paragraph in §10 or a new §12 noting that paclets with a single canonical Guide page may render `PacletOverviewGraph[graph]` near the top as a one-time addition. Phrased as guidance, not a per-function step.

### 3. QuantumFramework test run (validates the infrastructure)

This is the verification artifact, not part of the template itself, but the work to declare the template change "done":

- Run `build-function-graph.wl --paclet QuantumFramework --lookup-table OngoingProjects/Improving\ doc\ pages/per-function-runbook.md --docs OngoingProjects/Improving\ doc\ pages --exclude <§11 exclusion list> --rebuild`. Produces `OngoingProjects/Improving doc pages/function-graph.wl` (the QF override of the paclet-agnostic default).
- Pick `QuantumChannel.nb` (currently has an empty Properties & Relations stub). Apply the new §A.5.c and §A.5.h to it end-to-end.
- Add a one-cell `PacletOverviewGraph[...]` rendering near the top of `QuantumFramework/Documentation/English/Guides/WolframQuantumComputationFramework.nb`.
- Run the existing verification harness (`per-function-runbook.md:78-80`, §0.4) on `QuantumChannel.nb` — every Input cell must evaluate without messages.

If anything in the QF test run forces a change to the build script or runbook wording, that change is part of the deliverable. The QF artifacts themselves (the generated `function-graph.wl`, the modified `QuantumChannel.nb`, the Guide-page cell) are committed alongside the template change as the worked example, mirroring how `per-function-runbook.md:478-541` (§11) already serves as QF's worked example.

### Critical files to modify

Paclet-agnostic (the template change):

- `OngoingProjects/Improving doc pages/build-function-graph.wl` — new.
- `OngoingProjects/Improving doc pages/graph-helpers.wl` — new; exports `SymbolNeighborhoodGraph`, `PacletOverviewGraph`.
- `OngoingProjects/Improving doc pages/per-function-runbook.md` — add §0.6, edit §A.5.c, edit §B.6, add §A.5.h, add §B.6 step 14.
- `OngoingProjects/Improving doc pages/Doc writing skill Claude/refpage-template.md` — extend `:180-196` and `:216-218`.

QF-specific (the test bed, committed as worked example):

- `OngoingProjects/Improving doc pages/function-graph.wl` — generated. (QF's `$DOCS` is set to this directory per `per-function-runbook.md:528` §11 worked example, overriding the paclet-agnostic default of `$PACLET/scratch/docs-audit/`.)
- `QuantumFramework/Documentation/English/ReferencePages/Symbols/QuantumChannel.nb` — end-to-end test page.
- `QuantumFramework/Documentation/English/Guides/WolframQuantumComputationFramework.nb` — add overview Graph cell.

### Files to read (reuse, do not duplicate)

- `OngoingProjects/Improving doc pages/per-function-runbook.md:32-43`, `:118-120`, `:332-339`, `:498-541` — paclet-config discovery, grep fallback, lookup-table contract, QF worked example.
- `audit/function-index.md` — QF-only partial prose index; sanity check for the QF test run, not a source of truth.
- `QuantumFramework/Kernel/QuantumState/QuantumState.m`, `Kernel/QuantumOperator/QuantumOperator.m`, etc. — the kernel sources the generator reads during the QF test run.

## Verification

End-to-end checks before declaring this done. All paclet-agnostic except where marked.

1. **Generator self-check** — `wolframscript -file build-function-graph.wl --paclet $PACLET --rebuild` produces `$DOCS/function-graph.wl` whose `$Meta.KernelSHA` matches a fresh `shasum` over the union of inspected kernel files. Re-running without `--rebuild` prints "cache fresh."

2. **Coverage (QF)** — every symbol in QF's `PacletInfo.wl` `"Symbols"` list (minus the §11 exclusion list at `per-function-runbook.md:498-504`) has an entry. Spot-check 3 symbols (`QuantumState`, `QuantumOperator`, `QuantumChannel`) against `audit/function-index.md` for agreement on Forms; mismatches go into `$DOCS/issues-report.md` per `per-function-runbook.md:388-413` §7.

3. **Reverse-edge integrity** — for each symbol `S`, every `T` listed in `S["ConsumedBy"]` must have `S` appearing in `T["Accepts"]` or `T["OperatesOn"]`. Round-trip property; assert via `VerificationTest`.

4. **End-to-end doc-page test (QF)** — apply the new §A.5.c (code triples) and §A.5.h (Neat-Examples Graph cell) to `QuantumChannel.nb`. Run the verification harness at `$DOCS/verify-doc-page.wls`. Every Input cell evaluates without messages, and the rendered `Graph[]` visibly shows `QuantumChannel`'s neighborhood (`QuantumOperator` and `QuantumMeasurementOperator` as `Accepts` sources; `QuantumState` and `QuantumChannel` as `OperatesOn` targets).

5. **Guide-page overview (QF)** — open `Documentation/English/Guides/WolframQuantumComputationFramework.nb` after adding the overview cell. The `PacletOverviewGraph` rendering shows ≥ 25 nodes (one per public symbol) grouped by `"Kind"`, with no orphan nodes.

6. **Re-running on a changed kernel (SHA-drift trigger)** — touch one kernel file (e.g. add a no-op comment), re-run `--check`. Reports `STALE-BY-SHA` and does NOT auto-rebuild. Passing `--rebuild` regenerates.

6a. **Re-running after the age threshold (age trigger)** — force `$Meta.GeneratedAt` to a date > 30 days ago (e.g. via direct edit of the cache file), re-run `--check`. Reports `STALE-BY-AGE, auto-rebuilding` and writes a fresh graph. Passing `--no-auto-rebuild` makes the same run exit without writing.

6b. **API-level boundary** — grep the generated `function-graph.wl` for `SparseArray`, `Rule[`, and any internal-storage tokens. None should appear in the Accepts/Returns/OperatesOn/Properties slots. (`SparseArray` may legitimately appear as an `Accepts` `Head` for a symbol that takes sparse input — that's fine. What we're checking against is *content* of an Accepts entry containing `Rule[…SparseArray[…]…]` nesting, which would indicate the generator walked into FullForm.)

7. **Paclet-agnosticism smoke test** — run `build-function-graph.wl --paclet <other paclet>` on a second paclet (e.g. another paclet under `~/Documents/GitHub/`) with grep fallback (no lookup table). It should at minimum produce a non-empty `$Meta` and a partial graph without erroring. Genuine completeness for non-QF paclets is out of scope, but the script must not crash.

8. **No `Quiet`** — grep the generator and helper scripts for `Quiet[` after writing; must return zero matches (per `feedback_no_quiet.md`).
