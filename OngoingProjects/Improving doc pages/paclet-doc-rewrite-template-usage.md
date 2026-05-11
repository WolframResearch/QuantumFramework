# Using `paclet-doc-rewrite-template` for QuantumFramework doc rewrites

> Project-local cheat-sheet. For the full methodology, read the canonical template directly: [`~/.claude/prompts/paclet-doc-rewrite-template.md`](file:///Users/mohammadb/.claude/prompts/paclet-doc-rewrite-template.md).

## What this directory is

`OngoingProjects/Improving doc pages/` tracks ongoing work to clean up QuantumFramework's reference doc pages. The workflow is driven by the canonical template at `~/.claude/prompts/paclet-doc-rewrite-template.md` (paclet-agnostic, ~620 lines, 5 API regimes + bootstrap mode).

## Invoking for a QF page (the QF shortcut)

In a fresh Claude Code session, paste:

```
Use the template at ~/.claude/prompts/paclet-doc-rewrite-template.md.
Paclet config: use the QuantumFramework defaults from Appendix C.
<SYMBOL>       = QuantumOperator         (the symbol to document)
<SIBLING_PAGE> = QuantumState.nb         (or a closer cousin once it exists)
<PARENT>       = QuantumBasis            (look up in Appendix C parent-class map)
<SYMBOL-area>  = QuantumState            (kernel subdirectory)
Begin with Checkpoint 0.
```

The template reads Appendix C for QuantumFramework paclet-config and runs the 4-checkpoint workflow.

## What's done · what's next

| Symbol | Status | Notes |
|---|---|---|
| `QuantumState` | Partial — surgical edits applied (Options ExampleSection removed, XXXX placeholders cleared), spec for full rewrite locked in conversation history. Full mechanical rewrite handed off to a spawned task. | Backup at `QuantumState.nb.bak.20260510-182845` |
| `QuantumOperator`, `QuantumChannel`, `QuantumMeasurementOperator`, `QuantumCircuitOperator`, `QuantumEvolve`, others | Not started | Use the template; once `QuantumState.nb` is done, it becomes the `<SIBLING_PAGE>` for the rest |

## QF-specific quick references

| Item | Value |
|---|---|
| Paclet root | `/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework` |
| Docs dir | `<PACLET_ROOT>/Documentation/English/ReferencePages/Symbols` |
| Tests | `wolframscript -f Tests/RunTests.wls Stabilizer` |
| Regression baseline | `main@HEAD`; one pre-existing failure (`Phase2-PauliStabilizer-UsageMessage`) is unrelated to doc work |
| Empty-placeholder convention | `"XXXX"` (verified) |
| Most QF symbols | Regime (a) property-dispatch — except `QuantumEvolve` (regime b) and `QuantumCircuitOperator` (regime c) |

## When the canonical template needs updating

If you discover a new edge case during a rewrite (e.g. a sixth regime, a new kernel-source convention, a category-design rule that doesn't generalize), edit the canonical template at `~/.claude/prompts/paclet-doc-rewrite-template.md`. The QF shim at `~/.claude/prompts/qf-doc-rewrite-template.md` is a thin redirector — don't put new content there.

## Related files

- [Canonical template](file:///Users/mohammadb/.claude/prompts/paclet-doc-rewrite-template.md) — 620 lines, paclet-agnostic.
- [QF shim](file:///Users/mohammadb/.claude/prompts/qf-doc-rewrite-template.md) — 39-line redirect to the canonical, for muscle-memory continuity.
- [Pre-generalization backup](file:///Users/mohammadb/.claude/prompts/qf-doc-rewrite-template.md.bak.20260511-110943) — historical reference.
- `Doc writing skill Claude/` (sibling subdir) — older notes from the doc-improvement project, predates the canonical template.
- `per-function-runbook.md` (sibling file) — older per-function runbook, predates the canonical template.
