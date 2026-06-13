# CLAUDE.md: QuantumFramework working folder

## STRICT RULE: never read `Usage.m` for what a function does

`QuantumFramework/Kernel/Usage.m` is **generated automatically from the doc pages** (it is serialized `BoxData`/`Cell[...]` usage boxes, not human-authored source). It is not a source of truth for behavior, signatures, options, or return shapes, and it must never be cited as one.

- **Never** rely on `Usage.m` to learn or assert what a symbol does, its arguments, options, or output form.
- **Never** quote `Usage.m` as evidence in an audit, review, or answer.
- When usage/signature/behavior information is needed, **read the kernel definitions** for that symbol (the actual `*.m` implementation under `QuantumFramework/Kernel/`), or check the doc page directly, or test it with `wolframscript`.
- Treat any edit to `Usage.m` as a generated artifact: it is regenerated from doc pages, so hand edits there are pointless and will be overwritten.

This rule overrides any default instinct to grep `Usage.m` for a quick signature.
