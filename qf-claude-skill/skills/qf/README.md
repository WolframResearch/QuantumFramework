# `qf` skill: a deep working guide for the Wolfram QuantumFramework paclet

This is a Claude Code [Agent Skill](https://docs.claude.com/en/docs/claude-code/skills). When
you write, debug, or review Wolfram Language code that touches a `QuantumFramework` symbol,
Claude loads this skill and the bundled reference set so it codes against the *real* paclet
behavior instead of guessing.

## What's inside

```
qf/
├── SKILL.md                 # the skill: loading, top pitfalls, construction patterns, idioms
└── reference/
    ├── function-index.md        # every public function with file:line
    ├── mistakes.md              # ~90 documented pitfalls with citations
    ├── idioms-and-performance.md# fast paths, anti-patterns, notebook-confirmed workflows
    ├── architecture.md          # object map, property dispatch, dimension semantics
    ├── qf-tn-integration.md     # the TensorNetworks bridge + default circuit-apply path
    └── VERSION.md               # paclet version + commit the cites were verified against
```

The reference files are read on demand, so the skill costs almost nothing until it fires.

## Prerequisites

The skill is a *knowledge* layer; to actually run anything you need:

- **Wolfram Language / Mathematica** (the `wolframscript` CLI for verification), and
- the **`Wolfram/QuantumFramework` paclet** installed (`PacletInstall["Wolfram/QuantumFramework"]`)
  or a local checkout loaded via `PacletDirectoryLoad`.

The notes were verified against paclet **2.0.0**; see [reference/VERSION.md](reference/VERSION.md)
for the exact commit and how to re-check a cite on your own install (line numbers drift between
releases; behavior notes are stable).

## Install

**Per-user (all your projects):**

```bash
cp -r qf ~/.claude/skills/qf
```

**Per-project (everyone who clones the repo gets it):** commit this folder at
`.claude/skills/qf/` in the repo. Invoke as `/qf` or let Claude auto-load it on a QF symbol.

**As a plugin (versioned, shareable across teams):** this folder ships inside the
`qf-claude-skill/` plugin (one level up), which packages it for `claude plugin marketplace add`
/ `claude plugin install`.

## Companion skills (optional)

If you also have them installed, Claude will lean on `wl-style-guide` / `wl-native-style` for
idiomatic Wolfram Language and `nb-reader` for `.nb` doc pages. They are recommended, not
required: this skill stands on its own.

## Keeping it current

The reference set is a snapshot. When the paclet moves, regenerate the cites by re-reading the
touched `Kernel/*.m` files and bumping `reference/VERSION.md` to the new commit.
