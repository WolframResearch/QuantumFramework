# qf-claude-skill

A self-contained [Claude Code](https://docs.claude.com/en/docs/claude-code) plugin that bundles
the **`qf` skill**: a deep working guide for the Wolfram `QuantumFramework` paclet (top pitfalls,
construction patterns, performance idioms, a full function index, and the object-model
architecture). When you write or review Wolfram Language code that touches a `QuantumFramework`
symbol, Claude loads the skill and its bundled reference set so it codes against the real paclet
behavior instead of guessing.

This folder is both a plugin and a one-plugin marketplace, and it is the **single source** for
the skill (no duplicate copies to keep in sync).

## Layout

```
qf-claude-skill/
├── .claude-plugin/
│   ├── plugin.json        # the plugin manifest
│   └── marketplace.json   # advertises this folder as a marketplace hosting the plugin
└── skills/
    └── qf/                # the skill itself
        ├── SKILL.md
        ├── README.md
        └── reference/{function-index,mistakes,idioms-and-performance,architecture,qf-tn-integration,VERSION}.md
```

The reference files are read on demand, so the skill costs almost nothing until it fires.

## Install

**As a plugin** (versioned; from a local clone of the repo):

```bash
claude plugin marketplace add ./qf-claude-skill
claude plugin install qf@qf-skill
```

(The marketplace name is `qf-skill`, from `marketplace.json`; the plugin name is `qf`.) The
`marketplace add` form takes a path to the folder that holds `.claude-plugin/marketplace.json`,
so it works from a checkout. To enable a one-line `claude plugin marketplace add owner/repo`
straight from GitHub, the `marketplace.json` would need to sit at the repo root; this folder
keeps it self-contained instead.

**As a plain skill** (no plugin machinery): copy the skill folder into a skills directory.

```bash
cp -r qf-claude-skill/skills/qf ~/.claude/skills/qf     # personal, all projects
# or, per-project for collaborators:
cp -r qf-claude-skill/skills/qf <some-repo>/.claude/skills/qf
```

After either install the skill auto-loads on any `QuantumFramework` symbol, or invoke it with
`/qf`.

## Prerequisites

The skill is a knowledge layer. To actually run anything you also need **Wolfram Language /
Mathematica** (the `wolframscript` CLI) and the **`Wolfram/QuantumFramework` paclet** installed.
The cites were verified against paclet 2.0.0; see
[skills/qf/reference/VERSION.md](skills/qf/reference/VERSION.md) for the exact commit and how to
re-check a cite on your install (line numbers drift between releases; behavior notes are stable).

## Keeping it current

The reference set is a snapshot. When the paclet moves, regenerate the cites by re-reading the
touched `Kernel/*.m` files and bump `skills/qf/reference/VERSION.md` to the new commit.
