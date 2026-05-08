---
name: documentation-writing
description: Write Wolfram Language documentation notebooks — WFR function pages, paclet Guide pages, or Tech Note tutorials. Use whenever asked to document a function, write a reference page, create a guide page, write a tutorial, explain a concept in documentation form, create a tech note, or produce any .nb documentation notebook. Also trigger when documenting functions for the Wolfram Function Repository, writing paclet documentation, or creating narrative explanations of algorithms or workflows as notebooks.
---

# Documentation Writing

Write Wolfram Language documentation pages as `.nb` notebooks. Covers four page types:

- **WFR** — Wolfram Function Repository submissions (function definition + docs + tests + metadata)
- **Ref Page** — paclet symbol reference page (for paclets)
- **Guide** — index page listing symbols grouped by topic (for paclets)
- **Tech Note** — narrative tutorial explaining concepts with prose, math, and code (for paclets)

Use this skill whenever asked to document a function, write a reference page, create a guide page, or write a tech note / tutorial. One page per invocation.

## Required reading

1. **The appropriate template** in `references/`:
   - `references/wfr-template.md` — WFR function pages
   - `references/refpage-template.md` — Paclet symbol reference pages
   - `references/guide-template.md` — Guide pages
   - `references/technote-template.md` — Tech Note / tutorial pages
2. **The source material** — function definition, paclet source code, or concept to document
3. **An existing page for grounding** — read a shipped documentation page via `mcp__Wolfram__ReadNotebook` to see conventions in practice. System documentation is at `/Library/Wolfram/Documentation/<version>/en-us/Documentation/English/System/` under `ReferencePages/Symbols/`, `Guides/`, or `Tutorials/`. Paclet documentation is under `.../Paclets/<PacletName>/`.

### Recommended exemplar pages

Use these as grounding references (read via `mcp__Wolfram__ReadNotebook`). All paths below are for macOS (`/Library/Wolfram/Documentation/...`). On Windows the root is typically `C:\Program Files\Wolfram Research\Mathematica\<version>\Documentation\English\`. Adjust accordingly. Paths are relative to `.../Documentation/15.0/en-us/Documentation/English/`:

**WFR / Reference pages:**
- `System/ReferencePages/Symbols/Sin.nb` — good model for mathematical functions. Its structure (usage patterns, details, examples progression) is the template shared by all elementary and special function pages.

**Guide pages:**
- `System/Guides/GraphsAndNetworks.nb` — excellent example of a multi-level guide. It indexes many sub-areas, each linking to child guides. For a focused, single-topic guide, read one of its children instead (e.g., `System/Guides/GraphConstructionAndRepresentation.nb` or `System/Guides/GraphProperties.nb`).

**Tech Notes / Tutorials:**
- `Paclets/FEMDocumentation/Tutorials/FiniteElementOverview.nb` — masterclass in deep-dive tech notes: progressive structure, definition tables, interleaved math and code. The individual FEM tutorials (e.g., `SolvingPDEwithFEM.nb` — 23 sections, textbook-length) represent the extreme upper bound; they are essentially a textbook shipped as documentation. For most paclets, a deep-dive tech note targets 6–10 sections (~4,000–6,000 words) — beyond that, split into multiple focused tech notes. Use the FEM tutorials for structural patterns (definition tables, progressive complexity, math/code interleaving), not as a length target. Read selectively (first few sections) rather than loading the entire notebook.

---

## Process

### 1. Scope the page

**For Ref Page, WFR, and Guide pages**, skip to step 2 — scope is determined by the function or paclet.

**For Tech Notes**, start with a discovery conversation before doing any research. Ask the user:

1. **Audience** — who is this for? (beginner getting started, intermediate user wanting depth, advanced user needing internals)
2. **Goal** — what should the reader be able to do after reading this? (use the API, understand the algorithm, extend the framework)
3. **Depth** — pick a level:
   - **Overview** — what it is, when to use it, basic examples. 1–2 sections. Quick to write, quick to read.
   - **Practical** — how it works, all the options, common patterns, gotchas. 4–6 sections. The default.
   - **Deep dive** — theory, derivations, implementation details, advanced usage. 8+ sections. FEM-tutorial level. Requires significant back-and-forth.
4. **Theory** — how much math? (none / light notation / full derivations)
5. **Scope boundaries** — what's explicitly out of scope? (helps avoid scope creep)

These answers determine the section structure, code-to-prose ratio, and whether definition tables and mathematical derivations are needed. A "deep dive" tech note may need multiple writing passes and rounds of feedback — set that expectation upfront.

### 2. Gather context

**For all paclet page types (Ref Page, Guide, Tech Note)**, start by resolving the paclet identity. Evaluate via the MCP evaluator:

```wolfram
{$PublisherID, PacletObject["PacletName"]["Name"]}
```

or read the `PacletInfo.wl` in the paclet directory. This gives you the **PublisherID** and **PacletName** needed for all paclet-qualified links (`paclet:PublisherID/PacletName/ref/Symbol`) and the context line (`PublisherID`PacletName``). Use `paclet:ref/...` (short form) only for built-in Wolfram Language symbols.

**For Ref Page (paclet symbols):**
- All calling patterns — enumerate every valid argument types (raw pattern, lists, associations, pre-built object, etc.)
- Return type and structure (Association keys, wrapper heads, failure cases)
- Options and defaults (if any)
- Any `::usage` string
- Workflow neighbors — what functions produce inputs for this one, and what consumes its output
- Distinctive features — capabilities not available in the standard library (tracing, step execution, diagnostics). These justify the package and MUST be documented in Scope.
- Error/failure behavior — what happens with wrong argument count, wrong types

**For WFR functions:**
- All calling patterns (different arities, option variants)
- Option definitions and defaults
- Return types and edge cases
- Dependencies on other functions or packages
- Any `::usage` string

**For Guide pages:**
- Full list of public symbols in the paclet
- How symbols group by topic / functionality
- Which symbols are most important (featured) vs. secondary (bullet lists)
- Related guides and tech notes in the paclet

**For Tech Notes:**
- The concept or workflow to explain
- Which symbols are involved and how they relate
- Mathematical background if applicable (calibrated to the depth level from step 1)
- Good examples that build understanding progressively
- Related guides and other tech notes

### 3. Plan the page

Before writing any markdown, produce a short outline appropriate to the page type:

**Ref Page:**
- One-line description
- Usage lines — one per calling pattern. EVERY public calling pattern must appear as a separate line. This is the most visible part of the page and cannot be abbreviated.
- Details and Options bullet points
- Example plan: 2–3 basic examples, scope subsections (input forms, distinctive features, workflow context), application ideas
- See Also — both paclet symbols and system symbols. Use `WolframLanguageContext` via MCP evaluator for system symbols.
- Related Guides and Tech Notes

**WFR:**
- One-line description
- Usage rows (one per calling pattern)
- Which Details & Options points apply
- Example plan: 2–3 basic examples, scope subsections, application ideas
- Related Symbols — only documented system symbols, 5–10
- Keywords and Categories

**Guide:**
- Title and abstract (1–2 paragraphs)
- Topic groups with featured symbols (2–3 per group with descriptions) and secondary symbols (bullet list)
- Related Guides and Related Links

**Tech Note — Overview:** (~500–1,000 words of prose, 2–4 code blocks)
- Title and one-paragraph intro
- 1–2 `##` sections: what it is + basic usage
- Related Guides / Tech Notes

**Tech Note — Practical:** (~2,000–4,000 words of prose, 8–15 code blocks)
- Title and introductory paragraph
- 4–6 `##` sections with key concepts per section
- Definition tables where they help (one per stage/concept transition)
- Code examples plan per section
- Related Guides and Related Tech Notes

**Tech Note — Deep dive:** (~4,000–6,000 words of prose, 15–25 code blocks)
- Title and motivation (why this matters, what you'll learn)
- Section outline: `##` for major topics, `###` for subtopics
- Definition tables at each stage boundary (like the FEM tutorials: functions → descriptions before the code)
- Mathematical derivations plan: which formulas, how to build from simple to general
- Code examples plan: one per concept, progressing in complexity
- References section plan
- Related Guides and Related Tech Notes

> **Context limits.** A single LLM pass can reliably produce one page at a time — one reference page, one guide, or one tech note up to the Practical level. Deep-dive tech notes above ~6,000 words risk quality degradation; split them into 2–3 focused tech notes instead, or write them in multiple sessions (one section group per pass, then stitch together).
>
> **Documenting a whole paclet.** If the user asks you to document an entire package (guide + multiple reference pages + tech notes), plan the full page set first, then use subagents to write individual pages in parallel. Each subagent gets one page, the skill instructions, and the paclet source. Merge and cross-link afterward.

Wait for confirmation before proceeding.

### 4. Write the markdown

Follow the appropriate template exactly — it contains all sections and formatting rules. The template is the single source of truth for structure and style.

**Code block conventions:**
- Use `` ```wl `` fenced blocks.
- Show Input/Output pairs where appropriate: `In[n]:= ...` / `Out[n]= ...`
- Multiple outputs in one block: only the last expression omits the semicolon.
- Never use `Head[...]` in examples — show the expression directly.
- Prefer long variable names: `sampledData`, `fittedResult`, not abbreviations.
- Each `---` horizontal separator starts a fresh example with new variables. Plan splits so variables don't leak across examples.

**Plot conventions:**
- Add `Exclusions -> None` to `Plot`/`Plot3D` when plotting functions with potential discontinuities.
- `PlotLegends` inside the plot function, never in `Show`.

### 5. Validate examples

Run every code block through the MCP evaluator (`mcp__Wolfram__WolframLanguageEvaluator`). Do this as a batch after writing the full markdown — not interleaved with writing.

For WFR: evaluate the function definition first so it's available in the session.
For paclet pages: load the paclet first via `PacletDirectoryLoad` + `Get`.

For each code block:
- Run it. If it succeeds without errors, move on.
- If it errors, **do not fix it in the doc**. Flag it and continue.

Collect all issues (errors, unexpected results, warnings) into a list.

**Report validation results** to the user: total blocks run, blocks passed, blocks failed with error summary.

### 6. Write the notebook

Use the MCP tool directly:
```
mcp__Wolfram__WriteNotebook[
    file -> "<output-path>/<PageName>.nb",
    markdown -> <the full markdown string>
]
```

### 8. Report issues (if any)

If validation found issues or design tensions, present them to the user as a summary:
- Code that failed and the error message
- Any design tensions discovered during documentation (e.g. undocumented features, API gaps, missing edge case handling)

For non-trivial issues, create a separate `<PageName>_Issues.nb` notebook via `WriteNotebook` so issues are tracked alongside the doc page.

---

## Checklist — Ref Page

- [ ] Read paclet source code and `::usage` strings before writing anything
- [ ] Identify all calling patterns (different first-argument types)
- [ ] Identify distinctive features not in standard library
- [ ] Identify workflow neighbors (upstream producers, downstream consumers)
- [ ] Plan outline and get confirmation
- [ ] One-line description is clear and concise
- [ ] **Every** public calling pattern has its own usage line in the header
- [ ] Details and Options is thorough — return type, argument semantics, edge cases
- [ ] Every section from template is present (empty if nothing to say)
- [ ] Scope includes distinctive features with full working examples
- [ ] Scope includes at least one workflow-context example
- [ ] Applications section present (even if brief)
- [ ] Possible Issues is generous — all reasonable gotchas
- [ ] See Also has both paclet and system symbols
- [ ] Related Guides and Related Tech Notes listed
- [ ] All formatting follows `references/refpage-template.md`
- [ ] Every code block validated via MCP evaluator (batch, after writing)
- [ ] Validation results reported (total/passed/failed)
- [ ] Notebook written via `WriteNotebook` MCP tool

## Checklist — WFR

- [ ] Read and understand the function definition before writing anything
- [ ] Plan outline and get confirmation
- [ ] One-line description is clear and concise
- [ ] Usage section covers all calling patterns
- [ ] Details & Options is thorough
- [ ] Every section from template is present (empty if nothing to say)
- [ ] All formatting follows `references/wfr-template.md`
- [ ] Keywords are comprehensive
- [ ] Categories are selected
- [ ] Related Symbols has 5–10 documented system symbols
- [ ] Tests section has verification tests
- [ ] Compatibility section is filled in
- [ ] Every code block validated via MCP evaluator (batch, after writing)
- [ ] Validation results reported (total/passed/failed)
- [ ] Notebook written via `WriteNotebook` MCP tool

## Checklist — Guide

- [ ] Full inventory of public symbols before writing
- [ ] Plan topic groups and get confirmation
- [ ] Title and abstract are clear and informative
- [ ] Symbols grouped logically by topic
- [ ] Each group has 2–3 featured symbols with `[Symbol](paclet:ref/Symbol) -- description`
- [ ] Secondary symbols in bullet lists
- [ ] Related Guides listed
- [ ] Keywords present
- [ ] All formatting follows `references/guide-template.md`
- [ ] Notebook written via `WriteNotebook` MCP tool

## Checklist — Tech Note

- [ ] Scope the page: audience, goal, depth level, math level, boundaries
- [ ] Gather context based on depth level
- [ ] Plan section outline and get confirmation
- [ ] Title and intro paragraph set context
- [ ] Definition table if applicable
- [ ] Sections build understanding progressively
- [ ] Code examples are self-contained and runnable
- [ ] Math notation is clear (prose or code blocks)
- [ ] References cited where applicable
- [ ] Related Guides and Related Tech Notes listed
- [ ] Keywords present
- [ ] All formatting follows `references/technote-template.md`
- [ ] Every code block validated via MCP evaluator (batch, after writing)
- [ ] Validation results reported (total/passed/failed)
- [ ] Notebook written via `WriteNotebook` MCP tool
