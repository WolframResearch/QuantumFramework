# Template: Paclet Reference Page Structure

This is the canonical structure for a Wolfram paclet symbol reference page.
Every section must be present in the output. Leave the body empty if there is nothing to say.

This template is for **paclet symbols** — functions shipped in a `.paclet`, documented in `.nb` notebooks under `Documentation/English/ReferencePages/Symbols/`. It is NOT for Wolfram Function Repository submissions (use `wfr-template.md` for those).

## Formatting rules

- **Voice**: terse, Wolfram-documentation style. One-liner descriptions. No filler.
- **Code blocks**: `` ```wl `` fenced blocks. Show Input/Output pairs where appropriate.
- **Code captions**: plain text ending with a colon. One line only.
- **Tables**: all alignment headers must be explicit: `|:---|` or `|:---:|` or `|---:|`.
- **`---`** (horizontal rule) only use at the beginning of a new example.
- Multiple outputs in one block: only the last expression omits the semicolon.
- Never use `Head[...]` in examples — show the expression directly.
- Each `---` horizontal separator starts a fresh example with new variables. Plan splits so variables don't leak across examples.

### Paclet link format

Links must include the full paclet qualifier for non-system symbols. The skill resolves `PublisherID` and `PacletName` during the gather-context step.

- **System symbols**: `paclet:ref/Symbol`, `paclet:guide/GuideName`, `paclet:tutorial/TutorialName`
- **Paclet symbols**: `paclet:PublisherID/PacletName/ref/Symbol`
- **Paclet guides**: `paclet:PublisherID/PacletName/guide/GuideName`
- **Paclet tutorials**: `paclet:PublisherID/PacletName/tutorial/TechNoteName`

---

## Page structure

What follows is the page structure. Inline comments (in parentheses) are guidance — do not include them in the output.

## PublisherID\`PacletName\`

(Context line. Use the paclet's top-level context.)

# SymbolName

(One-line description: what the function does or represents.)

## Details and Options

(Bullet list covering behavior, argument semantics, return types, edge cases. Include everything a user needs before looking at examples.)

* When called with [form1], ``SymbolName`` does X.

* When called with [form2], ``SymbolName`` does Y.

* The function returns [description of return type / structure].

(If the return is structured — e.g. an Association with specific keys — use a table:)

|            |                                 |
| :--------- | :------------------------------ |
| "Key1"     | description of what Key1 holds  |
| "Key2"     | description of what Key2 holds  |

* [Detail about input types, constraints, or defaults.]

* [Detail about related functions for simpler/richer interfaces.]

* ``SymbolName`` returns a [Failure](paclet:ref/Failure) object if called with incorrect arguments.

(Options, if any:)

* The following options can be given:

| option   | default | description      |
| :------- | :------ | :--------------- |
| "Name"   | value   | what it controls |

## Examples (N)

(N = total number of example subsections that have content.)

```wl
In[1]:= Needs["PublisherID`PacletName`"]
```

### Basic Examples (N)

(2–3 simple cases giving immediate understanding. Each has a prose caption ending with a colon.)

Caption ending with colon:

```wl
In[1]:= SymbolName[args]

Out[1]= result
```

---

Another caption:

```wl
In[1]:= SymbolName[differentArgs]

Out[1]= result
```

### Scope (N)

(Systematic enumeration of all calling patterns, input types, and features. Organize by logical grouping — one `####` subsection per topic.)

#### Input Forms

(Show each valid first-argument type with examples.)

Caption:

```wl
In[1]:= SymbolName[form1, args]

Out[1]= result
```

---

Caption:

```wl
In[1]:= SymbolName[form2, args]

Out[1]= result
```

#### Distinctive Features

(IMPORTANT: If the package has capabilities not available in the standard library — tracing, step-by-step execution, bytecode inspection, specialized diagnostics — these MUST be documented here with full working examples. These features justify the package's existence and are what users cannot get elsewhere.)

Caption:

```wl
In[1]:= Code demonstrating distinctive feature

Out[1]= result
```

#### Workflow Context

(If the function participates in a pipeline with other package functions — e.g. compile -> create VM -> execute -> inspect — show at least one example demonstrating the function with its upstream or downstream neighbors.)

Caption:

```wl
In[1]:= upstreamResult = UpstreamFunction[args];
SymbolName[upstreamResult, expr]

Out[1]= result
```

### Options (N)

(One `####` per option with examples. Leave section present but empty body if no options.)

#### OptionName

Caption:

```wl
In[1]:= SymbolName[args, "OptionName" -> value]

Out[1]= result
```

### Applications (N)

(How the function is used in practice. Each application has a caption and tells a mini-story where appropriate.)

Caption:

```wl
In[1]:= Code

Out[1]= result
```

### Properties and Relations (N)

(How this function relates to others — both within the paclet and to system functions. Each entry is one fact demonstrated with code.)

(Source: the per-paclet function-relation graph at `$DOCS/function-graph.wl`. The graph's entry for this symbol has Accepts / Returns / OperatesOn / ConsumedBy / Properties — pick 3–6 of the highest-pedagogical-value edges and demonstrate each with `ExampleText` → `Input` → `Output`. See `per-function-runbook.md` §0.6 and §A.5.c for the build/check protocol.)

(Three canonical patterns must be included when applicable:)

(1) **Round-trip**: for any constructor/accessor pair `C`/`A`, the FIRST entry must be `A[C[x]] === x`. Origin: prior Phase 5c Y miss. See `per-function-runbook.md` §10.6.

(2) **Operator application**: when the graph's `OperatesOn` is non-empty, demonstrate at least one `F[args][x]` edge — e.g. `qo[QuantumState[...]]` produces a `QuantumState`.

(3) **Consumed-by spotlight**: when `ConsumedBy` is non-empty, demonstrate at least one downstream symbol that takes this one as input (e.g. `QuantumDistance[QuantumState[...], ...]` consuming a state).

(Additional patterns to look for:)
(- Equivalence to a system function under specific conditions)
(- Relationship between different calling patterns)
(- How return values feed into other functions)
(- Mathematical or logical properties)

Caption:

```wl
In[1]:= Code showing relationship

Out[1]= result
```

### Possible Issues (N)

(Gotchas that surprise users. Each entry shows what goes wrong and why. Be generous — include all reasonable gotchas and edge cases. Easier to trim later than to discover missing entries.)

(Patterns to look for:)
(- Argument order differences from similar system functions)
(- Input validation — what happens with wrong types or arity)
(- Edge cases — empty input, boundary values, degenerate cases)
(- Performance — when to prefer one calling pattern over another)

Caption:

```wl
In[1]:= Code that shows the issue

Out[1]= unexpected or error result
```

### Neat Examples (N)

(Two kinds of cells live here. The Graph cell is mandatory whenever the function-relation graph is meaningful for the symbol; visually-impressive demos are optional.)

(1) **Function-relation graph cell** — renders this symbol's immediate neighborhood by building the function-relation graph on demand via the `Wolfram/PacletFunctionGraph` paclet. There is no shipped `Resources/FunctionGraph.wl` to keep in sync; the build runs in a few seconds and the doc-page builder pre-evaluates this cell so the user opening the notebook sees the graph immediately:

```wl
(* hardcoded dev path until Wolfram/PacletFunctionGraph is published *)
PacletDirectoryLoad[
    "<repo>/OngoingProjects/Improving doc pages/PacletFunctionGraph"];
Needs["PublisherID`PacletName`"];
Needs["Wolfram`PacletFunctionGraph`"];

pacletRoot = ParentDirectory[NotebookDirectory[], 3];
graphData  = BuildFunctionGraph[pacletRoot, "WriteFile" -> False];
SymbolNeighborhoodGraph[graphData, "SymbolName"]
```

Returns a directed `Graph[]`. Edges point in the direction of data flow (`A → B` means "A flows into B" — A is one of B's accepted inputs, or A is what B produces on application). Vertex face color encodes `Kind` (Constructor / Operator / Transform / Predicate / Distance per `$FunctionGraphKindColors`); the focal symbol carries a heavier black border. No caption beyond a one-line ExampleText is needed.

(2) **Visually impressive demo** — optional. Only when genuinely impressive. Leave the cell empty otherwise.

## See Also

(Related symbols — both paclet symbols and system symbols. Paclet symbols use full qualified links, system symbols use short form.)

* [PacletSymbol](paclet:PublisherID/PacletName/ref/PacletSymbol)
* [SystemSymbol](paclet:ref/SystemSymbol)

## Related Guides

* [PacletGuide](paclet:PublisherID/PacletName/guide/GuideName)

## Related Tech Notes

(Optional — only if relevant tech notes exist.)

* [TechNote](paclet:PublisherID/PacletName/tutorial/TechNoteName)
