# Template: Guide Page Structure

This is the canonical structure for a Wolfram paclet Guide page.
A guide page is an index that lists all symbols in a paclet (or a functional area), grouped by topic. It has no code examples — just organized symbol listings with short descriptions.

## Formatting rules

- **Voice**: terse, Wolfram-documentation style. Descriptions are sentence fragments, not full sentences.
- **No code blocks** — guide pages list symbols, they don't demonstrate them.
- **Symbol entries** use the pattern: `[SymbolName](paclet:ref/SymbolName) -- short description`
- **Bullet lists** for secondary symbols (less prominent, no description): `* [Symbol](paclet:ref/Symbol)`
- **`---`** (horizontal rule) separates topic groups visually.

### Paclet link format

Links in documentation must include the full paclet qualifier for non-system symbols. The format depends on whether you're linking to system documentation or paclet documentation:

- **System symbols**: `paclet:ref/Symbol`, `paclet:guide/GuideName`, `paclet:tutorial/TutorialName`
- **Paclet symbols**: `paclet:PublisherID/PacletName/ref/Symbol`
- **Paclet guides**: `paclet:PublisherID/PacletName/guide/GuideName`
- **Paclet tutorials**: `paclet:PublisherID/PacletName/tutorial/TechNoteName`

The skill resolves `PublisherID` and `PacletName` during the gather-context step. Use `paclet:ref/...` (short form) only for built-in Wolfram Language symbols.

## Content patterns

There are two levels of symbol prominence in each topic group:

1. **Featured symbols** — shown with a description line. These are the 2–3 most important symbols in the group. Format: `[Symbol](paclet:ref/Symbol) -- what it does`
2. **Secondary symbols** — bullet list, no description. Related symbols that round out the group. Format: `* [Symbol](paclet:ref/Symbol)`

For groups with many symbols, a mix of both is ideal. For small groups (2–3 symbols), all can be featured.

## Guide levels

There are two common guide patterns:

### Top-level guide (multi-area index)

Covers a broad domain. Each topic group links to a child guide via its `###` heading. Groups are short — 1–2 featured symbols and 3–4 bullets, with a `[...](paclet:guide/ChildGuide)` pointer for more. The abstract is substantial (1–2 paragraphs) because it frames the whole area.

Example (from GraphsAndNetworks):

```markdown
### [Construction and Representation](paclet:guide/GraphConstructionAndRepresentation)

* [Graph](paclet:ref/Graph) -- graph object with vertex and edge properties

* [GraphData](paclet:ref/GraphData), [ExampleData](paclet:ref/ExampleData) -- curated collection of theoretical and empirical graphs

[CompleteGraph](paclet:ref/CompleteGraph) . [RandomGraph](paclet:ref/RandomGraph) . [AdjacencyGraph](paclet:ref/AdjacencyGraph) . [Import](paclet:ref/Import)
```

Key traits: group heading is a link to the child guide, `[...]` entries point readers deeper, groups stay concise.

### Focused guide (single-area)

Covers one functional area in depth. Groups are more granular — one per sub-topic, with more featured symbols and longer bullet lists. Some groups may have only featured symbols (no bullets) and some only bullets (no featured).

Example (from GraphConstructionAndRepresentation):

```markdown
* [Graph](paclet:ref/Graph) -- represent a general graph, or create it from vertices and edges

* [UndirectedEdge](paclet:ref/UndirectedEdge) -- an undirected edge (also entered as <->)

* [DirectedEdge](paclet:ref/DirectedEdge) -- a directed edge (also entered as ->)

### Basic Properties

* [VertexList](paclet:ref/VertexList), [EdgeList](paclet:ref/EdgeList) -- the list of vertices and edges in the graph

[AdjacencyList](paclet:ref/AdjacencyList) . [IncidenceList](paclet:ref/IncidenceList) . [EdgeRules](paclet:ref/EdgeRules)

### Graphs from Data

* [RelationGraph](paclet:ref/RelationGraph) -- generate a graph based on data and a binary relation

* [NearestNeighborGraph](paclet:ref/NearestNeighborGraph) -- generate the k-nearest neighbor graph

* [NestGraph](paclet:ref/NestGraph) -- generate a nested function graph

[CayleyGraph](paclet:ref/CayleyGraph) . [MoleculeGraph](paclet:ref/MoleculeGraph)
```

Key traits: top-level featured symbols before the first `###` (the "core" of the area), groups are detailed with more featured entries, child guide links in headings are optional (only when there's a further sub-guide).

---

## Page structure

What follows is the page structure. Inline comments (in parentheses) are guidance — do not include them in the output.

# Guide Title

(1–2 paragraph abstract. Sets context: what area this covers, why it matters, how the symbols relate to each other. Write in Wolfram documentation voice — informative, not promotional.)

---

### Topic Group Name

(Each group covers one functional area. The `###` heading names the topic. If the paclet has sub-guides for each topic, the heading can link to them.)

* [MainSymbol](paclet:ref/MainSymbol) -- what it does

* [SecondSymbol](paclet:ref/SecondSymbol), [ThirdSymbol](paclet:ref/ThirdSymbol) -- what these related symbols do together

[AdditionalSymbol1](paclet:ref/AdditionalSymbol1) . [AdditionalSymbol2](paclet:ref/AdditionalSymbol2) . [AdditionalSymbol3](paclet:ref/AdditionalSymbol3)

---

### Another Topic Group

(Repeat the pattern. Use `---` between groups.)

* [Symbol](paclet:ref/Symbol) -- description

[Symbol2](paclet:ref/Symbol2) . [Symbol3](paclet:ref/Symbol3)

## Related Guides

(Other guide pages in the paclet or related system guides.)

* [OtherGuide](paclet:guide/OtherGuide)

## Related Tech Notes

(Tutorials / tech notes that expand on topics covered by this guide.)

* [TechNoteName](paclet:tutorial/TechNoteName)

## Related Links

(External links — Wolfram U courses, blog posts, papers.)

* [Link text](url)
