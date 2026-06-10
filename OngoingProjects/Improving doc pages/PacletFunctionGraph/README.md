# PacletFunctionGraph

A paclet-agnostic tool that introspects any Wolfram Language paclet and produces
a **directed function-relation graph** showing how the paclet's head symbols
flow into one another.

The graph model is simple and self-contained:

- **Nodes** — only paclet head symbols (`QuantumState`, `QuantumOperator`, …).
  Built-in heads (`Real`, `Integer`, `List`, `String`, `Association`, …) are
  *never* nodes.
- **Edges** — directed, unlabeled. `A → B` means "data of head A flows into B"
  (either A is one of B's accepted inputs, or A produces something B consumes).
- **Sinks** — symbols that return non-paclet scalars (predicates, distances)
  appear as dead-end nodes; their `Kind` color signals what they are.

## Public API

| Symbol | Purpose |
|---|---|
| `BuildFunctionGraph[pacletPath, opts]` | Introspect the paclet, return the graph `Association`. Writes `Resources/FunctionGraph.wl` unless `"WriteFile" -> False`. |
| `LoadFunctionGraph[file]` | Load a previously-built graph. |
| `PacletFlowEdges[graph]` | List of `{fromName, toName}` string pairs — the edge primitive. |
| `SymbolNeighborhoodGraph[graph, symName]` | `Graph[]` showing all edges that touch `symName`. |
| `PacletOverviewGraph[graph]` | `Graph[]` of the whole paclet. Option `"IncludeIsolated" -> True` keeps disconnected symbols. |
| `$FunctionGraphKindColors` | `<|Kind -> RGBColor|>` palette used by the renderers. |

## Quick start

```wolfram
Needs["Wolfram`PacletFunctionGraph`"]

(* Build (writes <paclet>/Resources/FunctionGraph.wl) *)
g = BuildFunctionGraph["/path/to/MyPaclet"]

(* Or load a cached one *)
g = LoadFunctionGraph["/path/to/MyPaclet/Resources/FunctionGraph.wl"]

(* Render *)
PacletOverviewGraph[g]
SymbolNeighborhoodGraph[g, "MyConstructor"]
```

## Build options

```
BuildFunctionGraph[pacletPath,
    "OutputFile"     -> Automatic,        (* default: <paclet>/Resources/FunctionGraph.wl *)
    "WriteFile"      -> True,             (* False = return the Association without writing *)
    "Exclude"        -> {},               (* symbol names to skip *)
    "OnlySymbol"     -> None,             (* limit to one symbol *)
    "PredicateHeads" -> <||>,             (* extra predicate -> heads mappings, e.g. <|"QuantumStateQ" -> {"QuantumState"}|> *)
    "KindClassifier" -> Automatic         (* custom (sym :> "Kind") classifier; defaults to a name-based heuristic *)
]
```

### Extending the predicate table

The default `PredicateHeads` table covers built-ins (`NumericQ`, `MatrixQ`,
`AssociationQ`, …). When your paclet uses custom predicates such as
`QuantumStateQ`, pass them in:

```wolfram
g = BuildFunctionGraph["/path/to/QuantumFramework",
    "PredicateHeads" -> <|
        "QuantumStateQ"    -> {"QuantumState"},
        "QuantumOperatorQ" -> {"QuantumOperator"},
        "stateQ"           -> {"List", "SparseArray"}
    |>
]
```

### Custom Kind classification

If your paclet doesn't follow the default name conventions
(`*Distance` → `"Distance"`, `*Q` → `"Predicate"`, etc.), supply your own:

```wolfram
g = BuildFunctionGraph[..., "KindClassifier" -> (Function[s,
    Which[
        StringStartsQ[SymbolName[s], "Solve"], "Operator",
        StringEndsQ[SymbolName[s], "Plot"],    "Transform",
        True,                                  "Constructor"
    ]
])]
```

## Graph schema

`BuildFunctionGraph` returns an Association of the shape:

```
<|
    "$Meta" -> <|
        "Paclet"           -> "Publisher/Name",
        "PacletVersion"    -> "x.y.z",
        "PrimaryContext"   -> "Publisher`Name`",
        "GeneratedAt"      -> DateObject[...],
        "GenerationIssues" -> {...}
    |>,
    "Symbol1" -> <|
        "Kind"       -> "Constructor" | "Operator" | "Transform" | "Predicate" | "Distance" | ...,
        "Accepts"    -> {<|"Head" -> Sym, "Role" -> "Primary"|"Conversion"|"Mutation"|"NamedInstance"|"NamedArg"|"PatternMachinery", "Description" -> ...|>, ...},
        "OperatesOn" -> None  |  <|Head1 -> <|"Head" -> Head1, "Description" -> ...|>, ...|>,
        "ConsumedBy" -> {symbol references that consume this one}
    |>,
    "Symbol2" -> <|...|>,
    ...
|>
```

## Naming rationale

`PacletFunctionGraph` was chosen because:

- "Paclet" scopes the tool's target — distinguishing it from generic
  graph utilities or call-graph analyzers.
- "FunctionGraph" matches the shipped data file name (`FunctionGraph.wl`)
  and self-documents what the tool produces.

Other names considered: `PacletFlow` (catchier but less precise),
`PacletAtlas` (evocative but obscure).

## Status

`0.1.0` — initial release. Extracted from the QuantumFramework
documentation-rewrite project; tested only on QuantumFramework so far.
