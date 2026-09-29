---
Template: TechNote
Name: QuantumObjectComposition
Title: Quantum Object Composition
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/QuantumObjectComposition
Keywords: [composition, function relations, type taxonomy, dataflow, function graph]
RelatedGuides: [WolframQuantumComputationFramework]
RelatedTutorials: [GettingStarted, Quantumobjectabstraction]
---

The public objects of [QuantumFramework]() — states, operators, channels, circuits, measurements, transforms — compose with one another via typed transformations: a $QuantumOperator$ acting on a $QuantumState$ returns a $QuantumState$; a $QuantumDistance$ of two $QuantumState$s returns a real number. This tech note introduces a paclet-agnostic tool, Wolfram/PacletFunctionGraph, that introspects any paclet and produces a directed function-relation graph capturing exactly these compositions. After working through it, you will be able to read the type signature of every public symbol, distinguish type-preserving from type-changing morphisms, and render or query the composition graph for any neighborhood of interest.

## Loading the framework and building the graph

Load [QuantumFramework]() and the Wolfram`PacletFunctionGraph` helper paclet. The latter is not yet published, so we point PacletDirectoryLoad at the dev copy in this repository. The Remove call clears any Global\` stubs interned by the front-end before the paclet was loaded, which would otherwise shadow the package's exports:

```wl
#| eval: false
Needs["Wolfram`QuantumFramework`"];

PacletDirectoryLoad[
  "~/QuantumFramework/OngoingProjects/Improving doc \
pages/PacletFunctionGraph"];
Needs["Wolfram`PacletFunctionGraph`"];
```

Build the graph for QuantumFramework. BuildFunctionGraph loads the target paclet, walks its DownValues and SubValues to classify how every public symbol is called, and returns an Association keyed by symbol name. Predicate guards like _?QuantumStateQ or _?propQ are resolved automatically by walking each predicate's own DownValues — no manual configuration needed for any paclet:

```wl
pacletRoot = ParentDirectory[NotebookDirectory[], 3];
graphData = BuildFunctionGraph[pacletRoot, "WriteFile" -> False];
```

The result is an Association. A top-level `$Meta` entry records when it was generated, the paclet identity, and any issues from generation:

```wl
graphData["$Meta"]
```

Symbol entries (one per public head) populate the rest of the Association:

```wl
Keys[KeyDrop[graphData, "$Meta"]]
```

## The type taxonomy

Every public head is classified into one of five Kinds. The Kind controls how the head fits into the composition picture:

| Kind | Meaning | Examples |
|---|---|---|
| Constructor | Builds an object of its own type from raw data. | QuantumState, QuantumBasis, QuantumMeasurement |
| Operator | Object that also acts as a function on a target object, typically returning the same type. | QuantumOperator, QuantumChannel, QuantumCircuitOperator |
| Transform | Function that takes one object and returns the same or a related type. | QuantumPartialTrace, QuantumTensorProduct, QuantumWignerTransform |
| Predicate | Function returning a Boolean (Symbol True or False). | QuantumEntangledQ |
| Distance | Function returning a scalar (Real) from two states. | QuantumDistance, QuantumSimilarity |

Inspect the Kind distribution across the paclet:

```wl
KeyValueMap[{#1, #2["Kind"]} &, 
  KeyDrop[graphData, "$Meta"]] // TableForm
```

The Kind tag drives vertex coloring in the rendered graphs (next section) and is the first filter when querying for symbols of a given role. A "Predicate" head is, for example, never a source of any outgoing flow edge — it returns a scalar, not a paclet object.

## Reading an entry

Each symbol's entry records four fields: Kind, Accepts, OperatesOn, and ConsumedBy. Inspect QuantumState as an example:

```wl
entry = graphData["QuantumState"];
Keys[entry]
```

**Accepts** is a list of input forms the head's constructor recognizes, classified by the WL $Head$ of the user-typed argument and tagged with a Role: "Primary" (the canonical input), "Mutation" (constructor takes its own head to re-wrap), "Conversion" (takes another paclet head and reroutes), "NamedInstance" (a String like "Bell" selecting a catalog entry), "NamedArg" (a Rule like "Eigenvalues" -> v), or "PatternMachinery" (internal). Tally the heads to see the dominant input forms:

```wl
Tally[#["Head"] & /@ entry["Accepts"]]
```

Filter to just the entries whose Role is Primary or Conversion (the ones that produce graph edges):

```wl
Select[entry["Accepts"], 
 MemberQ[{"Primary", "Mutation", "Conversion"}, #["Role"]] &]
```

**OperatesOn** is non-$None$ only for heads whose application has been recorded as F[args][target] in the kernel. $QuantumState$ is a Constructor, not an Operator, so this is $None$:

```wl
entry["OperatesOn"]
```

**ConsumedBy** is the reverse-edge list: every other paclet symbol that takes $QuantumState$ as input. This answers “what can I pass a state into?”:

```wl
entry["ConsumedBy"]
```

## Worked compositions

The graph captures two structurally distinct kinds of composition. We walk through one example of each, then show a property-dispatch case that the graph deliberately does not record.

### Type-preserving: operator application

An $Operator$-kind head applied to a target of type $T$ typically returns a result of the same type $T$. In the flow graph this is a self-loop on $T$. Build a state, apply Pauli-X, and check the head:

```wl
state = QuantumState["0"];
operator = QuantumOperator["X"];
result = operator[state];
Head[result]
```

Same head as the target. That preservation is precisely the OperatesOn field of an Operator-kind head, recording “what target type, what output type” for each application form:

```wl
graphData["QuantumOperator"]["OperatesOn"]
```

### Type-changing: distance returns a Real

A $Distance$-kind head consumes objects and produces a scalar. The composition is type-changing — from $QuantumState$s out to $Real$. Because the graph drops scalars (Real / Integer / Boolean) from the node set, a Distance head appears as a vertex with only incoming edges (a sink):

```wl
d = QuantumDistance[QuantumState["0"], QuantumState["Plus"], 
   "Fidelity"];
Head[N[d]]
```

QuantumDistance's Accepts list confirms it consumes QuantumStates and a metric-name String:

```wl
graphData["QuantumDistance"]["Accepts"]
```

### Property dispatch (not captured by the graph)

Every paclet object responds to property-string queries: qs["DensityMatrix"] returns a matrix, qs["NormalizedState"] returns a fresh $QuantumState$, and so on. These are runtime accessors, not constructor calls, so they don't appear as edges in the function-graph by design — the graph is about *object-to-object* composition, not property lookup. Still, the writer of a ref page documents them; here is a quick survey on a Bell state:

```wl
qs = QuantumState["PhiPlus"];
{Head[qs["DensityMatrix"]], Head[qs["NormalizedState"]], 
 Head[qs["VonNeumannEntropy"]]}
```

## Rendering the graph

Three rendering helpers ship with Wolfram`PacletFunctionGraph`: PacletFlowEdges (the edge primitive), SymbolNeighborhoodGraph (one symbol's flow neighborhood), and PacletOverviewGraph (the whole paclet). All three accept the graph Association and return either an edge list or a $Graph$.

The flow neighborhood of QuantumState: every incoming arrow is a head that flows into QuantumState (as constructor input or as the output of an operator that produces a state), and every outgoing arrow is a head that consumes a QuantumState. The bold black border marks the focal vertex; node color encodes Kind.

```wl
SymbolNeighborhoodGraph[graphData, "QuantumState"]
```

The paclet-wide flow view. Operator-like heads (green) cluster in the upper rows with self-loops indicating they preserve their target type; Transforms (yellow) sit one layer below; Distance / Predicate / Similarity (purple / red) are pure sinks with only incoming arrows:

```wl
PacletOverviewGraph[graphData]
```

If you only want the raw edge data, PacletFlowEdges returns a deduplicated list of {from, to} pairs. This is useful when feeding the graph into other tooling (e.g. computing centrality or partitioning):

```wl
edges = PacletFlowEdges[graphData];
{Length[edges], Take[edges, 5]}
```

## Querying the graph

The graph is data — queryable with standard Association/list patterns.

### What produces a QuantumState?

Find the symbols with an outgoing edge into $QuantumState$:

```wl
producers = 
 Cases[PacletFlowEdges[graphData], {f_, "QuantumState"} :> f]
```

### What consumes a QuantumState?

The mirror query: outgoing edges from QuantumState. Equivalent to (but easier to read than) the symmetric look-up in ConsumedBy:

```wl
consumers = 
 Cases[PacletFlowEdges[graphData], {"QuantumState", t_} :> t]
```

### Which heads are closed under their own action?

An operator-like head is closed if its OperatesOn maps some paclet head $K$ back to itself — in the flow graph, a self-loop or a parallel pair $K$ -> $F$ -> $K$.

```wl
selfLoops[g_] := Cases[PacletFlowEdges[g], {x_, x_}];
selfLoops[graphData]
```

### Which heads are sinks?

A sink has incoming edges but no outgoing ones — symbols whose output is a scalar ($Real$, $Boolean$, $Association$) rather than a paclet object. $QuantumDistance$ and $QuantumEntangledQ$ are typical sinks:

```wl
sinks[g_] := Module[{e = PacletFlowEdges[g]},
     Complement[e[[All, 2]], e[[All, 1]]] // Union
   ];
sinks[graphData]
```

## Where this leaves us

Three takeaways:

- • Every public QuantumFramework head has a typed signature recorded in the function-graph. The composition rules are not buried in kernel files — one call to BuildFunctionGraph produces a queryable Association.

- • Operator-like heads are characterized by self-loops in the flow graph (their OperatesOn maps their target type back to itself). $QuantumOperator$, $QuantumChannel$, and $QuantumCircuitOperator$ all preserve $QuantumState$ under application.

- • Type-changing heads ($QuantumDistance$, $QuantumEntangledQ$) return scalars and therefore appear as sink vertices: incoming edges only.

The same workflow applies to any paclet: Needs["Wolfram`PacletFunctionGraph`"], BuildFunctionGraph["/path/to/SomePaclet"], then render or query.
