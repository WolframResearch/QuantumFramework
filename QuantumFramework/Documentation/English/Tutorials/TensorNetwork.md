---
Template: TechNote
Name: TensorNetwork
Title: Tensor Network
Context: Wolfram`QuantumFramework`
CellContext: Global`
Paclet: Wolfram/QuantumFramework
URI: Wolfram/QuantumFramework/tutorial/TensorNetwork
RelatedGuides: [WolframQuantumComputationFramework]
---

<!-- #| style: DefinitionBox -->
|   |   |
|---|---|
| [QuantumCircuitOperator]()[…]["TensorNetwork"] | returns the tensor network of a circuit |
| [TensorNetworkIndexGraph]()[…] | transforms a tensor network into a new graph with indices as vertices |
| [TensorNetworkFreeIndices]()[…] | returns free indices in a tensor network index graph |
| [ContractTensorNetwork]()[…] | contracts indices in a tensor network |

<!-- #| style: MathCaption -->
Create a quantum circuit:

```wl
circuit = 
  QuantumCircuitOperator[{"S", "H" -> 2, "X" -> 3, "CNOT", 
    "SWAP" -> {2, 3}, {1}, {2}, {3}}];
circuit["Diagram"]
```

In the Wolfram Quantum Framework, measurement outcomes are recorded in an ancillary quantum register (serving as a detector), so the complete circuit diagram includes both the system and this detector subsystem.

```wl
circuit["Diagram", "ShowExtraQudits" -> True]
```

The first measurement result is stored in wire index 0, with subsequent results assigned to decreasing (negative) wire indices.

Returns the tensor network representation corresponding to the given quantum circuit:

```wl
net = circuit["TensorNetwork", 
  GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}]
```

The tensor network is a graph annotated with tensors and contraction indices.

```wl
GraphQ[net] && TensorNetworkQ[net]
```

For a graph to be a tensor network, two main conditions must be met, which is what [TensorNetworkQ]() does. First, every vertex in the graph must have properly formatted indices consisting only of Subscript or Superscript expressions, and each vertex cannot contain any duplicate indices - all indices at any single vertex must be unique from each other.

Second, the tensor rank compatibility condition must be satisfied. This can happen in two ways: either each vertex's tensor rank must exactly match the number of indices assigned to that vertex, and when counting all index occurrences throughout the entire graph, every distinct index must appear exactly twice across all vertices (ensuring proper pairing for tensor contraction), or alternatively, all tensors in the network must have rank zero representing scalar tensors with no indices.

Lists all annotation keys available for the tensor network:

```wl
AnnotationKeys[{net, 0}]
```

Let’s look at tensor, vertex labels, and corresponding indexes in the tensor network:

```wl
list = Developer`FromPackedArray@VertexList[net];
TableForm[
 Transpose[
  Prepend[AnnotationValue[{net, list}, #] & /@ {"Tensor", "Index"}, 
   Sort@AnnotationValue[net, VertexLabels]]], TableDepth -> 2, 
 TableHeadings -> {None, (Style[#1, 
       Bold] &) /@ {"Vertex index ->\n Vertex label", "Tensor", 
     "Leg indices"}}]
```

Tensors in the tensor network are mixed type, meaning they consist of so-called “contravariant” (upper) indices and “covariant” (lower) indices. For example, the 2nd measurement (7th vertex), acts on qubit-2 (denoted by contravariant and covariant indices 2) and its result is saved on wire, denoted by the index “-1”.

Vertices correspond to circuit’s operators/gate indices, in addition to “Initial” tensor with index 0 (for the initial state):

```wl
VertexList[net]
```

```wl
Length@circuit["Flatten"]["Operators"]
```

```wl
% == Max[%%]
```

Note that each edge represents (i.e., is tagged by) a contraction:

```wl
EdgeList[net]
```

Corresponding contraction:

```wl
EdgeTags[net]
```

Perform the contraction:

```wl
finalTensor = ContractTensorNetwork[net]
```

Confirm that the result is the same as default circuit application:

```wl
circuit[]["Tensor"] == finalTensor
```

Another tensor network representation uses indices as graph vertices with tensors as cliques:

```wl
indexNet = 
 TensorNetworkIndexGraph[net, 
  GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}]
```

In the above graph, the directed edges imply tensor contraction; also tensors are cliques in above graph:

```wl
HighlightGraph[indexNet, 
 Subgraph[indexNet, #] & /@ FindClique[indexNet, Infinity, All]]
```

Free indices are the ones that left after contraction:

```wl
TensorNetworkFreeIndices[net]
```

Highlight free indices:

```wl
HighlightGraph[indexNet, TensorNetworkFreeIndices[net]]
```

Free indices can be extracted as vertices with zero in- and out- degree:

```wl
Pick[VertexList[indexNet], 
 VertexInDegree[indexNet] + VertexOutDegree[indexNet], 0]
```

```wl
ContainsExactly[TensorNetworkFreeIndices[net], 
 TensorNetworkFreeIndices[net]]
```

### Contraction and Einstein Summation

Many useful information of a tensor network can be extracted using TensorNetworkData:

```wl
TensorNetworkData[net] // Dataset
```

Also, one can get the useful data using graph functionalities, too:

```wl
indices = 
 AnnotationValue[{net, list}, "Index"] /. Rule @@@ EdgeTags[net]
out = TensorNetworkFreeIndices[net]
tensors = AnnotationValue[{net, list}, "Tensor"];
```

Show ContractTensorNetwork is the same as EinsteinSummation:

```wl
ContractTensorNetwork[net] === 
 ActivateTensor@EinsteinSummation[indices -> out, tensors]
```

One can compare the performance, on how the relevant computation is done

Perform the contraction in the order of network’s EdgeList:

```wl
ContractTensorNetwork[net, Method -> "Naive"]; // RepeatedTiming
```

Optimize the order for contraction, using EinsteinSummation and symbolic tensors package:

```wl
ContractTensorNetwork[net]; // RepeatedTiming
```

### Initial state different from ground state

Note that the initial tensor in the tensor network we studied here was a registered state. Additionally, one can start from any initial state.

Generate a random state:

```wl
\[Psi]0 = QuantumState[{"RandomPure", 3}];
```

Initialize the tensor network from above state:

```wl
net2 = QuantumCircuitOperator[{\[Psi]0, circuit}]["TensorNetwork"];
```

See supplement info, for package-scoped symbols.

Show the tensor contraction is the same as transformation of state by the circuit:

```wl
ContractTensorNetwork[net2] == circuit[\[Psi]0]["Tensor"]
```

## Supplement info

### FromTensorNetwork

Any directed graph can be turned into a tensor network, even if the graph is not annotated.

```wl
SeedRandom[1];
graph = DirectedGraph[RandomGraph[{10, 10}], "Acyclic", 
  VertexLabels -> v_ :> HoldForm[v]]
```

Convert the graph into a tensor network such that tensors for each vertex is randomly generated in the dimension $\{2,2,\ldots \}$ where the overall length of dimensions is determined by the tensor rank (overall legs/edges):

```wl
qcTN = GraphTensorNetwork@graph
```

Get the tensor network data:

```wl
TensorNetworkData@qcTN
```

One can assign symbolic tensor into vertices too:

```wl
qcTNSymbolic = GraphTensorNetwork[graph, Method -> "Symbolic"];
TensorNetworkData@qcTNSymbolic
```

### Package-scoped symbols

In the Wolfram quantum framework, there are some package-scoped symbols that live in the `` Wolfram`QuantumFramework`PackageScope` `` context.

```wl
?Wolfram`QuantumFramework`PackageScope`$*
```
