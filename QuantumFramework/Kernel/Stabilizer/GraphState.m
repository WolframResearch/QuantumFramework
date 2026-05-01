Package["Wolfram`QuantumFramework`"]

PackageExport[GraphState]
PackageExport[LocalComplement]
PackageScope[GraphStateQ]



(* ============================================================================ *)
(* Phase 5 \[Dash] Graph-state representation (Anders & Briegel 2005).          *)
(*                                                                              *)
(* Every stabilizer state is local-Clifford-equivalent to a graph state         *)
(* (AndBri05 \[Section]2; cite van den Nest, Dehaene, De Moor 2004).            *)
(* The graph-state representation has memory cost O(N \[CenterDot]\[OverBar]d)  *)
(* with \[OverBar]d = avg degree, beating the O(N^2) tableau for sparse codes.  *)
(*                                                                              *)
(* Internal representation:                                                     *)
(*   GraphState[<|"Graph" -> g_Graph, "VOPs" -> vops:{Integer..}|>]             *)
(*                                                                              *)
(* The VOPs are integer indices into the 24-element single-qubit Clifford       *)
(* group; index 0 = identity, others enumerated per AndBri05 \[Section]2        *)
(* footnote (the "And05" reference). For Phase 5 v1, only index 0 (identity)    *)
(* is supported; richer VOPs deferred (synthesis \[Section]2.2).                *)
(* ============================================================================ *)


(* Predicate *)
GraphStateQ[GraphState[KeyValuePattern[{
    "Graph" -> _Graph,
    "VOPs" -> _List
}]]] := True

GraphStateQ[_] := False


(* Property contract *)
_GraphState["Properties"] = {
    "Graph", "VOPs", "Vertices", "Edges", "VertexCount", "EdgeCount",
    "AdjacencyMatrix", "Stabilizers", "Qubits", "PauliStabilizer"
}


(* ============================================================================ *)
(* Constructors                                                                 *)
(* ============================================================================ *)

(* From a Graph: VOPs default to all-identity (0) *)
GraphState[g_Graph] := GraphState[<|
    "Graph" -> g,
    "VOPs" -> ConstantArray[0, VertexCount[g]]
|>]


(* From a PauliStabilizer: build the canonical graph state for a "graph-form"
   stabilizer state. For Phase 5 v1, only handles the case where the input is
   already in graph-state form (stabilizers of the form X_i Z_{neighbors of i}). *)
GraphState[ps_PauliStabilizer ? PauliStabilizerQ] := Module[{n, stabs, g},
    n = ps["Qubits"];
    stabs = ps["Stabilizers"];
    (* Try to extract graph from stabilizers of the form X_i ... Z_{neighbors of i} ... I ...  *)
    Module[{edges},
        edges = Catenate @ Table[
            With[{stabStr = stabs[[i]]},
                Module[{chars = Characters @ If[StringStartsQ[stabStr, "-"], StringDrop[stabStr, 1], stabStr]},
                    (* Position i should have "X"; other positions with "Z" are neighbors *)
                    If[chars[[i]] === "X",
                        Position[chars, "Z"][[All, 1]] /. j_Integer :> If[i < j, UndirectedEdge[i, j], Nothing],
                        {}
                    ]
                ]
            ],
            {i, n}
        ];
        g = Graph[Range[n], DeleteDuplicates @ edges];
        GraphState[<|"Graph" -> g, "VOPs" -> ConstantArray[0, n]|>]
    ]
]


(* ============================================================================ *)
(* Direct property handlers                                                     *)
(* ============================================================================ *)

GraphState[assoc_Association][prop_String] /; KeyExistsQ[assoc, prop] := assoc[prop]

gs_GraphState["Vertices"] := VertexList @ gs["Graph"]
gs_GraphState["Edges"] := EdgeList @ gs["Graph"]
gs_GraphState["VertexCount"] := VertexCount @ gs["Graph"]
gs_GraphState["EdgeCount"] := EdgeCount @ gs["Graph"]
gs_GraphState["Qubits"] := VertexCount @ gs["Graph"]
gs_GraphState["AdjacencyMatrix"] := Normal @ AdjacencyMatrix @ gs["Graph"]


(* ============================================================================ *)
(* Stabilizer extraction                                                        *)
(*                                                                              *)
(* For a graph state with all-identity VOPs, the stabilizer at vertex i is      *)
(*   K_i = X_i \[CircleTimes] Product_{j in neighbors(i)} Z_j                   *)
(* (AndBri05 \[Section]2 Eq 1).                                                 *)
(* ============================================================================ *)

gs_GraphState["Stabilizers"] := Module[{n = gs["VertexCount"], adj = gs["AdjacencyMatrix"]},
    Table[
        StringJoin @ Table[
            Which[
                i == j, "X",
                adj[[i, j]] == 1, "Z",
                True, "I"
            ],
            {j, n}
        ],
        {i, n}
    ]
]


(* Convert to PauliStabilizer (for n <= some threshold, e.g. 30) *)
gs_GraphState["PauliStabilizer"] := PauliStabilizer @ gs["Stabilizers"]


(* ============================================================================ *)
(* LocalComplement[g, v]: AndBri05 Definition 1                                 *)
(*                                                                              *)
(* Complement all edges among the neighbors of vertex v (i.e. for each pair     *)
(* (a, b) of distinct neighbors of v, toggle the edge between a and b).         *)
(*                                                                              *)
(* Theorem (AndBri05 Thm 1): the resulting graph state differs from |G> by      *)
(* a known local unitary U \[Proportional] sqrt(K_G^(v)), captured by VOP       *)
(* updates per Eq (8) of AndBri05.                                              *)
(*                                                                              *)
(* Phase 5 v1: returns the new Graph (with edges toggled) but does not update   *)
(* VOPs. For VOP-tracked LC, defer to a richer Phase 5 follow-up.               *)
(* ============================================================================ *)

LocalComplement[g_Graph, v_] := Module[{neighbors, neighborPairs},
    neighbors = AdjacencyList[g, v];
    neighborPairs = Subsets[neighbors, {2}];
    Fold[
        If[EdgeQ[#1, UndirectedEdge @@ #2],
            EdgeDelete[#1, UndirectedEdge @@ #2],
            EdgeAdd[#1, UndirectedEdge @@ #2]
        ] &,
        g,
        neighborPairs
    ]
]


LocalComplement[gs_GraphState, v_] :=
    GraphState[<|
        "Graph" -> LocalComplement[gs["Graph"], v],
        "VOPs" -> gs["VOPs"]    (* TODO Phase 5+: track VOP updates per AndBri05 Eq 8 *)
    |>]


(* ============================================================================ *)
(* Formatting                                                                   *)
(* ============================================================================ *)

MakeBoxes[gs_GraphState ? GraphStateQ, form_] ^:= With[{n = gs["VertexCount"], m = gs["EdgeCount"]},
    BoxForm`ArrangeSummaryBox["GraphState",
        gs,
        Framed["\[ScriptCapitalG]"],
        {{BoxForm`SummaryItem[{"Vertices: ", n}]},
         {BoxForm`SummaryItem[{"Edges: ", m}]}},
        {{BoxForm`SummaryItem[{"Graph: ", gs["Graph"]}]}},
        form
    ]
]
