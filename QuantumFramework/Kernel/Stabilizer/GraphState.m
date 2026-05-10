Package["Wolfram`QuantumFramework`"]

PackageExport[GraphState]
PackageExport[LocalComplement]
PackageScope[GraphStateQ]



(* ============================================================================ *)
(* Graph-state representation (Anders & Briegel 2005).                          *)
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
(* footnote. Currently only index 0 (identity) is supported; richer VOPs are    *)
(* on the roadmap.                                                              *)
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


(* GraphState::nongraph emitted when the input PauliStabilizer is not in     *)
(* graph-state form. Only graph-form stabilizers (X_i Z_{N(i)}) are handled. *)
GraphState::nongraph = "PauliStabilizer `1` is not in graph-state form (each stabilizer must have X at exactly one position and Z elsewhere). Use a local-Clifford conversion (AndBri05 Lemma 1) before calling GraphState."

(* From a PauliStabilizer: build the canonical graph state for a "graph-form" *)
(* stabilizer state. Emits ::nongraph and returns $Failed when the input is  *)
(* not graph-form.                                                            *)
GraphState[ps_PauliStabilizer ? PauliStabilizerQ] := Module[{n, stabs, charSets, isGraphForm, edges, g},
    n = ps["Qubits"];
    stabs = ps["Stabilizers"];
    charSets = Characters @ If[StringStartsQ[#, "-"], StringDrop[#, 1], #] & /@ stabs;

    (* Graph-form check: there must be n stabilizer rows, the i-th row must
       have "X" at position i, and every other position must be "Z" or "I" (no
       "Y" anywhere, no "X" off the diagonal). *)
    isGraphForm = Length[charSets] == n &&
        AllTrue[
            Range[n],
            charSets[[#, #]] === "X" &&
            AllTrue[Drop[charSets[[#]], {#}], MatchQ[# /. "Z" -> "I", "I"] &] &
        ];

    If[!isGraphForm,
        Message[GraphState::nongraph, ps];
        Return[$Failed]
    ];

    edges = Catenate @ Table[
        Position[charSets[[i]], "Z"][[All, 1]] /. j_Integer :> If[i < j, UndirectedEdge[i, j], Nothing],
        {i, n}
    ];
    g = Graph[Range[n], DeleteDuplicates @ edges];
    GraphState[<|"Graph" -> g, "VOPs" -> ConstantArray[0, n]|>]
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
(* Currently returns the new Graph (with edges toggled) but does not update    *)
(* VOPs. VOP-tracked LC is on the roadmap.                                     *)
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
        "VOPs" -> gs["VOPs"]    (* TODO: track VOP updates per AndBri05 Eq 8 *)
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
