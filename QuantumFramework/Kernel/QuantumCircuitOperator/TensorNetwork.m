Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`Cotengra`"]

PackageExport["TensorNetworkIndexGraph"]
PackageExport["GraphTensorNetwork"]
PackageExport["TensorNetworkQuantumCircuit"]
PackageExport["TensorNetworkQ"]
PackageExport["ContractTensorNetwork"]
PackageExport["TensorNetworkIndices"]
PackageExport["TensorNetworkTensors"]
PackageExport["TensorNetworkData"]
PackageExport["TensorNetworkIndexDimensions"]
PackageExport["TensorNetworkFreeIndices"]
PackageExport["TensorNetworkAdd"]
PackageExport["RemoveTensorNetworkCycles"]
PackageExport["TensorNetworkContractionPath"]
PackageExport["TensorNetworkContractPath"]
PackageExport["TensorNetworkContraction"]
PackageExport["TensorNetworkNetGraph"]
PackageExport["TensorNetworkIndexReplace"]

PackageScope["InitializeTensorNetwork"]
PackageScope["TensorNetworkApply"]
PackageScope["TensorNetworkCompile"]
PackageScope["QuantumTensorNetwork"]
PackageScope["QuantumCircuitHypergraph"]



TensorNetworkQ::msg1 = "Not all vertices are annotated with tensors."
TensorNetworkQ::msg2 = "Not all indices are duplicate free lists of Subscripts and Superscripts."
TensorNetworkQ::msg3 = "Tensor ranks and indices are not compatible."
TensorNetworkQ::msg4 = "Not all edges are tagged with indices."

TensorNetworkQ[net_Graph, verbose : _ ? BooleanQ : False] := Module[{
    tensors, indices
},
    {tensors, indices} = AssociationThread[VertexList[net] -> #] & /@
        (AnnotationValue[{net, Developer`FromPackedArray[VertexList[net]]}, #] & /@ {"Tensor", "Index"});
    (* (
        AllTrue[tensors, TensorQ] ||
        (If[verbose, Message[TensorNetworkQ::msg1]]; False)
    ) && *)
    (
        AllTrue[indices, MatchQ[#, {(Superscript | Subscript)[_, _] ...}] && DuplicateFreeQ[#] &] ||
        (If[verbose, Message[TensorNetworkQ::msg2]]; False)
     ) &&
    With[{ranks = tensorRank /@ Values[tensors]},
        (And @@ Thread[ranks == Length /@ Values[indices]] && Total[ranks] == CountDistinct[Catenate[Values[indices]]]) ||
        Total[ranks] == 0 ||
        (If[verbose, Message[TensorNetworkQ::msg3]]; False)
    ] (* &&
    (
        AllTrue[
            EdgeList[net],
            MatchQ[
                _[from_, to_, {i_, j_}] /;
                MemberQ[indices[from], i] && MemberQ[indices[to], j]
            ]
        ] ||
        (If[verbose, Message[TensorNetworkQ::msg4]]; False)
    ) *)
]

TensorNetworkQ[verbose : _ ? BooleanQ : False][net_Graph] := TensorNetworkQ[net, verbose]


NeighborhoodEdges[g_, vs_List] := Catenate[EdgeList[g, _[#, __] | _[_, #, ___]] & /@ vs]
NeighborhoodEdges[g_, v_] := NeighborhoodEdges[g, {v}]

(* temporary due to EdgeContract bug *)
edgeContract[g_, edge_] := With[{edges = NeighborhoodEdges[g, edge[[1]]]},
	EdgeAdd[EdgeDelete[g, edges], Replace[DeleteCases[edges, edge], {
		head_[edge[[1]], edge[[1]], rest___] :> head[edge[[2]], edge[[2]], rest],
		head_[edge[[1]], rest___] :> head[edge[[2]], rest],
		head_[v_, edge[[1]], rest___] :> head[v, edge[[2]], rest]
	}, {1}]]
]

edgeContractWithIndex[g_, edge : _[from_, to_, {i_, j_}]] := Annotate[
	{edgeContract[g, edge], to},
	"Index" -> If[
		from === to,
		DeleteCases[AnnotationValue[{g, to}, "Index"], i | j],
		DeleteCases[Join[AnnotationValue[{g, from}, "Index"], AnnotationValue[{g, to}, "Index"]], i | j]
	]
]

ContractEdge[g_, edge : _[from_, to_, {i_, j_}]] := Enclose @ Module[{
	tensors = ConfirmMatch[AnnotationValue[{g, {from, to}}, "Tensor"], {__ ? TensorQ}],
	indices = Confirm[AnnotationValue[{g, {from, to}}, "Index"]],
	rank
},
	rank = tensorRank[tensors[[1]]];
	Annotate[
		{edgeContractWithIndex[g, edge], to},
		"Tensor" -> TensorContract[
			TensorProduct[tensors[[1]], tensors[[2]]],
			{Confirm[FirstPosition[indices[[1]], i]][[1]], rank + Confirm[FirstPosition[indices[[2]], j]][[1]]}
		]
	]
]

ContractEdge[g_, edge : _[v_, v_, {i_, j_}]] := Enclose @ Module[{
	tensor = ConfirmBy[AnnotationValue[{g, v}, "Tensor"], TensorQ],
	index = Confirm[AnnotationValue[{g, v}, "Index"]]
},
	Annotate[{edgeContractWithIndex[g, edge], v}, "Tensor" -> TensorContract[tensor, Catenate @ Position[index, i | j]]]
]


NaiveContractTensorNetwork[net_Graph] := Enclose @ Module[{g, edges},
	{g, {edges}} = Reap @ NestWhile[Confirm @ ContractEdge[#, Sow @ First[EdgeList[#]]] &, net, EdgeCount[#] > 0 &];
	Transpose[
        AnnotationValue[{g, edges[[-1, 2]]}, "Tensor"],
        Ordering @ OrderingBy[AnnotationValue[{g, edges[[-1, 2]]}, "Index"], Replace[{Subscript[_, x_] :> {x}, Superscript[_, x_] :> x}]]
    ]
]

FastContractTensorNetwork[net_Graph] := Enclose[
    Block[{indices, outIndices, tensors, scalarPositions, scalars},
        indices = TensorNetworkIndices[net] /. Rule @@@ EdgeTags[net];
        tensors =  TensorNetworkTensors[net];
        outIndices = TensorNetworkFreeIndices[net];
        If[MemberQ[tensors, {}], Return[ArrayReshape[{}, Append[Table[1, Length[outIndices] - 1], 0]]]];
        scalarPositions = Position[indices, {}, {1}, Heads -> False];
        scalars = Extract[tensors, scalarPositions];
        indices = Delete[indices, scalarPositions];
        tensors = Delete[tensors, scalarPositions];
        Times @@ scalars * ActivateTensor @ Confirm[EinsteinSummation[indices -> outIndices, tensors]]
    ],
    (ReleaseHold[#["HeldMessageCall"]]; #) &
]

Options[ContractTensorNetwork] = {Method -> Automatic}

ContractTensorNetwork[net_Graph ? (TensorNetworkQ[True]), OptionsPattern[]] := Switch[
    OptionValue[Method],
    "Naive",
    NaiveContractTensorNetwork[net],
    _,
    FastContractTensorNetwork[net]
]


TensorNetworkIndices[net_Graph ? TensorNetworkQ] := AnnotationValue[{net, Developer`FromPackedArray[VertexList[net]]}, "Index"]

TensorNetworkTensors[net_Graph ? TensorNetworkQ] := AnnotationValue[{net, Developer`FromPackedArray[VertexList[net]]}, "Tensor"]


TensorNetworkData[net_Graph ? TensorNetworkQ] := With[{vs = Developer`FromPackedArray[VertexList[net]]}, {
    tensors = AnnotationValue[{net, vs}, "Tensor"],
    indices = AnnotationValue[{net, vs}, "Index"],
    edges = EdgeList[net]
}, {
    dimensions = TensorNetworkIndexDimensions[indices, tensors],
    tags = edges[[All, 3]]
}, {
    contractions = Replace[indices, Catenate[Thread[# -> #, List, 1] & /@ tags], {2}]
},
    <|
        "Tensors" -> tensors,
        "Dimensions" -> TensorDimensions /@ tensors,
        "Indices" -> indices,
        "Vertices" -> vs,
        "FreeIndices" -> TensorNetworkFreeIndices[indices, tags],
        "Bonds" -> (#1 -> {##2} & @@@ MapAt[Lookup[dimensions, #[[1]]] &, edges, {All, 3}]),
        "Contractions" -> contractions,
        "ContractionIndices" -> Replace[contractions, {i_, _} :> i, {2}]
    |>
]


TensorNetworkFreeIndices[net_Graph ? TensorNetworkQ] :=
    TensorNetworkFreeIndices[TensorNetworkIndices[net], EdgeTags[net]]

TensorNetworkFreeIndices[indices_List, tags_List] :=
    SortBy[Replace[{Superscript[_, x_] :> {0, x}, Subscript[_, x_] :> {1, x}}]] @ DeleteElements[Catenate[indices], Catenate[tags]]


TensorNetworkIndexDimensions[net_Graph ? TensorNetworkQ] := TensorNetworkIndexDimensions @@ Lookup[TensorNetworkData[net], {"Indices", "Tensors"}]

TensorNetworkIndexDimensions[indices_List, tensors_List] := Catenate @ MapThread[Thread[#1 -> TensorDimensions[#2]] &, {indices, tensors}]


TensorNetworkIndexReplace[net_ ? TensorNetworkQ, rules_] :=
    Graph[net, AnnotationRules -> MapThread[#1 -> {"Index" -> #2} &, {VertexList[net], Replace[TensorNetworkIndices[net], rules, {2}]}]]


InitializeTensorNetwork[net_Graph ? TensorNetworkQ, tensor_, index_List : Automatic] := Annotate[
    {
        EdgeAdd[
            VertexDelete[net, _ ? NonPositive],
            MapIndexed[
                Replace[#1, DirectedEdge[_, i_, {_, to_}] :> DirectedEdge[0, i, {Superscript[0, #2[[1]]], to}]] &,
                EdgeList[net, DirectedEdge[_ ? NonPositive, __]]
            ]
        ],
        0
    },
    {
        "Tensor" -> tensor,
        "Index" -> Replace[index, Automatic :> (Superscript[0, #] & /@ Range[tensorRank[tensor]])],
        VertexLabels -> "Initial"
    }
]

TensorNetworkAdd[net_Graph ? TensorNetworkQ, Labeled[tensor_, label_ : None], autoIndex : _List | Automatic : Automatic] := Enclose @ With[{
    newVertex = Max[VertexList[net]] + 1,
    toIndex = Replace[autoIndex, Automatic :> Take[SortBy[TensorNetworkFreeIndices[net], Replace[{Superscript[_, x_] :> {1, x}, Subscript[_, x_] :> {0, x}}]], UpTo[tensorRank[tensor]]]]
},
{
    index = Join[
        Replace[toIndex, {Superscript[_, q_] :> Subscript[newVertex, q], Subscript[_, q_] :> Superscript[newVertex, q]}, {1}],
        Subscript[newVertex, #] & /@ Range[tensorRank[tensor] - Length[toIndex]]
    ]
},
    ConfirmAssert[tensorRank[tensor] == Length[index]];
    Annotate[
        {
            EdgeAdd[
                net,
                MapThread[
                    If[MatchQ[#1, _Superscript], DirectedEdge[newVertex, First[#2], {#1, #2}], DirectedEdge[First[#2], newVertex, {#2, #1}]] &,
                    {Take[index, UpTo[Length[toIndex]]], toIndex}
                ]
            ],
            newVertex
        },
        {
            "Tensor" -> tensor,
            "Index" -> index,
            VertexLabels -> label
        }
    ]
]

TensorNetworkAdd[net_Graph ? TensorNetworkQ, tensor_, index : _List | Automatic : Automatic] := TensorNetworkAdd[net, Labeled[tensor, None], index]

VertexCompleteGraph[vs_List] := With[{n = Length[vs]}, AdjacencyGraph[vs, SparseArray[Band[{1, 1}] -> 0, {n, n}, 1]]]

TensorNetworkIndexGraph[net_Graph ? (TensorNetworkQ[True]), opts : OptionsPattern[Graph]] := GraphUnion[
    DirectedEdge @@@ EdgeTags[net],
    Sequence @@ (VertexCompleteGraph /@ TensorNetworkIndices[net]),
    opts,
    VertexLabels -> Automatic
]

Options[QuantumTensorNetwork] = Join[{"PrependInitial" -> True, "Computational" -> True}, Options[Graph]]

QuantumTensorNetwork[qco_QuantumCircuitOperator, opts : OptionsPattern[]] := Enclose @ Block[{
    circuit = qco["Sort"], width, min, ops, orders, arity, vertices, edges, tensors
},
	ConfirmAssert[AllTrue[circuit["Operators"], #["Order"] === #["FullOrder"] &]];
    width = circuit["Width"];
    min = circuit["Min"];
    ops = circuit["NormalOperators"];
    If[ TrueQ[OptionValue["Computational"]],
        ops = Splice[{
            If[#["InputQudits"] > 0 && ! #["Input"]["ComputationalQ"], QuantumOperator[MatrixInverse[#["Input"]["ReducedMatrix"]], {#, #} & @ #["InputOrder"], QuantumBasis[#, #] & @ #["InputDimensions"], "Label" -> "I"[#["Label"]]], Nothing],
            #,
            If[#["OutputQudits"] > 0 && ! #["Output"]["ComputationalQ"], QuantumOperator[#["Output"]["ReducedMatrix"], {#, #} & @ #["OutputOrder"], QuantumBasis[#, #] & @ #["OutputDimensions"], "Label" -> "I"[#["Label"]]], Nothing]
        }] & /@ ops
    ];
    arity = circuit["Arity"];
    MapThread[
        PrependTo[ops, QuantumOperator[QuantumState[{1}, #2, "Label" -> "0"], {#1}]] &,
        Reverse /@ {circuit["InputOrder"], PadLeft[circuit["InputDimensions"], arity, 2]}
    ];
	orders = #["Order"] & /@ ops;
    vertices = Range[Length[ops]] - arity;
	edges = Catenate @ FoldPairList[
		{nprev, order} |-> Block[{output, input, n, prev, next, indices},
            {n, prev} = nprev;
            n += 1;
            next = prev;
			{output, input} = order;
            next[[ output - min + 1 ]] = n;
			next[[ Complement[input, output] - min + 1 ]] = None[n];
            indices = {Superscript[prev[[# - min + 1]], #], Subscript[next[[# - min + 1]], #]} & /@ input;
			{
                Replace[
                    Thread[DirectedEdge[prev[[ input - min + 1 ]], next[[ input - min + 1]], indices]],
                    {
                        DirectedEdge[None[_], ___] :> Nothing,
                        DirectedEdge[from_, None[to_], tag_] :> DirectedEdge[from, to, tag /. None[i_] :> i]
                    },
                    {1}
                ],
                {n, next}
            }
		],
		{1 - arity, Table[1 - arity, width]},
		Rest[orders]
	];
	tensors = If[#["MatrixQ"], #["Double"], #]["Tensor"] & /@ ops;
	ConfirmBy[
        If[TrueQ[OptionValue["PrependInitial"]] && circuit["Arity"] > 0, Identity, VertexDelete[#, _ ? NonPositive] &] @ Graph[
            vertices,
            edges,
            FilterRules[{opts}, Options[Graph]],
            AnnotationRules ->
                MapThread[#1 -> {
                        "Tensor" -> #2,
                        "Index" -> Join[OperatorApplied[Superscript, 2][#1] /@ Sort[#3[[1]]], OperatorApplied[Subscript, 2][#1] /@ Sort @ #3[[2]]]
                    } &,
                    {vertices, tensors, orders}
                ],
            VertexLabels ->
                Thread[vertices -> (Replace[#["Label"], {
                    label : Subscript["C", cop_][__] :> Interpretation[Row[{"C", cop}], label],
                    label : Subscript["R", rops__][angle_] :> Interpretation[Subscript["R", rops][angle], label]
                }] & /@ ops)],
            GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}
        ],
        TensorNetworkQ
    ]
]

QuantumCircuitHypergraph[qc_ ? QuantumCircuitOperatorQ, opts : OptionsPattern[]] := Enclose @ Block[{
    net = QuantumTensorNetwork[qc["Flatten"], FilterRules[{opts}, Except[Options[Graph], Options[QuantumTensorNetwork]]]], vs, indices, labels, edges
},
	Confirm[Needs["WolframInstitute`Hypergraph`" -> "H`"], "Hypergraph paclet is not installed."];
	vs = Developer`FromPackedArray @ VertexList[net];
	indices = TensorNetworkIndices[net];
	labels = AnnotationValue[{net, vs}, VertexLabels];
	edges = Replace[indices, Rule @@@ Reverse /@ EdgeTags[net], {2}];
	H`Hypergraph[Union @@ edges, edges, FilterRules[{opts}, Options[H`Hypergraph]], EdgeLabels -> Thread[edges -> labels]]
]


Options[TensorNetworkCompile] = Join[{"ReturnCircuit" -> False, "ReturnTensorNetwork" -> False, "Trace" -> True}, Options[QuantumTensorNetwork], Options[ContractTensorNetwork]]

TensorNetworkCompile[qco_QuantumCircuitOperator, opts : OptionsPattern[]] := Enclose @ Block[{
    circuit = qco["Normal"], width, net, phaseSpaceQ, bendQ, order, res,
    traceOrder, eigenOrder, computationalQ, basis
},
    width = circuit["Width"];
    basis = Confirm @ circuit["TensorNetworkBasis"];
    phaseSpaceQ = basis["Picture"] === "PhaseSpace";
    computationalQ = TrueQ[OptionValue["Computational"]] || ! phaseSpaceQ &&
        TrueQ[Or @@ Unequal @@@
            ConfirmBy[
                With[{info = Confirm @ circuit["TensorNetworkInfo"]},
                    Partition[Lookup[info["QuditBases"], Catenate[info["ContractionIndices"]]], 2]
                ],
                AllTrue[Equal @@ Through[#["Dimension"]] &]
            ]
        ];
    traceOrder = circuit["TraceOrder"];
    eigenOrder = circuit["Eigenorder"];
    order = Sort /@ circuit["Order"];
    bendQ = (AnyTrue[circuit["Operators"], #["MatrixQ"] &] || circuit["TraceQudits"] > 0) && ! phaseSpaceQ;
    If[ bendQ,
        (* TODO: handle this case somehow *)
        (* ConfirmAssert[AllTrue[Join[DeleteElements[order[[1]], Join[traceOrder, eigenOrder]], order[[2]]], Positive]]; *)
        circuit = circuit["Bend"];
        If[ TrueQ[OptionValue["Trace"]] && circuit["TraceQudits"] > 0,
            circuit = circuit /* MapThread[{"Cap", #2} -> {#1, #1 + width} &, {circuit["TraceOrder"], circuit["TraceBasis"]["Decompose"]}];
        ]
    ];
    If[TrueQ[OptionValue["ReturnCircuit"]], Return[circuit]];
    net = ConfirmBy[QuantumTensorNetwork[circuit, "Computational" -> computationalQ, FilterRules[{opts}, Options[QuantumTensorNetwork]], "PrependInitial" -> False], TensorNetworkQ];
    If[TrueQ[OptionValue["ReturnTensorNetwork"]], Return[net]];
    res = Confirm @ ContractTensorNetwork[net, FilterRules[{opts}, Options[ContractTensorNetwork]]];
    res = With[{basis = Confirm @ circuit["TensorNetworkBasis"]},
        QuantumState[
            SparseArrayFlatten[res],
            If[computationalQ, QuantumBasis[QuditBasis[basis["OutputDimensions"]], QuditBasis[basis["InputDimensions"]]], basis]
        ]
    ];
    If[ TrueQ[OptionValue["Trace"]] && traceOrder =!= {},
        If[ bendQ,
            basis = QuantumPartialTrace[basis, traceOrder - qco["Min"] + 1],

            res = QuantumPartialTrace[res, traceOrder - qco["Min"] + 1]
        ];
        order = {DeleteElements[order[[1]], traceOrder], order[[2]]};
    ];
    If[ bendQ,
        If[ eigenOrder =!= {},
            order = {Join[Take[Select[order[[1]], NonPositive], - Length[eigenOrder]], Select[order[[1]], Positive]], order[[2]]};
        ];
        res = QuantumState[If[computationalQ, res, res["State"]], QuantumBasis[{basis, basis["Conjugate"]}]]["Unbend"]
        ,
        res = QuantumState[If[computationalQ, res, res["State"]], basis];
    ];
    res = Which[
        eigenOrder =!= {},
        QuantumMeasurementOperator[QuantumOperator[res, order], qco["Targets"]],
        ! TrueQ[OptionValue["Trace"]] && traceOrder =!= {},
        QuantumChannel[QuantumOperator[res, order]],
        True,
        QuantumOperator[res, order]
    ];
    res
]

Options[TensorNetworkApply] = Options[TensorNetworkCompile]

TensorNetworkApply[qco_QuantumCircuitOperator, qs_QuantumState, opts : OptionsPattern[]] := Block[{
    circuit = qco["Sort"], res
},
    If[ qs["Qudits"] > 0,
        circuit = QuantumCircuitOperator[{qs -> circuit["FullInputOrder"]} /* circuit, "Label" -> None]
    ];
    res = TensorNetworkCompile[circuit, opts];
    Which[
        QuantumMeasurementOperatorQ[res],
        If[ContainsAll[res["OutputOrder"], res["TargetOrder"]], QuantumMeasurement[res], res],
        QuantumOperatorQ[res],
        res["State"],
        True,
        res
    ]
]

Options[GraphTensorNetwork] = {Method -> "Random"}

GraphTensorNetwork[g_ /; DirectedGraphQ[g] && AcyclicGraphQ[g], OptionsPattern[]] := Enclose @ Block[{
	vs = Developer`FromPackedArray[TopologicalSort[g]], edges = EdgeList[g],
	labels,
	inputs, outputs, inputOrder, outputOrder,
    orders, output,
	indices, outOrders, inOrders
},
	orders = <||>;
	output = {};
	Do[
		inputs = VertexInComponent[g, v, {1}];
		outputs = VertexOutComponent[g, v, {1}];
		inputOrder = Lookup[First /@ PositionIndex[output], inputs, Nothing];
		inputOrder = Take[Join[inputOrder, Length[output] + Range[Length[inputs] - Length[inputOrder]]], UpTo[Length[inputs]]];
		If[ Length[outputs] <= Length[inputs],
			outputOrder = Take[inputOrder, Length[outputs]],
			outputOrder = Take[Join[inputOrder, Lookup[PositionIndex[output], None, {}], Length[output] + Range[Length[outputs]]], Length[outputs]]
		];
		output = PadRight[output, Max[Length[output], outputOrder], None];
		output[[Complement[inputOrder, outputOrder]]] = None;
		Scan[(output[[#]] = v) &, outputOrder];
		orders[v] = {outputOrder, inputOrder};
		,
		{v, vs}
	];
	indices = AssociationThread[vs, AnnotationValue[{g, vs}, "Index"]];
	outOrders = KeyValueMap[{v, i} |->
		Replace[i, {
			$Failed :> orders[v][[1]],
			is_List :> Cases[is, Superscript[_, o_] :> o]
		}],
		indices
	];
	inOrders = KeyValueMap[{v, i} |->
		Replace[i, {
			$Failed :> orders[v][[2]],
			is_List :> Cases[is, Subscript[_, o_] :> o]
		}],
		indices
	];
    indices = Association @ MapThread[
        {v, idx, in, out} |-> v -> Replace[idx, $Failed :> Join[Superscript[v, #] & /@ Sort[out], Subscript[v, #] & /@ Sort[in]]],
        {vs, Values[indices], inOrders, outOrders}
    ];
	labels = KeyValueMap[
        Replace[#2, {$Failed | Automatic :> #1, Interpretation[_, label_] :> label}] &,
        AssociationThread[vs, AnnotationValue[{g, vs}, VertexLabels]]
    ];
    Block[{in = Cases[_Subscript] /@ indices, out = Cases[_Superscript] /@ indices, outIdx, inIdx},
        edges = Table[
            {outIdx, inIdx} = FirstCase[#, {_[_, x_], _[_, x_]}, First[#]] & @ Tuples[{Lookup[out, edge[[1]]], Lookup[in, edge[[2]]]}];
            out[edge[[1]]] //= DeleteCases[outIdx];
            in[edge[[2]]] //= DeleteCases[inIdx];
            Append[edge[[;; 2]], {outIdx, inIdx}]
            ,
            {edge, Catenate @ Lookup[GroupBy[edges, First], vs, {}]}   
        ]
    ];
	Graph[
        vs, edges,
        AnnotationRules -> MapThread[
            {v, t, i} |-> With[{label = labels[[i]], index = indices[[i]]},
                v -> {
                    VertexLabels -> label,
                    "Tensor" -> Replace[t, {
                        $Failed :> With[{qubits = Length[index]}, Switch[
                            OptionValue[Method],
                            "Zero", QuantumState[{"Register", qubits}]["Tensor"],
                            "Symbolic", ArraySymbol[label, ConstantArray[2, qubits]],
                            "Random", QuantumState[{"RandomPure", qubits}]["Tensor"]
                        ]]
                    }],
                    "Index" -> index
                }
            ],
            {vs, AnnotationValue[{g, vs}, "Tensor"], Range[Length[vs]]}
        ],
        Options[g]
    ]
]

GraphTensorNetwork[g_ ? DirectedGraphQ, opts : OptionsPattern[]] :=
    Enclose @ GraphTensorNetwork[ConfirmBy[RemoveTensorNetworkCycles[g], AcyclicGraphQ], opts]

GraphTensorNetwork[g_ ? GraphQ, opts : OptionsPattern[]] := GraphTensorNetwork[DirectedGraph[g, "Acyclic"], opts]

TensorNetworkQuantumCircuit[tn_ ? TensorNetworkQ] := QuantumCircuitOperator @ MapThread[{label, tensor, indices} |->
    With[{
        order = Lookup[GroupBy[indices, Head -> Last], {Superscript, Subscript}, {}]
    },
        QuantumOperator[
            QuantumState[Flatten[tensor]]["SplitDual", Length[order[[1]]]],
            order,
            "Label" -> label
        ]
    ] // Replace[label, {Subscript["Measurement", target___] :> (QuantumMeasurementOperator[#, {target}] &), _ -> Identity}],
    AnnotationValue[{tn, Developer`FromPackedArray[VertexList[tn]]}, #] & /@ {VertexLabels, "Tensor", "Index"}
]

RemoveTensorNetworkCycles[inputNet_ ? DirectedGraphQ, opts : OptionsPattern[Graph]] := Enclose @ Block[{
    net = IndexGraph[inputNet], cycles, id, q, r, edge, tag, cup, cap, cupIndex, capIndex, dim
},
	id = Max[VertexList[net], 0] + 1;
	{q, r} = MinMax[EdgeTags[net][[All, All, 2]]];
    q = Min[q, 1];
    r = Max[r, 1];


	While[
        Length[cycles = FindCycle[net, Infinity, 1]] > 0,

        q--;
		edge = cycles[[1, -1]];
        tag = If[Length[edge] == 3 && MatchQ[edge[[3]], {Superscript[_Integer, _Integer], Subscript[_Integer, _Integer]}],
            edge[[3]],
            None
        ];
		cup = id++;
		cap = id++;
		net = EdgeDelete[net, edge];
		net = VertexAdd[net, {cap, cup}];
        If[ tag =!= None,

            cupIndex = {Superscript[cup, q], Superscript[cup, tag[[2, 2]]]};
            capIndex = {Subscript[cap, q], Subscript[cap, tag[[1, 2]]]};
            net = EdgeAdd[net, {
                DirectedEdge[cup, cap, {cupIndex[[1]], capIndex[[1]]}],
                DirectedEdge[cup, edge[[2]], {cupIndex[[2]], tag[[2]]}],
                DirectedEdge[edge[[1]], cap, {tag[[1]], capIndex[[2]]}]
            }];
            dim = Enclose[
                ConfirmBy[
                    Dimensions[Confirm[AnnotationValue[{net, edge[[1]]}, "Tensor"]]][[ Confirm @ Lookup[PositionIndex[AnnotationValue[{net, edge[[1]]}, "Index"]], tag[[1]], $Failed, First] ]],
                    IntegerQ
                ],
                2 &
            ];
            net = Annotate[{net, cup}, "Index" -> cupIndex];
		    net = Annotate[{net, cap}, "Index" -> capIndex];

            ,

            net = EdgeAdd[net, {
                DirectedEdge[cup, cap],
                DirectedEdge[cup, edge[[2]]],
                DirectedEdge[edge[[1]], cap]
            }];
            dim = Enclose[First[Dimensions[Confirm[AnnotationValue[{net, edge[[1]]}, "Tensor"]]], 1], 2 &]
        ];

		net = Annotate[{net, cup}, {"Tensor" -> QuantumOperator["Cup"[dim]]["Tensor"], VertexLabels -> "Cup"}];
		net = Annotate[{net, cap}, {"Tensor" -> QuantumOperator["Cap"[dim]]["Tensor"], VertexLabels -> "Cap"}];
	];
	Graph[net, opts]
]



Options[TensorNetworkContractionPath] = {"ReturnParameters" -> False, "Optimal" -> False}

TensorNetworkContractionPath[net_ ? TensorNetworkQ, OptionsPattern[]] := Enclose @ Block[{
	tensors, indices, pairs, freeIndices, normalIndices, input, output, dimensions
},
	tensors = TensorNetworkTensors[net];
    indices = TensorNetworkIndices[net];
	pairs = EdgeTags[net];
	freeIndices = TensorNetworkFreeIndices[net];
	dimensions = AssociationThread[Catenate[indices], Catenate[Dimensions /@ tensors]];
	ConfirmAssert[AllTrue[Partition[Lookup[dimensions, Catenate[pairs]], 2], Apply[Equal]]];
	pairs = Rule @@@ pairs;
	dimensions = KeyMap[Replace[pairs], dimensions];
	indices = Replace[indices, pairs, {2}];
	normalIndices = Thread[# -> Range[Length[#]]] & [Union @@ indices];
	input = Replace[indices, normalIndices, {2}];
	output = Replace[freeIndices, normalIndices, {1}];
	dimensions = KeyMap[Replace[normalIndices], dimensions];
	If[TrueQ[OptionValue["ReturnParameters"]], Return[{input, output, dimensions}]];
	If[ TrueQ[OptionValue["Optimal"]],
        OptimalPath[input, output, dimensions, "size"],
        GreedyPath[input, output, dimensions]
    ]
]


einsum[{i_, j_} -> k_, a_, b_] := Block[{c = Complement[Join[i, j], k], adim = Dimensions[a], bdim = Dimensions[b], al, br, ac, bc, ad, bd},
	ac = Catenate @ Lookup[PositionIndex[i], c];
	bc = Catenate @ Lookup[PositionIndex[j], c];
	al = Complement[Range[Length[i]], ac];
	br = Complement[Range[Length[j]], bc];
	ad = Times @@ adim[[ac]];
	bd = Times @@ bdim[[bc]];
	Transpose[
		ArrayReshape[
			ArrayReshape[Transpose[a, FindPermutation[Join[al, ac]]], {Times @@ adim / ad, ad}] . ArrayReshape[Transpose[b, FindPermutation[Join[bc, br]]], {bd, Times @@ bdim / bd}],
			Join[adim[[al]], bdim[[br]]]
		],
		FindPermutation[Join[i[[al]], j[[br]]], k]
	]
]


TensorNetworkContractPath[
	KeyValuePattern[{
		"Tensors" -> initTensors_,
		"ContractionIndices" -> initIndices_,
		"FreeIndices" -> freeIndices_
	}],
	path_
] := Enclose @ Block[{tensors = initTensors, indices = initIndices},
    Do[
		Replace[p, {
			{i_} :> (
				{tensors, indices} = Append[Delete[#, {i}], #[[i]]] & /@ {tensors, indices}
			),
			{i_, j_} :> Block[{out, tensor},
				out = SymmetricDifference @@ Extract[indices, {{i}, {j}}];
				tensor = ActivateTensor @ EinsteinSummation[indices[[{i, j}]] -> out, {tensors[[i]], tensors[[j]]}];
				(* tensor = einsum[indices[[{i, j}]] -> out, tensors[[i]], tensors[[j]]]; *)

				tensors = Append[Delete[tensors, {{i}, {j}}], tensor];
				indices = Append[Delete[indices, {{i}, {j}}], out]
			]
		}],
		{p, path}
	];
	ConfirmAssert[Length[tensors] == Length[indices] == 1];
	ConfirmAssert[ContainsAll[indices[[1]], freeIndices]];
	Transpose[tensors[[1]], FindPermutation[indices[[1]], freeIndices]]
]

TensorNetworkContractPath[net_ ? TensorNetworkQ, path_] := TensorNetworkContractPath[TensorNetworkData[net], path]


einsumArrayDot[{i_, j_} -> out_, a_, b_, inactiveQ : _ ? BooleanQ : False] := Block[{
	c = DeleteElements[DeleteDuplicates @ Join[i, j], Replace[out, Automatic :> SymmetricDifference[i, j]]],
	k, perm,
	al, br, x, y,
	aIndex = First /@ PositionIndex[i], bIndex = First /@ PositionIndex[j],
	inactive = If[inactiveQ, Function[f, Inactive[f][##] &], Identity]
},
	al = DeleteElements[i, c];
	br = DeleteElements[j, c];
	k = Length[c];
	If[ k == 0
		,
		x = inactive[TensorProduct][a, b]
		,
		x = inactive[ArrayDot][a, b, Thread[{Lookup[aIndex, c], Lookup[bIndex, c]}]];
	];
	If[ out === Automatic,
		{x, Join[al, br]}
		,
		perm = FindPermutation[Join[al, br], out];
		If[perm === Cycles[{}], x, inactive[Transpose][x, perm]]
	]
]

einsumArrayDotTranspose[{i_, j_} -> out_, a_, b_, inactiveQ : _ ? BooleanQ : False] := Block[{
	c = DeleteElements[DeleteDuplicates @ Join[i, j], Replace[out, Automatic :> SymmetricDifference[i, j]]],
	k, perm,
	al, br, x, y,
	inactive = If[inactiveQ, Function[f, Inactive[f][##] &], Identity]
},
	al = DeleteElements[i, c];
	br = DeleteElements[j, c];
	k = Length[c];
	If[ k == 0
		,
		x = inactive[TensorProduct][a, b]
		,
		perm = FindPermutation[i, Join[al, c]];
		x = If[perm === Cycles[{}], a, inactive[Transpose][a, perm]];
		perm = FindPermutation[j, Join[c, br]];
		y = If[perm === Cycles[{}], b, inactive[Transpose][b, perm]];
		x = inactive[ArrayDot][x, y, k];
	];
	If[ out === Automatic,
		{x, Join[al, br]}
		,
		perm = FindPermutation[Join[al, br], out];
		If[perm === Cycles[{}], x, inactive[Transpose][x, perm]]
	]
]

einsumTensorContract[{i_, j_} -> out_, a_, b_, inactiveQ : _ ? BooleanQ : False] := Block[{
	c = DeleteElements[DeleteDuplicates @ Join[i, j], Replace[out, Automatic :> SymmetricDifference[i, j]]],
	k, perm,
	al, br, x,
	aIndex = First /@ PositionIndex[i], bIndex = First /@ PositionIndex[j],
	inactive = If[inactiveQ, Function[f, Inactive[f][##] &], Identity]
},
	al = DeleteElements[i, c];
	br = DeleteElements[j, c];
	k = Length[c];
	If[ k == 0
		,
		x = inactive[TensorProduct][a, b]
		,
		x = inactive[TensorContract][
			inactive[TensorProduct][a, b],
			MapThread[{#1, #2 + Length[i]} &,
				{Lookup[aIndex, c], Lookup[bIndex, c]}
			]
		]
	];
	If[ out === Automatic,
		{x, Join[al, br]}
		,
		perm = FindPermutation[Join[al, br], out];
		If[perm === Cycles[{}], x, inactive[Transpose][x, perm]]
	]
]

einsumDot[{i_, j_} -> out_, a_, b_, inactiveQ : _ ? BooleanQ : False] := Block[{
	c = DeleteElements[DeleteDuplicates @ Join[i, j], Replace[out, Automatic :> SymmetricDifference[i, j]]],
	k, perm,
	al, br, x, y,
	inactive = If[inactiveQ, Function[f, Inactive[f][##] &], Identity],
	aIndex = PositionIndex[i], bIndex = PositionIndex[j],
	aDim = symbolicTensorDimensions[a], bDim = symbolicTensorDimensions[b],
	reshape
},
	reshape[t_, newShape_] := If[symbolicTensorDimensions[t] === newShape, t, inactive[ArrayReshape][t, newShape]];
	al = DeleteElements[i, c];
	br = DeleteElements[j, c];
	k = Length[c];
	If[ k == 0,
		x = inactive[TensorProduct][a, b]
		,
		perm = FindPermutation[i, Join[al, c]];
		x = If[perm === Cycles[{}], a, inactive[Transpose][a, perm]];
		perm = FindPermutation[j, Join[c, br]];
		y = If[perm === Cycles[{}], b, inactive[Transpose][b, perm]];
		If[ k == 1,
			x = inactive[Dot][x, y]
			,
			x = inactive[Dot][
				reshape[x, Catenate[MapAt[{Times @@ #} &, {2}] @ (Extract[aDim, Lookup[aIndex, #]] & /@ {al, c})]],
				reshape[y, Catenate[MapAt[{Times @@ #} &, {1}] @ (Extract[bDim, Lookup[bIndex, #]] & /@ {c, br})]]
			];
			x = reshape[x, Join[Extract[aDim, Lookup[aIndex, al]], Extract[bDim, Lookup[bIndex, br]]]]
		];
	];
	If[ out === Automatic,
		{x, Join[al, br]}
		,
		perm = FindPermutation[Join[al, br], out];
		If[perm === Cycles[{}], x, inactive[Transpose][x, perm]]
	]
]

einsumTableSum[{i_, j_} -> out_, a_, b_, inactiveQ : _ ? BooleanQ : False] := Block[{
	c = DeleteElements[DeleteDuplicates @ Join[i, j], Replace[out, Automatic :> SymmetricDifference[i, j]]],
	k, perm,
	al, br, x, y,
	inactive = If[inactiveQ, Function[f, Inactive[f][##] &], Identity],
	aIndex = PositionIndex[i], bIndex = PositionIndex[j],
	aDim = symbolicTensorDimensions[a], bDim = symbolicTensorDimensions[b]
},
	
	al = DeleteElements[i, c];
	br = DeleteElements[j, c];
	k = Length[c];
	
	If[ k == 0
		,
		x = With[{
			p1 = Symbol["\[FormalI]" <> ToString[#]] & /@ Range[Length[al]],
			p2 = Symbol["\[FormalJ]" <> ToString[#]] & /@ Range[Length[br]]
		},
			inactive[With][{inactive[Set][\[FormalCapitalA], a], inactive[Set][\[FormalCapitalB], b]},
				inactive[Table][(inactive[Part][\[FormalCapitalA], ##] & @@ p1) * (inactive[Part][\[FormalCapitalB], ##] & @@ p2), ##] & @@ Join[
					MapIndexed[{Symbol["\[FormalI]" <> ToString[#2[[1]]]], #1} &, aDim],
					MapIndexed[{Symbol["\[FormalJ]" <> ToString[#2[[1]]]], #1} &, bDim]
				]
			]
		]
		,
		x = With[{
			p1 = ReplacePart[i, Join[
				Thread[Lookup[aIndex, al] -> (Symbol["\[FormalI]" <> ToString[#]] & /@ Range[Length[al]])],
				Thread[Lookup[aIndex, c] -> (Symbol["\[FormalC]" <> ToString[#]] & /@ Range[Length[c]])]
			]],
			p2 = ReplacePart[j, Join[
				Thread[Lookup[bIndex, br] -> (Symbol["\[FormalJ]" <> ToString[#]] & /@ Range[Length[br]])],
				Thread[Lookup[bIndex, c] -> (Symbol["\[FormalC]" <> ToString[#]] & /@ Range[Length[c]])]
			]],
			cs = MapIndexed[{Symbol["\[FormalC]" <> ToString[#2[[1]]]], #1} &, Extract[aDim, Lookup[aIndex, c]]]
		},
			inactive[With][{inactive[Set][\[FormalCapitalA], a], inactive[Set][\[FormalCapitalB], b]},
				inactive[Table][
					inactive[Sum][(inactive[Part][\[FormalCapitalA], ##] & @@ p1) * (inactive[Part][\[FormalCapitalB], ##] & @@ p2), ##] & @@ cs,
					## 
				] & @@ Join[
					MapIndexed[{Symbol["\[FormalI]" <> ToString[#2[[1]]]], #1} &, Extract[aDim, Lookup[aIndex, al]]],
					MapIndexed[{Symbol["\[FormalJ]" <> ToString[#2[[1]]]], #1} &, Extract[bDim, Lookup[bIndex, br]]]
				]
			]	
		]
	];
	If[ out === Automatic,
		{x, Join[al, br]}
		,
		perm = FindPermutation[Join[al, br], out];
		If[perm === Cycles[{}], x, inactive[Transpose][x, perm]]
	]
]


Options[contractTensorPair] = {"EinsumFunction" -> "ArrayDot", "Inactivate" -> True}

contractTensorPair[{\[FormalCapitalT][tensor1_, indices1_, flop1_], \[FormalCapitalT][tensor2_, indices2_, flop2_]}, opts : OptionsPattern[]] :=
	\[FormalCapitalT] @@ Append[0] @ Switch[
		OptionValue["EinsumFunction"],
			"ArrayDotTranspose", einsumArrayDotTranspose,
			"ArrayDot", einsumArrayDot,
			"Dot", einsumDot,
			"TensorContract", einsumTensorContract,
			"TableSum", einsumTableSum
		][{indices1, indices2} -> Automatic, tensor1, tensor2, TrueQ[OptionValue["Inactivate"]]]


Options[TensorNetworkContraction] = Join[Options[contractTensorPair], {"TransposeFunction" -> Transpose}]

TensorNetworkContraction[net_Graph ? TensorNetworkQ, path_, opts : OptionsPattern[]] :=
    TensorNetworkContraction[TensorNetworkData[net], path, opts]

TensorNetworkContraction[netData : KeyValuePattern["Vertices" -> vertices_], path : {{_Integer, _Integer} ...}, opts : OptionsPattern[]] := 
    TensorNetworkContraction[netData, PathToTreePath[path, vertices], opts]

TensorNetworkContraction[
    KeyValuePattern[{
		"Vertices" -> vertices_,
		"Tensors" -> tensors_,
		"Contractions" -> contractions_,
        "FreeIndices" -> freeIndices_
	}],
    treePath_,
    opts : OptionsPattern[]
] := With[{
    contractOpts = FilterRules[{opts}, Options[contractTensorPair]]
}, {
    tensorPath = NestWhile[
        ReplaceAll[tensorPair : {_\[FormalCapitalT], _\[FormalCapitalT]} :> contractTensorPair[tensorPair, contractOpts]],
        Replace[treePath, MapThread[{#1} -> \[FormalCapitalT][#2, #3, 0] &, {vertices, tensors, contractions}], {-2}],
        Not @* MatchQ[_\[FormalCapitalT]]
    ],
    transposeFunction = OptionValue["TransposeFunction"]
},
    With[{perm = FindPermutation[tensorPath[[2]], freeIndices]},
		If[ perm === Cycles[{}],
			tensorPath[[1]],
			If[TrueQ[OptionValue["Inactivate"]], Inactive, Identity][transposeFunction][
				tensorPath[[1]],
				FindPermutation[tensorPath[[2]], freeIndices]
			]
		]
    ]
]


Options[ToSymbolicTensor] = {"ArrayDotExpand" -> False}

ToSymbolicTensor[t_, OptionsPattern[]] := Activate[t /. {
	c_Cycles :> c,
	IgnoringInactive[Verbatim[With][vars_, body_]] :> ArraySymbol["T", symbolicTensorDimensions[body(* /. Rule @@@ vars*)]],
	IgnoringInactive[Verbatim[ArrayReshape][a_, d_]] :> ArraySymbol[ToSymbolicTensor[a], d],
	IgnoringInactive[Verbatim[ArrayDot][a_, b_, k_Integer]] :> If[TrueQ[OptionValue["ArrayDotExpand"]],
		With[{d1 = symbolicTensorDimensions[a], d2 = symbolicTensorDimensions[b]},
			Inactive[Dot][
				ArraySymbol[ToSymbolicTensor[a], Append[d1[[;; - k - 1]], Times @@ d1[[- k ;;]]]],
				ArraySymbol[ToSymbolicTensor[b], Prepend[d2[[k + 1 ;;]], Times @@ d2[[;; k]]]]
			]
		],
		Inactive[ArrayDot][ToSymbolicTensor[a], ToSymbolicTensor[b], k]
	],
	IgnoringInactive[Verbatim[ArrayDot][a_, b_, indices : {{_Integer, _Integer} ...}]] :> If[TrueQ[OptionValue["ArrayDotExpand"]],
		With[{d1 = symbolicTensorDimensions[a], d2 = symbolicTensorDimensions[b]},
			Inactive[Dot][
				ArraySymbol[
					Inactive[Transpose][ToSymbolicTensor[a], FindPermutation[Range[Length[d1]], Join[Complement[Range[Length[d1]], #], #] & @ indices[[All, 1]]]],
					Append[Delete[d1, List /@ indices[[All, 1]]], Times @@ d1[[indices[[All, 1]]]]]
				],
				ArraySymbol[
					Inactive[Transpose][ToSymbolicTensor[b], FindPermutation[Range[Length[d2]], Join[#, Complement[Range[Length[d2]], #]] & @ indices[[All, 2]]]],
					Append[Times @@ d2[[indices[[All, 2]]]], Delete[d2, List /@ indices[[All, 2]]]]
				]
			]
		],
		Inactive[ArrayDot][ToSymbolicTensor[a], ToSymbolicTensor[b], indices]
	],
	IgnoringInactive[Verbatim[Table][_, iter : {_, _Integer} ..]] :> ArraySymbol["T", {iter}[[All, 2]]],
	IgnoringInactive[(h : TensorContract | Transpose)[body_, args__]] :> Inactive[h][ToSymbolicTensor[body], args],
	tt_ ? TensorQ :> ArraySymbol["T", TensorDimensions[tt]]
}, Except[TensorContract | Dot | ArrayDot | TensorProduct]]

symbolicTensorDimensions[t_] := Replace[TensorDimensions[ToSymbolicTensor[t]], Except[_List] -> {}]

symbolicTensorRank[t_] := Length[symbolicTensorDimensions[t]]


TensorNetworkNetGraph[net_ ? TensorNetworkQ] := TensorNetworkNetGraph[net, TensorNetworkContractionPath[net]]

TensorNetworkNetGraph[net_ ? TensorNetworkQ, path_] := Enclose @ Block[{tensors, indices, freeIndices, g, tensorQueue, addEinsumLayer, n},
    tensors = Chop @* FullSimplify @* N @* Normal /@ TensorNetworkTensors[net];
    indices = TensorNetworkIndices[net];
    freeIndices = TensorNetworkFreeIndices[net];
    indices = Replace[indices, Rule @@@ EdgeTags[net], {2}];
    n = Length[tensors];
    g = NetGraph[
        Block[{body, symbols},
            symbols = Union @@ Reap[body = Replace[#, s_Symbol /; ! NumericQ[s] :> Slot @@ {Sow[ToString[s]]}, Infinity]][[2]];
            If[ symbols === {},
                Confirm @ NetArrayLayer["Array" -> NumericArray[#]],
                Confirm @ FunctionLayer[Function[Evaluate[body]], Sequence @@ (# -> {} & /@ symbols)]
            ]
        ] & /@ tensors,
        # -> NetPort["T" <> ToString[#]] & /@ Range[n]
    ];
    tensorQueue = "T" <> ToString[#] & /@ Range[n];
    addEinsumLayer[{i_, j_} -> k_, a_, b_] :=
		With[{
			c = Complement[Join[i, j], k],
			adim = Flatten[{NetExtract[g, a]}],
			bdim = Flatten[{NetExtract[g, b]}]
	    },
		With[{
			ac = Catenate @ Lookup[PositionIndex[i], c],
			bc = Catenate @ Lookup[PositionIndex[j], c]
		},
		With[{
			al = Complement[Range[Length[i]], ac],
			br = Complement[Range[Length[j]], bc],
			ad = Times @@ adim[[ac]],
			bd = Times @@ bdim[[bc]]
		},
		With[{
			reshapea = {Times @@ adim / ad, ad},
			reshapeb = {bd, Times @@ bdim / bd},
			perma = PermutationList @ FindPermutation[Join[al, ac]],
			permb = PermutationList @ FindPermutation[Join[bc, br]],
			reshape = Join[adim[[al]], bdim[[br]]],
			perm = PermutationList @ FindPermutation[Join[i[[al]], j[[br]]], k]
		},
		n++;
		g = Confirm @ NetGraph[
			{
				g,
				Confirm @ NetGraph[<|
                    (*ToString[n] -> FunctionLayer[Apply[
                        Transpose[
                            ArrayReshape[
                                ArrayReshape[Transpose[#1, perma], reshapea] . ArrayReshape[Transpose[#2, permb], reshapeb],
                                reshape
                            ],
                            perm
                        ] &
                    ]]*)
					"a" -> NetChain[{If[perma === {}, Nothing, TransposeLayer[perma]], ReshapeLayer[reshapea]}],
					"b" -> NetChain[{If[permb === {}, Nothing, TransposeLayer[permb]], ReshapeLayer[reshapeb]}],
					"dot" -> DotLayer[],
					"post" -> NetChain[{ReshapeLayer[reshape], If[perm === {}, Nothing, TransposeLayer[perm]]}]
				|>, {NetPort[a] -> "a", NetPort[b] -> "b", {"a", "b"} -> "dot" -> "post" -> NetPort["T" <> ToString[n]]}]
			},
			{NetPort[{1, a}] -> NetPort[{2, a}], NetPort[{1, b}] -> NetPort[{2, b}]}
		];
	]]]];
    Do[
		Replace[p, {
			{i_} :> (
				{tensorQueue, indices} = Append[Delete[#, {i}], #[[i]]] & /@ {tensorQueue, indices}
			),
			{i_, j_} :> With[{out = SymmetricDifference @@ Extract[indices, {{i}, {j}}]},
				addEinsumLayer[indices[[{i, j}]] -> out, tensorQueue[[i]], tensorQueue[[j]]];
				tensorQueue = Append[Delete[tensorQueue, {{i}, {j}}], "T" <> ToString[n]];
				indices = Append[Delete[indices, {{i}, {j}}], out];
			]
		}],
		{p, path}
	];
	ConfirmAssert[Length[tensorQueue] == Length[indices] == 1];
	ConfirmAssert[ContainsAll[indices[[1]], freeIndices]];
	With[{perm = PermutationList @ FindPermutation[indices[[1]], freeIndices]},
		If[perm === {}, g, NetAppend[g, TransposeLayer[perm]]]
	]
]

