Package["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`TensorNetworks`"]

PackageExport["TensorNetworkQuantumCircuit"]

PackageScope[tensorNetworkIndexSort]

PackageScope["TensorNetworkApply"]
PackageScope["TensorNetworkCompile"]
PackageScope["QuantumTensorNetworkGraph"]
PackageScope["QuantumTensorNetwork"]
PackageScope["QuantumCircuitHypergraph"]


tensorNetworkIndexSort[indices_List] := Catenate[Lookup[GroupBy[indices, MatchQ[_Superscript], SortBy[Last]], {True, False}]]

Options[QuantumTensorNetwork] = {"PrependInitial" -> True, "Computational" -> True, "ReturnIndices" -> False}

Options[QuantumTensorNetworkGraph] = Join[Options[QuantumTensorNetwork], Options[Graph]]

QuantumTensorNetworkGraph[qco_QuantumCircuitOperator, opts : OptionsPattern[]] := Enclose @ Block[{
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
        TensorNetworkGraphQ
    ]
]

(* TODO: refactor with above *)
QuantumTensorNetwork[qco_QuantumCircuitOperator, OptionsPattern[]] := Enclose @ Block[{
    circuit = qco["Sort"], width, min, ops, orders, arity, vertices, rules, tensors, indices
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
	rules = Catenate @ FoldPairList[
		{nprev, order} |-> Block[{output, input, n, prev, next, indexRules},
            {n, prev} = nprev;
            n += 1;
            next = prev;
			{output, input} = order;
            next[[ output - min + 1 ]] = n;
			next[[ Complement[input, output] - min + 1 ]] = None[n];
            indexRules = Superscript[prev[[# - min + 1]], #] -> Subscript[next[[# - min + 1]], #] & /@ input;
			{
                Replace[
                    Thread[DirectedEdge[prev[[ input - min + 1 ]], next[[ input - min + 1]], indexRules]],
                    {
                        DirectedEdge[None[_], ___] :> Nothing,
                        DirectedEdge[_, None[_], tag_] :> (tag /. None[i_] :> i),
                        DirectedEdge[_, _, tag_] :> tag
                    },
                    1
                ],
                {n, next}
            }
            
		],
		{1 - arity, Table[1 - arity, width]},
		Rest[orders]
	];
    If[ ! TrueQ[OptionValue["PrependInitial"]],
        vertices = Drop[vertices, arity];
        orders = Drop[orders, arity];
        ops = Drop[ops, arity];
    ];
    indices = Replace[
        MapThread[
            Join[OperatorApplied[Superscript, 2][#1] /@ Sort[#2[[1]]], OperatorApplied[Subscript, 2][#1] /@ Sort @ #2[[2]]] &,
            {vertices, orders}
        ],
        rules,
        {2}
    ];
    If[ TrueQ[OptionValue["ReturnIndices"]],
        Return[indices];
    ];
	tensors = If[#["MatrixQ"], #["Double"], #]["Tensor"] & /@ ops;
	ConfirmBy[
        TensorNetwork[
            tensors, indices,
            With[{free = Keys[Select[Counts[Catenate[indices]], # == 1 &]]},
                FindPermutation[
                    free, tensorNetworkIndexSort[free]
                ]
            ]
        ],
        TensorNetworkQ
    ]
]


QuantumCircuitHypergraph[qc_ ? QuantumCircuitOperatorQ, opts : OptionsPattern[]] := Enclose @ Block[{
    net = QuantumTensorNetworkGraph[qc["Flatten"], FilterRules[{opts}, Except[Options[Graph], Options[QuantumTensorNetworkGraph]]]], vs, indices, labels, edges
},
	Confirm[Needs["WolframInstitute`Hypergraph`" -> "H`"], "Hypergraph paclet is not installed."];
	vs = Developer`FromPackedArray @ VertexList[net];
	indices = TensorNetworkIndices[net];
	labels = AnnotationValue[{net, vs}, VertexLabels];
	edges = Replace[indices, Rule @@@ Reverse /@ EdgeTags[net], {2}];
	H`Hypergraph[Union @@ edges, edges, FilterRules[{opts}, Options[H`Hypergraph]], EdgeLabels -> Thread[edges -> labels]]
]


Options[TensorNetworkCompile] = Join[{"ReturnCircuit" -> False, "ReturnTensorNetwork" -> False, "Trace" -> True}, Options[QuantumTensorNetworkGraph], Options[TensorNetworkContract]]

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
    res = Confirm @ TensorNetworkContract[net, FilterRules[{opts}, Options[TensorNetworkContract]]];
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


TensorNetworkQuantumCircuit[tn_ ? TensorNetworkGraphQ] := QuantumCircuitOperator @ MapThread[{label, tensor, indices} |->
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
