Package["Wolfram`QuantumFramework`Experimental`"]

PackageImport["Wolfram`QuantumFramework`"]


PackageExport[DecomposedQuantumStateProbabilities]
PackageExport[QuantumBeamSearch]
PackageExport[QuantumDiagramProcess]

PackageExport[QuantumMPS]
PackageExport[QuantumMPO]
PackageExport[QuantumMPSApply]

PackageExport[FeynmanBacktracking]



DecomposedQuantumStateProbabilities[states : {{__QuantumState}..}] :=
	Normalize[Abs[Total[Flatten[Outer[Times, Sequence @@ #, 1], Length[#]] & /@ Map[#["Computational"]["StateVector"] &, states, {2}]]] ^ 2, Total]

DecomposedQuantumStateProbabilities[states : {(_ -> {__QuantumState})..}] :=
	Normalize[Abs[Total[states[[All, 1]] (Flatten[Outer[Times, Sequence @@ #, 1], Length[#]] & /@ Map[#["Computational"]["StateVector"] &, states[[All, 2]], {2}])]] ^ 2, Total]

BeamBranch[prob_ -> states_, op_] := With[{
	decompose = op["State"][
		QuantumTensorProduct[states[[op["InputOrder"]]]]
	]["DecomposeWithProbabilities"]
},
	prob #1 -> ReplacePart[states, Thread[op["InputOrder"] -> #2]] & @@@ decompose
]


Options[QuantumBeamSearch] = {"Width" -> 8, "Deterministic" -> False, "Shots" -> 1};

QuantumBeamSearch[states_List, ops_List, OptionsPattern[]] := Module[{
	width = OptionValue["Width"],
	random = !TrueQ[OptionValue["Deterministic"]],
	shots = OptionValue["Shots"]
},
	Normalize[#, Total] & @ Total @ Table[
		DecomposedQuantumStateProbabilities @ Fold[
			{beam, op} |-> SubsetMap[Normalize[#, Total] &, {All, 1}] @
				With[{candidates = Catenate[BeamBranch[#, op] & /@ beam]},
					If[ random,
						RandomSample[candidates[[All, 1]] -> candidates, UpTo[width]],
						TakeLargestBy[candidates, First, UpTo[width]]
					]
				],
			{1 -> states},
			N @ ops
		],
		shots
	]
]


DiagramProcess := DiagramProcess = ResourceFunction["https://www.wolframcloud.com/obj/murzin.nikolay/DeployedResources/Function/DiagramProcess"]

QuantumDiagramProcess[qco_QuantumCircuitOperator] := With[{
    ops = qco["Operators"], net = qco["TensorNetworkGraph", "PrependInitial" -> False], n = qco["GateCount"]
},
    With[{
        map = GroupBy[EdgeTags[net], #[[2]] &, #[[1, 1]] &],
        freeIndices = TensorNetworkFreeIndices[net]
    },
        DiagramProcess[
            Subsuperscript[
                With[{mat = ops[[#]]["Computational"]["Tensor"]}, Labeled[Part[mat, ##] &, ops[[#]]["Label"]]],
                Sequence @@ Reverse @ Replace[
                    MapAt[ReplaceAll[map], TakeDrop[HoldForm /@ AnnotationValue[{net, #}, "Index"], ops[[#]]["OutputQudits"]], 2],
                    With[{
                        outs = Alternatives @@ Cases[freeIndices, _Superscript],
                        ins = Alternatives @@ Cases[freeIndices, _Subscript]
                    },
                        {in : HoldForm[ins] :> Overscript[in, ˘], out : HoldForm[outs] :> Overscript[out, \[DownBreve]]}
                    ],
                    {2}
                ]
            ] & /@ Range[n]
        ]
    ]
]


Options[QuantumMPS] = {"Ordered" -> False, "Sides" -> True}

QuantumMPS[qs_QuantumState, m : _Integer | Infinity : Infinity, OptionsPattern[]] := Block[{
	decompose = If[VectorQ[#[[All, 1]], NumericQ], TakeLargestBy[#, First, UpTo[m]], #] & @ qs["DecomposeWithAmplitudes", qs["Dimensions"]],
	dimensions = qs["Dimensions"],
	proba, n, rowVector, colVector, matrices, result
},
	n = Length[decompose];
	proba = Keys[decompose];

	matrices = If[n > 1,
		colVector = QuantumOperator[{Table[{1}, n]}, {{0}, {}}, QuantumBasis[{n}, {}, "Label" -> TraditionalForm[Bra[{" "}]]]];
		rowVector = QuantumOperator[{proba}, {{}, {0}}, QuantumBasis[{}, {n}, "Label" -> TraditionalForm[Ket[{" "}]]]];
		matrices = MapIndexed[
			QuantumOperator[
				Transpose[
					ReplacePart[ConstantArray[Table[0, #1[[2]]], {n, n}], Thread[{#, #} & /@ Range[n] -> #1[[1]], List, 2]],
					2 <-> 3
				],
				If[#2[[1]] <= qs["OutputQudits"], {Prepend[#2, 0], {0}}, {{0}, Append[#2, 0]}],
				QuantumBasis[{n, #1[[2]]}, {n}]
			] &,
			Thread[{Transpose @ Map[#["Computational"]["StateVector"] &, Values[decompose], {2}], dimensions}]
		],

		colVector = Nothing;
		rowVector = QuantumOperator[{proba}, {{}, {}}, QuantumBasis[{}, {}, "Label" -> TraditionalForm[First[proba]]]];
		matrices = MapIndexed[
			QuantumOperator[{{#1[[1]]}}, {#2, {}}, QuantumBasis[{#1[[2]]}, {}]] &,
			Thread[{Transpose @ Map[#["Computational"]["StateVector"] &, Values[decompose], {2}], dimensions}]
		]
	];
	result = If[TrueQ[OptionValue["Sides"]],
		QuantumCircuitOperator[{colVector, Splice @ matrices, rowVector}],
		QuantumCircuitOperator[MapAt[rowVector, -1] @ MapAt[If[colVector === Nothing, #, #[colVector]] &, 1] @ matrices]
	];
	If[	TrueQ[OptionValue["Ordered"]],
		result = Reverse[MapIndexed["I"[#1] -> #2 -> qs["OutputQudits"] + #2 &, qs["InputDimensions"]]] /* result
	];
	QuantumCircuitOperator[result, "MPS"]
]

QuantumMPS[qo_QuantumOperator, m : _Integer | Infinity : Infinity, opts : OptionsPattern[]] :=
	With[{range = Range[Length[qo["InputOrder"]]]},
		{"Permutation"[qo["InputDimensions"], range] -> qo["InputOrder"] -> range + Length[qo["OutputOrder"]]}
	] /*
	QuantumMPS[qo["State"], m, "Ordered" -> False, opts] /*
	With[{range = Range[Length[qo["OutputOrder"]]]},
		{"Permutation"[qo["OutputDimensions"], range] -> range -> qo["OutputOrder"]}
	]


Options[QuantumMPO] = {"Ordered" -> True}

QuantumMPO[qo_QuantumOperator, m : _Integer | Infinity : Infinity, OptionsPattern[]] := Block[{
	mps = QuantumMPS[qo["ReverseInput"]["State"], m, "Sides" -> False],
	top, bot, split,
	result
},
	{top, bot} = TakeDrop[mps["Operators"], qo["OutputQudits"]];
	split = Length[top] - Length[bot];
	bot = #["Transpose", {{0, -1}, {-1, 0}}] & /@ Reverse[bot];
	result = QuantumCircuitOperator[
		If[ split > 0,
			Join[top[[;; split]], MapThread[Construct, {top[[split + 1 ;;]], bot}]],
			Join[bot[[;; - split]], MapThread[Construct, {top, bot[[- split + 1 ;;]]}]]
		] //
		MapAt[QuantumOperator["Cap"[#["OutputDimensions"][[1]]], {-1, 0}][#] &, -1] //
		MapAt[
			If[
				ContainsAll[#["OutputOrder"], {-1 ,0}],
				QuantumOperator["Curry"[Replace[{-1, 0}, #["OutputOrderDimensions"], {1}]], {-1, 0} -> {0}][#],
				QuantumOperator[#["State"], {#["OutputOrder"] /. -1 -> 0, #["InputOrder"]}]
			] &,
			{;; -2}
		] //
		MapAt[
			If[
				ContainsAll[#["InputOrder"], {-1 ,0}],
				#[QuantumOperator["Uncurry"[Replace[{-1, 0}, #["InputOrderDimensions"], {1}]], {0} -> {-1, 0}]],
				QuantumOperator[#["State"], {#["OutputOrder"], #["InputOrder"] /. -1 -> 0}]
			] &,
		 	{2 ;;}
		] //
		Map[QuantumOperator[
			#["State"],
			{#["OutputOrder"], Replace[#["InputOrder"], i_ /; i > qo["OutputQudits"] :> qo["Qudits"] - i + Max[split, 0] + 1, {1}]},
			"Label" -> None] &
		]
	];
	If[	TrueQ[OptionValue["Ordered"]],
		result =
			{"Permutation"[qo["InputDimensions"]] -> qo["InputOrder"] -> Reverse[Range[Length[qo["InputOrder"]]] + Max[split, 0]]} /*
			result /*
			{"Permutation"[qo["OutputDimensions"]] -> Range[Length[qo["OutputOrder"]]] -> qo["OutputOrder"]}
	];
	QuantumCircuitOperator[result, "MPO"]
]

(* An MPS here is a circuit (the QuantumMPS of a state, or a QuantumMPSApply result) whose
   elementary operators are site tensors, one per physical wire, chained by a bond carried on
   wire 0: a site reads the bond of the site before it on input wire 0 and passes a bond on to
   the next site on output wire 0. A bond of dimension 1 is no wire at all. A vector on wire 0
   may open the chain, a covector on wire 0 may close it, and 0-qudit scalars may scale it.
   mpsChain reads the chain, in circuit order, as {wires, tensors} with each tensor a rank-3
   array T[bondIn, physical, bondOut], the boundaries absorbed so the outer bonds have
   dimension 1. *)

mpsSiteQ[op_] := MatchQ[op["InputOrder"], {} | {0}] && MatchQ[Sort[op["OutputOrder"]], {_Integer ? Positive} | {0, _Integer ? Positive}]

(* "MatrixRepresentation" of the sorted operator has rows (bondOut, physical), bondOut slowest,
   and columns bondIn *)
mpsSiteTensor[op_] := With[{s = op["Sort"]},
	With[{
		bondOut = If[MemberQ[s["OutputOrder"], 0], First[s["OutputDimensions"]], 1],
		d = Last[s["OutputDimensions"]],
		bondIn = If[s["InputOrder"] === {0}, First[s["InputDimensions"]], 1]
	},
		Transpose[ArrayReshape[Normal[s["MatrixRepresentation"]], {bondOut, d, bondIn}], {3, 2, 1}]
	]
]

mpsElement[op_] := Which[
	mpsSiteQ[op], "Site" -> {First[DeleteCases[op["OutputOrder"], 0]], mpsSiteTensor[op]},
	op["OutputOrder"] === {0} && op["InputOrder"] === {}, "Left" -> Flatten[Normal[op["MatrixRepresentation"]]],
	op["OutputOrder"] === {} && op["InputOrder"] === {0}, "Right" -> Flatten[Normal[op["MatrixRepresentation"]]],
	op["OutputOrder"] === {} && op["InputOrder"] === {}, "Scalar" -> First[Flatten[Normal[op["MatrixRepresentation"]]]],
	True, "Other" -> op
]

mpsChain[mps_QuantumCircuitOperator] := Enclose @ Block[{elements, wires, tensors, left, right},
	elements = GroupBy[mpsElement /@ mps["Flatten"]["Operators"], First -> Last];
	ConfirmAssert[! KeyExistsQ[elements, "Other"],
		"QuantumMPSApply expects the MPS of a state: site tensors on positive wires, chained by a bond on wire 0"];
	{wires, tensors} = Transpose @ ConfirmBy[Lookup[elements, "Site", {}], # =!= {} &, "The MPS has no site tensors"];
	left = Lookup[elements, "Left", {}];
	right = Lookup[elements, "Right", {}];
	ConfirmAssert[DuplicateFreeQ[wires] && Length[left] <= 1 && Length[right] <= 1, "The MPS chain is malformed"];
	ConfirmAssert[left === {} || Length[First[left]] == First[Dimensions[First[tensors]]], "The opening bond vector does not fit the first site"];
	ConfirmAssert[right === {} || Length[First[right]] == Last[Dimensions[Last[tensors]]], "The closing bond covector does not fit the last site"];
	tensors = MapAt[(Times @@ Lookup[elements, "Scalar", {}]) If[left === {}, #, {First[left] . #}] &, tensors, 1];
	tensors = MapAt[If[right === {}, #, Map[List, # . First[right], {2}]] &, tensors, -1];
	ConfirmAssert[
		First[Dimensions[First[tensors]]] == 1 && Last[Dimensions[Last[tensors]]] == 1 &&
			And @@ MapThread[Last[Dimensions[#1]] == First[Dimensions[#2]] &, {Most[tensors], Rest[tensors]}],
		"The MPS bond dimensions do not chain"
	];
	{wires, tensors}
]

(* Factor a matrix as m == left . right through its rank r, the new bond dimension. Approximate
   numbers: the singular value decomposition, left = u and right = sigma . v^dagger, keeping every
   nonzero singular value. Exact numbers whose radicals are all of rational numbers (1/Sqrt[2],
   Exp[I Pi/4], ...): the exact rank factorization by row reduction, left = the pivot columns of m
   and right = the nonzero rows of RowReduce[m]; an exact singular value decomposition would
   build the algebraic singular values of m, which swell the amplitudes and can exhaust the
   kernel. Nested radicals or Root objects (QuantumMPS of a generic exact state has them) make
   row reduction slow and the contracted amplitudes huge, so such a matrix is split numerically,
   with a message. A zero matrix keeps a bond of dimension 1. *)
QuantumMPSApply::numsplit = "The exact site tensors hold nested radicals or Root objects; the updated sites are split numerically.";

mpsSimpleExactQ[m_] := Precision[m] === Infinity && FreeQ[m, Power[Except[_Integer | _Rational], _Rational] | _Root]

mpsFactor[m_] /; Precision[m] === Infinity && ! mpsSimpleExactQ[m] := (Message[QuantumMPSApply::numsplit]; mpsFactor[N[m]])

mpsFactor[m_] /; mpsSimpleExactQ[m] := With[{rows = DeleteCases[RowReduce[m], {0 ..}]},
	If[ rows === {},
		{ConstantArray[0, {Length[m], 1}], ConstantArray[0, {1, Last[Dimensions[m]]}]},
		{m[[All, First[FirstPosition[#, Except[0], {1}, {1}, Heads -> False]] & /@ rows]], rows}
	]
]

mpsFactor[m_] := With[{usv = SingularValueDecomposition[m, Max[1, Length[SingularValueList[m]]]]},
	{usv[[1]], usv[[2]] . ConjugateTranspose[usv[[3]]]}
]

(* Split a block B[bondIn, d_1, ..., d_k, bondOut] back into k site tensors from the left, one
   factorization per cut. Nothing is truncated, so the block is reproduced exactly. A new bond
   has the rank of the block across its cut, which equals the Schmidt rank of the state only
   when the sites outside the block are isometries; QuantumMPS builds a sum of product states,
   which are not, so the bonds can exceed the Schmidt ranks. *)
mpsSplitStep[{sites_, rest_}] := With[{dims = Dimensions[rest]},
	With[{lr = mpsFactor[ArrayReshape[rest, {Times @@ Take[dims, 2], Times @@ Drop[dims, 2]}]]},
		{
			Append[sites, ArrayReshape[lr[[1]], {dims[[1]], dims[[2]], Last[Dimensions[lr[[1]]]]}]],
			ArrayReshape[lr[[2]], Prepend[Drop[dims, 2], First[Dimensions[lr[[2]]]]]]
		}
	]
]

mpsSplit[block_] := Append @@ Nest[mpsSplitStep, {{}, block}, ArrayDepth[block] - 3]

mpsSiteOperator[wire_, t_] := With[{bondIn = Dimensions[t][[1]], d = Dimensions[t][[2]], bondOut = Dimensions[t][[3]]},
	QuantumOperator[
		QuantumState[
			Flatten[Transpose[t, {3, 2, 1}]],
			QuantumBasis[If[bondOut > 1, {bondOut, d}, {d}], If[bondIn > 1, {bondIn}, {}]]
		],
		{If[bondOut > 1, {0, wire}, {wire}], If[bondIn > 1, {0}, {}]}
	]
]

(* Apply op to the sites on its wires: contract the run of sites from the first to the last of
   them into one block, act with op's matrix (sorted-wire order, so a reversed order such as
   CNOT on {2, 1} keeps its control on wire 2) on those physical legs, and split the block back
   into sites. Sites strictly between non-adjacent wires are part of the block, so the result
   is the exact MPS of op applied to the state. A wire of op that carries no site joins the
   end of the chain in |0>, the state a circuit wire starts in. *)
QuantumMPSApply[mps_QuantumCircuitOperator, op_QuantumOperator] /; op["InputOrder"] == op["OutputOrder"] := Enclose @ Block[{
	wires, tensors, gateWires, newWires, sites, first, last, block, legs, order, moved, gateDim
},
	{wires, tensors} = Confirm @ mpsChain[mps];
	gateWires = Sort[op["InputOrder"]];
	newWires = Complement[gateWires, wires];
	tensors = Join[tensors, {List /@ UnitVector[#, 1]} & /@ Lookup[AssociationThread[gateWires, op["Sort"]["InputDimensions"]], newWires]];
	wires = Join[wires, newWires];
	sites = Flatten[Position[wires, #, {1}, 1] & /@ gateWires];
	ConfirmAssert[op["Sort"]["InputDimensions"] === op["Sort"]["OutputDimensions"] === (Dimensions[#][[2]] & /@ tensors[[sites]]),
		"The operator dimensions do not match the MPS sites"];
	{first, last} = MinMax[sites];
	block = Fold[Dot, tensors[[first ;; last]]];
	(* move the gate's physical legs to the front, in sorted-wire order; leg 1 is the opening bond *)
	legs = sites - first + 2;
	order = Join[legs, Complement[Range[ArrayDepth[block]], legs]];
	moved = Transpose[block, Ordering[order]];
	gateDim = Times @@ Take[Dimensions[moved], Length[legs]];
	moved = ArrayReshape[
		Normal[op["MatrixRepresentation"]] . ArrayReshape[moved, {gateDim, (Times @@ Dimensions[moved]) / gateDim}],
		Dimensions[moved]
	];
	ConfirmAssert[first == last || ArrayQ[moved, _, NumericQ], "Splitting the updated sites needs numeric tensors"];
	QuantumCircuitOperator[
		MapThread[mpsSiteOperator, {wires, Join[tensors[[;; first - 1]], mpsSplit[Transpose[moved, order]], tensors[[last + 1 ;;]]]}],
		mps["Label"] /* op["Label"]
	]
]


$FeynmanBacktrackingCache = <||>;

iFeynmanBacktracking[ops : {___QuantumOperator}, state : {___Integer ? Positive} | Automatic, width_Integer, implicitInputOrder_List : {}] := Enclose @ Block[{
	key, amplitude, op, outputPos, inputPos, padState, inputs
},
	key = {Length[ops], state};
	If[ key[[1]] == 0,
		amplitude = Boole @ MatchQ[Delete[state, List /@ implicitInputOrder], {1 ...}];
		(*$FeynmanBacktrackingCache[key] = {amplitude, 1};*)
		Return[amplitude]
	];
	amplitude = Lookup[$FeynmanBacktrackingCache, Key[key]];
	If[! MissingQ[amplitude], $FeynmanBacktrackingCache[[Key[key], 2]]++; Return[amplitude[[1]]]];

	op = ops[[-1]];
	{outputPos, inputPos} = op["Order"];
	padState = PadRight[Replace[state, Automatic :> ConstantArray[1, Length[width]]], width, 1];
	inputs = Extract[op["Tensor"], padState[[outputPos]]];
	amplitude = If[ ArrayQ[inputs],
		inputs = SparseArray[inputs];
		Total @ MapThread[#2 iFeynmanBacktracking[Most[ops], ReplacePart[padState, Thread[inputPos -> #1]], width, implicitInputOrder] &, {inputs["ExplicitPositions"], inputs["ExplicitValues"]}]
		,
		(* no input - operator is a state *)
		If[inputs == 0, 0, inputs iFeynmanBacktracking[Most[ops], padState, width, Join[implicitInputOrder, outputPos]]]
	];
	$FeynmanBacktrackingCache[key] = {amplitude, 1};
	amplitude
]

FeynmanBacktracking[circuit_QuantumCircuitOperator, state : {___Integer ? Positive} : Automatic] := (
	$FeynmanBacktrackingCache = <||>;
	iFeynmanBacktracking[circuit["Flatten"]["NormalOperators"], state, circuit["Width"]]
)

