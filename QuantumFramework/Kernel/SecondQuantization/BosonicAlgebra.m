(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicRelations"]

PackageExport["BosonicNormalOrder"]

PackageExport["BosonicBCH"]

PackageExport["BosonicZassenhaus"]


BosonicRelations::usage =
"\!\(BosonicRelations[vars]\) Returns the list of bosonic commutation relations for the operators in vars.";

BosonicRelations[vars_List] :=
Block[{pairs},
    pairs = Subsets[vars, {2}];
    Map[
        If[ MatchQ[#, {x_, SuperDagger[x_]} | {SuperDagger[x_], x_}],
            Commutator[#[[1]], #[[2]]] - 1,
            Commutator[#[[1]], #[[2]]]
        ] &,
        pairs
    ]
]


grobnerNormalOrder[expr_, vars_List, scalars_List] :=
Block[{
    orderedVars = OrderVariables[vars],
    rels,
    alg
},
    rels = BosonicRelations[vars];
    If[ Length[scalars] > 0,
        alg = NonCommutativeAlgebra["ScalarVariables" -> scalars];
        NonCommutativePolynomialReduce[expr, rels, orderedVars, alg][[2]],
        NonCommutativePolynomialReduce[expr, rels, orderedVars][[2]]
    ]
]


BosonicNormalOrder::usage =
"\!\(BosonicNormalOrder[expr, vars]\) Brings expr into normal order using bosonic commutation relations.
\!\(BosonicNormalOrder[expr, vars, Method -> m]\) Specifies the reduction method: \"GrobnerBasis\" (default) or \"Blasiak\".
\!\(BosonicNormalOrder[expr, vars, \"Scalars\" -> syms]\) Treats syms as commuting scalars during reduction.";

BosonicNormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicNormalOrder] = {Method -> "GrobnerBasis", "Scalars" -> {}}

BosonicNormalOrder[expr_, vars_List, opts : OptionsPattern[]] :=
Block[{
    method = OptionValue[Method],
    scalars = OptionValue["Scalars"]
},
    If[ MatchQ[expr, _Plus],
        Return[Total[ParallelMap[
            BosonicNormalOrder[#, vars, "Scalars" -> scalars, Method -> method] &,
            List @@ expr
        ]]]
    ];
    If[ method === "GrobnerBasis", Return[grobnerNormalOrder[expr, vars, scalars]] ];
    If[ method === "Blasiak",      Return[MultiModeBlasiakOrder[expr, vars, scalars]] ];
    Message[BosonicNormalOrder::unknownMethod, method];
    $Failed
]


bosonicReduce[poly_, ncVars_List] :=
Block[{rels, scalars, alg},
    rels = BosonicRelations[ncVars];
    scalars = DeleteDuplicates @ Cases[
        {poly},
        s_Symbol /; !FormalSymbolQ[s] && !MemberQ[ncVars, s] && !NumericQ[s],
        Infinity
    ];
    alg = NonCommutativeAlgebra["ScalarVariables" -> scalars];
    NonCommutativePolynomialReduce[poly, rels, ncVars, alg][[2]]
]


BosonicBCH::usage =
"\!\(BosonicBCH[ops, n]\) Computes the BCH series log(e\!\(\*SuperscriptBox[\()\), \(X1\)]\) e\!\(\*SuperscriptBox[\()\), \(X2\)]\)\[Ellipsis]) truncated at order n and reduces the result using bosonic commutation relations.
\!\(BosonicBCH[ops, n, ncVars]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicBCH[ops_List, n_Integer] := BosonicBCH[ops, n, ExtractNCVars[ops]]

BosonicBCH[ops_List, n_Integer, ncVars_List] :=
Block[{bchTerm},
    bchTerm = ResourceFunction["BakerCampbellHausdorffTerms"][ops, n];
    bosonicReduce[NonCommutativeExpand[bchTerm], ncVars]
]


BosonicZassenhaus::usage =
"\!\(BosonicZassenhaus[ops, n]\) Computes the Zassenhaus factorisation e\!\(\*SuperscriptBox[\()\), \(X1+X2+\[Ellipsis]\)]\) = e\!\(\*SuperscriptBox[\()\), \(X1\)]\) e\!\(\*SuperscriptBox[\()\), \(X2\)]\)\[Ellipsis] truncated at order n and reduces each factor using bosonic commutation relations.
\!\(BosonicZassenhaus[ops, n, ncVars]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicZassenhaus[ops_List, n_Integer] := BosonicZassenhaus[ops, n, ExtractNCVars[ops]]

BosonicZassenhaus[ops_List, n_Integer, ncVars_List] :=
Block[{zassTerm},
    zassTerm = ResourceFunction["ZassenhausTerms"][ops, n];
    bosonicReduce[NonCommutativeExpand[zassTerm], ncVars]
]
