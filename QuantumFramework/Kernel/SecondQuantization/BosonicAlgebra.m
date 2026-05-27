(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicRelations"]

PackageExport["BosonicNormalOrder"]

PackageExport["BosonicBCH"]

PackageExport["BosonicZassenhaus"]

PackageExport["BosonicBCHExact"]


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
    NonCommutativePolynomialReduce[poly, rels, OrderVariables[ncVars], alg][[2]]
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


$nonuls  = {0. -> 0, 0. I -> 0, Complex[0.,0.] -> 0, Complex[x_, 0.]->x,
          Complex[0.,y_] -> I y};
          
BosonicBCHExact::usage = "\!\(BosonicBCHExact[expr, vars]\) Exactly evaluates a product of exponentials or an exponential conjugation when the algebra closes, using BCH and adjoint action identities.
\!\(BosonicBCHExact[expr, vars, \"Scalars\" -> syms]\) Treats syms as commuting scalars.";

Options[BosonicBCHExact] = {"Scalars" -> {}};

BosonicBCHExact[expr_] :=
    Block[{vars, scalars},
        vars = ExtractNCVars[{expr}];
        scalars = DeleteDuplicates @ Cases[{expr}, s_Symbol /; !FormalSymbolQ[
            s] && !NumericQ[s], Infinity];
        BosonicBCHExact[expr, vars, "Scalars" -> scalars]
    ]
        

BosonicBCHExact[Exp[b1_] ** Exp[b2_], vars_List, opts : OptionsPattern[
    ]] :=
    Block[{scalars = OptionValue["Scalars"], cmt},
        cmt = Simplify @ BosonicNormalOrder[Commutator[b1, b2], vars,
             "Scalars" -> scalars];
        Which[
            cmt === 0,
                Exp[b1 + b2]
            ,
            FreeQ[cmt /. $nonuls, Alternatives @@ vars],
                Exp[cmt / 2] * Exp[b1 + b2]
            ,
            True,
                Exp[b1] ** Exp[b2]
        ]
    ]

BosonicBCHExact[Exp[b1_] ** ops__ ** Exp[b3_], vars_List, opts : OptionsPattern[]] /;
    b3 === -b1 && Length[{ops}] > 1 :=
    NonCommutativeMultiply @@ (BosonicBCHExact[Exp[b1] ** # ** Exp[b3], vars, opts] & /@ {ops})

BosonicBCHExact[Exp[b1_] ** op_ ** Exp[b3_], vars_List, opts : OptionsPattern[
    ]] /; b3 === -b1 :=
    Block[{scalars = OptionValue["Scalars"], cmt, cmt2, lam},
        cmt = Simplify @ BosonicNormalOrder[Commutator[b1, op], vars,
             "Scalars" -> scalars];
        cmt2 = Simplify @ BosonicNormalOrder[Commutator[b1, cmt], vars,
             "Scalars" -> scalars];
        Which[
            cmt === 0,
                op
            ,
            FreeQ[cmt /. $nonuls, Alternatives @@ vars],
                op + cmt
            ,
            (cmt2 /. $nonuls) === 0,
                op + cmt
            ,
            FreeQ[Simplify[cmt2 / op] /. $nonuls, Alternatives @@ vars
                ],
                Block[{},
                    lam = Simplify[cmt2 / op] /. $nonuls;
                    Simplify[Cos[I Sqrt[lam]] * op + Sin[I Sqrt[lam]]
                         / (I Sqrt[lam]) * cmt] /. $nonuls
                ]
            ,
            True,
                Exp[b1] ** op ** Exp[b3]
        ]
    ]

BosonicBCHExact[op_ ** Exp[b3_], vars_List, opts : OptionsPattern[]] :=
    Block[{scalars = OptionValue["Scalars"], cmt, cmt2},
        cmt = Simplify @ BosonicNormalOrder[Commutator[b3, op], vars,
             "Scalars" -> scalars];
        cmt2 = Simplify @ BosonicNormalOrder[Commutator[b3, cmt], vars,
             "Scalars" -> scalars];
        Which[
            cmt === 0,
                Exp[b3] ** op
            ,
            (cmt2 /. $nonuls) === 0,
                Exp[b3] ** op - Exp[b3] ** cmt
            ,
            True,
                op ** Exp[b3]
        ]
    ]
    
BosonicBCHExact[Exp[b1_] ** Exp[b2_] ** Exp[b3_], vars_List, opts : OptionsPattern[]] /;
    b3 === -b1 :=
	Block[{scalars = OptionValue["Scalars"], newb2},
    newb2 = BosonicBCHExact[Exp[b1] ** b2 ** Exp[b3], vars, opts];
    If[FreeQ[newb2, BosonicBCHExact],
        Exp[newb2],
        Exp[b1] ** Exp[b2] ** Exp[b3]
    ]
]
