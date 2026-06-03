(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicRelations"]

PackageExport["BosonicNormalOrder"]

PackageExport["BosonicBCHTerms"]

PackageExport["BosonicZassenhausTerms"]

PackageExport["BosonicBCHExact"]

PackageExport["BosonVEV"]

PackageExport["BosonicMatrixElement"]


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


BosonicBCHTerms::usage =
"\!\(BosonicBCHTerms[ops, n]\) Computes the nth-order term of the BCH series log(\!\(\*SuperscriptBox[\(e\), \(op1\)]\) \!\(\*SuperscriptBox[\(e\), \(op2\)]\)\[Ellipsis]), reduced using bosonic commutation relations.
\!\(BosonicBCHTerms[ops, n, ncVars]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicBCHTerms[ops_List, n_Integer] := BosonicBCH[ops, n, ExtractNCVars[ops]]

BosonicBCHTerms[ops_List, n_Integer, ncVars_List] :=
Block[{bchTerm},
    bchTerm = ResourceFunction["BakerCampbellHausdorffTerms"][ops, n];
    bosonicReduce[NonCommutativeExpand[bchTerm], ncVars]
]


BosonicZassenhausTerms::usage =
"\!\(BosonicZassenhausTerms[ops, n]\) Computes the nth term of the Zassenhaus formula and reduces the expression using bosonic commutation relations.
\!\(BosonicZassenhausTerms[ops, n, ncVars]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicZassenhausTerms[ops_List, n_Integer] := BosonicZassenhaus[ops, n, ExtractNCVars[ops]]

BosonicZassenhausTerms[ops_List, n_Integer, ncVars_List] :=
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


BosonVEV::usage =
"\!\(BosonVEV[expr, vars]\) Computes the vacuum expectation value \[LeftAngleBracket]0\[VerticalSeparator]expr\[VerticalSeparator]0\[RightAngleBracket] by normal-ordering expr and extracting the scalar (c-number) part.
\!\(BosonVEV[expr]\) Auto-detects non-commutative variables from expr.
\!\(BosonVEV[expr, vars, Method -> m]\) Passes method options to BosonicNormalOrder.";

Options[BosonVEV] = Options[BosonicNormalOrder];

BosonVEV[expr_, opts : OptionsPattern[]] :=
    Block[{vars = ExtractNCVars[{expr}]},
        VacuumExpectation[expr, vars, opts]
    ]

BosonVEV[expr_, vars_List, opts : OptionsPattern[]] :=
    Block[{normalOrdered, terms},
        normalOrdered = BosonicNormalOrder[expr, vars, opts];
        terms = If[Head[normalOrdered] === Plus, List @@ normalOrdered, {normalOrdered}];
        Total @ Select[terms, FreeQ[#, Alternatives @@ vars] &]
    ]


BosonicMatrixElement::usage =
"\!\(BosonicMatrixElement[{m, n}, DisplacementOperator[\[Alpha]]]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket] in closed form via associated Laguerre polynomials.
\!\(BosonicMatrixElement[{m, n}, SqueezeOperator[\[Xi]]]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket] in closed form (zero when m+n is odd).";

SetAttributes[BosonicMatrixElement, HoldRest]

BosonicMatrixElement[{m_Integer?NonNegative, n_Integer?NonNegative}, DisplacementOperator[\[Alpha]_]] :=
    With[{a2 = Abs[\[Alpha]]^2},
        If[m >= n,
            E^(-a2/2) Sqrt[n!/m!] \[Alpha]^(m - n) LaguerreL[n, m - n, a2],
            E^(-a2/2) Sqrt[m!/n!] (-Conjugate[\[Alpha]])^(n - m) LaguerreL[m, n - m, a2]
        ]
    ]


BosonicMatrixElement[{m_Integer?NonNegative, n_Integer?NonNegative}, SqueezeOperator[\[Xi]_]] :=
    Module[{r = Abs[\[Xi]], \[Phi] = Arg[\[Xi]], \[Tau], s},
        If[OddQ[m + n], Return[0]];
        \[Tau] = Exp[I \[Phi]] Tanh[r];
        If[m >= n,
            s = (m - n)/2;
            (-\[Tau]/2)^s / s! Sqrt[m!/n!] Sech[r]^(n + 1/2) *
                Hypergeometric2F1[-n/2, (1 - n)/2, s + 1, -Sinh[r]^2],
            s = (n - m)/2;
            (Conjugate[\[Tau]]/2)^s / s! Sqrt[n!/m!] Sech[r]^(m + 1/2) *
                Hypergeometric2F1[-m/2, (1 - m)/2, s + 1, -Sinh[r]^2]
        ]
    ]
