(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicRelations"]

PackageExport["BosonicNormalOrder"]

PackageExport["BosonicAntinormalOrder"]

PackageExport["BosonicBCHTerms"]

PackageExport["BosonicZassenhausTerms"]

PackageExport["BosonicVEV"]

PackageExport["BosonicMatrixElement"]


BosonicRelations::usage =
"\!\(\*RowBox[{\"BosonicRelations\", \"[\", RowBox[{StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Returns the list of bosonic commutation relations for the operators in vars.";

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


grobnerNormalOrder[expr_, vars_List, scalars_List, direction_String : "Normal"] :=
Block[{
    orderedVars = OrderVariables[vars, direction],
    rels,
    alg
},
    rels = BosonicRelations[vars];
    If[ Length[scalars] > 0,
        alg = NonCommutativeAlgebra["ScalarVariables" -> scalars];
        NonCommutativePolynomialReduce[expr, rels, orderedVars, alg][[2]],
        NonCommutativePolynomialReduce[expr, rels, NonCommutativeAlgebra[<|"Generators"-> orderedVars|>]][[2]]
    ]
]


BosonicNormalOrder::usage =
"\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings expr into normal order using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the reduction method: \"GrobnerBasis\" (default) or \"Blasiak\".\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) Treats syms as commuting scalars during reduction.";

BosonicNormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicNormalOrder] = {Method -> "GrobnerBasis", "Scalars" -> {}}

(* With no field variables there is nothing to reduce: a c-number is already in
   normal order, so return it unchanged. *)
BosonicNormalOrder[expr_, {}, OptionsPattern[]] := expr

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


(* a -> a\[Dagger], a\[Dagger] -> -a is an automorphism ([a\[Dagger], -a] = 1) carrying a normal-ordered
   monomial to an anti-ordered one, so conjugating by it turns any normal-ordering
   method into an anti-normal one.  "Inverse" maps the result back. *)
ladderSwap[expr_, vars_List, direction_String : "Forward"] :=
Block[{
    swap = Catenate @ Map[
        If[ direction === "Inverse",
            {SuperDagger[#] -> #, # -> -SuperDagger[#]} &,
            {SuperDagger[#] -> -#, # -> SuperDagger[#]} &
        ],
        Select[vars, FreeQ[#, SuperDagger] &]
    ]
},
    LiftNCScalars[expr /. swap, vars]
]


BosonicAntinormalOrder::usage =
"\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings the polynomial expr into anti-normal order, with every annihilation operator to the left of every creation operator, using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the reduction method: \"GrobnerBasis\" (default) or \"Blasiak\".\n\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) Treats syms as commuting scalars during reduction.";

BosonicAntinormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicAntinormalOrder] = Options[BosonicNormalOrder]

(* With no field variables there is nothing to reorder: a c-number is already in
   anti-normal order, so return it unchanged. *)
BosonicAntinormalOrder[expr_, {}, OptionsPattern[]] := expr

BosonicAntinormalOrder[expr_, vars_List, opts : OptionsPattern[]] :=
With[{
    method = OptionValue[Method],
    scalars = OptionValue["Scalars"]
},
    If[ MatchQ[expr, _Plus],
        Return[Total[ParallelMap[
            BosonicAntinormalOrder[#, vars, "Scalars" -> scalars, Method -> method] &,
            List @@ expr
        ]]]
    ];
    Replace[method, {
        "GrobnerBasis" :> grobnerNormalOrder[expr, vars, scalars, "Antinormal"],
        (* normal-order the swapped monomial, then swap back *)
        "Blasiak" :> ladderSwap[
            MultiModeBlasiakOrder[ladderSwap[expr, vars], vars, scalars],
            vars,
            "Inverse"
        ],
        _ :> (Message[BosonicAntinormalOrder::unknownMethod, method]; $Failed)
    }]
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
"\!\(\*RowBox[{\"BosonicBCHTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"]\"}]\) Computes the nth-order term of the BCH series log(\!\(\*SuperscriptBox[\(e\), \(op1\)]\) \!\(\*SuperscriptBox[\(e\), \(op2\)]\)\[Ellipsis]), reduced using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicBCHTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"ncVars\", \"TI\"]}], \"]\"}]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicBCHTerms[ops_List, n_Integer] := BosonicBCHTerms[ops, n, ExtractNCVars[ops]]

BosonicBCHTerms[ops_List, n_Integer, ncVars_List] :=
Block[{bchTerm},
    bchTerm = ResourceFunction["BakerCampbellHausdorffTerms"][ops, n];
    bosonicReduce[NonCommutativeExpand[bchTerm], ncVars]
]


BosonicZassenhausTerms::usage =
"\!\(\*RowBox[{\"BosonicZassenhausTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"]\"}]\) Computes the nth term of the Zassenhaus formula and reduces the expression using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicZassenhausTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"ncVars\", \"TI\"]}], \"]\"}]\) Uses the explicitly supplied list of non-commutative variables ncVars.";

BosonicZassenhausTerms[ops_List, n_Integer] := BosonicZassenhausTerms[ops, n, ExtractNCVars[ops]]

BosonicZassenhausTerms[ops_List, n_Integer, ncVars_List] :=
Block[{zassTerm},
    zassTerm = ResourceFunction["ZassenhausTerms"][ops, n];
    bosonicReduce[NonCommutativeExpand[zassTerm], ncVars]
]


BosonicVEV::usage =
"\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Computes the vacuum expectation value \[LeftAngleBracket]0\[VerticalSeparator]expr\[VerticalSeparator]0\[RightAngleBracket] by normal-ordering expr and extracting the scalar (c-number) part.\n\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"]}], \"]\"}]\) Auto-detects non-commutative variables from expr.\n\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Passes method options to BosonicNormalOrder.";

Options[BosonicVEV] = Options[BosonicNormalOrder];

BosonicVEV[expr_, opts : OptionsPattern[]] :=
    Block[{vars = ExtractNCVars[{expr}]},
        (* A c-number carries no field variables; it is its own vacuum expectation.
           Short-circuit before re-dispatching: BosonicVEV[expr, {}] would re-match
           this same rule ({} is absorbed as empty options), so calling it here would
           recurse without bound. *)
        If[vars === {}, expr, BosonicVEV[expr, vars, opts]]
    ]

BosonicVEV[expr_, vars_List, opts : OptionsPattern[]] :=
    Block[{normalOrdered, terms},
        normalOrdered = BosonicNormalOrder[expr, vars, opts];
        terms = If[Head[normalOrdered] === Plus, List @@ normalOrdered, {normalOrdered}];
        Total @ Select[terms, FreeQ[#, Alternatives @@ vars] &]
    ]


BosonicMatrixElement::usage =
"\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket] in closed form via associated Laguerre polynomials.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket] in closed form (zero when m+n is odd).";

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
