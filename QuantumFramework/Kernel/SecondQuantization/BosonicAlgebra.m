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
    pairs = Subsets[LexicographicSort@vars, {2}];
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
"\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings expr into normal order using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the reduction method: Automatic (default), \"GrobnerBasis\", or \"Blasiak\".\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) Treats syms as commuting scalars during reduction.";

BosonicNormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicNormalOrder] = {Method -> Automatic, "Scalars" -> {}}

(* Ordering is linear: a sum maps over its summands, and the map stays on
   this kernel (parallel sub-kernels do not load the paclet).  A summand the
   monomial reducer rejects fails the whole sum rather than leaving a
   $Failed inside a Plus. *)
linearOverSums[f_][expr_] :=
    If[ MatchQ[expr, _Plus],
        Replace[Map[f, expr], p_ /; ! FreeQ[p, $Failed] :> $Failed],
        f[expr]
    ]

(* With no field variables there is nothing to reduce: a c-number is already in
   normal order, so return it unchanged. *)
BosonicNormalOrder[expr_, {}, OptionsPattern[]] := expr

BosonicNormalOrder[expr_, vars_List, opts : OptionsPattern[]] :=
With[{
    method = OptionValue[Method],
    scalars = OptionValue["Scalars"]
},
    Replace[method, {
        (* NonCommutativePolynomialReduce reduces a whole Plus itself. *)
        Automatic | "GrobnerBasis" :> grobnerNormalOrder[expr, vars, scalars],
        (* MultiModeBlasiakOrder takes one monomial at a time. *)
        "Blasiak" :> linearOverSums[MultiModeBlasiakOrder[#, vars, scalars] &][expr],
        _ :> (Message[BosonicNormalOrder::unknownMethod, method]; $Failed)
    }]
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
"\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings the polynomial expr into anti-normal order, with every annihilation operator to the left of every creation operator, using bosonic commutation relations.\n\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the reduction method: Automatic (default), \"GrobnerBasis\", or \"Blasiak\".\n\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) Treats syms as commuting scalars during reduction.";

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
    Replace[method, {
        (* NonCommutativePolynomialReduce reduces a whole Plus itself. *)
        Automatic | "GrobnerBasis" :> grobnerNormalOrder[expr, vars, scalars, "Antinormal"],
        (* Normal-order each swapped monomial, then swap back. *)
        "Blasiak" :> linearOverSums[
            ladderSwap[
                MultiModeBlasiakOrder[ladderSwap[#, vars], vars, scalars],
                vars,
                "Inverse"
            ] &
        ][expr],
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

BosonicVEV::novars =
"No field variables were recognized in `1`. Ladder operators must be Formal symbols (see FieldVariables), or be given explicitly as BosonicVEV[expr, vars]."

Options[BosonicVEV] = Options[BosonicNormalOrder];

BosonicVEV[expr_, opts : OptionsPattern[]] :=
    With[{vars = ExtractNCVars[{expr}]},
        (* Non-formal symbols are scalars, so with no field variables the only
           operator structure left is a ** joining c-numbers; SuperDagger instead
           means ladder operators were meant but not written as Formal symbols.
           Short-circuit before re-dispatching: BosonicVEV[expr, {}] would re-match
           this same rule ({} is absorbed as empty options), so calling it here would
           recurse without bound. *)
        With[{value = Which[
            vars =!= {},
                BosonicVEV[expr, vars, opts],
            FreeQ[expr, SuperDagger],
                vevCollapse[expr],
            True,
                Message[BosonicVEV::novars, expr]; $Failed
        ]},
            (* An unevaluated inner call means there is nothing to read off; let
               this call stay unevaluated too rather than echo the vars back. *)
            value /; Head[value] =!= BosonicVEV
        ]
    ]

vevVars[vars_List] :=
    DeleteDuplicates @
        Catenate[{#, SuperDagger[#]} & /@ DeleteDuplicates[vars /. SuperDagger[v_] :> v]]

vevCollapse[expr_] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    NonCommutativeMultiply -> Times
}

vevPolynomialQ[expr_, vars_] :=
    With[{syms = Table[Unique["vev"], Length[vars]]},
        PolynomialQ[vevCollapse[expr] /. Thread[vars -> syms], syms]
    ]

(* The condition sits on the left-hand side so a non-polynomial never reaches
   BosonicNormalOrder: the reducer would both fail and leak its own messages. *)
BosonicVEV[expr_, vars_List, opts : OptionsPattern[]] /;
        vevPolynomialQ[expr, vevVars[vars]] :=
    With[{all = vevVars[vars]},
        vevCollapse[BosonicNormalOrder[expr, all, opts] /. $nonuls] /.
            (Alternatives @@ all -> 0)
    ]


BosonicMatrixElement::usage =
"\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket] in closed form via associated Laguerre polynomials.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket] in closed form (zero when m+n is odd).\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"poly\", \"TI\"]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]poly\[VerticalSeparator]n\[RightAngleBracket] for a polynomial in the field variables of a single mode, by normal ordering and the closed form for \[LeftAngleBracket]m\[VerticalSeparator]\!\(\*SuperscriptBox[\"a\", RowBox[{\"\[Dagger]\", \"p\"}]]\)\!\(\*SuperscriptBox[\"a\", \"q\"]\)\[VerticalSeparator]n\[RightAngleBracket]. The Fock indices may be symbolic.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"\[Beta]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"poly\", \"TI\"], \",\", \"\\\"Basis\\\"->\\\"Coherent\\\"\"}], \"]\"}]\) Returns \[LeftAngleBracket]\[Alpha]\[VerticalSeparator]poly\[VerticalSeparator]\[Beta]\[RightAngleBracket] between coherent states of amplitude \[Alpha] and \[Beta], by normal ordering and \[LeftAngleBracket]\[Alpha]\[VerticalSeparator]\!\(\*SuperscriptBox[\"a\", RowBox[{\"\[Dagger]\", \"p\"}]]\)\!\(\*SuperscriptBox[\"a\", \"q\"]\)\[VerticalSeparator]\[Beta]\[RightAngleBracket] = \!\(\*SuperscriptBox[OverscriptBox[\"\[Alpha]\", \"_\"], \"p\"]\)\!\(\*SuperscriptBox[\"\[Beta]\", \"q\"]\)\[LeftAngleBracket]\[Alpha]\[VerticalBar]\[Beta]\[RightAngleBracket].";

Options[BosonicMatrixElement] = {"Basis" -> "Fock", "Generators" -> Automatic}

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


fockLadderElement[m_, n_, p_, q_] :=
    KroneckerDelta[m - p, n - q] Sqrt[FactorialPower[m, p]] Sqrt[FactorialPower[n, q]]

coherentOverlap[\[Alpha]_, \[Beta]_] :=
    Exp[Conjugate[\[Alpha]] \[Beta] - (Abs[\[Alpha]]^2 + Abs[\[Beta]]^2)/2]

ladderSum[eval_, expr_, v_] :=
    Module[{x, y},
        Total[(#[[2]] eval @@ #[[1]]) & /@
            CoefficientRules[
                Expand[
                    vevCollapse[BosonicNormalOrder[expr, {v, SuperDagger[v]}]] /.
                        {SuperDagger[v] -> x, v -> y}
                ],
                {x, y}
            ]
        ]
    ]

coherentBasisQ[opts___] :=
    OptionValue[BosonicMatrixElement, {opts}, "Basis"] === "Coherent"

(* Any formal symbol may name a mode, so a formal scalar is told apart only by declaring
   the generators; left Automatic the variables are inferred as before. *)
fieldVariable[expr_, opts___] :=
    Replace[
        DeleteDuplicates[
            Replace[OptionValue[BosonicMatrixElement, {opts}, "Generators"],
                Automatic :> ExtractNCVars[{expr}]] /. SuperDagger[w_] :> w
        ],
        {{u_} :> u, _ :> None}
    ]

BosonicMatrixElement[{m_, n_}, c_ ? NumericQ, opts : OptionsPattern[]] :=
    c If[coherentBasisQ[opts], coherentOverlap[m, n], KroneckerDelta[m, n]]

BosonicMatrixElement[{m_, n_}, expr_, opts : OptionsPattern[]] :=
    With[{v = fieldVariable[expr, opts]},
        If[ coherentBasisQ[opts],
            coherentOverlap[m, n] ladderSum[Conjugate[m]^#1 n^#2 &, expr, v],
            ladderSum[fockLadderElement[m, n, ##] &, expr, v]
        ] /;
            v =!= None && FreeQ[expr, _QuantumOperator | _QuantumState] &&
                vevPolynomialQ[expr, {v, SuperDagger[v]}]
    ]


(* Normal ordering changes the function: :Exp[c n]: is (1 + c)^n on the diagonal. *)
numberDiagonal[expr_, v_, k_] :=
    With[{nn = SuperDagger[v] ** v},
        With[{f = Replace[expr, {
                    NormalOrdered[Exp[c_. w_]] /; w === nn :> (1 + c)^k,
                    e_ :> ReplaceAll[e, nn -> k]
                }]},
            If[FreeQ[f, v | SuperDagger[v] | NormalOrdered], f, $Failed]
        ]
    ]

(* A wing is an exponential of one operator alone, so its element is a Taylor coefficient
   of Exp[P]: a displacement's linear wing and a squeeze's quadratic one share a formula,
   and the latter's parity rule is the vanishing of the odd coefficients. *)
wing[m_, k_, p_, x_] := SeriesCoefficient[Exp[p], {x, 0, m - k}] Sqrt[m!/k!]

wingPolynomial[s_, target_, v_, x_] :=
    With[{p = s /. target -> x /. NonCommutativeMultiply -> Times},
        If[FreeQ[p, v | _SuperDagger] && PolynomialQ[p, x] && TrueQ[(p /. x -> 0) == 0],
            p, $Failed]
    ]

normalFactor[fac_, v_, k_, x_] :=
    With[{s = Replace[fac, {Exp[e_] :> e, _ :> None}]},
        With[{
            up = If[s === None, $Failed, wingPolynomial[s, SuperDagger[v], v, x]],
            dn = If[s === None, $Failed, wingPolynomial[s, v, v, x]]
        },
            Which[
                up =!= $Failed, "Up" -> up,
                dn =!= $Failed, "Down" -> dn,
                True, Replace[numberDiagonal[fac, v, k], f : Except[$Failed] :> "Diagonal" -> f]
            ]
        ]
    ]

(* scalar Exp[P[ad]] f[n] Exp[Q[a]], any factor absent.  Wings of a kind commute so their
   exponents add; between Fock states the middle sum runs only to Min[m, n]. *)
normalProductElement[expr_, v_, m_, n_, j_] :=
    Module[{x, terms, scalar, tagged, u, w, f, r},
        terms = If[Head[expr] === Times, List @@ expr, {expr}];
        scalar = Times @@ Select[terms, FreeQ[#, v | SuperDagger[v]] &];
        tagged = Replace[Times @@ Select[terms, ! FreeQ[#, v | SuperDagger[v]] &],
            {q_NonCommutativeMultiply :> List @@ q, q_ :> {q}}];
        tagged = normalFactor[#, v, j, x] & /@ tagged;
        If[MemberQ[tagged, $Failed], Return[$Failed, Module]];
        u = Total[Cases[tagged, ("Up" -> c_) :> c]];
        w = Total[Cases[tagged, ("Down" -> c_) :> c]];
        f = Times @@ Cases[tagged, ("Diagonal" -> g_) :> g];
        r = Which[
            u === 0 && w === 0, scalar (f /. j -> n) KroneckerDelta[m, n],
            w === 0, scalar wing[m, n, u, x] (f /. j -> n),
            u === 0, scalar (f /. j -> m) wing[n, m, w, x],
            IntegerQ[m] && IntegerQ[n],
                scalar Sum[wing[m, j, u, x] f wing[n, j, w, x], {j, 0, Min[m, n]}],
            True, $Failed
        ];
        If[FreeQ[r, SeriesCoefficient], r, $Failed]
    ]

BosonicMatrixElement[{m_, n_}, expr_, opts : OptionsPattern[]] :=
    Module[{v = fieldVariable[expr, opts], j, r},
        r = Which[
            v === None, $Failed,
            coherentBasisQ[opts],
                Replace[numberDiagonal[expr, v, j],
                    f : Except[$Failed] :>
                        Exp[-(Abs[m]^2 + Abs[n]^2)/2] Sum[f (Conjugate[m] n)^j/j!, {j, 0, Infinity}]],
            True,
                (* A mixed exponential is not yet a product of wings; disentangle and retry. *)
                Replace[normalProductElement[expr, v, m, n, j],
                    $Failed /; ! FreeQ[expr, E] :>
                        normalProductElement[BosonicExpOrder[expr], v, m, n, j]]
        ];
        r /; r =!= $Failed && FreeQ[expr, _QuantumOperator | _QuantumState]
    ]
