(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicExpOrder"]
PackageExport["NormalOrdered"]
PackageExport["AntinormalOrdered"]
PackageExport["BosonicExpSimplify"]


algA[v_]   := algA[v] = NonCommutativeAlgebra[<|"Generators" -> {v, SuperDagger[v]}|>];


getNCPower[var_, var_] := 1;

getNCPower[GeneralizedPower[NonCommutativeMultiply, var_, n_], var_] := n;

getNCPower[_, _] := 0;


thm44RHS[L_, R_, lambda_, v_] :=
    With[{e = L + R - 1},
        With[{adagE = If[e == 1,
                SuperDagger[v],
                GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], e]
            ]},
            (1 - e lambda adagE)^(-(R/e)) **
                Exp[((1 - e lambda adagE)^(-1/e) - 1) ** SuperDagger[v] ** v]
        ]
    ]


thm44ANRHS[L_, R_, lambda_, v_] :=
    With[{e = L + R - 1},
        With[{aE = If[e == 1,
                v,
                GeneralizedPower[NonCommutativeMultiply, v, e]
            ]},
            (1 + e lambda aE)^(-(R/e)) **
                Exp[(1 - (1 + e lambda aE)^(-1/e)) ** v ** SuperDagger[v]]
        ]
    ]


rulesNO = {
    Exp[alpha_ SuperDagger[v_?FormalSymbolQ] + beta_ v_?FormalSymbolQ] :>
        Exp[(alpha beta)/2] Exp[alpha SuperDagger[v]] ** Exp[beta v],

    Exp[alpha_ v_?FormalSymbolQ] ** Exp[beta_ SuperDagger[v_?FormalSymbolQ]] :>
        Exp[alpha beta] Exp[beta SuperDagger[v]] ** Exp[alpha v],

    Exp[lambda_ SuperDagger[v_?FormalSymbolQ] ** v_?FormalSymbolQ] :>
        NormalOrdered[Exp[(E^lambda - 1) SuperDagger[v] ** v]],

    Exp[lambda_ v_?FormalSymbolQ ** SuperDagger[v_?FormalSymbolQ]] :>
        E^lambda NormalOrdered[Exp[(E^lambda - 1) SuperDagger[v] ** v]],

    Exp[alpha_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2])] **
            Exp[beta_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] :>
        Exp[(beta GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/(1 - 4 alpha beta)] **
            Exp[-Log[1 - 4 alpha beta] SuperDagger[v] ** v] **
            Exp[(alpha GeneralizedPower[NonCommutativeMultiply, v, 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[alpha_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2]) +
            beta_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] :>
        With[{omega = 2 Sqrt[alpha beta]},
            Exp[((beta Tan[omega]) GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/omega] **
                Exp[Log[Sec[omega]] (1/2 + SuperDagger[v] ** v)] **
                Exp[((alpha Tan[omega]) GeneralizedPower[NonCommutativeMultiply, v, 2])/omega]
        ],

    Exp[alpha_ v_?FormalSymbolQ] ** expr_ /; NonCommutativePolynomialQ[expr, algA[v]] :>
        BosonicNormalOrder[expr /. {SuperDagger[v] -> SuperDagger[v] + alpha}, {v, SuperDagger[v]}] ** Exp[alpha v],

    expr_ ** Exp[beta_ SuperDagger[v_?FormalSymbolQ]] /; NonCommutativePolynomialQ[expr, algA[v]] :>
        Exp[beta SuperDagger[v]] ** BosonicNormalOrder[expr /. {SuperDagger[v]-> SuperDagger[v], v -> v + beta}, {v, SuperDagger[v]}],

    Exp[lambda_*(adagL : (SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _])) **
            v_?FormalSymbolQ ** (adagR : (SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _]))] :>
        NormalOrdered @ thm44RHS[getNCPower[adagL, SuperDagger[v]], getNCPower[adagR, SuperDagger[v]], lambda, v],

    Exp[lambda_*v_?FormalSymbolQ ** (adagR : GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _])] :>
        NormalOrdered @ thm44RHS[0, getNCPower[adagR, SuperDagger[v]], lambda, v],

    Exp[lambda_*(adagL : GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _]) ** v_?FormalSymbolQ] :>
        NormalOrdered @ thm44RHS[getNCPower[adagL, SuperDagger[v]], 0, lambda, v],

    Exp[alpha_ SuperDagger[v1_?FormalSymbolQ] ** v2_?FormalSymbolQ + beta_ v1_?FormalSymbolQ ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(alpha Tanh[Omega]/Omega) SuperDagger[v1] ** v2] **
                (Exp[-Log[Cosh[Omega]] SuperDagger[v1] ** v1]
                Exp[Log[Cosh[Omega]] SuperDagger[v2] ** v2]) **
                Exp[(beta Tanh[Omega]/Omega) v1 ** SuperDagger[v2]]
        ],

    Exp[alpha_ v1_?FormalSymbolQ ** v2_?FormalSymbolQ + beta_ SuperDagger[v1_?FormalSymbolQ] ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(beta Tan[Omega]/Omega) SuperDagger[v1] ** SuperDagger[v2]] **
                Sec[Omega] Exp[Log[Sec[Omega]] (SuperDagger[v1] ** v1 + SuperDagger[v2] ** v2)] **
                Exp[(alpha Tan[Omega]/Omega) v1 ** v2]
        ]
};



rulesAN = {
    Exp[alpha_ SuperDagger[v_?FormalSymbolQ] + beta_ v_?FormalSymbolQ] :>
        Exp[-(alpha beta)/2] Exp[beta v] ** Exp[alpha SuperDagger[v]],

    Exp[beta_ SuperDagger[v_?FormalSymbolQ]] ** Exp[alpha_ v_?FormalSymbolQ] :>
        Exp[-alpha beta] Exp[alpha v] ** Exp[beta SuperDagger[v]],

    Exp[lambda_ v_?FormalSymbolQ ** SuperDagger[v_?FormalSymbolQ]] :>
        AntinormalOrdered[Exp[(1 - E^(-lambda)) v ** SuperDagger[v]]],

    Exp[lambda_ SuperDagger[v_?FormalSymbolQ] ** v_?FormalSymbolQ] :>
        E^(-lambda) AntinormalOrdered[Exp[(1 - E^(-lambda)) v ** SuperDagger[v]]],

    Exp[beta_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] **
            Exp[alpha_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2])] :>
        Exp[(alpha GeneralizedPower[NonCommutativeMultiply, v, 2])/(1 - 4 alpha beta)] **
            Exp[Log[1 - 4 alpha beta] v ** SuperDagger[v]] **
            Exp[(beta GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[alpha_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2]) +
            beta_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] :>
        With[{omega = 2 Sqrt[alpha beta]},
            Exp[((alpha Tan[omega]) GeneralizedPower[NonCommutativeMultiply, v, 2])/omega] **
                Exp[Log[Sec[omega]] (1/2 - v ** SuperDagger[v])] **
                Exp[((beta Tan[omega]) GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/omega]
        ],

    expr_ ** Exp[alpha_ v_?FormalSymbolQ] /; NonCommutativePolynomialQ[expr, algA[v]] :>
        Exp[alpha v] ** BosonicAntinormalOrder[expr /. {SuperDagger[v] -> SuperDagger[v] - alpha}, {v, SuperDagger[v]}],

    Exp[beta_ SuperDagger[v_?FormalSymbolQ]] ** expr_ /; NonCommutativePolynomialQ[expr, algA[v]] :>
        BosonicAntinormalOrder[expr /. {SuperDagger[v] -> SuperDagger[v], v -> v - beta}, {v, SuperDagger[v]}] ** Exp[beta SuperDagger[v]],

    Exp[lambda_*(aL : (v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, _])) **
            SuperDagger[v_?FormalSymbolQ] ** (aR : (v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, _]))] :>
        AntinormalOrdered @ thm44ANRHS[getNCPower[aL, v], getNCPower[aR, v], lambda, v],

    Exp[lambda_*SuperDagger[v_?FormalSymbolQ] ** (aR : GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, _])] :>
        AntinormalOrdered @ thm44ANRHS[0, getNCPower[aR, v], lambda, v],

    Exp[lambda_*(aL : GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, _]) ** SuperDagger[v_?FormalSymbolQ]] :>
        AntinormalOrdered @ thm44ANRHS[getNCPower[aL, v], 0, lambda, v],

    Exp[alpha_ SuperDagger[v1_?FormalSymbolQ] ** v2_?FormalSymbolQ + beta_ v1_?FormalSymbolQ ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(beta Tanh[Omega]/Omega) v1 ** SuperDagger[v2]] **
                (Exp[Log[Cosh[Omega]] v1 ** SuperDagger[v1]]
                Exp[-Log[Cosh[Omega]] v2 ** SuperDagger[v2]]) **
                Exp[(alpha Tanh[Omega]/Omega) SuperDagger[v1] ** v2]
        ],

    Exp[alpha_ v1_?FormalSymbolQ ** v2_?FormalSymbolQ + beta_ SuperDagger[v1_?FormalSymbolQ] ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(alpha Tan[Omega]/Omega) v1 ** v2] **
                Sec[Omega] Exp[-Log[Sec[Omega]] (v1 ** SuperDagger[v1] + v2 ** SuperDagger[v2])] **
                Exp[(beta Tan[Omega]/Omega) SuperDagger[v1] ** SuperDagger[v2]]
        ]
};


canonicalizeModeOrder[expr_] :=
    expr //. NonCommutativeMultiply[
            a___,
            x : (mx_?FormalSymbolQ | SuperDagger[mx_?FormalSymbolQ]),
            y : (my_?FormalSymbolQ | SuperDagger[my_?FormalSymbolQ]),
            b___
        ] /; mx =!= my && !OrderedQ[{mx, my}] :>
        NonCommutativeMultiply[a, y, x, b]




ladderPower[op_, 0] := 1
ladderPower[op_, p_] := GeneralizedPower[NonCommutativeMultiply, op, p]

ladderMonomial[adag_, 0, a_, 0] := 1
ladderMonomial[adag_, m_, a_, 0] := ladderPower[adag, m]
ladderMonomial[adag_, 0, a_, k_] := ladderPower[a, k]
ladderMonomial[adag_, m_, a_, k_] := ladderPower[adag, m] ** ladderPower[a, k]

orderedMonomial[NormalOrdered, v_, {m_, k_}] := ladderMonomial[SuperDagger[v], m, v, k]
orderedMonomial[AntinormalOrdered, v_, {m_, k_}] := ladderMonomial[v, k, SuperDagger[v], m]

expSeriesTerm[ord_, coeff_, v_] := With[{k = If[FreeQ[coeff, Global`k], Global`k, K[1]]},
    Inactivate[Sum[
        (coeff^k/k!) orderedMonomial[ord, v, {k, k}],
        {k, 0, Infinity}
    ], Sum]
]



NormalOrdered /: Times[c_, NormalOrdered[expr_]] /; FreeQ[c, _?FormalSymbolQ] := NormalOrdered[c expr]

AntinormalOrdered /: Times[c_, AntinormalOrdered[expr_]] /; FreeQ[c, _?FormalSymbolQ] := AntinormalOrdered[c expr]


NormalOrdered[expr_]["Series"] := expr /. {
    Exp[c_. + coeff_ SuperDagger[v_] ** v_] /; FreeQ[c, _?FormalSymbolQ] :>
        Exp[c] expSeriesTerm[NormalOrdered, coeff, v]
}

AntinormalOrdered[expr_]["Series"] := expr /. {
    Exp[c_. + coeff_ v_?FormalSymbolQ ** SuperDagger[v_?FormalSymbolQ]] /; FreeQ[c, _?FormalSymbolQ] :>
        Exp[c] expSeriesTerm[AntinormalOrdered, coeff, v]
}


dequantizeLadder[expr_, v_, x_, y_] := expr /. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    SuperDagger[v] -> x,
    v -> y,
    NonCommutativeMultiply -> Times
}


ladderMode[expr_] := First[Cases[expr, SuperDagger[w_] :> w, Infinity], None]


orderedSeries[ord_, expr_, lambda_] := Module[{v, x, y, n, coefN},
    v = ladderMode[expr];
    If[v === None, Return[$Failed, Module]];
    n = If[lambda =!= Global`n && FreeQ[expr, Global`n], Global`n, K[1]];
    coefN = SeriesCoefficient[dequantizeLadder[expr, v, x, y], {lambda, 0, n}];

    coefN = coefN /. Piecewise[{{val_, _}}, _] :> val;
    ord @ Inactivate[
        Sum[lambda^n (coefN /. {x -> SuperDagger[v], y -> v}), {n, 0, Infinity}],
        Sum
    ]
]

NormalOrdered[expr_]["Series", lambda_Symbol] := orderedSeries[NormalOrdered, expr, lambda]

AntinormalOrdered[expr_]["Series", lambda_Symbol] := orderedSeries[AntinormalOrdered, expr, lambda]


orderedSeriesOrder[ord_, expr_, lambda_, order_] := Module[{v, x, y},
    v = ladderMode[expr];
    If[v === None, Return[$Failed, Module]];
    ord[
        Normal @ Series[dequantizeLadder[expr, v, x, y], {lambda, 0, order}] /.
            {x -> SuperDagger[v], y -> v}
    ]
]

NormalOrdered[expr_]["Series", lambda_Symbol, order_Integer] :=
    orderedSeriesOrder[NormalOrdered, expr, lambda, order]

AntinormalOrdered[expr_]["Series", lambda_Symbol, order_Integer] :=
    orderedSeriesOrder[AntinormalOrdered, expr, lambda, order]



ncMonomials[ord_, expr_, v_] := Module[{x, y},
    Total[
        (#[[2]] orderedMonomial[ord, v, #[[1]]]) & /@
        CoefficientRules[Expand[dequantizeLadder[expr, v, x, y]], {x, y}]
    ]
]


ladderPolynomialQ[expr_] := Module[{v = ladderMode[expr], x, y},
    v =!= None && PolynomialQ[dequantizeLadder[expr, v, x, y], {x, y}]
]

NormalOrdered[expr_] /; ladderPolynomialQ[expr] :=
    ncMonomials[NormalOrdered, expr, ladderMode[expr]]

AntinormalOrdered[expr_] /; ladderPolynomialQ[expr] :=
    ncMonomials[AntinormalOrdered, expr, ladderMode[expr]]


expApplyRules[expr_, assum_, rules_] :=
    Assuming[assum,
        Simplify[canonicalizeModeOrder[expr] /. rules,
            ExcludedForms -> _GeneralizedPower
        ]
    ]


BosonicExpOrder::usage =
"\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) brings an expression built from exponentials of a formal field variable (e.g. \[FormalA]) and its \!\(\*RowBox[{StyleBox[\"var\", \"TI\"], SuperscriptBox[\"\[Dagger]\", \"\"]}]\) into normal order, using the boson disentangling identities.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", \"\\\"Ordering\\\"->\", \"\\\"Antinormal\\\"\"}], \"]\"}]\) brings expr into anti-normal order, with every annihilation operator to the left of every creation operator.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"assum\", \"TI\"]}], \"]\"}]\) uses \!\(\*StyleBox[\"assum\", \"TI\"]\) as additional simplifying assumptions.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{\[Ellipsis], \",\", \"TimeConstraint->\", StyleBox[\"t\", \"TI\"]}], \"]\"}]\) limits simplification to \!\(\*StyleBox[\"t\", \"TI\"]\) seconds (default 30).";

BosonicExpOrder::ordering = "`1` is not a valid setting for \"Ordering\"; use \"Normal\" or \"Antinormal\"."

Options[BosonicExpOrder] = {
    "Ordering" -> "Normal",
    TimeConstraint -> 30,
    Assumptions -> $Assumptions
};

BosonicExpOrder[expr_, opts : OptionsPattern[]] := BosonicExpOrder[expr, True, opts]

BosonicExpOrder[expr_, assum_, opts : OptionsPattern[]] :=
    Block[{tmax = OptionValue[TimeConstraint], fullAssum},
        fullAssum = Union[Flatten[{assum}], Flatten[{OptionValue[Assumptions]}]];
        Replace[OptionValue["Ordering"], {
            "Normal" :>
                TimeConstrained[expApplyRules[expr, fullAssum, rulesNO], tmax, expr],
            "Antinormal" :>
                TimeConstrained[expApplyRules[expr, fullAssum, rulesAN], tmax, expr],
            other_ :> (Message[BosonicExpOrder::ordering, other]; $Failed)
        }]
    ]


$nonuls = {0. -> 0, 0. I -> 0, Complex[0., 0.] -> 0, Complex[x_, 0.] -> x,
    Complex[0., y_] -> I y};


BosonicExpSimplify::usage =
"\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) simplifies a product of exponentials, or a conjugation by an exponential, of the bosonic operators vars, returning a closed form when the algebra closes and expr unchanged otherwise.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) detects the non-commutative variables and scalars appearing in expr.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) treats syms as commuting scalars.";

Options[BosonicExpSimplify] = {"Scalars" -> {}};

BosonicExpSimplify[expr_] :=
    Block[{vars, scalars},
        vars = ExtractNCVars[{expr}];
        scalars = DeleteDuplicates @ Cases[{expr}, s_Symbol /; !FormalSymbolQ[
            s] && !NumericQ[s], Infinity];
        BosonicExpSimplify[expr, vars, "Scalars" -> scalars]
    ]


BosonicExpSimplify[Exp[b1_] ** Exp[b2_], vars_List, opts : OptionsPattern[
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

BosonicExpSimplify[Exp[b1_] ** ops__ ** Exp[b3_], vars_List, opts : OptionsPattern[]] /;
    Simplify[b3 + b1] === 0 && Length[{ops}] > 1 :=
    NonCommutativeMultiply @@ (BosonicExpSimplify[Exp[b1] ** # ** Exp[b3], vars, opts] & /@ {ops})

BosonicExpSimplify[Exp[b1_] ** op_ ** Exp[b3_], vars_List, opts : OptionsPattern[
    ]] /; Simplify[b3 + b1] === 0 :=
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

BosonicExpSimplify[op_ ** Exp[b3_], vars_List, opts : OptionsPattern[]] :=
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

BosonicExpSimplify[Exp[b1_] ** Exp[b2_] ** Exp[b3_], vars_List, opts : OptionsPattern[]] /;
    Simplify[b3 + b1] === 0 :=
    Block[{scalars = OptionValue["Scalars"], newb2},
        newb2 = BosonicExpSimplify[Exp[b1] ** b2 ** Exp[b3], vars, opts];
        If[newb2 === Exp[b1] ** b2 ** Exp[b3],
            Exp[b1] ** Exp[b2] ** Exp[b3],
            Exp[newb2]
        ]
    ]
