(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicExpOrder"]

PackageExport["NormalOrdered"]

PackageExport["AntinormalOrdered"]

PackageExport["BosonicExpSimplify"]

PackageExport["BosonicExpBraid"]


algA[v_]   := algA[v] = NonCommutativeAlgebra[<|"Generators" -> {v, SuperDagger[v]}|>];

algA[vars_List] := algA[vars] = NonCommutativeAlgebra[<|"Generators" -> vars|>];


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


BosonicExpSimplify::usage =
"\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) simplifies a product of exponentials, or a conjugation by an exponential, of the bosonic operators vars, returning a closed form when the algebra closes and expr unchanged otherwise.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) detects the non-commutative variables and scalars appearing in expr.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Scalars\\\"->\", StyleBox[\"syms\", \"TI\"]}], \"]\"}]\) treats syms as commuting scalars.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{\[Ellipsis], \",\", \"Assumptions->\", StyleBox[\"assum\", \"TI\"]}], \"]\"}]\) uses \!\(\*StyleBox[\"assum\", \"TI\"]\) as simplifying assumptions while reducing.";

Options[BosonicExpSimplify] = {"Scalars" -> {}, Assumptions :> $Assumptions};

ncPolyQ[exprs_List, vars_] :=
    AllTrue[exprs, NonCommutativePolynomialQ[#, algA[vars]] &]

ncJoin[xs___] := ncJoinList @ DeleteCases[{xs}, 1]

ncJoinList[ys_] /; MemberQ[ys, 0] := 0
ncJoinList[{}] := 1
ncJoinList[{x_}] := x
ncJoinList[ys_] := NonCommutativeMultiply @@ ys

zeroQ[expr_] := (expr /. $nonuls) === 0

inverseExponentsQ[b_, c_] := TrueQ[Simplify[b + c] == 0]

ncCommutator[a_, b_, vars_, scalars_] :=
    ncCommutator[a, b, vars, scalars, $Assumptions]

ncCommutator[a_, b_, vars_, scalars_, assum_] :=
    ncCommutator[a, b, vars, scalars, assum] =
        Simplify[
            BosonicNormalOrder[Commutator[a, b], vars, "Scalars" -> scalars],
            assum
        ]

adRotate[op_, cmt_, lam_, vars_] /; cNumberQ[lam, vars] :=
    With[{w = I Sqrt[lam]}, Simplify[Cos[w] op + Sin[w] / w cmt] /. $nonuls]

adRotate[___] := $Failed

adAction[b_, op_, vars_, scalars_] :=
    With[{cmt = ncCommutator[b, op, vars, scalars]},
        With[{cmt2 = ncCommutator[b, cmt, vars, scalars]},
            Which[
                cmt === 0,                          op,
                cNumberQ[cmt, vars] || zeroQ[cmt2], op + cmt,
                True,                               adRotate[op, cmt, Simplify[cmt2 / op], vars]
            ]
        ]
    ]

expOrFailed[$Failed] := $Failed
expOrFailed[h_] := Exp[h]

expConjugate[b_, Exp[g_], vars_, scalars_] /; ncPolyQ[{b, g}, vars] :=
    expOrFailed @ adAction[b, g, vars, scalars]

expConjugate[b_, op_, vars_, scalars_] /; ncPolyQ[{b, op}, vars] :=
    adAction[b, op, vars, scalars]

expConjugate[___] := $Failed

expReduceRules[vars_, scalars_] := {
    NonCommutativeMultiply[l___, 1, r___] :> ncJoin[l, r],

    NonCommutativeMultiply[l___, Exp[b_], mid_, Exp[c_], r___] /; inverseExponentsQ[b, c] :>
        With[{g = expConjugate[b, mid, vars, scalars]},
            ncJoin[l, g, r] /; g =!= $Failed
        ],

    NonCommutativeMultiply[l___, Exp[b_], mid__, Exp[c_], r___] /;
            Length[{mid}] > 1 && inverseExponentsQ[b, c] :>
        With[{gs = expConjugate[b, #, vars, scalars] & /@ {mid}},
            ncJoin[l, ##, r] & @@ gs /; FreeQ[gs, $Failed]
        ],

    NonCommutativeMultiply[l___, Exp[b1_], Exp[b2_], r___] /; ncPolyQ[{b1, b2}, vars] :>
        With[{c = ncCommutator[b1, b2, vars, scalars]},
            ncJoin[l, Exp[c / 2] Exp[b1 + b2], r] /; cNumberQ[c, vars]
        ],

    NonCommutativeMultiply[l___, op_, Exp[b_], r___] /; ncPolyQ[{op, b}, vars] :>
        With[{cmt = ncCommutator[b, op, vars, scalars]},
            With[{cmt2 = ncCommutator[b, cmt, vars, scalars]},
                ncJoin[l, Exp[b], op - cmt, r] /; zeroQ[cmt2]
            ]
        ]
}

BosonicExpSimplify[expr_, opts : OptionsPattern[]] :=
    With[{
        vars = ExtractNCVars[{expr}],
        scalars = DeleteDuplicates @ Cases[{expr},
            s_Symbol /; ! FormalSymbolQ[s] && ! NumericQ[s], Infinity]
    },
        BosonicExpSimplify[expr, vars, opts, "Scalars" -> scalars]
    ]

BosonicExpSimplify[expr_, vars_List, opts : OptionsPattern[]] :=
    With[{assum = OptionValue[Assumptions], scalars = OptionValue["Scalars"]},
        Assuming[assum, expr //. expReduceRules[vars, scalars]]
    ]


braidDequantize[gen_, v_, x_, y_] := gen //. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    SuperDagger[v] -> x,
    v -> y,
    NonCommutativeMultiply -> Times
}

adFamilyRules[cr_, v_, vars_] /; ! FreeQ[Values[cr], Alternatives @@ vars] := $Failed

adFamilyRules[cr_, v_, vars_] /; Keys[cr] === {{1, 1}} :=
    With[{lam = Lookup[cr, Key[{1, 1}]]},
        {SuperDagger[v] -> Exp[lam] SuperDagger[v], v -> Exp[-lam] v}
    ]

adFamilyRules[cr_, v_, vars_] /; SubsetQ[{{1, 0}, {0, 1}}, Keys[cr]] :=
    With[{mu = Lookup[cr, Key[{1, 0}], 0], nu = Lookup[cr, Key[{0, 1}], 0]},
        {SuperDagger[v] -> SuperDagger[v] + nu, v -> v - mu}
    ]

adFamilyRules[cr_, v_, vars_] /; SubsetQ[{{2, 0}, {0, 2}}, Keys[cr]] :=
    With[{p = Lookup[cr, Key[{2, 0}], 0], q = Lookup[cr, Key[{0, 2}], 0]},
        With[{lam = 2 Sqrt[p q]},
            {SuperDagger[v] -> Cos[lam] SuperDagger[v] + 2 q Sinc[lam] v,
             v -> Cos[lam] v - 2 p Sinc[lam] SuperDagger[v]}
        ]
    ]

adFamilyRules[___] := $Failed

modeAdRules[gen_, v_, vars_] :=
    Module[{x, y},
        adFamilyRules[
            KeyDrop[
                Association @ CoefficientRules[
                    Expand @ braidDequantize[gen, v, x, y], {x, y}],
                Key[{0, 0}]
            ],
            v, vars
        ]
    ]

adRules[gen_, vars_] :=
    With[{perMode = modeAdRules[gen, #, vars] & /@ Select[vars, FreeQ[#, SuperDagger] &]},
        Catenate[perMode] /; FreeQ[perMode, $Failed]
    ]

adRules[___] := $Failed

adApply[rules_, expr_, vars_] := LiftNCScalars[expr /. rules, vars]

rulesBraidRight[vars_] := {
    NonCommutativeMultiply[l___, Exp[bA_], Exp[bB_], r___] :>
        With[{ar = adRules[bA, vars]},
            ncJoin[l, Exp[adApply[ar, bB, vars]], Exp[bA], r] /; ar =!= $Failed
        ]
}

rulesBraidLeft[vars_] := {
    NonCommutativeMultiply[l___, Exp[bA_], Exp[bB_], r___] :>
        With[{br = adRules[-bB, vars]},
            ncJoin[l, Exp[bB], Exp[adApply[br, bA, vars]], r] /; br =!= $Failed
        ]
}

BosonicExpBraid::usage =
"\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) moves the left factor of an adjacent pair of exponentials past the right one, transporting the right generator by the adjoint action of the left: \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) becomes \!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{\"A\", \",\", \"B\"}], \"]\"}]]\)\!\(\*SuperscriptBox[\"e\", \"A\"]\). It reorders only, performing no ordering, merging or cancellation.\n\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", \"\\\"Direction\\\"->\", \"\\\"Left\\\"\"}], \"]\"}]\) instead moves the right factor left, as \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) becomes \!\(\*SuperscriptBox[\"e\", \"B\"]\)\!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{RowBox[{\"-\", \"B\"}], \",\", \"A\"}], \"]\"}]]\).\n\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) uses vars as the field variables.";

BosonicExpBraid::direction = "`1` is not a valid setting for \"Direction\"; use \"Right\" or \"Left\"."

Options[BosonicExpBraid] = {"Direction" -> "Right"};

BosonicExpBraid[expr_, opts : OptionsPattern[]] :=
    BosonicExpBraid[expr, ExtractNCVars[{expr}], opts]

BosonicExpBraid[expr_, vars_List, opts : OptionsPattern[]] :=
    Replace[OptionValue["Direction"], {
        "Right" :> expr /. rulesBraidRight[vars],
        "Left" :> expr /. rulesBraidLeft[vars],
        other_ :> (Message[BosonicExpBraid::direction, other]; $Failed)
    }]
