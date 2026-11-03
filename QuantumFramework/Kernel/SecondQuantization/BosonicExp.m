(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicExpOrder"]

PackageExport["NormalOrdered"]

PackageExport["AntinormalOrdered"]

PackageExport["BosonicExpSimplify"]

PackageExport["BosonicExpBraid"]


ladderVars[lv_] := {lv["Creators"], lv["Annihilators"]}

ladderPairs[lv_] := Transpose[ladderVars[lv]]

shiftPolynomialQ[expr_, lv_] :=
    NonCommutativePolynomialQ[expr, lv["Algebra"]] && FreeQ[expr, Power[E, _?(! FreeQ[#, Alternatives @@ lv["Generators"]] &)]];

ladderRank[lv_][f_] := Which[
    FreeQ[f, Alternatives @@ lv["Creators"]], 3,
    FreeQ[f /. Thread[lv["Creators"] -> 1], Alternatives @@ lv["Annihilators"]], 1,
    True, 2
];

ladderModes[lv_, f_] := Union @ Lookup[lv["Modes"], Cases[f, Alternatives @@ lv["Generators"], {0, Infinity}]]

modeSwap[lv_, rank_] := NonCommutativeMultiply[l___,
        x : (_Power | _NormalOrdered | _AntinormalOrdered), y : (_Power | _NormalOrdered | _AntinormalOrdered), r___] /;
    ! FreeQ[x, Alternatives @@ lv["Generators"]] && ! FreeQ[y, Alternatives @@ lv["Generators"]] && rank[x] > rank[y] &&
        DisjointQ[ladderModes[lv, x], ladderModes[lv, y]] :>
    NonCommutativeMultiply[l, y, x, r];

modeMerge[lv_] := NonCommutativeMultiply[l___, Power[E, x_], Power[E, y_], r___] /;
    ! FreeQ[x, Alternatives @@ lv["Generators"]] && ladderRank[lv][x] === ladderRank[lv][y] && ladderRank[lv][x] =!= 2 &&
        ladderModes[lv, x] === ladderModes[lv, y] :>
    NonCommutativeMultiply[l, Exp[x + y], r];


getNCPower[var_, var_] := 1;

getNCPower[GeneralizedPower[NonCommutativeMultiply, var_, n_], var_] := n;

getNCPower[_, _] := 0;


thm44RHS[L_, R_, lambda_, {cr_, an_}] :=
    With[{e = L + R - 1},
        With[{adagE = If[e == 1, cr, GeneralizedPower[NonCommutativeMultiply, cr, e]]},
            (1 - e lambda adagE)^(-(R/e)) **
                Exp[((1 - e lambda adagE)^(-1/e) - 1) ** cr ** an]
        ]
    ]


thm44ANRHS[L_, R_, lambda_, {cr_, an_}] :=
    With[{e = L + R - 1},
        With[{aE = If[e == 1, an, GeneralizedPower[NonCommutativeMultiply, an, e]]},
            (1 + e lambda aE)^(-(R/e)) **
                Exp[(1 - (1 + e lambda aE)^(-1/e)) ** an ** cr]
        ]
    ]


modeRulesNO[{cr_, an_}, vv_, lv_] := {
    Exp[alpha_. cr + beta_. an] :>
        Exp[(alpha beta)/2] Exp[alpha cr] ** Exp[beta an],

    Exp[alpha_. an] ** Exp[beta_. cr] :>
        Exp[alpha beta] Exp[beta cr] ** Exp[alpha an],

    Exp[lambda_. cr ** an] :>
        NormalOrdered[Exp[(E^lambda - 1) cr ** an], vv],

    Exp[lambda_. an ** cr] :>
        E^lambda NormalOrdered[Exp[(E^lambda - 1) cr ** an], vv],

    Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2])] **
            Exp[beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] :>
        Exp[(beta GeneralizedPower[NonCommutativeMultiply, cr, 2])/(1 - 4 alpha beta)] **
            Exp[-Log[1 - 4 alpha beta] cr ** an] **
            Exp[(alpha GeneralizedPower[NonCommutativeMultiply, an, 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2]) +
            beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] :>
        With[{omega = 2 Sqrt[alpha beta]},
            Exp[((beta Tan[omega]) GeneralizedPower[NonCommutativeMultiply, cr, 2])/omega] **
                Exp[Log[Sec[omega]] (1/2 + cr ** an)] **
                Exp[((alpha Tan[omega]) GeneralizedPower[NonCommutativeMultiply, an, 2])/omega]
        ],

    Exp[alpha_. an] ** expr_ /; shiftPolynomialQ[expr, lv] :>
        BosonicNormalOrder[expr /. {cr -> cr + alpha}, vv] ** Exp[alpha an],

    expr_ ** Exp[beta_. cr] /; shiftPolynomialQ[expr, lv] :>
        Exp[beta cr] ** BosonicNormalOrder[expr /. {cr -> cr, an -> an + beta}, vv],

    Exp[lambda_*(adagL : (cr | GeneralizedPower[NonCommutativeMultiply, cr, _])) **
            an ** (adagR : (cr | GeneralizedPower[NonCommutativeMultiply, cr, _]))] :>
        NormalOrdered[thm44RHS[getNCPower[adagL, cr], getNCPower[adagR, cr], lambda, {cr, an}], vv],

    Exp[lambda_*an ** (adagR : GeneralizedPower[NonCommutativeMultiply, cr, _])] :>
        NormalOrdered[thm44RHS[0, getNCPower[adagR, cr], lambda, {cr, an}], vv],

    Exp[lambda_*(adagL : GeneralizedPower[NonCommutativeMultiply, cr, _]) ** an] :>
        NormalOrdered[thm44RHS[getNCPower[adagL, cr], 0, lambda, {cr, an}], vv]
}

pairRulesNO[{cr1_, an1_}, {cr2_, an2_}] := {
    Exp[alpha_. cr1 ** an2 + beta_. an1 ** cr2] :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(alpha Tanh[Omega]/Omega) cr1 ** an2] **
                (Exp[-Log[Cosh[Omega]] cr1 ** an1]
                Exp[Log[Cosh[Omega]] cr2 ** an2]) **
                Exp[(beta Tanh[Omega]/Omega) an1 ** cr2]
        ],

    Exp[alpha_. an1 ** an2 + beta_. cr1 ** cr2] :>
        With[{Omega = Sqrt[alpha beta]},
            Sec[Omega] (Exp[(beta Tan[Omega]/Omega) cr1 ** cr2] **
                Exp[Log[Sec[Omega]] (cr1 ** an1 + cr2 ** an2)] **
                Exp[(alpha Tan[Omega]/Omega) an1 ** an2])
        ]
}


modeRulesAN[{cr_, an_}, vv_, lv_] := {
    Exp[alpha_. cr + beta_. an] :>
        Exp[-(alpha beta)/2] Exp[beta an] ** Exp[alpha cr],

    Exp[beta_. cr] ** Exp[alpha_. an] :>
        Exp[-alpha beta] Exp[alpha an] ** Exp[beta cr],

    Exp[lambda_. an ** cr] :>
        AntinormalOrdered[Exp[(1 - E^(-lambda)) an ** cr], vv],

    Exp[lambda_. cr ** an] :>
        E^(-lambda) AntinormalOrdered[Exp[(1 - E^(-lambda)) an ** cr], vv],

    Exp[beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] **
            Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2])] :>
        Exp[(alpha GeneralizedPower[NonCommutativeMultiply, an, 2])/(1 - 4 alpha beta)] **
            Exp[Log[1 - 4 alpha beta] an ** cr] **
            Exp[(beta GeneralizedPower[NonCommutativeMultiply, cr, 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2]) +
            beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] :>
        With[{omega = 2 Sqrt[alpha beta]},
            Exp[((alpha Tan[omega]) GeneralizedPower[NonCommutativeMultiply, an, 2])/omega] **
                Exp[Log[Sec[omega]] (1/2 - an ** cr)] **
                Exp[((beta Tan[omega]) GeneralizedPower[NonCommutativeMultiply, cr, 2])/omega]
        ],

    expr_ ** Exp[alpha_. an] /; shiftPolynomialQ[expr, lv] :>
        Exp[alpha an] ** BosonicAntinormalOrder[expr /. {cr -> cr - alpha}, vv],

    Exp[beta_. cr] ** expr_ /; shiftPolynomialQ[expr, lv] :>
        BosonicAntinormalOrder[expr /. {cr -> cr, an -> an - beta}, vv] ** Exp[beta cr],

    Exp[lambda_*(aL : (an | GeneralizedPower[NonCommutativeMultiply, an, _])) **
            cr ** (aR : (an | GeneralizedPower[NonCommutativeMultiply, an, _]))] :>
        AntinormalOrdered[thm44ANRHS[getNCPower[aL, an], getNCPower[aR, an], lambda, {cr, an}], vv],

    Exp[lambda_*cr ** (aR : GeneralizedPower[NonCommutativeMultiply, an, _])] :>
        AntinormalOrdered[thm44ANRHS[0, getNCPower[aR, an], lambda, {cr, an}], vv],

    Exp[lambda_*(aL : GeneralizedPower[NonCommutativeMultiply, an, _]) ** cr] :>
        AntinormalOrdered[thm44ANRHS[getNCPower[aL, an], 0, lambda, {cr, an}], vv]
}

pairRulesAN[{cr1_, an1_}, {cr2_, an2_}] := {
    Exp[alpha_. cr1 ** an2 + beta_. an1 ** cr2] :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(beta Tanh[Omega]/Omega) an1 ** cr2] **
                (Exp[Log[Cosh[Omega]] an1 ** cr1]
                Exp[-Log[Cosh[Omega]] an2 ** cr2]) **
                Exp[(alpha Tanh[Omega]/Omega) cr1 ** an2]
        ],

    Exp[alpha_. an1 ** an2 + beta_. cr1 ** cr2] :>
        With[{Omega = Sqrt[alpha beta]},
            Sec[Omega] (Exp[(alpha Tan[Omega]/Omega) an1 ** an2] **
                Exp[-Log[Sec[Omega]] (an1 ** cr1 + an2 ** cr2)] **
                Exp[(beta Tan[Omega]/Omega) cr1 ** cr2])
        ]
}


expRules[lv_, "Normal"] := expRules[lv, "Normal"] = Join[
    Flatten @ Transpose[modeRulesNO[#, ladderVars[lv], lv] & /@ ladderPairs[lv]],
    Flatten[pairRulesNO @@@ Permutations[ladderPairs[lv], {2}]],
    {modeSwap[lv, ladderRank[lv]], modeMerge[lv]}
]

expRules[lv_, "Antinormal"] := expRules[lv, "Antinormal"] = Join[
    Flatten @ Transpose[modeRulesAN[#, ladderVars[lv], lv] & /@ ladderPairs[lv]],
    Flatten[pairRulesAN @@@ Permutations[ladderPairs[lv], {2}]],
    {modeSwap[lv, 4 - ladderRank[lv][#] &], modeMerge[lv]}
]


canonicalizeGenerators[expr_, lv_] :=
    With[{g = Alternatives @@ lv["Generators"]},
        (expr /. Exp[e_] :> Exp[Expand[e]]) //.
            NonCommutativeMultiply[l___, x : g, y : g, r___] /; lv["Modes"][x] > lv["Modes"][y] :>
                NonCommutativeMultiply[l, y, x, r]
    ]




ladderPower[op_, 0] := 1
ladderPower[op_, p_] := GeneralizedPower[NonCommutativeMultiply, op, p]

ladderMonomial[adag_, 0, a_, 0] := 1
ladderMonomial[adag_, m_, a_, 0] := ladderPower[adag, m]
ladderMonomial[adag_, 0, a_, k_] := ladderPower[a, k]
ladderMonomial[adag_, m_, a_, k_] := ladderPower[adag, m] ** ladderPower[a, k]

orderedMonomial[NormalOrdered, {cr_, an_}, {m_, k_}] := ladderMonomial[cr, m, an, k]
orderedMonomial[AntinormalOrdered, {cr_, an_}, {m_, k_}] := ladderMonomial[an, k, cr, m]

MakeBoxes[SeriesIndex[s_String], StandardForm] ^:= InterpretationBox[s, SeriesIndex[s]]

MakeBoxes[SeriesIndex[s_String], TraditionalForm] ^:= InterpretationBox[s, SeriesIndex[s]]

expSeriesTerm[ord_, coeff_, pair_] := With[{k = SeriesIndex["k"]},
    Inactivate[Sum[
        (coeff^k/k!) orderedMonomial[ord, pair, {k, k}],
        {k, 0, Infinity}
    ], Sum]
]



NormalOrdered /: Times[s_, NormalOrdered[expr_, vars_]] :=
    With[{lv = ResolveLadderVariables[NormalOrdered, vars]},
        NormalOrdered[s expr, vars] /; AssociationQ[lv] && FreeQ[s, Alternatives @@ lv["Generators"]]
    ]

AntinormalOrdered /: Times[s_, AntinormalOrdered[expr_, vars_]] :=
    With[{lv = ResolveLadderVariables[AntinormalOrdered, vars]},
        AntinormalOrdered[s expr, vars] /; AssociationQ[lv] && FreeQ[s, Alternatives @@ lv["Generators"]]
    ]


NormalOrdered[expr_, vars_]["Series"] :=
    With[{lv = ResolveLadderVariables[NormalOrdered, vars]},
        (expr /. ((Exp[s_. + coeff_ #1 ** #2] /; FreeQ[s, Alternatives @@ lv["Generators"]] :>
            Exp[s] expSeriesTerm[NormalOrdered, coeff, {#1, #2}]) & @@@ ladderPairs[lv])) /; AssociationQ[lv]
    ]

AntinormalOrdered[expr_, vars_]["Series"] :=
    With[{lv = ResolveLadderVariables[AntinormalOrdered, vars]},
        (expr /. ((Exp[s_. + coeff_ #2 ** #1] /; FreeQ[s, Alternatives @@ lv["Generators"]] :>
            Exp[s] expSeriesTerm[AntinormalOrdered, coeff, {#1, #2}]) & @@@ ladderPairs[lv])) /; AssociationQ[lv]
    ]


dequantizeLadder[expr_, {cr_, an_}, x_, y_] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    cr -> x,
    an -> y,
    NonCommutativeMultiply -> Times
}


ladderMode[expr_, lv_] :=
    Replace[Select[ladderPairs[lv], ! FreeQ[expr, Alternatives @@ #] &], {{p_} :> p, _ :> None}]


orderedSeries[ord_, expr_, lambda_, lv_] := Module[{v, x, y, n, coefN},
    v = ladderMode[expr, lv];
    If[v === None, Return[$Failed, Module]];
    n = SeriesIndex["n"];
    coefN = SeriesCoefficient[dequantizeLadder[expr, v, x, y], {lambda, 0, n}];

    coefN = coefN /. Piecewise[{{val_, _}}, _] :> val;
    ord[
        Inactivate[Sum[lambda^n (coefN /. {x -> First[v], y -> Last[v]}), {n, 0, Infinity}], Sum],
        ladderVars[lv]
    ]
]

NormalOrdered[expr_, vars_]["Series", lambda_Symbol] :=
    With[{lv = ResolveLadderVariables[NormalOrdered, vars]}, orderedSeries[NormalOrdered, expr, lambda, lv] /; AssociationQ[lv]]

AntinormalOrdered[expr_, vars_]["Series", lambda_Symbol] :=
    With[{lv = ResolveLadderVariables[AntinormalOrdered, vars]}, orderedSeries[AntinormalOrdered, expr, lambda, lv] /; AssociationQ[lv]]


orderedSeriesOrder[ord_, expr_, lambda_, order_, lv_] := Module[{v, x, y},
    v = ladderMode[expr, lv];
    If[v === None, Return[$Failed, Module]];
    ord[
        Normal @ Series[dequantizeLadder[expr, v, x, y], {lambda, 0, order}] /. {x -> First[v], y -> Last[v]},
        ladderVars[lv]
    ]
]

NormalOrdered[expr_, vars_]["Series", lambda_Symbol, order_Integer] :=
    With[{lv = ResolveLadderVariables[NormalOrdered, vars]}, orderedSeriesOrder[NormalOrdered, expr, lambda, order, lv] /; AssociationQ[lv]]

AntinormalOrdered[expr_, vars_]["Series", lambda_Symbol, order_Integer] :=
    With[{lv = ResolveLadderVariables[AntinormalOrdered, vars]}, orderedSeriesOrder[AntinormalOrdered, expr, lambda, order, lv] /; AssociationQ[lv]]



ncMonomials[ord_, expr_, pair_] := Module[{x, y},
    Total[
        (#[[2]] orderedMonomial[ord, pair, #[[1]]]) & /@
        CoefficientRules[Expand[dequantizeLadder[expr, pair, x, y]], {x, y}]
    ]
]


ladderPolynomialQ[expr_, pair_] :=
    Module[{x, y},
        With[{d = dequantizeLadder[expr, pair, x, y]},
            PolynomialQ[d, {x, y}] && AllTrue[Variables[d], MatchQ[#, x | y] || FreeQ[#, x | y] &]
        ]
    ]

NormalOrdered[expr_, vars_] :=
    With[{p = Replace[Quiet @ ResolveLadderVariables[NormalOrdered, vars], {lv_Association :> ladderMode[expr, lv], _ -> None}]},
        ncMonomials[NormalOrdered, expr, p] /; p =!= None && ladderPolynomialQ[expr, p]
    ]

AntinormalOrdered[expr_, vars_] :=
    With[{p = Replace[Quiet @ ResolveLadderVariables[AntinormalOrdered, vars], {lv_Association :> ladderMode[expr, lv], _ -> None}]},
        ncMonomials[AntinormalOrdered, expr, p] /; p =!= None && ladderPolynomialQ[expr, p]
    ]


expApplyRules[expr_, assum_, rules_, lv_] :=
    Assuming[assum,
        Simplify[
            FixedPoint[
                (# /. Prepend[rules, n : (_NormalOrdered | _AntinormalOrdered) :> n]) //.
                    (NonCommutativeMultiply[l___, Times[s_, m_NonCommutativeMultiply], r___] /; FreeQ[s, Alternatives @@ lv["Generators"]] :>
                        s NonCommutativeMultiply[l, m, r]) &,
                canonicalizeGenerators[expr, lv]
            ],
            ExcludedForms -> _GeneralizedPower
        ]
    ]


BosonicExpOrder::usage =
"\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) brings an expression built from exponentials of the ladder operators vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra into normal order, using the boson disentangling identities.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Ordering\\\"->\", \"\\\"Antinormal\\\"\"}], \"]\"}]\) brings expr into anti-normal order, with every annihilator to the left of every creator.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{\[Ellipsis], \",\", \"TimeConstraint->\", StyleBox[\"t\", \"TI\"]}], \"]\"}]\) limits simplification to \!\(\*StyleBox[\"t\", \"TI\"]\) seconds (default 30).";

BosonicExpOrder::ordering = "`1` is not a valid setting for \"Ordering\"; use \"Normal\" or \"Antinormal\"."

BosonicExpOrder::assum = "`1` reads as assumptions rather than ladder operators; give the ladder operators second and the assumptions as Assumptions -> `1`."

Options[BosonicExpOrder] = {
    "Ordering" -> "Normal",
    TimeConstraint -> 30,
    Assumptions -> $Assumptions
};

BosonicExpOrder[expr_, opts : OptionsPattern[]] := (Message[BosonicExpOrder::argr, BosonicExpOrder, 2]; $Failed)

With[{p = True | False | _Greater | _GreaterEqual | _Less | _LessEqual | _Equal | _Unequal | _Element | _NotElement | _And | _Or},
    BosonicExpOrder[expr_, assum : (p | {p ..}), OptionsPattern[]] := (Message[BosonicExpOrder::assum, assum]; $Failed)
]

BosonicExpOrder[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicExpOrder, vars]},
        If[ lv === $Failed,
            $Failed,
            Replace[OptionValue["Ordering"], {
                ord : "Normal" | "Antinormal" :>
                    TimeConstrained[expApplyRules[expr, OptionValue[Assumptions], expRules[lv, ord], lv], OptionValue[TimeConstraint], expr],
                other_ :> (Message[BosonicExpOrder::ordering, other]; $Failed)
            }]
        ]
    ]


BosonicExpSimplify::usage =
"\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) simplifies a product of exponentials, or a conjugation by an exponential, of the ladder operators vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra, returning a closed form when the algebra closes and expr unchanged otherwise.\n\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{\[Ellipsis], \",\", \"Assumptions->\", StyleBox[\"assum\", \"TI\"]}], \"]\"}]\) uses \!\(\*StyleBox[\"assum\", \"TI\"]\) as simplifying assumptions while reducing.";

Options[BosonicExpSimplify] = {Assumptions :> $Assumptions};

ncPolyQ[exprs_List, lv_] :=
    AllTrue[exprs, NonCommutativePolynomialQ[#, lv["Algebra"]] && vevPolynomialQ[#, lv["Generators"]] &]

ncJoin[xs___] := ncJoinList @ DeleteCases[{xs}, 1]

ncJoinList[ys_] /; MemberQ[ys, 0] := 0
ncJoinList[{}] := 1
ncJoinList[{x_}] := x
ncJoinList[ys_] := NonCommutativeMultiply @@ ys

zeroQ[expr_] := (expr /. $nonuls) === 0

inverseExponentsQ[b_, c_] := TrueQ[Simplify[b + c] == 0]

ncCommutator[a_, b_, lv_] :=
    ncCommutator[a, b, lv, $Assumptions]

ncCommutator[a_, b_, lv_, assum_] :=
    ncCommutator[a, b, lv, assum] =
        Simplify[BosonicNormalOrder[Commutator[a, b], ladderVars[lv]], assum]

adRotate[op_, cmt_, lam_, lv_] /; cNumberQ[lam, lv["Generators"]] :=
    With[{w = I Sqrt[lam]}, Simplify[Cos[w] op + Sin[w] / w cmt] /. $nonuls]

adRotate[___] := $Failed

adAction[b_, op_, lv_] :=
    With[{cmt = ncCommutator[b, op, lv]},
        With[{cmt2 = ncCommutator[b, cmt, lv]},
            Which[
                cmt === 0,                                      op,
                cNumberQ[cmt, lv["Generators"]] || zeroQ[cmt2], op + cmt,
                True,                                           adRotate[op, cmt, Simplify[cmt2 / op], lv]
            ]
        ]
    ]

expOrFailed[$Failed] := $Failed
expOrFailed[h_] := Exp[h]

expConjugate[b_, Exp[g_], lv_] /; ncPolyQ[{b, g}, lv] :=
    expOrFailed @ adAction[b, g, lv]

expConjugate[b_, op_, lv_] /; ncPolyQ[{b, op}, lv] :=
    adAction[b, op, lv]

expConjugate[___] := $Failed

expReduceRules[lv_] := {
    NonCommutativeMultiply[l___, 1, r___] :> ncJoin[l, r],

    NonCommutativeMultiply[l___, Exp[b_], mid_, Exp[c_], r___] /; inverseExponentsQ[b, c] :>
        With[{g = expConjugate[b, mid, lv]},
            ncJoin[l, g, r] /; g =!= $Failed
        ],

    NonCommutativeMultiply[l___, Exp[b_], mid__, Exp[c_], r___] /;
            Length[{mid}] > 1 && inverseExponentsQ[b, c] :>
        With[{gs = expConjugate[b, #, lv] & /@ {mid}},
            ncJoin[l, ##, r] & @@ gs /; FreeQ[gs, $Failed]
        ],

    NonCommutativeMultiply[l___, Exp[b1_], Exp[b2_], r___] /; ncPolyQ[{b1, b2}, lv] :>
        With[{c = ncCommutator[b1, b2, lv]},
            ncJoin[l, Exp[c / 2] Exp[b1 + b2], r] /; cNumberQ[c, lv["Generators"]]
        ],

    NonCommutativeMultiply[l___, op_, Exp[b_], r___] /; ncPolyQ[{op, b}, lv] :>
        With[{cmt = ncCommutator[b, op, lv]},
            With[{cmt2 = ncCommutator[b, cmt, lv]},
                ncJoin[l, Exp[b], op - cmt, r] /; zeroQ[cmt2]
            ]
        ]
}

BosonicExpSimplify[expr_, opts : OptionsPattern[]] := (Message[BosonicExpSimplify::argr, BosonicExpSimplify, 2]; $Failed)

BosonicExpSimplify[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicExpSimplify, vars]},
        If[lv === $Failed, $Failed, Assuming[OptionValue[Assumptions], expr //. expReduceRules[lv]]]
    ]


adFamilyRules[cr_, _, gens_] /; ! FreeQ[Values[cr], Alternatives @@ gens] := $Failed

adFamilyRules[cr_, {cre_, ann_}, _] /; Keys[cr] === {{1, 1}} :=
    With[{lam = Lookup[cr, Key[{1, 1}]]},
        {cre -> Exp[lam] cre, ann -> Exp[-lam] ann}
    ]

adFamilyRules[cr_, {cre_, ann_}, _] /; SubsetQ[{{1, 0}, {0, 1}}, Keys[cr]] :=
    With[{mu = Lookup[cr, Key[{1, 0}], 0], nu = Lookup[cr, Key[{0, 1}], 0]},
        {cre -> cre + nu, ann -> ann - mu}
    ]

adFamilyRules[cr_, {cre_, ann_}, _] /; SubsetQ[{{2, 0}, {0, 2}}, Keys[cr]] :=
    With[{p = Lookup[cr, Key[{2, 0}], 0], q = Lookup[cr, Key[{0, 2}], 0]},
        With[{lam = 2 Sqrt[p q]},
            {cre -> Cos[lam] cre + 2 q Sinc[lam] ann,
             ann -> Cos[lam] ann - 2 p Sinc[lam] cre}
        ]
    ]

adFamilyRules[___] := $Failed

modeAdRules[gen_, pair_, gens_] :=
    Module[{x, y},
        adFamilyRules[
            KeyDrop[
                Association @ CoefficientRules[
                    Expand @ dequantizeLadder[gen, pair, x, y], {x, y}],
                Key[{0, 0}]
            ],
            pair, gens
        ]
    ]

adRules[gen_, lv_] :=
    With[{perMode = modeAdRules[gen, #, lv["Generators"]] & /@ ladderPairs[lv]},
        Catenate[perMode] /; FreeQ[perMode, $Failed]
    ]

adRules[___] := $Failed

adApply[rules_, expr_, lv_] := LiftNCScalars[expr /. rules, lv["Generators"]]

rulesBraidRight[lv_] := {
    NonCommutativeMultiply[l___, Exp[bA_], Exp[bB_], r___] :>
        With[{ar = adRules[bA, lv]},
            ncJoin[l, Exp[adApply[ar, bB, lv]], Exp[bA], r] /; ar =!= $Failed
        ]
}

rulesBraidLeft[lv_] := {
    NonCommutativeMultiply[l___, Exp[bA_], Exp[bB_], r___] :>
        With[{br = adRules[-bB, lv]},
            ncJoin[l, Exp[bB], Exp[adApply[br, bA, lv]], r] /; br =!= $Failed
        ]
}

BosonicExpBraid::usage =
"\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) moves the left factor of an adjacent pair of exponentials of the ladder operators vars past the right one, transporting the right generator by the adjoint action of the left: \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) becomes \!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{\"A\", \",\", \"B\"}], \"]\"}]]\)\!\(\*SuperscriptBox[\"e\", \"A\"]\). It reorders only, performing no ordering, merging or cancellation.\n\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Direction\\\"->\", \"\\\"Left\\\"\"}], \"]\"}]\) instead moves the right factor left, as \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) becomes \!\(\*SuperscriptBox[\"e\", \"B\"]\)\!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{RowBox[{\"-\", \"B\"}], \",\", \"A\"}], \"]\"}]]\).";

BosonicExpBraid::direction = "`1` is not a valid setting for \"Direction\"; use \"Right\" or \"Left\"."

Options[BosonicExpBraid] = {"Direction" -> "Right"};

BosonicExpBraid[expr_, opts : OptionsPattern[]] := (Message[BosonicExpBraid::argr, BosonicExpBraid, 2]; $Failed)

BosonicExpBraid[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicExpBraid, vars]},
        If[ lv === $Failed,
            $Failed,
            Replace[OptionValue["Direction"], {
                "Right" :> expr /. rulesBraidRight[lv],
                "Left" :> expr /. rulesBraidLeft[lv],
                other_ :> (Message[BosonicExpBraid::direction, other]; $Failed)
            }]
        ]
    ]
