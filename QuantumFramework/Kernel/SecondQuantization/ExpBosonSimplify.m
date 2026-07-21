(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["ExpBosonNormalOrder"]
PackageExport["NormalOrdered"]


algA[v_]    := algA[v]    = NonCommutativeAlgebra[<|"Generators" -> {v}|>];
algAdag[v_] := algAdag[v] = NonCommutativeAlgebra[<|"Generators" -> {SuperDagger[v]}|>];


getNCPower[var_, var_] := 1;

getNCPower[GeneralizedPower[NonCommutativeMultiply, var_, n_], var_] := n;

getNCPower[_, _] := 0;


thm44RHS[L_, R_, \[Lambda]_, v_] :=
    With[{e = L + R - 1},
        With[{adagE = If[e == 1,
                SuperDagger[v],
                GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], e]
            ]},
            (1 - e \[Lambda] adagE)^(-(R/e)) **
                Exp[((1 - e \[Lambda] adagE)^(-1/e) - 1) ** SuperDagger[v] ** v]
        ]
    ]


rulesNO = {
    Exp[\[Alpha]_ SuperDagger[v_?FormalSymbolQ] + \[Beta]_ v_?FormalSymbolQ] :>
        Exp[(\[Alpha] \[Beta])/2] Exp[\[Alpha] SuperDagger[v]] ** Exp[\[Beta] v],

    Exp[\[Alpha]_ v_?FormalSymbolQ] ** Exp[\[Beta]_ SuperDagger[v_?FormalSymbolQ]] :>
        Exp[\[Alpha] \[Beta]] Exp[\[Beta] SuperDagger[v]] ** Exp[\[Alpha] v],

    Exp[\[Lambda]_ SuperDagger[v_?FormalSymbolQ] ** v_?FormalSymbolQ] :>
        NormalOrdered[Exp[(E^\[Lambda] - 1) SuperDagger[v] ** v]],

    Exp[\[Lambda]_ v_?FormalSymbolQ ** SuperDagger[v_?FormalSymbolQ]] :>
        E^\[Lambda] NormalOrdered[Exp[(E^\[Lambda] - 1) SuperDagger[v] ** v]],

    Exp[\[Alpha]_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2])] **
            Exp[\[Beta]_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] :>
        Exp[(\[Beta] GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/(1 - 4 \[Alpha] \[Beta])] **
            Exp[-Log[1 - 4 \[Alpha] \[Beta]] SuperDagger[v] ** v] **
            Exp[(\[Alpha] GeneralizedPower[NonCommutativeMultiply, v, 2])/(1 - 4 \[Alpha] \[Beta])] / Sqrt[1 - 4 \[Alpha] \[Beta]],

    Exp[\[Alpha]_ (v_?FormalSymbolQ ** v_?FormalSymbolQ | GeneralizedPower[NonCommutativeMultiply, v_?FormalSymbolQ, 2]) +
            \[Beta]_ (SuperDagger[v_?FormalSymbolQ] ** SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], 2])] :>
        With[{\[Omega] = 2 Sqrt[\[Alpha] \[Beta]]},
            Exp[((\[Beta] Tan[\[Omega]]) GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], 2])/\[Omega]] **
                Sqrt[Sec[\[Omega]]] Exp[Log[Sec[\[Omega]]] SuperDagger[v] ** v] **
                Exp[((\[Alpha] Tan[\[Omega]]) GeneralizedPower[NonCommutativeMultiply, v, 2])/\[Omega]]
        ],

    Exp[\[Alpha]_ v_?FormalSymbolQ] ** expr_ /; NonCommutativePolynomialQ[expr, algAdag[v]] :>
        NonCommutativeExpand[expr /. {SuperDagger[v] -> SuperDagger[v] + \[Alpha]}, algAdag[v]] ** Exp[\[Alpha] v],

    expr_ ** Exp[\[Beta]_ SuperDagger[v_?FormalSymbolQ]] /; NonCommutativePolynomialQ[expr, v, algA[v]] :>
        Exp[\[Beta] SuperDagger[v]] ** NonCommutativeExpand[expr /. {v -> v + \[Beta]}, algA[v]],

    Exp[\[Lambda]_*(adagL : (SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _])) **
            v_?FormalSymbolQ ** (adagR : (SuperDagger[v_?FormalSymbolQ] | GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _]))] :>
        NormalOrdered @ thm44RHS[getNCPower[adagL, SuperDagger[v]], getNCPower[adagR, SuperDagger[v]], \[Lambda], v],

    Exp[\[Lambda]_*v_?FormalSymbolQ ** (adagR : GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _])] :>
        NormalOrdered @ thm44RHS[0, getNCPower[adagR, SuperDagger[v]], \[Lambda], v],

    Exp[\[Lambda]_*(adagL : GeneralizedPower[NonCommutativeMultiply, SuperDagger[v_?FormalSymbolQ], _]) ** v_?FormalSymbolQ] :>
        NormalOrdered @ thm44RHS[getNCPower[adagL, SuperDagger[v]], 0, \[Lambda], v],

    Exp[\[Alpha]_ SuperDagger[v1_?FormalSymbolQ] ** v2_?FormalSymbolQ + \[Beta]_ v1_?FormalSymbolQ ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{\[CapitalOmega] = Sqrt[\[Alpha] \[Beta]]},
            Exp[(\[Alpha] Tanh[\[CapitalOmega]]/\[CapitalOmega]) SuperDagger[v1] ** v2] **
                Exp[-Log[Cosh[\[CapitalOmega]]] SuperDagger[v1] ** v1] **
                Exp[Log[Cosh[\[CapitalOmega]]] SuperDagger[v2] ** v2] **
                Exp[(\[Beta] Tanh[\[CapitalOmega]]/\[CapitalOmega]) v1 ** SuperDagger[v2]]
        ],

    Exp[\[Alpha]_ v1_?FormalSymbolQ ** v2_?FormalSymbolQ + \[Beta]_ SuperDagger[v1_?FormalSymbolQ] ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{\[CapitalOmega] = Sqrt[\[Alpha] \[Beta]]},
            Exp[(\[Beta] Tan[\[CapitalOmega]]/\[CapitalOmega]) SuperDagger[v1] ** SuperDagger[v2]] **
                Sec[\[CapitalOmega]] Exp[Log[Sec[\[CapitalOmega]]] (SuperDagger[v1] ** v1 + SuperDagger[v2] ** v2)] **
                Exp[(\[Alpha] Tan[\[CapitalOmega]]/\[CapitalOmega]) v1 ** v2]
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


NormalOrdered[pref_. Exp[\[Mu]_ SuperDagger[v_?FormalSymbolQ] ** v_?FormalSymbolQ]]["Series"] :=
    With[{k = Symbol["Global`k"]},
        HoldForm[pref Sum[
            (\[Mu]^k/k!) GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], k] **
                GeneralizedPower[NonCommutativeMultiply, v, k],
            {k, 0, Infinity}
        ]]
    ]


expApplyRules[expr_, assum_] :=
    Assuming[assum,
        Simplify[canonicalizeModeOrder[expr] /. rulesNO,
            ExcludedForms -> _GeneralizedPower
        ]
    ]


ExpBosonNormalOrder::usage =
"\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) brings an expression built from exponentials of a formal field variable (e.g. \[FormalA]) and its \!\(\*RowBox[{StyleBox[\"var\", \"TI\"], SuperscriptBox[\"\[Dagger]\", \"\"]}]\) into normal order, using the boson disentangling identities.\n\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"assum\", \"TI\"]}], \"]\"}]\) uses \!\(\*StyleBox[\"assum\", \"TI\"]\) as additional simplifying assumptions.\n\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", RowBox[{\[Ellipsis], \",\", \"TimeConstraint->\", StyleBox[\"t\", \"TI\"]}], \"]\"}]\) limits simplification to \!\(\*StyleBox[\"t\", \"TI\"]\) seconds (default 30).";

Options[ExpBosonNormalOrder] = {
    TimeConstraint -> 30,
    Assumptions -> $Assumptions
};

ExpBosonNormalOrder[expr_, opts : OptionsPattern[]] := ExpBosonNormalOrder[expr, True, opts]

ExpBosonNormalOrder[expr_, assum_, opts : OptionsPattern[]] :=
    Block[{tmax = OptionValue[TimeConstraint], fullAssum},
        fullAssum = Union[Flatten[{assum}], Flatten[{OptionValue[Assumptions]}]];
        TimeConstrained[expApplyRules[expr, fullAssum], tmax, expr]
    ]
