(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["ExpBosonOrder"]
PackageExport["NormalOrdered"]


algA[v_]    := algA[v]    = NonCommutativeAlgebra[<|"Generators" -> {v}|>];
algAdag[v_] := algAdag[v] = NonCommutativeAlgebra[<|"Generators" -> {SuperDagger[v]}|>];


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
                Sqrt[Sec[omega]] Exp[Log[Sec[omega]] SuperDagger[v] ** v] **
                Exp[((alpha Tan[omega]) GeneralizedPower[NonCommutativeMultiply, v, 2])/omega]
        ],

    Exp[alpha_ v_?FormalSymbolQ] ** expr_ /; NonCommutativePolynomialQ[expr, algAdag[v]] :>
        NonCommutativeExpand[expr /. {SuperDagger[v] -> SuperDagger[v] + alpha}, algAdag[v]] ** Exp[alpha v],

    expr_ ** Exp[beta_ SuperDagger[v_?FormalSymbolQ]] /; NonCommutativePolynomialQ[expr, v, algA[v]] :>
        Exp[beta SuperDagger[v]] ** NonCommutativeExpand[expr /. {v -> v + beta}, algA[v]],

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
                Exp[-Log[Cosh[Omega]] SuperDagger[v1] ** v1] **
                Exp[Log[Cosh[Omega]] SuperDagger[v2] ** v2] **
                Exp[(beta Tanh[Omega]/Omega) v1 ** SuperDagger[v2]]
        ],

    Exp[alpha_ v1_?FormalSymbolQ ** v2_?FormalSymbolQ + beta_ SuperDagger[v1_?FormalSymbolQ] ** SuperDagger[v2_?FormalSymbolQ]] /; v1 =!= v2 :>
        With[{Omega = Sqrt[alpha beta]},
            Exp[(beta Tan[Omega]/Omega) SuperDagger[v1] ** SuperDagger[v2]] **
                Sec[Omega] Exp[Log[Sec[Omega]] (SuperDagger[v1] ** v1 + SuperDagger[v2] ** v2)] **
                Exp[(alpha Tan[Omega]/Omega) v1 ** v2]
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


NormalOrdered[expr_]["Series"] := expr /. {Exp[l_ SuperDagger[v_] ** v_] :> With[{k = Global`k}, Inactivate[Sum[
            (l^k/k!) GeneralizedPower[NonCommutativeMultiply, SuperDagger[v], k] **
                GeneralizedPower[NonCommutativeMultiply, v, k],
            {k, 0, Infinity}
        ], Sum]] }


expApplyRules[expr_, assum_] :=
    Assuming[assum,
        Simplify[canonicalizeModeOrder[expr] /. rulesNO,
            ExcludedForms -> _GeneralizedPower
        ]
    ]


ExpBosonOrder::usage =
"\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", StyleBox[\"expr\", \"TI\"], \"]\"}]\) brings an expression built from exponentials of a formal field variable (e.g. \[FormalA]) and its \!\(\*RowBox[{StyleBox[\"var\", \"TI\"], SuperscriptBox[\"\[Dagger]\", \"\"]}]\) into normal order, using the boson disentangling identities.\n\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"assum\", \"TI\"]}], \"]\"}]\) uses \!\(\*StyleBox[\"assum\", \"TI\"]\) as additional simplifying assumptions.\n\!\(\*RowBox[{\"ExpBosonNormalOrder\", \"[\", RowBox[{\[Ellipsis], \",\", \"TimeConstraint->\", StyleBox[\"t\", \"TI\"]}], \"]\"}]\) limits simplification to \!\(\*StyleBox[\"t\", \"TI\"]\) seconds (default 30).";

Options[ExpBosonOrder] = {
    TimeConstraint -> 30,
    Assumptions -> $Assumptions
};

ExpBosonOrder[expr_, opts : OptionsPattern[]] := ExpBosonOrder[expr, True, opts]

ExpBosonOrder[expr_, assum_, opts : OptionsPattern[]] :=
    Block[{tmax = OptionValue[TimeConstraint], fullAssum},
        fullAssum = Union[Flatten[{assum}], Flatten[{OptionValue[Assumptions]}]];
        TimeConstrained[expApplyRules[expr, fullAssum], tmax, expr]
    ]
