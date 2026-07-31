(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageScope["BlasiakNormalOrder"]

PackageScope["MultiModeBlasiakOrder"]

PackageScope["ParseBlasiakMonomial"]



ncPower[_, 0] := 1
ncPower[op_, 1] := op
ncPower[op_, p_] := GeneralizedPower[NonCommutativeMultiply, op, p]

blasiakBuildTerm[coeff_, 0, 0, b_, bdag_] := coeff
blasiakBuildTerm[coeff_, 0, q_, b_, bdag_] := coeff ncPower[b, q]
blasiakBuildTerm[coeff_, p_, 0, b_, bdag_] := coeff ncPower[bdag, p]
blasiakBuildTerm[coeff_, p_, q_, b_, bdag_] := coeff ncPower[bdag, p] ** ncPower[b, q]



blasiakExponentsQ[v_] := v =!= {} && VectorQ[v, IntegerQ[#] && NonNegative[#] &]


blasiakCoefficient[cumVec_, powers_, k_] :=
    (1 / k!) Sum[
        Binomial[k, j] (-1)^(k - j) Times @@ MapThread[FactorialPower, {Most[cumVec] + j, powers}],
        {j, 0, k}
    ]


BlasiakNormalOrder::badexponents =
    "expected equal-length vectors of non-negative integers, got `1` and `2`"

BlasiakNormalOrder[rList_ ? blasiakExponentsQ, sList_ ? blasiakExponentsQ, var_] /;
    Length[rList] === Length[sList] :=
With[{
    excess = Total[rList - sList],
    b = var,
    bdag = SuperDagger[var]
},
    If[ excess >= 0,
        With[{cumulative = FoldList[Plus, 0, rList - sList]},
            Sum[
                blasiakBuildTerm[
                    blasiakCoefficient[cumulative, sList, k],
                    excess + k, k, b, bdag
                ],
                {k, First[sList], Total[sList]}
            ]
        ],

        With[{
            cumulative = FoldList[Plus, 0, Reverse[sList] - Reverse[rList]],
            reversedPowers = Reverse[rList]
        },
            Sum[
                blasiakBuildTerm[
                    blasiakCoefficient[cumulative, reversedPowers, k],
                    k, k - excess, b, bdag
                ],
                {k, Last[rList], Total[rList]}
            ]
        ]
    ]
]

BlasiakNormalOrder[rList_, sList_, var_] :=
    (Message[BlasiakNormalOrder::badexponents, rList, sList]; $Failed)



ParseBlasiakMonomial::badfactor = "unsupported non-ladder factor: `1`"

blasiakFactorType[f_, var_, bdag_] := Replace[f, {
    x_ /; x === var  :> {"b", 1},
    x_ /; x === bdag :> {"bdag", 1},
    GeneralizedPower[NonCommutativeMultiply, x_, n_] /; x === var  :> {"b", n},
    GeneralizedPower[NonCommutativeMultiply, x_, n_] /; x === bdag :> {"bdag", n},
    Power[x_, n_] /; x === var  :> {"b", n},
    Power[x_, n_] /; x === bdag :> {"bdag", n},
    _ :> {"unknown", 0}
}]

ParseBlasiakMonomial[expr_, var_] := Module[{factors, typed, bad, blocks},

    factors = Replace[expr, {ncm_NonCommutativeMultiply :> List @@ ncm, e_ :> {e}}];
    typed = blasiakFactorType[#, var, SuperDagger[var]] & /@ factors;

    bad = FirstPosition[typed, {"unknown", _}, None, {1}];
    If[ bad =!= None,
        Message[ParseBlasiakMonomial::badfactor, factors[[First[bad]]]];
        Return[$Failed, Module]
    ];

    blocks = {#[[1, 1]], Total[#[[All, 2]]]} & /@ Split[typed, First[#1] === First[#2] &];

    (* padding to a creator-led, annihilator-closed word makes blocks alternate with even
       length, so odd and even positions pair up as r with s *)
    blocks = Join[
        If[blocks =!= {} && blocks[[1, 1]] === "b", {{"bdag", 0}}, {}],
        blocks,
        If[blocks =!= {} && blocks[[-1, 1]] === "bdag", {{"b", 0}}, {}]
    ];

    {Reverse[blocks[[1 ;; ;; 2, 2]]], Reverse[blocks[[2 ;; ;; 2, 2]]]}
]



blasiakModePolynomial[factors_List, {b_, bdag_}] :=
With[{modeFactors = Select[factors, ! FreeQ[#, b | bdag] &]},
    If[ modeFactors === {},
        1,
        Replace[
            ParseBlasiakMonomial[NonCommutativeMultiply @@ modeFactors, b],
            {
                {rList_, sList_} :> BlasiakNormalOrder[rList, sList, b],
                _ :> $Failed
            }
        ]
    ]
]


MultiModeBlasiakOrder[exprIn_, vars_List, scalars_List] :=
With[{altVars = Alternatives @@ vars},
Module[{expr, factors, scalarFactors, modePolys, result},

    expr = exprIn //. {
        GeneralizedPower[NonCommutativeMultiply,
            NonCommutativeMultiply[args__], n_Integer ? Positive] :>
                NonCommutativeMultiply @@ Flatten[ConstantArray[{args}, n]],
        Power[NonCommutativeMultiply[args__], n_Integer ? Positive] :>
            NonCommutativeMultiply @@ Flatten[ConstantArray[{args}, n]],
        NonCommutativeMultiply[a___, NonCommutativeMultiply[b__], c___] :>
            NonCommutativeMultiply[a, b, c]
    };

    expr = LiftNCScalars[expr, vars];

    factors = Catenate @ Map[
        Replace[{ncm_NonCommutativeMultiply :> List @@ ncm, f_ :> {f}}],
        Replace[expr, {p : (_Times | _NonCommutativeMultiply) :> List @@ p, e_ :> {e}}]
    ];

    {scalarFactors, factors} =
        Lookup[GroupBy[factors, FreeQ[#, altVars] &], {True, False}, {}];

    modePolys = blasiakModePolynomial[factors, #] & /@ Partition[vars, 2];

    If[ MemberQ[modePolys, $Failed], Return[$Failed, Module] ];

    result = NonCommutativeMultiply @@ DeleteCases[modePolys, 1];

    result = result //. {
        NonCommutativeMultiply[a___, b_ + c_, d___] :>
            NonCommutativeMultiply[a, b, d] + NonCommutativeMultiply[a, c, d],
        NonCommutativeMultiply[a___, Times[c_, b__], d___] /; FreeQ[c, altVars] :>
            c NonCommutativeMultiply[a, b, d],
        NonCommutativeMultiply[a___, c_, d___] /; FreeQ[c, altVars] :>
            c NonCommutativeMultiply[a, d]
    };

    result = result /. {
        NonCommutativeMultiply[x_] :> x,
        NonCommutativeMultiply[] -> 1
    };

    result = result /. NonCommutativeMultiply[args__] :>
        NonCommutativeMultiply @@ SortBy[{args}, {
            FreeQ[#, SuperDagger] &,
            First @ FirstPosition[vars, v_ /; ! FreeQ[#, v], {Infinity}, {1}] &
        }];

    (Times @@ scalarFactors) result
]]
