(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageScope["MultiModeBlasiakOrder"]


General::ladpoly = "`1` is not a polynomial in the ladder operators."

(* A normal-ordered polynomial is <|{p1, \[Ellipsis], pn, q1, \[Ellipsis], qn} -> coefficient|>, and each mode of a product of two
   normal-ordered monomials reorders by a^q a\[Dagger]^r = \[Sum]_k Binomial[q, k] FactorialPower[r, k] a\[Dagger]^(r-k) a^(q-k). *)
blasiakProduct[n_][kp_ -> cp_, kq_ -> cq_] :=
    (Join[#[[All, 1]], #[[All, 2]]] -> cp cq Times @@ #[[All, 3]]) & /@ Tuples @ Table[
        With[{p = kp[[i]], q = kp[[n + i]], r = kq[[i]], s = kq[[n + i]]},
            Table[{p + r - k, q + s - k, Binomial[q, k] FactorialPower[r, k]}, {k, 0, Min[q, r]}]
        ],
        {i, n}
    ]

blasiakTimes[n_][pp_, qq_] := DeleteCases[Merge[Flatten @ Outer[blasiakProduct[n], Normal[pp], Normal[qq], 1], Total], 0]

blasiakPolynomial[lv_, e_] :=
    With[{gens = Alternatives @@ lv["Generators"], zero = ConstantArray[0, Length[lv["Generators"]]], n = Length[lv["Creators"]]},
        Replace[e, {
            s_ /; FreeQ[s, gens] :> <|zero -> s|>,
            g : gens :> <|Boole[SameQ[g, #]] & /@ lv["Generators"] -> 1|>,
            (GeneralizedPower[NonCommutativeMultiply, g : gens, k_Integer ? NonNegative] | Power[g : gens, k_Integer ? NonNegative]) :>
                <|k (Boole[SameQ[g, #]] & /@ lv["Generators"]) -> 1|>,
            (GeneralizedPower[NonCommutativeMultiply, x_, k_Integer ? NonNegative] | Power[x_, k_Integer ? NonNegative]) :>
                With[{px = blasiakPolynomial[lv, x]}, Nest[blasiakTimes[n][#, px] &, <|zero -> 1|>, k]],
            s_Plus :> Merge[blasiakPolynomial[lv, #] & /@ List @@ s, Total],
            m_NonCommutativeMultiply :> Fold[blasiakTimes[n], blasiakPolynomial[lv, #] & /@ List @@ m],
            t_Times :> With[{op = Select[List @@ t, ! FreeQ[#, gens] &]},
                If[Length[op] === 1, t / First[op] blasiakPolynomial[lv, First[op]], Throw[t, "ladpoly"]]
            ],
            f_ :> Throw[f, "ladpoly"]
        }]
    ]

blasiakWord[lv_, e_] :=
    Replace[
        MapThread[Replace[{#1, #2}, {{_, 0} -> Nothing, {g_, 1} :> g, {g_, k_} :> GeneralizedPower[NonCommutativeMultiply, g, k]}] &,
            {lv["Generators"], e}],
        {{} -> 1, {x_} :> x, w_ :> NonCommutativeMultiply @@ w}
    ]

MultiModeBlasiakOrder[expr_, lv_Association, f_] :=
    Catch[
        Total @ KeyValueMap[#2 blasiakWord[lv, #1] &, DeleteCases[Expand /@ blasiakPolynomial[lv, expr], 0]],
        "ladpoly",
        (Message[f::ladpoly, #1]; $Failed) &
    ]
