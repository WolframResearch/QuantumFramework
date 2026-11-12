(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicExpOrder"]

PackageExport["NormalOrdered"]

PackageExport["AntinormalOrdered"]

PackageExport["BosonicExpSimplify"]

PackageExport["BosonicExpBraid"]

PackageScope["ncPolynomial"]


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


(* e^(w N + p ad^2 + q a^2 + mu ad + nu a) from the corner D = Cos[l] - w Sinc[l], l = Sqrt[4 p q - w^2], of
   exp[{{w, 2 p}, {-2 q, -w}}]; every coefficient but the root is regular in l, so l = 0 is its Series limit.
   The root D^(-1/2) winds in oscillating evolutions, so it is e^(-f/2) (e^(-f) D)^(-1/2) with the phase f of rootPhase.
   Antinormal order is the normal one under a -> ad, ad -> -a: (w, p, q, mu, nu) -> (-w, q, p, -nu, mu). *)
quadraticExp[g_, {cr_, an_}, lv_, s_] :=
    Module[{x, y, l, f, coeffs, c, w, p, q, mu, nu, t, m, u, v, e, r},
        coeffs = Association @ CoefficientRules[Expand @ dequantize[#, Thread[{cr, an} -> {x, y}]], {x, y}] &;
        c = coeffs[g];
        If[ ! SubsetQ[{{0, 0}, {1, 1}, {2, 0}, {0, 2}, {1, 0}, {0, 1}}, Keys[c]] || ! FreeQ[Values[c], Alternatives @@ lv["Generators"]], Return[$Failed, Module]];
        c = coeffs[BosonicNormalOrder[g, {{cr}, {an}}] /. $nonuls];
        If[ AnyTrue[{{{1, 1}}, {{2, 0}, {1, 0}}, {{0, 2}, {0, 1}}, {{1, 0}, {0, 1}}}, SubsetQ[Append[#, {0, 0}], Keys[c]] &], Return[$Failed, Module]];
        {w, p, q, mu, nu} = {s, 1, 1, s, 1} Lookup[c, Key /@ If[s == 1, Identity, Map[Reverse]] @ {{1, 1}, {2, 0}, {0, 2}, {1, 0}, {0, 1}}, 0];
        {t, m, u, v, e, r} = Which[
            p === q === 0, {0, w, mu (E^w - 1)/w, nu (E^w - 1)/w, mu nu (E^w - 1 - w)/w^2, w/2},
            PossibleZeroQ[4 p q - w^2], Coefficient[Normal @ Series[# /. f -> 0, {l, 0, 0}], l, 0],
            True, # /. {l -> Sqrt[4 p q - w^2], f -> With[{x = {{w, 2 p}, {-2 q, -w}}},
                If[MatrixQ[x / SelectFirst[Flatten[x], ! PossibleZeroQ[#] &], NumericQ], rootPhase[x],
                    Piecewise[{{I Sign[Re[I w]] Sqrt[4 p q - w^2], Im[4 p q - w^2] == 0 && Re[4 p q - w^2] > 0}}, 0]]]} /. $nonuls
        ] & @ Join[
            {Tan[l]/(l - w Tan[l]), Log[Sec[l] l/(l - w Tan[l])]},
            {mu Sinc[l] - (w mu - 2 p nu) (1 - Cos[l])/l^2,
             nu Sinc[l] - (w nu - 2 q mu) (1 - Cos[l])/l^2,
             ((q mu^2 + p nu^2) (-2 w - (l^2 - 2 w) Cos[l] + l (1 + w) Sin[l]) + mu nu (l^2 + 2 w^2 + (l^2 (w - 1) - 2 w^2) Cos[l] - l w (1 + w) Sin[l]))/l^4
            } / (Cos[l] - w Sinc[l]),
            {(-f - Log[E^-f (Cos[l] - w Sinc[l])])/2}];
        Exp[Lookup[c, Key[{0, 0}], 0] - s w/2 + e] (ncJoin[
            Exp[p t GeneralizedPower[NonCommutativeMultiply, #1, 2] + s u #1],
            Exp[If[w === 0, m (1/2 + s #1 ** #2), r + m s #1 ** #2]],
            Exp[q t GeneralizedPower[NonCommutativeMultiply, #2, 2] + v #2]
        ] & @@ If[s == 1, {cr, an}, {an, cr}])
    ]

ncPolynomial[poly_, vars_, ops_] := Total[#2 ncJoin @@ MapThread[ladderPower, {ops, #1}] & @@@ CoefficientRules[poly, vars]]

splitModes[g_, lv_] :=
    With[{parts = KeySort @ GroupBy[List @@ g, ladderModes[lv, #] &, Total]},
        If[ AllTrue[Keys[parts], Length[#] < 2 &] && Count[Keys[parts], {_}] > 1,
            Exp[Lookup[parts, Key[{}], 0]] ncJoin @@ Exp /@ Values @ KeySelect[parts, # =!= {} &],
            $Failed
        ]
    ]

(* Det[D] = Sum c_S e^(tau lambda_S) over n-subsets S of the eigenvalues of x = tau x0 (Cauchy-Binet).  If the largest |c_S|
   outweighs the rest and the others only rotate, pulling tau lambda_S out of the root keeps its principal branch continuous. *)
rootPhase[x_] := rootPhase[x, SelectFirst[Flatten[x], ! PossibleZeroQ[#] &]]

rootPhase[x_, tau_] /; MatrixQ[x / tau, NumericQ] :=
    Module[{n = Length[x] / 2, lam, v, terms, best},
        {lam, v} = MapAt[N @* Transpose, Eigensystem[x / tau], 2];
        (
            terms = {#[[1, 1]], Abs @ Total[#[[All, 2]]]} & /@ GatherBy[
                With[{vi = Inverse[v]}, {Total[lam[[#]]], Det[v[[n + 1 ;;, #]]] Det[vi[[#, n + 1 ;;]]]} & /@ Subsets[Range[2 n], {n}]],
                Round[N @ First[#], 10^-8] &
            ];
            best = First @ MaximalBy[terms, Last];
            If[ 2 Last[best] > Total[terms[[All, 2]]] && AllTrue[terms[[All, 1]] - First[best], PossibleZeroQ @ Chop @ ComplexExpand @ Re[tau #] &],
                tau First[best],
                0
            ]
        ) /; MatrixRank[v] == 2 n
    ]

rootPhase[__] := 0

adjointMatrix[w_, p_, q_] := ArrayFlatten[{{w, -p}, {q, -Transpose[w]}}]

affineShift[m_, ell_, n_] /; AllTrue[ell, PossibleZeroQ] := {ConstantArray[0, n], ConstantArray[0, n], 0}

affineShift[m_, ell_, n_] :=
    With[{aug = MatrixExp[ArrayFlatten[{
            {Transpose[m], ArrayFlatten[{{0, -IdentityMatrix[n]}, {IdentityMatrix[n], 0}}], ConstantArray[0, {2 n, 1}]},
            {ConstantArray[0, {2 n, 2 n}], -m, Transpose[{ell}]},
            {ConstantArray[0, {1, 2 n}], ConstantArray[0, {1, 2 n}], {{0}}}}]]},
        Simplify @ Append[TakeDrop[aug[[2 n + 1 ;; 4 n, -1]], n], (aug[[2 n + 1 ;; 4 n, 2 n + 1 ;; 4 n]] . ell) . aug[[;; 2 n, -1]] / 2]
    ]

(* e^(ad P ad/2 + ad W a + a Q a/2 + mu ad + nu a) over coupled modes, from the blocks {{A, B}, {C, D}} of exp[{{W, P}, {-Q, -W^T}}]:
   e^(-Tr W/2) Det[D]^(-1/2) e^(ad B D^-1 ad/2) :e^(ad (D^-T - 1) a): e^(-a D^-1 C a/2), which for P = Q = 0 is
   :e^(ad (e^W - 1) a):; linear terms come in through affineShift and braid to the wings.  Antinormal order is the same
   under a -> ad, ad -> -a.  An exact numeric matrix is exponentiated as z m/tau, z -> tau: MatrixExp of the number form
   expands nested radicals twentyfold, the symbolic one keeps them in RootSum. *)
multimodeQuadraticExp[g_, lv_, s_] :=
    Module[{n = Length[ladderModes[lv, g]], x, y, z, xs, ys, cr, an, op, middle, mexp, poly, w, p, q, lam, rho, kap, m, k, e, d, yy, l},
        {xs, ys} = {Array[x, n], Array[y, n]};
        mexp[mm_] := If[ MatrixQ[mm, NumericQ] && Precision[mm] === Infinity,
            With[{tau = SelectFirst[Flatten[mm], ! PossibleZeroQ[#] &, 1]}, MatrixExp[z mm / tau] /. z -> tau],
            MatrixExp[mm]
        ];
        {cr, an} = If[s == 1, Identity, Reverse] @ ladderVars[lv][[All, ladderModes[lv, g]]];
        op = ncPolynomial[Expand[# /. Thread[ys -> s ys]], Join[xs, ys], Join[cr, an]] &;
        middle = If[ DiagonalMatrixQ[#],
            Exp[op[Log[Diagonal[#]] . (xs ys)]],
            If[s == 1, NormalOrdered, AntinormalOrdered][Exp[op[xs . (# - IdentityMatrix[n]) . ys]], ladderVars[lv]]
        ] &;
        poly = Expand @ dequantize[If[s == 1, BosonicNormalOrder, BosonicAntinormalOrder][g, ladderVars[lv]] /. $nonuls, Thread[Join[cr, an] -> Join[xs, s ys]]];
        If[ n < 2 || FreeQ[poly, Alternatives @@ xs] || FreeQ[poly, Alternatives @@ ys] ||
            ! PolynomialQ[poly, Join[xs, ys]] || ! SubsetQ[{0, 1, 2}, Total /@ Keys @ CoefficientRules[poly, Join[xs, ys]]], Return[$Failed, Module]];
        {w, p, q} = {D[poly, {xs}, {ys}], D[poly, {xs, 2}], D[poly, {ys, 2}]};
        m = ArrayFlatten[{{w, p}, {-q, -Transpose[w]}}];
        {lam, rho, kap} = affineShift[adjointMatrix[w, p, q], Flatten[{D[poly, {xs}], D[poly, {ys}]} /. Thread[Join[xs, ys] -> 0]], n];
        Exp[(poly /. Thread[Join[xs, ys] -> 0]) + kap + lam . rho / 2] If[ AllTrue[Flatten[{p, q}], PossibleZeroQ],
            k = mexp[w];
            ncJoin[Exp[op[Simplify[k . lam] . xs]], middle[k], Exp[op[rho . ys]]],
            e = mexp[m];
            d = Inverse[e[[n + 1 ;;, n + 1 ;;]]];
            yy = -d . e[[n + 1 ;;, ;; n]];
            l = rootPhase[m];
            E^(Simplify[lam . yy . lam] / 2 - (Tr[w] + l) / 2) (E^-l Det[e[[n + 1 ;;, n + 1 ;;]]])^(-1/2) ncJoin[
                Exp[op[xs . e[[;; n, n + 1 ;;]] . d . xs / 2 + Simplify[Transpose[d] . lam] . xs]],
                middle[Transpose[d]],
                Exp[op[ys . yy . ys / 2 + Simplify[yy . lam + rho] . ys]]
            ]
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

    Exp[alpha_. an] ** NormalOrdered[Exp[k_. + c_. cr ** an], vv] /; FreeQ[k, cr | an] :>
        NormalOrdered[Exp[k + c cr ** an], vv] ** Exp[alpha (1 + c) an],

    NormalOrdered[Exp[k_. + c_. cr ** an], vv] ** Exp[beta_. cr] /; FreeQ[k, cr | an] :>
        Exp[beta (1 + c) cr] ** NormalOrdered[Exp[k + c cr ** an], vv],

    Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2])] **
            Exp[beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] :>
        Exp[(beta GeneralizedPower[NonCommutativeMultiply, cr, 2])/(1 - 4 alpha beta)] **
            Exp[-Log[1 - 4 alpha beta] cr ** an] **
            Exp[(alpha GeneralizedPower[NonCommutativeMultiply, an, 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[g_Plus] :> With[{d = quadraticExp[g, {cr, an}, lv, 1]}, d /; d =!= $Failed],

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
        ],

    Exp[w1_. cr1 ** an1 + w2_. cr2 ** an2 + alpha_. cr1 ** cr2 + beta_. an1 ** an2 + c_.] /;
            FreeQ[{w1, w2, alpha, beta, c}, cr1 | an1 | cr2 | an2] && ! (PossibleZeroQ[alpha] && PossibleZeroQ[beta]) :>
        With[{s = (w1 + w2)/2, d = (w1 - w2)/2, mu = Sqrt[alpha beta - (w1 + w2)^2/4]},
            With[{e = Cos[mu] - s Sinc[mu]},
                E^(c - s)/e (Exp[(alpha Sinc[mu]/e) cr1 ** cr2] **
                    Exp[(d - Log[e]) cr1 ** an1 - (d + Log[e]) cr2 ** an2] **
                    Exp[(beta Sinc[mu]/e) an1 ** an2])
            ]
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

    Exp[beta_. cr] ** AntinormalOrdered[Exp[k_. + c_. an ** cr], vv] /; FreeQ[k, cr | an] :>
        AntinormalOrdered[Exp[k + c an ** cr], vv] ** Exp[beta (1 - c) cr],

    AntinormalOrdered[Exp[k_. + c_. an ** cr], vv] ** Exp[alpha_. an] /; FreeQ[k, cr | an] :>
        Exp[alpha (1 - c) an] ** AntinormalOrdered[Exp[k + c an ** cr], vv],

    Exp[beta_. (cr ** cr | GeneralizedPower[NonCommutativeMultiply, cr, 2])] **
            Exp[alpha_. (an ** an | GeneralizedPower[NonCommutativeMultiply, an, 2])] :>
        Exp[(alpha GeneralizedPower[NonCommutativeMultiply, an, 2])/(1 - 4 alpha beta)] **
            Exp[Log[1 - 4 alpha beta] an ** cr] **
            Exp[(beta GeneralizedPower[NonCommutativeMultiply, cr, 2])/(1 - 4 alpha beta)] / Sqrt[1 - 4 alpha beta],

    Exp[g_Plus] :> With[{d = quadraticExp[g, {cr, an}, lv, -1]}, d /; d =!= $Failed],

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
        ],

    Exp[w1_. cr1 ** an1 + w2_. cr2 ** an2 + alpha_. cr1 ** cr2 + beta_. an1 ** an2 + c_.] /;
            FreeQ[{w1, w2, alpha, beta, c}, cr1 | an1 | cr2 | an2] && ! (PossibleZeroQ[alpha] && PossibleZeroQ[beta]) :>
        With[{s = (w1 + w2)/2, d = (w1 - w2)/2, mu = Sqrt[alpha beta - (w1 + w2)^2/4]},
            With[{e = Cos[mu] + s Sinc[mu]},
                E^(c - s)/e (Exp[(beta Sinc[mu]/e) an1 ** an2] **
                    Exp[(d + Log[e]) an1 ** cr1 + (Log[e] - d) an2 ** cr2] **
                    Exp[(alpha Sinc[mu]/e) cr1 ** cr2])
            ]
        ]
}


expRules[lv_, "Normal"] := expRules[lv, "Normal"] = Join[
    Flatten @ Transpose[modeRulesNO[#, ladderVars[lv], lv] & /@ ladderPairs[lv]],
    Flatten[pairRulesNO @@@ Permutations[ladderPairs[lv], {2}]],
    {Exp[g_Plus] :> With[{d = splitModes[g, lv]}, d /; d =!= $Failed],
     Exp[g_Plus] :> With[{d = multimodeQuadraticExp[g, lv, 1]}, d /; d =!= $Failed]},
    {modeSwap[lv, ladderRank[lv]], modeMerge[lv]}
]

expRules[lv_, "Antinormal"] := expRules[lv, "Antinormal"] = Join[
    Flatten @ Transpose[modeRulesAN[#, ladderVars[lv], lv] & /@ ladderPairs[lv]],
    Flatten[pairRulesAN @@@ Permutations[ladderPairs[lv], {2}]],
    {Exp[g_Plus] :> With[{d = splitModes[g, lv]}, d /; d =!= $Failed],
     Exp[g_Plus] :> With[{d = multimodeQuadraticExp[g, lv, -1]}, d /; d =!= $Failed]},
    {modeSwap[lv, 4 - ladderRank[lv][#] &], modeMerge[lv]}
]


exponentTerm[t_, g_] :=
    With[{f = Replace[t, {u_Times :> List @@ u, u_ :> {u}}]},
        Times @@ Select[f, ! FreeQ[#, g] &] -> Times @@ Select[f, FreeQ[#, g] &]
    ]

canonicalizeGenerators[expr_, lv_] :=
    With[{g = Alternatives @@ lv["Generators"]},
        (expr /. Exp[e_] :> Exp[Total @ KeyValueMap[#1 #2 &,
            GroupBy[exponentTerm[#, g] & /@ Replace[Expand[e], {s_Plus :> List @@ s, s_ :> {s}}], First -> Last, Total]]]) //.
            NonCommutativeMultiply[l___, x : g, y : g, r___] /; lv["Modes"][x] > lv["Modes"][y] :>
                NonCommutativeMultiply[l, y, x, r]
    ]


ladderPower[op_, 0] := 1
ladderPower[op_, p_] := GeneralizedPower[NonCommutativeMultiply, op, p]

ladderMonomial[adag_, m_, a_, k_] := ncJoin[ladderPower[adag, m], ladderPower[a, k]]

orderedMonomial[NormalOrdered, {cr_, an_}, {m_, k_}] := ladderMonomial[cr, m, an, k]
orderedMonomial[AntinormalOrdered, {cr_, an_}, {m_, k_}] := ladderMonomial[an, k, cr, m]

MakeBoxes[SeriesIndex[s_String], StandardForm | TraditionalForm] ^:= InterpretationBox[s, SeriesIndex[s]]

expSeriesTerm[ord_, coeff_, pair_] := With[{k = SeriesIndex["k"]},
    Inactivate[Sum[
        (coeff^k/k!) orderedMonomial[ord, pair, {k, k}],
        {k, 0, Infinity}
    ], Sum]
]



NormalOrdered::usage =
"\!\(\*RowBox[{\"NormalOrdered\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) represents the normal-ordered product :expr: in the ladder operators vars.";

AntinormalOrdered::usage =
"\!\(\*RowBox[{\"AntinormalOrdered\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) represents the anti-normal-ordered product of expr in the ladder operators vars.";

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


ladderMode[expr_, lv_] :=
    Replace[Select[ladderPairs[lv], ! FreeQ[expr, Alternatives @@ #] &], {{p_} :> p, _ :> None}]


orderedSeries[ord_, expr_, lambda_, lv_] := Module[{v, x, y, n, coefN},
    v = ladderMode[expr, lv];
    If[v === None, Return[$Failed, Module]];
    n = SeriesIndex["n"];
    coefN = SeriesCoefficient[dequantize[expr, Thread[v -> {x, y}]], {lambda, 0, n}];

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
        Normal @ Series[dequantize[expr, Thread[v -> {x, y}]], {lambda, 0, order}] /. {x -> First[v], y -> Last[v]},
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
        CoefficientRules[Expand[dequantize[expr, Thread[pair -> {x, y}]]], {x, y}]
    ]
]


ladderPolynomialQ[expr_, pair_] :=
    Module[{x, y},
        With[{d = dequantize[expr, Thread[pair -> {x, y}]]},
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


expApplyRules[expr_, assum_, rules_, lv_, tc_] :=
    Assuming[assum,
        Replace[
            TimeConstrained[
                FixedPoint[
                    (# /. Prepend[rules, n : (_NormalOrdered | _AntinormalOrdered) :> n]) //.
                        (NonCommutativeMultiply[l___, Times[s_, m_NonCommutativeMultiply], r___] /; FreeQ[s, Alternatives @@ lv["Generators"]] :>
                            s NonCommutativeMultiply[l, m, r]) &,
                    canonicalizeGenerators[expr, lv]
                ],
                tc
            ],
            r : Except[$Aborted] :> TimeConstrained[Simplify[r, ExcludedForms -> _GeneralizedPower | _Piecewise, Trig -> False], tc, r]
        ]
    ]


BosonicExpOrder::usage =
"\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) brings a product of exponentials of the ladder operators vars into normal order.\n\!\(\*RowBox[{\"BosonicExpOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Ordering\\\"->\", \"\\\"Antinormal\\\"\"}], \"]\"}]\) brings it into anti-normal order.";

BosonicExpOrder::ordering = "`1` is not a valid setting for \"Ordering\"; use \"Normal\" or \"Antinormal\"."

BosonicExpOrder::time = "Ordering stopped after the TimeConstraint of `1` seconds; the input is returned unordered."

BosonicExpOrder::assum ="`1` reads as assumptions rather than ladder operators; give the ladder operators second and the assumptions as Assumptions -> `1`."

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
                    Replace[expApplyRules[expr, OptionValue[Assumptions], expRules[lv, ord], lv, OptionValue[TimeConstraint]],
                        $Aborted :> (Message[BosonicExpOrder::time, OptionValue[TimeConstraint]]; expr)],
                other_ :> (Message[BosonicExpOrder::ordering, other]; $Failed)
            }]
        ]
    ]


BosonicExpSimplify::usage =
"\!\(\*RowBox[{\"BosonicExpSimplify\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) simplifies products of exponentials and functions of number operators in the ladder operators vars.";

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
                True,                                           Replace[adRotate[op, cmt, Simplify[cmt2 / op], lv], $Failed :> linearConjugate[b, op, lv]]
            ]
        ]
    ]

linearConjugate[b_, op_, lv_] :=
    Module[{modes = ladderModes[lv, {b, op}], n, x, y, xs, ys, cr, an, zero, pb, po, w, p, q, e},
        n = Length[modes];
        {xs, ys} = {Array[x, n], Array[y, n]};
        {cr, an} = ladderVars[lv][[All, modes]];
        zero = Thread[Join[xs, ys] -> 0];
        {pb, po} = Expand @ dequantize[BosonicNormalOrder[#, ladderVars[lv]] /. $nonuls, Thread[Join[cr, an] -> Join[xs, ys]]] & /@ {b, op};
        If[ ! AllTrue[{pb, po}, PolynomialQ[#, Join[xs, ys]] &] || Max[Total /@ Keys @ CoefficientRules[pb, Join[xs, ys]]] > 2 ||
            Max[Total /@ Keys @ CoefficientRules[po, Join[xs, ys]]] > 1, Return[$Failed, Module]];
        {w, p, q} = {D[pb, {xs}, {ys}], D[pb, {xs, 2}], D[pb, {ys, 2}]};
        e = If[ AllTrue[Flatten[{w, {D[pb, {xs}], D[pb, {ys}]} /. zero}], PossibleZeroQ],
            bogoliubovExp[p, q],
            $Failed
        ];
        If[ FreeQ[e, $Failed | Indeterminate | ComplexInfinity | MatrixFunction],
            Join[cr, an] . (e . (Join[D[po, {xs}], D[po, {ys}]] /. zero)) + (po /. zero),
            e = MatrixExp[ArrayFlatten[{
                    {adjointMatrix[w, p, q], ConstantArray[0, {2 n, 1}]},
                    {{Join[D[pb, {ys}], -D[pb, {xs}]] /. zero}, {{0}}}}]] .
                Append[Join[D[po, {xs}], D[po, {ys}]] /. zero, po /. zero];
            Simplify[Join[cr, an] . Most[e] + Last[e]]
        ]
    ]

bogoliubovExp[p_, q_] :=
    With[{c = MatrixFunction[Cosh[Sqrt[#]] &, #] &, s = MatrixFunction[If[PossibleZeroQ[#], 1, Sinh[Sqrt[#]] / Sqrt[#]] &, #] &},
        ArrayFlatten[{{c[-p . q], -p . s[-q . p]}, {q . s[-p . q], c[-q . p]}}]
    ]

expConjugate[b_, Exp[g_], lv_] /; ncPolyQ[{b, g}, lv] :=
    Replace[expConjugate[b, g, lv], h : Except[$Failed] :> Exp[h]]

expConjugate[b_, op_, lv_] /; ncPolyQ[{b, op}, lv] :=
    Replace[adAction[b, op, lv], $Failed :> conjugateParts[b, op, lv]]

expConjugate[___] := $Failed

conjugateParts[b_, op_Plus, lv_] :=
    With[{gs = expConjugate[b, #, lv] & /@ List @@ op}, If[MemberQ[gs, $Failed], $Failed, Total[gs]]]

conjugateParts[b_, op_Times, lv_] :=
    With[{s = Select[op, FreeQ[#, Alternatives @@ lv["Generators"]] &]},
        If[s === 1, $Failed, Replace[expConjugate[b, op / s, lv], g : Except[$Failed] :> s g]]
    ]

conjugateParts[b_, op_NonCommutativeMultiply, lv_] :=
    With[{gs = expConjugate[b, #, lv] & /@ List @@ op}, If[MemberQ[gs, $Failed], $Failed, NonCommutativeMultiply @@ gs]]

conjugateParts[b_, GeneralizedPower[NonCommutativeMultiply, x_, n_Integer ? Positive], lv_] :=
    Replace[expConjugate[b, x, lv], g : Except[$Failed] :> GeneralizedPower[NonCommutativeMultiply, g, n]]

conjugateParts[___] := $Failed

numberFunctionQ[f_, {cr_, an_}] := Module[{k}, ! FreeQ[f, cr ** an] && FreeQ[f /. cr ** an -> k, cr | an]]

balancedWordQ[w_List, {cr_, an_}] :=
    MatchQ[w, {(cr | an | GeneralizedPower[NonCommutativeMultiply, cr | an, _Integer ? Positive]) ..}] &&
        Total[getNCPower[#, cr] - getNCPower[#, an] & /@ w] == 0

numberPolynomial[w_List, {cr_, an_}] :=
    BosonicNormalOrder[NonCommutativeMultiply @@ w, {{cr}, {an}}] /.
        GeneralizedPower[NonCommutativeMultiply, cr, k_] ** GeneralizedPower[NonCommutativeMultiply, an, k_] :>
            Product[cr ** an - j, {j, 0, k - 1}]

(* a^q f(N) = f(N + q) a^q and f(N) ad^p = ad^p f(N + p), and functions of one N commute.
   A word with as many ad as a beside a function of N is itself a polynomial in N, ad^k a^k
   = N (N - 1) ... (N - k + 1). Ladder powers only move outward, never into ad f(N) a = N f(N - 1),
   whose cancellation would lose the vacuum. *)
shiftRules[{cr_, an_}] := {
    NonCommutativeMultiply[l___, f_, w__, r___] /; numberFunctionQ[f, {cr, an}] && balancedWordQ[{w}, {cr, an}] :>
        ncJoin[l, Together[f numberPolynomial[{w}, {cr, an}]], r],
    NonCommutativeMultiply[l___, w__, f_, r___] /; numberFunctionQ[f, {cr, an}] && balancedWordQ[{w}, {cr, an}] :>
        ncJoin[l, Together[numberPolynomial[{w}, {cr, an}] f], r],
    NonCommutativeMultiply[l___, f_, g_, r___] /; numberFunctionQ[f, {cr, an}] && numberFunctionQ[g, {cr, an}] :>
        ncJoin[l, Together[f g], r],
    NonCommutativeMultiply[l___, x : (an | GeneralizedPower[NonCommutativeMultiply, an, _Integer ? Positive]), f_, r___] /;
            numberFunctionQ[f, {cr, an}] :>
        ncJoin[l, f /. cr ** an -> cr ** an + getNCPower[x, an], x, r],
    NonCommutativeMultiply[l___, f_, x : (cr | GeneralizedPower[NonCommutativeMultiply, cr, _Integer ? Positive]), r___] /;
            numberFunctionQ[f, {cr, an}] :>
        ncJoin[l, x, f /. cr ** an -> cr ** an + getNCPower[x, cr], r]
}

expReduceRules[lv_] := Join[{
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
            ncJoin[l, Exp[c / 2] Exp[Expand[b1 + b2]], r] /; cNumberQ[c, lv["Generators"]]
        ],

    NonCommutativeMultiply[l___, op_, Exp[b_], r___] /; ncPolyQ[{op, b}, lv] :>
        With[{cmt = ncCommutator[b, op, lv]},
            With[{cmt2 = ncCommutator[b, cmt, lv]},
                ncJoin[l, Exp[b], op - cmt, r] /; zeroQ[cmt2]
            ]
        ]
}, Catenate[shiftRules /@ ladderPairs[lv]]]

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

adFamilyRules[cr_, {cre_, ann_}, _] /; SubsetQ[{{2, 0}, {0, 2}, {1, 1}, {1, 0}, {0, 1}}, Keys[cr]] :=
    With[{w = Lookup[cr, Key[{1, 1}], 0], p = Lookup[cr, Key[{2, 0}], 0], q = Lookup[cr, Key[{0, 2}], 0],
            mu = Lookup[cr, Key[{1, 0}], 0], nu = Lookup[cr, Key[{0, 1}], 0]},
        With[{lam = Sqrt[4 p q - w^2]},
            {cre -> Cos[lam] cre + Sinc[lam] (w cre + 2 q ann + nu) + (w nu - 2 q mu) Sinc[lam/2]^2/2,
             ann -> Cos[lam] ann - Sinc[lam] (2 p cre + w ann + mu) + (w mu - 2 p nu) Sinc[lam/2]^2/2}
        ]
    ]

adFamilyRules[___] := $Failed

modeAdRules[gen_, pair_, gens_] :=
    Module[{x, y},
        adFamilyRules[
            KeyDrop[
                Association @ CoefficientRules[
                    Expand @ dequantize[gen, Thread[pair -> {x, y}]], {x, y}],
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
"\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) moves the left exponential of an adjacent pair past the right one, \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) \[Rule] \!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{\"A\", \",\", \"B\"}], \"]\"}]]\)\!\(\*SuperscriptBox[\"e\", \"A\"]\).\n\!\(\*RowBox[{\"BosonicExpBraid\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Direction\\\"->\", \"\\\"Left\\\"\"}], \"]\"}]\) moves the right one left, \!\(\*SuperscriptBox[\"e\", \"A\"]\)\!\(\*SuperscriptBox[\"e\", \"B\"]\) \[Rule] \!\(\*SuperscriptBox[\"e\", \"B\"]\)\!\(\*SuperscriptBox[\"e\", RowBox[{\"Ad\", \"[\", RowBox[{RowBox[{\"-\", \"B\"}], \",\", \"A\"}], \"]\"}]]\).";

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
