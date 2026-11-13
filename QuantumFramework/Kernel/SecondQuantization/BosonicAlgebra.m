(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageExport["BosonicNormalOrder"]

PackageExport["BosonicAntinormalOrder"]

PackageExport["BosonicBCHTerms"]

PackageExport["BosonicZassenhausTerms"]

PackageExport["BosonicVEV"]

PackageExport["BosonicMatrixElement"]

PackageScope["vevPolynomialQ"]


BosonicNormalOrder::usage =
"\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"a\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"]\"}]\) brings expr into normal order in the annihilators \!\(\*SubscriptBox[StyleBox[\"a\", \"TI\"], \"i\"]\) and creators \!\(\*SubsuperscriptBox[StyleBox[\"a\", \"TI\"], \"i\", \"\[Dagger]\"]\).\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"c\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"a\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"}\"}]}], \"]\"}]\) takes the creators \!\(\*SubscriptBox[StyleBox[\"c\", \"TI\"], \"i\"]\).\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"alg\", \"TI\"]}], \"]\"}]\) uses the WeylAlgebra alg.";

BosonicNormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicNormalOrder] = {Method -> Automatic}

expandCoefficients[expr_, lv_] :=
    With[{mono = _NonCommutativeMultiply | GeneralizedPower[NonCommutativeMultiply, Alternatives @@ lv["Generators"], _] | Alternatives @@ lv["Generators"]},
        Total @ KeyValueMap[Expand[Total[#2]] #1 &,
            GroupBy[
                Replace[Replace[expr, {s_Plus :> List @@ s, s_ :> {s}}],
                    {m : mono :> {m, 1}, Times[c__, m : mono] :> {m, Times[c]}, c_ :> {1, c}}, {1}],
                First -> Last
            ]
        ]
    ]

(* NonCommutativeExpand refuses a power of a sum above "ExpandMaxExponent"; such a power is
   built one factor at a time, normal-ordering after each product.  Multiplying by the short
   base keeps every step cheap, where squaring would multiply two long polynomials. *)
stepPower[w_, k_, alg_] := Fold[NonCommutativeExpand[#1 ** w, alg] &, w, Range[k - 1]]

expandLargePowers[expr_, lv_] :=
    With[{cap = Lookup[Last @ First @ SystemOptions["NonCommutativeAlgebraOptions" -> "ExpandMaxExponent"], "ExpandMaxExponent"]},
        expr /. GeneralizedPower[NonCommutativeMultiply, w : Except[Alternatives @@ lv["Generators"]], k_Integer /; k > cap] :>
            stepPower[NonCommutativeExpand[expandLargePowers[w, lv], lv["Algebra"]], k, lv["Algebra"]]
    ]

ladderOrder[lv_, Automatic, _] :=
    expandCoefficients[
        NonCommutativeExpand[
            expandLargePowers[
                # //. GeneralizedPower[NonCommutativeMultiply, w_NonCommutativeMultiply, k_Integer ? Positive] :>
                    NonCommutativeMultiply @@ Catenate[ConstantArray[List @@ w, k]],
                lv
            ],
            lv["Algebra"]
        ],
        lv
    ] &

ladderOrder[lv_, "Blasiak", f_] := MultiModeBlasiakOrder[#, lv, f] &

BosonicNormalOrder[expr_, vars_, OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicNormalOrder, vars], method = OptionValue[Method]},
        Which[
            lv === $Failed, $Failed,
            ! MatchQ[method, Automatic | "Blasiak"], Message[BosonicNormalOrder::unknownMethod, method]; $Failed,
            True, ladderOrder[lv, method, BosonicNormalOrder][expr]
        ]
    ]


(* a -> a\[Dagger], a\[Dagger] -> -a is an automorphism ([a\[Dagger], -a] = 1) carrying a normal-ordered
   monomial to an anti-ordered one, so conjugating by it turns any normal-ordering
   method into an anti-normal one.  "Inverse" maps the result back. *)
ladderSwap[expr_, lv_Association, direction_String : "Forward"] :=
    LiftNCScalars[
        expr /. If[ direction === "Inverse",
            Join[Thread[lv["Creators"] -> lv["Annihilators"]], Thread[lv["Annihilators"] -> -lv["Creators"]]],
            Join[Thread[lv["Creators"] -> -lv["Annihilators"]], Thread[lv["Annihilators"] -> lv["Creators"]]]
        ],
        lv["Generators"]
    ]


BosonicAntinormalOrder::usage =
"\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) brings expr into anti-normal order in the ladder operators vars.";

BosonicAntinormalOrder::unknownMethod = "Unknown Method: `1`."

Options[BosonicAntinormalOrder] = Options[BosonicNormalOrder]

BosonicAntinormalOrder[expr_, vars_, OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicAntinormalOrder, vars], method = OptionValue[Method]},
        Which[
            lv === $Failed, $Failed,
            ! MatchQ[method, Automatic | "Blasiak"], Message[BosonicAntinormalOrder::unknownMethod, method]; $Failed,
            True, expandCoefficients[ladderSwap[ladderOrder[lv, method, BosonicAntinormalOrder][ladderSwap[expr, lv]], lv, "Inverse"], lv]
        ]
    ]


BosonicBCHTerms::usage =
"\!\(\*RowBox[{\"BosonicBCHTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives the nth-order term of log(\!\(\*SuperscriptBox[\(e\), \(op1\)]\) \!\(\*SuperscriptBox[\(e\), \(op2\)]\)\[Ellipsis]) in normal order.";

BosonicBCHTerms[ops_List, n_Integer, vars_] :=
    Replace[ResolveLadderVariables[BosonicBCHTerms, vars],
        lv_Association :> NonCommutativeExpand[ResourceFunction["BakerCampbellHausdorffTerms"][ops, n], lv["Algebra"]]
    ]


BosonicZassenhausTerms::usage =
"\!\(\*RowBox[{\"BosonicZassenhausTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives the nth term of the Zassenhaus formula in normal order.";

BosonicZassenhausTerms[ops_List, n_Integer, vars_] :=
    Replace[ResolveLadderVariables[BosonicZassenhausTerms, vars],
        lv_Association :> NonCommutativeExpand[ResourceFunction["ZassenhausTerms"][ops, n], lv["Algebra"]]
    ]


BosonicVEV::usage =
"\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives the vacuum expectation value \[LeftAngleBracket]0\[VerticalSeparator]expr\[VerticalSeparator]0\[RightAngleBracket].";

Options[BosonicVEV] = {Assumptions :> $Assumptions};

vevPolynomialQ[expr_, vars_] :=
    With[{syms = Table[Unique["vev"], Length[vars]]},
        PolynomialQ[dequantize[expr, Thread[vars -> syms]], syms]
    ]

BosonicVEV[exprs_List, vars_, opts : OptionsPattern[]] := BosonicVEV[#, vars, opts] & /@ exprs

BosonicVEV[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicVEV, vars]},
        With[{r = If[AssociationQ[lv], BosonicMatrixElement[ConstantArray[0, {2, Length[lv["Annihilators"]]}], expr, vars, opts], $Failed]},
            r /; Head[r] =!= BosonicMatrixElement
        ]
    ]


BosonicMatrixElement::usage =
"\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket].\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket].\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]expr\[VerticalSeparator]n\[RightAngleBracket] in the ladder operators vars.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"m\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"n\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"}\"}], \",\", StyleBox[\"poly\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives the element of poly between multimode Fock states.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"\[Beta]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Basis\\\"->\\\"Coherent\\\"\"}], \"]\"}]\) gives \[LeftAngleBracket]\[Alpha]\[VerticalSeparator]expr\[VerticalSeparator]\[Beta]\[RightAngleBracket] between coherent states.";

Options[BosonicMatrixElement] = {"Basis" -> "Fock", Assumptions :> $Assumptions}

SetAttributes[BosonicMatrixElement, HoldRest]

fockIndexQ[x_] := ! ListQ[x] && (! NumericQ[x] || IntegerQ[x] && x >= 0)

SetAttributes[ladderBranch, HoldAll]

ladderBranch[up_, down_, {m_, n_}, d_, c_] := With[{a = {m, n} \[Element] NonNegativeIntegers},
    Which[
        NumericQ[d], If[d >= 0, up, down],
        TrueQ[Refine[d >= 0, a]], Refine[Piecewise[{{up, c}}, 0], a],
        TrueQ[Refine[d <= 0, a]], Refine[Piecewise[{{down, c}}, 0], a],
        True, Refine[Piecewise[{{up, d >= 0 && c}, {down, d < 0 && c}}, 0], a]
    ]
]

BosonicMatrixElement[{m_ ? fockIndexQ, n_ ? fockIndexQ}, (DisplacementOperator | SqueezeOperator)[_ ? PossibleZeroQ | {_ ? PossibleZeroQ, _}]] :=
    KroneckerDelta[m, n]

BosonicMatrixElement[{m_ ? fockIndexQ, n_ ? fockIndexQ}, DisplacementOperator[\[Alpha]_]] :=
    With[{a2 = Conjugate[\[Alpha]] \[Alpha], d = Simplify[m - n]},
        E^(-a2/2) ladderBranch[
            Sqrt[n!/m!] \[Alpha]^d LaguerreL[n, d, a2],
            Sqrt[m!/n!] (-Conjugate[\[Alpha]])^(-d) LaguerreL[m, -d, a2],
            {m, n}, d, True
        ]
    ]

BosonicMatrixElement[{m_ ? fockIndexQ, n_ ? fockIndexQ}, SqueezeOperator[\[Xi]_]] :=
    With[{r = If[ListQ[\[Xi]], First[\[Xi]], Abs[\[Xi]]], \[Phi] = If[ListQ[\[Xi]], Last[\[Xi]], Arg[\[Xi]]], d = Simplify[m - n]},
        If[ NumericQ[d] && OddQ[d], 0,
            ladderBranch[
                (-E^(I \[Phi]) Tanh[r]/2)^(d/2)/(d/2)! Sqrt[m!/n!] Sech[r]^(n + 1/2) Hypergeometric2F1[-n/2, (1 - n)/2, d/2 + 1, -Sinh[r]^2],
                (E^(-I \[Phi]) Tanh[r]/2)^(-d/2)/(-d/2)! Sqrt[n!/m!] Sech[r]^(m + 1/2) Hypergeometric2F1[-m/2, (1 - m)/2, 1 - d/2, -Sinh[r]^2],
                {m, n}, d, Mod[d, 2] == 0
            ]
        ]
    ]


fockLadderElement[m_, n_, p_, q_] :=
    KroneckerDelta[m - p, n - q] Sqrt[FactorialPower[m, p]] Sqrt[FactorialPower[n, q]]

coherentOverlap[\[Alpha]_, \[Beta]_] :=
    Exp[Conjugate[\[Alpha]] \[Beta] - (Conjugate[\[Alpha]] \[Alpha] + Conjugate[\[Beta]] \[Beta])/2]

ladderSum[element_, {ms_, ns_}, expr_, lv_Association] :=
    Module[{x, y, k = Length[lv["Annihilators"]], vars, p, r, e},
        vars = Join[Array[x, k], Array[y, k]];
        p = Expand[dequantize[BosonicNormalOrder[expr, {lv["Creators"], lv["Annihilators"]}] /. $nonuls, Thread[lv["Generators"] -> vars]]];
        If[ p === 0, 0,
            r = CoefficientRules[p, vars];
            e = Transpose[Keys[r]];
            Values[r] . Times @@ Table[MapThread[element[ms[[i]], ns[[i]], ##] &, {e[[i]], e[[k + i]]}], {i, k}]
        ]
    ]

coherentBasisQ[opts___] :=
    OptionValue[BosonicMatrixElement, {opts}, "Basis"] === "Coherent"

BosonicMatrixElement[{{m_}, {n_}}, rest__] := BosonicMatrixElement[{m, n}, rest]

BosonicMatrixElement[{m_, n_}, c_ ? NumericQ, opts : OptionsPattern[]] /;
        Length[Flatten[{m}]] === Length[Flatten[{n}]] :=
    c Times @@ MapThread[
        If[coherentBasisQ[opts], coherentOverlap[##], KroneckerDelta[##]] &,
        {Flatten[{m}], Flatten[{n}]}
    ]

BosonicMatrixElement[{_, _}, Except[_DisplacementOperator | _SqueezeOperator | _ ? NumericQ], OptionsPattern[]] :=
    (Message[BosonicMatrixElement::argrx, BosonicMatrixElement, 2, 3]; $Failed)

BosonicMatrixElement[{m_, n_}, expr_, vars : Except[_Rule | _RuleDelayed], opts : OptionsPattern[]] :=
    Module[{lv = ResolveLadderVariables[BosonicMatrixElement, vars], ms = Flatten[{m}], ns = Flatten[{n}], pair, j, d, r},
        pair = If[AssociationQ[lv], {First[lv["Creators"]], First[lv["Annihilators"]]}];
        r = Which[
            lv === $Failed || ! FreeQ[expr, _QuantumOperator | _QuantumState], $Failed,
            Length[ms] === Length[ns] === Length[lv["Annihilators"]] && vevPolynomialQ[expr, lv["Generators"]],
                If[ coherentBasisQ[opts],
                    Times @@ MapThread[coherentOverlap, {ms, ns}] ladderSum[ZeroLimitPower[Conjugate[#1], #3] ZeroLimitPower[#2, #4] &, {ms, ns}, expr, lv],
                    ladderSum[fockLadderElement, {ms, ns}, expr, lv]
                ],
            Length[ms] === Length[ns] === Length[lv["Annihilators"]] && Length[ms] > 1 && ! coherentBasisQ[opts] &&
                    FreeQ[d = dequantize[expr, Thread[MapThread[NonCommutativeMultiply, {lv["Creators"], lv["Annihilators"]}] -> Array[j, Length[ns]]]],
                        Alternatives @@ lv["Generators"] | _NormalOrdered | _AntinormalOrdered],
                (d /. Thread[Array[j, Length[ns]] -> ns]) Times @@ MapThread[KroneckerDelta, {ms, ns}],
            If[ Length[lv["Annihilators"]] === 1, ListQ[m] || ListQ[n],
                Length[ms] =!= Length[lv["Annihilators"]] || ! coherentBasisQ[opts] && ! MatrixQ[{ms, ns}, IntegerQ]], $Failed,
            True,
                With[{el = If[Length[lv["Annihilators"]] === 1, normalProductElement[#, pair, m, n, j, coherentBasisQ[opts]] &,
                        symbolElement[#, lv, ms, ns, coherentBasisQ[opts]] &]},
                    Replace[sumProductElement[el, expr, lv["Generators"]],
                        $Failed /; ! FreeQ[expr, E] :>
                            sumProductElement[el, BosonicExpOrder[expr, {lv["Creators"], lv["Annihilators"]}, Assumptions -> OptionValue[Assumptions]], lv["Generators"]]]
                ]
        ];
        r /; lv === $Failed || r =!= $Failed
    ]


(* Normal ordering changes the function: :Exp[c n]: is (1 + c)^n on the diagonal.  A base
   that vanishes leaves only the vacuum, so 0^n with n >= 0 is KroneckerDelta[n, 0] rather
   than 0^0 once n is set to 0; the factor 1 + b keeps a machine zero inexact. *)
numberDiagonal[expr_, {cr_, an_}, k_] :=
    With[{nn = cr ** an},
        With[{f = Replace[expr, {
                    NormalOrdered[Exp[c_. w_], _] /; w === nn :> (1 + c)^k,
                    e_ :> ReplaceAll[e, nn -> k]
                }] /. {
                    Power[b_ /; TrueQ[b == 0], p_ /; TrueQ[Simplify[p >= 0, k >= 0]]] :> (1 + b) KroneckerDelta[p, 0],
                    Power[b_ /; TrueQ[b == 0], p_ /; ! FreeQ[p, k]] :> Piecewise[{{1 + b, p == 0}}, b^p]
                }},
            If[FreeQ[f, cr | an | NormalOrdered], f, $Failed]
        ]
    ]

(* A wing is a function of one operator alone, x^p Exp[P], so its element is a Taylor
   coefficient: a displacement's linear wing and a squeeze's quadratic one share a formula,
   and the latter's parity rule is the vanishing of the odd coefficients. *)
wing[m_, k_, g_, x_] := SeriesCoefficient[g, {x, 0, m - k}] Sqrt[m!/k!]

wingPolynomial[s_, pair_, side_, x_] :=
    Module[{o},
        With[{p = dequantize[s, Thread[pair -> If[side === "Up", {x, o}, {o, x}]]]},
            If[FreeQ[p, o] && PolynomialQ[p, x] && TrueQ[(p /. x -> 0) == 0], p, $Failed]
        ]
    ]

normalFactor[cr_, {cr_, _}, _, _] := "Creators" -> 1

normalFactor[GeneralizedPower[NonCommutativeMultiply, cr_, p_Integer ? Positive], {cr_, _}, _, _] := "Creators" -> p

normalFactor[an_, {_, an_}, _, _] := "Annihilators" -> 1

normalFactor[GeneralizedPower[NonCommutativeMultiply, an_, q_Integer ? Positive], {_, an_}, _, _] := "Annihilators" -> q

normalFactor[fac_, pair_, k_, x_] :=
    With[{s = Replace[fac, {Exp[e_] :> e, _ :> None}]},
        With[{
            up = If[s === None, $Failed, wingPolynomial[s, pair, "Up", x]],
            dn = If[s === None, $Failed, wingPolynomial[s, pair, "Down", x]]
        },
            Which[
                up =!= $Failed, "Up" -> up,
                dn =!= $Failed, "Down" -> dn,
                True, Replace[numberDiagonal[fac, pair, k], f : Except[$Failed] :> "Diagonal" -> f]
            ]
        ]
    ]

(* scalar Exp[P[ad]] ad^p f[n] a^q Exp[Q[a]], any factor absent: the up side is x^p Exp[P], the down side x^q Exp[Q],
   and the middle sum runs to Min[m, n].  The closed form vanishes for n < q before f[n - q] is evaluated, so
   a^dagger f[N] a keeps the vacuum when f[n - 1] is singular at n = 0; coherent elements sum the middle before the
   amplitudes go in, so a zero amplitude meets no 0^j.  Factors must arrive in normal order: an antinormal product
   would be silently reordered, wrong by the exponential of the commutator. *)
normalProductElement[expr_, pair_, m_, n_, j_, coherent_] :=
    Module[{x, t, scalar, tagged, kinds, u, w, p, q, fu, fd, f, mono, r},
        {scalar, tagged} = Times @@@ Lookup[GroupBy[Replace[expr, {e_Times :> List @@ e, e_ :> {e}}], FreeQ[Alternatives @@ pair]], {True, False}, {}];
        tagged = normalFactor[#, pair, j, x] & /@ Replace[tagged, {e_NonCommutativeMultiply :> List @@ e, e_ :> {e}}];
        If[ MemberQ[tagged, $Failed] || ! OrderedQ[Keys[tagged] /. {"Up" | "Creators" -> 1, "Diagonal" -> 2, "Down" | "Annihilators" -> 3}],
            Return[$Failed, Module]
        ];
        kinds = Merge[tagged, Identity];
        {u, w, p, q} = Total /@ Lookup[kinds, {"Up", "Down", "Creators", "Annihilators"}, {}];
        fu = x^p Exp[u];
        fd = x^q Exp[w];
        f = Times @@ Lookup[kinds, "Diagonal", {}];
        mono = scalar KroneckerDelta[m - p, n - q] Sqrt[FactorialPower[m, p]] Sqrt[FactorialPower[n, q]] (f /. j -> n - q) &;
        r = Which[
            coherent,
                scalar (fu /. x -> Conjugate[m]) (fd /. x -> n) Exp[-(Conjugate[m] m + Conjugate[n] n)/2] *
                    (Sum[f t^j/j!, {j, 0, Infinity}] /. t -> Conjugate[m] n),
            fu === 1 && fd === 1, scalar (f /. j -> n) KroneckerDelta[m, n],
            u === 0 && w === 0,
                If[ IntegerQ[m] && IntegerQ[n],
                    If[m >= p && n >= q, mono[], 0],
                    Piecewise[{{mono[], m >= p && n >= q}}, 0]
                ],
            fd === 1, scalar wing[m, n, fu, x] (f /. j -> n),
            fu === 1, scalar (f /. j -> m) wing[n, m, fd, x],
            IntegerQ[m] && IntegerQ[n],
                scalar Sum[wing[m, j, fu, x] f wing[n, j, fd, x], {j, 0, Min[m, n]}],
            True, $Failed
        ];
        If[FreeQ[r, SeriesCoefficient], r, $Failed]
    ]

sumProductElement[f_, expr_, gens_] :=
    With[{t = Expand[expr /. p_NonCommutativeMultiply :> Distribute[p]] //.
            NonCommutativeMultiply[l___, s_ x_, r___] /; FreeQ[s, Alternatives @@ gens] :> s NonCommutativeMultiply[l, x, r]},
        If[ Head[t] === Plus,
            With[{e = f /@ List @@ t}, If[MemberQ[e, $Failed], $Failed, Total[e]]],
            f[t]
        ]
    ]


symbolElement[expr_, lv_, ms_, ns_, coherent_] :=
    Module[{x, y, k, val, xs, ys, ks, toXY, factors, num, middles},
        {xs, ys, ks} = Array[#, Length[ms]] & /@ {x, y, k};
        toXY = dequantize[#, Join[Thread[MapThread[NonCommutativeMultiply, {lv["Creators"], lv["Annihilators"]}] -> ks],
            Thread[lv["Creators"] -> xs], Thread[lv["Annihilators"] -> ys]]] &;
        factors = Replace[FixedPoint[Flatten[Replace[#, p : _Times | _NonCommutativeMultiply :> List @@ p, {1}]] &, {expr}], {
            NormalOrdered[f_, ___] :> {2, toXY[f] /. Thread[ks -> xs ys]},
            f_ :> With[{c = toXY[f]}, {Replace[Union[Cases[c, (h : x | k | y)[_] :> h, {0, Infinity}]], {{} -> 0, {x} -> 1, {k} -> 2, {y} -> 3, _ -> 4}], c}]
        }, {1}];
        num = Times @@ Select[factors[[All, 2]], ! FreeQ[#, k] &];
        middles = Append[Cases[factors, {2, f_} /; FreeQ[f, k] :> f], num];
        If[ MemberQ[factors[[All, 1]], 4] || ! OrderedQ[DeleteCases[factors[[All, 1]], 0]] ||
                ! DuplicateFreeQ[Catenate[Union[Cases[#, (x | y | k)[i_] :> i, {0, Infinity}]] & /@ middles]],
            $Failed,
            With[{
                s = Times @@ Select[factors[[All, 2]], FreeQ[#, k] &] Exp[-xs . ys] *
                    Sum[val @@ ks Times @@ ((xs ys)^ks / ks!), Evaluate[Sequence @@ Thread[{ks, 0, If[coherent, Infinity, MapThread[Min, {ms, ns}]]}]]],
                vals = val[kk__] :> (num /. Thread[ks -> {kk}])},
                If[ coherent,
                    Times @@ MapThread[coherentOverlap, {ms, ns}] (s /. vals /. Join[Thread[xs -> Conjugate[ms]], Thread[ys -> ns]]),
                    (D[Exp[xs . ys] s, Sequence @@ Thread[{Join[xs, ys], Join[ms, ns]}]] /. Thread[Join[xs, ys] -> 0] /. vals) / Sqrt[Times @@ (ms! ns!)]
                ]
            ]
        ]
    ]
