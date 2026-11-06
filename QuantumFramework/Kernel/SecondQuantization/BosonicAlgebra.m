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

ladderOrder[lv_, Automatic, _] :=
    expandCoefficients[
        NonCommutativeExpand[
            # //. GeneralizedPower[NonCommutativeMultiply, w_NonCommutativeMultiply, k_Integer ? Positive] :>
                NonCommutativeMultiply @@ Catenate[ConstantArray[List @@ w, k]],
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

Options[BosonicVEV] = Options[BosonicNormalOrder];

vevCollapse[expr_] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    NonCommutativeMultiply -> Times
}

vevPolynomialQ[expr_, vars_] :=
    With[{syms = Table[Unique["vev"], Length[vars]]},
        PolynomialQ[vevCollapse[expr] /. Thread[vars -> syms], syms]
    ]

BosonicVEV[exprs_List, vars_, opts : OptionsPattern[]] := BosonicVEV[#, vars, opts] & /@ exprs

BosonicVEV[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicVEV, vars]},
        If[ lv === $Failed,
            $Failed,
            vevCollapse[BosonicNormalOrder[expr, vars, opts] /. $nonuls] /. Thread[lv["Generators"] -> 0]
        ] /; lv === $Failed || vevPolynomialQ[expr, lv["Generators"]]
    ]

BosonicVEV[expr_, vars_, opts : OptionsPattern[]] :=
    With[{r = BosonicMatrixElement[ConstantArray[0, {2, Length[ResolveLadderVariables[BosonicVEV, vars]["Annihilators"]]}], expr, vars, ##] & @@ FilterRules[{opts}, Assumptions]},
        r /; Head[r] =!= BosonicMatrixElement
    ]


BosonicMatrixElement::usage =
"\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket].\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket].\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives \[LeftAngleBracket]m\[VerticalSeparator]expr\[VerticalSeparator]n\[RightAngleBracket] in the ladder operators vars.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"m\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"n\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"}\"}], \",\", StyleBox[\"poly\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives the element of poly between multimode Fock states.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"\[Beta]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Basis\\\"->\\\"Coherent\\\"\"}], \"]\"}]\) gives \[LeftAngleBracket]\[Alpha]\[VerticalSeparator]expr\[VerticalSeparator]\[Beta]\[RightAngleBracket] between coherent states.";

Options[BosonicMatrixElement] = {"Basis" -> "Fock", Assumptions :> $Assumptions}

SetAttributes[BosonicMatrixElement, HoldRest]

BosonicMatrixElement[{m_Integer?NonNegative, n_Integer?NonNegative}, (DisplacementOperator | SqueezeOperator)[_?PossibleZeroQ]] :=
    KroneckerDelta[m, n]

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

ladderSum[element_, {ms_, ns_}, expr_, lv_Association] :=
    Module[{x, y, k = Length[lv["Annihilators"]], vars, p, r, e},
        vars = Join[Array[x, k], Array[y, k]];
        p = Expand[vevCollapse[BosonicNormalOrder[expr, {lv["Creators"], lv["Annihilators"]}]] /. Thread[lv["Generators"] -> vars]];
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
                    FreeQ[d = expr //. Thread[MapThread[NonCommutativeMultiply, {lv["Creators"], lv["Annihilators"]}] -> ns],
                        Alternatives @@ lv["Generators"] | _NormalOrdered | _AntinormalOrdered],
                (d /. NonCommutativeMultiply -> Times) Times @@ MapThread[KroneckerDelta, {ms, ns}],
            If[ Length[lv["Annihilators"]] === 1, ListQ[m] || ListQ[n],
                ! MatrixQ[{ms, ns}, IntegerQ] || Length[ms] =!= Length[lv["Annihilators"]] || coherentBasisQ[opts]], $Failed,
            True,
                With[{el = If[Length[lv["Annihilators"]] === 1, normalProductElement[#, pair, m, n, j, coherentBasisQ[opts]] &, multimodeProductElement[#, lv, ms, ns] &]},
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
    Module[{c, a},
        With[{p = s /. Thread[pair -> {c, a}] /. If[side === "Up", c, a] -> x /. NonCommutativeMultiply -> Times},
            If[FreeQ[p, c | a] && PolynomialQ[p, x] && TrueQ[(p /. x -> 0) == 0], p, $Failed]
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

(* scalar Exp[P[ad]] ad^p f[n] a^q Exp[Q[a]], any factor absent.  Factors of a kind commute,
   so the up side is x^p Exp[P] and the down side x^q Exp[Q]; between Fock states the middle
   sum runs only to Min[m, n].  Without exponentials the element is read in closed form, and
   it vanishes for n < q before f[n - q] is evaluated, so a^dagger f[N] a keeps the vacuum
   even when f[n - 1] is singular at n = 0.  Between coherent states both sides are numbers
   and the middle is summed in t before the amplitudes go in, so a zero amplitude meets no 0^j.
   The factors must arrive in normal order: an antinormal product would otherwise be silently
   reordered and keep its own scalar, wrong by the exponential of the commutator. *)
normalProductElement[expr_, pair_, m_, n_, j_, coherent_] :=
    Module[{x, t, terms, scalar, tagged, u, w, p, q, fu, fd, f, mono, r},
        terms = If[Head[expr] === Times, List @@ expr, {expr}];
        scalar = Times @@ Select[terms, FreeQ[#, Alternatives @@ pair] &];
        tagged = Replace[Times @@ Select[terms, ! FreeQ[#, Alternatives @@ pair] &],
            {q_NonCommutativeMultiply :> List @@ q, q_ :> {q}}];
        tagged = normalFactor[#, pair, j, x] & /@ tagged;
        If[MemberQ[tagged, $Failed], Return[$Failed, Module]];
        If[ ! OrderedQ[First /@ tagged /. {"Up" | "Creators" -> 1, "Diagonal" -> 2, "Down" | "Annihilators" -> 3}],
            Return[$Failed, Module]
        ];
        u = Total[Cases[tagged, ("Up" -> c_) :> c]];
        w = Total[Cases[tagged, ("Down" -> c_) :> c]];
        p = Total[Cases[tagged, ("Creators" -> c_) :> c]];
        q = Total[Cases[tagged, ("Annihilators" -> c_) :> c]];
        fu = x^p Exp[u];
        fd = x^q Exp[w];
        f = Times @@ Cases[tagged, ("Diagonal" -> g_) :> g];
        mono = scalar KroneckerDelta[m - p, n - q] Sqrt[FactorialPower[m, p]] Sqrt[FactorialPower[n, q]] (f /. j -> n - q) &;
        r = Which[
            coherent,
                scalar (fu /. x -> Conjugate[m]) (fd /. x -> n) Exp[-(Abs[m]^2 + Abs[n]^2)/2] *
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

(* Creators commute among themselves, as do annihilators, so each multimode wing is a function of
   commuting variables and its Fock element a multivariate Taylor coefficient. *)
multimodeProductElement[expr_, lv_, ms_, ns_] :=
    Module[{x, y, k, xs, ys, ks, fs, side, up, f, dn},
        {xs, ys, ks} = Array[#, Length[ms]] & /@ {x, y, k};
        fs = Flatten[Replace[Replace[expr, {t_Times :> List @@ t, t_ :> {t}}], q_NonCommutativeMultiply :> List @@ q, {1}]] //.
            Join[Thread[MapThread[NonCommutativeMultiply, {lv["Creators"], lv["Annihilators"]}] -> ks], Thread[lv["Creators"] -> xs], Thread[lv["Annihilators"] -> ys]] /.
            {GeneralizedPower[NonCommutativeMultiply, b_, p_] :> b^p, NonCommutativeMultiply -> Times};
        side = Boole[Outer[! FreeQ[#1, #2] &, fs, Alternatives @@@ {xs, ks, ys}, 1]];
        If[ Max[Total[side, {2}]] > 1 || ! OrderedQ[DeleteCases[side = side . {1, 2, 3}, 0]], Return[$Failed, Module]];
        {up, f, dn} = Times @@ Pick[fs, side, #] & /@ {1, 2, 3};
        Times @@ Pick[fs, side, 0] Times @@ Sqrt[ms! ns!] Total[
            SeriesCoefficient[up dn, Sequence @@ Thread[{Join[xs, ys], 0, Join[ms - #, ns - #]}]] (f /. Thread[ks -> #]) / Times @@ (#!) & /@
                Tuples[Range[0, MapThread[Min, {ms, ns}]]]]
    ]
