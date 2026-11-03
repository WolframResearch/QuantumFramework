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
"\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings expr into normal order, where vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra declares the creators and annihilators, with [\!\(\*SubscriptBox[\"a\", \"i\"]\), \!\(\*SubscriptBox[\"c\", \"i\"]\)] = 1.\n\!\(\*RowBox[{\"BosonicNormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the method: Automatic (default) or \"Blasiak\".";

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
"\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Brings the polynomial expr into anti-normal order, with every annihilator to the left of every creator, where vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra declares the ladder operators.\n\!\(\*RowBox[{\"BosonicAntinormalOrder\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Specifies the method: Automatic (default) or \"Blasiak\".";

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
"\!\(\*RowBox[{\"BosonicBCHTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Computes the nth-order term of the BCH series log(\!\(\*SuperscriptBox[\(e\), \(op1\)]\) \!\(\*SuperscriptBox[\(e\), \(op2\)]\)\[Ellipsis]) in normal order, for the ladder operators vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra.";

BosonicBCHTerms[ops_List, n_Integer, vars_] :=
    Replace[ResolveLadderVariables[BosonicBCHTerms, vars],
        lv_Association :> NonCommutativeExpand[ResourceFunction["BakerCampbellHausdorffTerms"][ops, n], lv["Algebra"]]
    ]


BosonicZassenhausTerms::usage =
"\!\(\*RowBox[{\"BosonicZassenhausTerms\", \"[\", RowBox[{StyleBox[\"ops\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Computes the nth term of the Zassenhaus formula in normal order, for the ladder operators vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra.";

BosonicZassenhausTerms[ops_List, n_Integer, vars_] :=
    Replace[ResolveLadderVariables[BosonicZassenhausTerms, vars],
        lv_Association :> NonCommutativeExpand[ResourceFunction["ZassenhausTerms"][ops, n], lv["Algebra"]]
    ]


BosonicVEV::usage =
"\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Computes the vacuum expectation value \[LeftAngleBracket]0\[VerticalSeparator]expr\[VerticalSeparator]0\[RightAngleBracket] for the ladder operators vars = {{\!\(\*SubscriptBox[\"c\", \"1\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[\"a\", \"1\"]\), \[Ellipsis]}} or a WeylAlgebra, by normal ordering expr and keeping its c-number part.\n\!\(\*RowBox[{\"BosonicVEV\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"Method->\", StyleBox[\"m\", \"TI\"]}], \"]\"}]\) Passes the method to BosonicNormalOrder.";

Options[BosonicVEV] = Options[BosonicNormalOrder];

vevCollapse[expr_] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, op_, p_] :> op^p,
    NonCommutativeMultiply -> Times
}

vevPolynomialQ[expr_, vars_] :=
    With[{syms = Table[Unique["vev"], Length[vars]]},
        PolynomialQ[vevCollapse[expr] /. Thread[vars -> syms], syms]
    ]

BosonicVEV[expr_, vars_, opts : OptionsPattern[]] :=
    With[{lv = ResolveLadderVariables[BosonicVEV, vars]},
        If[ lv === $Failed,
            $Failed,
            vevCollapse[BosonicNormalOrder[expr, vars, opts] /. $nonuls] /. Thread[lv["Generators"] -> 0]
        ] /; lv === $Failed || vevPolynomialQ[expr, lv["Generators"]]
    ]


BosonicMatrixElement::usage =
"\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]D(\[Alpha])\[VerticalSeparator]n\[RightAngleBracket] in closed form via associated Laguerre polynomials.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]S(\[Xi])\[VerticalSeparator]n\[RightAngleBracket] in closed form (zero when m+n is odd).\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"m\", \"TI\"], \",\", StyleBox[\"n\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Returns \[LeftAngleBracket]m\[VerticalSeparator]expr\[VerticalSeparator]n\[RightAngleBracket] for expr built from the ladder operators vars = {{c}, {a}} or a WeylAlgebra: a polynomial through its normal order and the closed form for \[LeftAngleBracket]m\[VerticalSeparator]\!\(\*SuperscriptBox[\"c\", \"p\"]\)\!\(\*SuperscriptBox[\"a\", \"q\"]\)\[VerticalSeparator]n\[RightAngleBracket], a product of exponentials through its disentangled form. The Fock indices may be symbolic.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"m\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"n\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"}\"}]}], \"}\"}], \",\", StyleBox[\"poly\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) Returns the element between multimode Fock states for a polynomial in the modes of vars, index k pairing with mode k.\n\!\(\*RowBox[{\"BosonicMatrixElement\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"\[Beta]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", \"\\\"Basis\\\"->\\\"Coherent\\\"\"}], \"]\"}]\) Returns \[LeftAngleBracket]\[Alpha]\[VerticalSeparator]expr\[VerticalSeparator]\[Beta]\[RightAngleBracket] between coherent states, the indices then being amplitudes.";

Options[BosonicMatrixElement] = {"Basis" -> "Fock"}

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
    Module[{x, y, k = Length[lv["Annihilators"]], vars, p},
        vars = Join[Array[x, k], Array[y, k]];
        p = Expand[vevCollapse[BosonicNormalOrder[expr, {lv["Creators"], lv["Annihilators"]}]] /. Thread[lv["Generators"] -> vars]];
        Which[
            p === 0, 0,
            k == 1,
                With[{t = Replace[p, {s_Plus :> List @@ s, s_ :> {s}}]},
                    (t /. Thread[vars -> 1]) . MapThread[element[First[ms], First[ns], #1, #2] &, Exponent[t, #] & /@ vars]
                ],
            True,
                With[{r = CoefficientRules[p, vars]},
                    With[{e = Transpose[Keys[r]]},
                        Values[r] . Times @@ Table[
                            With[{pq = Transpose[{e[[i]], e[[k + i]]}]},
                                With[{u = DeleteDuplicates[pq]},
                                    Lookup[AssociationThread[u, element[ms[[i]], ns[[i]], ##] & @@@ u], pq]
                                ]
                            ],
                            {i, k}
                        ]
                    ]
                ]
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
    Module[{lv = ResolveLadderVariables[BosonicMatrixElement, vars], ms = Flatten[{m}], ns = Flatten[{n}], pair, j, r},
        pair = If[AssociationQ[lv], {First[lv["Creators"]], First[lv["Annihilators"]]}];
        r = Which[
            lv === $Failed || ! FreeQ[expr, _QuantumOperator | _QuantumState], $Failed,
            Length[ms] === Length[ns] === Length[lv["Annihilators"]] && vevPolynomialQ[expr, lv["Generators"]],
                If[ coherentBasisQ[opts],
                    Times @@ MapThread[coherentOverlap, {ms, ns}] ladderSum[Conjugate[#1]^#3 #2^#4 &, {ms, ns}, expr, lv],
                    ladderSum[fockLadderElement, {ms, ns}, expr, lv]
                ],
            Length[lv["Annihilators"]] =!= 1 || ListQ[m] || ListQ[n], $Failed,
            coherentBasisQ[opts],
                Replace[numberDiagonal[expr, pair, j],
                    f : Except[$Failed] :> Exp[-(Abs[m]^2 + Abs[n]^2)/2] Sum[f (Conjugate[m] n)^j/j!, {j, 0, Infinity}]],
            True,
                Replace[normalProductElement[expr, pair, m, n, j],
                    $Failed /; ! FreeQ[expr, E] :>
                        normalProductElement[BosonicExpOrder[expr, {lv["Creators"], lv["Annihilators"]}], pair, m, n, j]]
        ];
        r /; lv === $Failed || r =!= $Failed
    ]


(* Normal ordering changes the function: :Exp[c n]: is (1 + c)^n on the diagonal. *)
numberDiagonal[expr_, {cr_, an_}, k_] :=
    With[{nn = cr ** an},
        With[{f = Replace[expr, {
                    NormalOrdered[Exp[c_. w_], _] /; w === nn :> (1 + c)^k,
                    e_ :> ReplaceAll[e, nn -> k]
                }]},
            If[FreeQ[f, cr | an | NormalOrdered], f, $Failed]
        ]
    ]

(* A wing is an exponential of one operator alone, so its element is a Taylor coefficient
   of Exp[P]: a displacement's linear wing and a squeeze's quadratic one share a formula,
   and the latter's parity rule is the vanishing of the odd coefficients. *)
wing[m_, k_, p_, x_] := SeriesCoefficient[Exp[p], {x, 0, m - k}] Sqrt[m!/k!]

wingPolynomial[s_, pair_, side_, x_] :=
    Module[{c, a},
        With[{p = s /. Thread[pair -> {c, a}] /. If[side === "Up", c, a] -> x /. NonCommutativeMultiply -> Times},
            If[FreeQ[p, c | a] && PolynomialQ[p, x] && TrueQ[(p /. x -> 0) == 0], p, $Failed]
        ]
    ]

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

(* scalar Exp[P[ad]] f[n] Exp[Q[a]], any factor absent.  Wings of a kind commute so their
   exponents add; between Fock states the middle sum runs only to Min[m, n]. *)
normalProductElement[expr_, pair_, m_, n_, j_] :=
    Module[{x, terms, scalar, tagged, u, w, f, r},
        terms = If[Head[expr] === Times, List @@ expr, {expr}];
        scalar = Times @@ Select[terms, FreeQ[#, Alternatives @@ pair] &];
        tagged = Replace[Times @@ Select[terms, ! FreeQ[#, Alternatives @@ pair] &],
            {q_NonCommutativeMultiply :> List @@ q, q_ :> {q}}];
        tagged = normalFactor[#, pair, j, x] & /@ tagged;
        If[MemberQ[tagged, $Failed], Return[$Failed, Module]];
        (* The factors are assembled as one normally ordered product, so they must arrive
           in that order: an antinormal product would otherwise be silently reordered and
           keep its own scalar, which is wrong by the exponential of the commutator. *)
        If[ ! OrderedQ[First /@ tagged /. {"Up" -> 1, "Diagonal" -> 2, "Down" -> 3}],
            Return[$Failed, Module]
        ];
        u = Total[Cases[tagged, ("Up" -> c_) :> c]];
        w = Total[Cases[tagged, ("Down" -> c_) :> c]];
        f = Times @@ Cases[tagged, ("Diagonal" -> g_) :> g];
        r = Which[
            u === 0 && w === 0, scalar (f /. j -> n) KroneckerDelta[m, n],
            w === 0, scalar wing[m, n, u, x] (f /. j -> n),
            u === 0, scalar (f /. j -> m) wing[n, m, w, x],
            IntegerQ[m] && IntegerQ[n],
                scalar Sum[wing[m, j, u, x] f wing[n, j, w, x], {j, 0, Min[m, n]}],
            True, $Failed
        ];
        If[FreeQ[r, SeriesCoefficient], r, $Failed]
    ]
