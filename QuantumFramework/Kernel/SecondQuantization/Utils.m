(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["$FockSize"]

PackageExport["SetFockSpaceSize"]

PackageExport["OperatorVariance"]

PackageExport["FieldVariables"]

PackageScope["ResolveLadderVariables"]

PackageExport["ToBosonicOperator"]

PackageScope["LiftNCScalars"]

PackageScope["$nonuls"]

PackageScope["cNumberQ"]

PackageScope["dequantize"]

PackageExport["G2Coherence"]

PackageExport["G1Correlation"]

PackageExport["CovarianceMatrix"]


$FockSize::usage = "Global variable holding the current Fock space truncation size (default 16).";


SetFockSpaceSize::usage = "\!\(\*RowBox[{\"SetFockSpaceSize\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) sets $FockSize to \!\(\*StyleBox[\"size\", \"TI\"]\). Default is 16.";

SetFockSpaceSize[size:_Integer?Positive:16]:= $FockSize = size;


OperatorVariance::usage =
"\!\(\*RowBox[{\"OperatorVariance\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", StyleBox[\"op\", \"TI\"]}], \"]\"}]\) computes \[LeftAngleBracket]\!\(\*SuperscriptBox[\"O\", \"2\"]\)\[RightAngleBracket] - \!\(\*SuperscriptBox[\"\[LeftAngleBracket]O\[RightAngleBracket]\", \"2\"]\) of \!\(\*StyleBox[\"op\", \"TI\"]\) in \!\(\*StyleBox[\"state\", \"TI\"]\).";

OperatorVariance[state_QuantumState, op_QuantumOperator] :=
	If[ state["VectorQ"],
		With[{opState = op @ state},
			(state["Dagger"] @ op @ opState)["Scalar"] - (state["Dagger"] @ opState)["Scalar"]^2
		],
		Tr[state["Operator"]@(op^2)]-Tr[state["Operator"]@op]^2
	]


General::ladvars = "`1` does not declare ladder operators; give {{c1, \[Ellipsis], cn}, {a1, \[Ellipsis], an}}, creators first, or a WeylAlgebra.";
General::ladflat = "The list `1` mixes creators and annihilators; give the annihilators alone, `2`, or {{creators}, {annihilators}}.";
General::ladpairs = "`1` lists one {creator, annihilator} pair per mode; write `2`, the creators of every mode first.";
General::laddagger = "The annihilator `1` carries a dagger; give {{c1, \[Ellipsis], cn}, {a1, \[Ellipsis], an}}, creators first.";
General::ladpair = "The creator `1` is paired with `2` instead of `3`.";
General::ladalg = "`1` is not a Weyl algebra.";
General::ladone = "`1` declares more than one mode; this operator acts on a single mode.";

daggerVariables[vars_] := With[{as = DeleteDuplicates[Flatten[vars] /. SuperDagger[v_] :> v]}, {SuperDagger /@ as, as}]

ResolveLadderVariables[f_, alg_NonCommutativeAlgebra] :=
    Replace[alg["Structure"], {
        {"Weyl", vars_} :> ResolveLadderVariables[f, vars],
        _ :> (Message[f::ladalg, alg]; $Failed)
    }]

ResolveLadderVariables[f_, vars : {Except[_List] ..}] :=
    If[ FreeQ[vars, SuperDagger],
        ResolveLadderVariables[f, {SuperDagger /@ vars, vars}],
        Message[f::ladflat, vars, First @ daggerVariables[vars] /. SuperDagger[v_] :> v]; $Failed
    ]

ResolveLadderVariables[f_, vars : {cs : {__}, as : {__}}] :=
    With[{cross = If[Length[cs] === Length[as],
            FirstCase[Transpose[{cs, as}], {SuperDagger[x_], a_} /; x =!= a && MemberQ[as, x] :> {SuperDagger[x], a, x}, None], None]},
        Which[
            ! FreeQ[vars, {___, x_, ___, SuperDagger[x_], ___} | {___, SuperDagger[x_], ___, x_, ___}],
                Message[f::ladpairs, vars, daggerVariables[vars]]; $Failed,
            MemberQ[as, _SuperDagger], Message[f::laddagger, FirstCase[as, _SuperDagger]]; $Failed,
            cross =!= None, Message[f::ladpair, Sequence @@ cross]; $Failed,
            True, Replace[Quiet[Check[WeylAlgebra[vars], $Failed]], {
                $Failed :> (Message[f::ladvars, vars]; $Failed),
                alg_ :> <|"Creators" -> cs, "Annihilators" -> as, "Generators" -> Join[cs, as], "Algebra" -> alg,
                    "Modes" -> AssociationThread[Join[cs, as], Join[#, #] & @ Range[Length[cs]]]|>
            }]
        ]
    ]

ResolveLadderVariables[f_, vars_] := (Message[f::ladvars, vars]; $Failed)


$nonuls = {0. -> 0, 0. I -> 0, Complex[0., 0.] -> 0, Complex[x_, 0.] :> x,
    Complex[0., y_] :> I y};

cNumberQ[expr_, vars_] := FreeQ[expr /. $nonuls, Alternatives @@ vars]

dequantize[expr_, rules_] := expr //. Join[rules, {GeneralizedPower[NonCommutativeMultiply, b_, p_] :> b^p, NonCommutativeMultiply -> Times}]


FieldVariables::usage =
"\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{}], \"]\"}]\) gives the ladder variables {{\!\(\*SuperscriptBox[\"\[FormalA]\", \"\[Dagger]\"]\)}, {\[FormalA]}}.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"var\", \"TI\"]}], \"]\"}]\) gives {{\!\(\*SuperscriptBox[StyleBox[\"var\", \"TI\"], \"\[Dagger]\"]\)}, {var}}.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"var\", \"TI\"], \"1\"], \",\", \"\[Ellipsis]\"}], \"}\"}]}], \"]\"}]\) gives {{\!\(\*SubsuperscriptBox[StyleBox[\"var\", \"TI\"], \"1\", \"\[Dagger]\"]\), \[Ellipsis]}, {\!\(\*SubscriptBox[StyleBox[\"var\", \"TI\"], \"1\"]\), \[Ellipsis]}}.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"var\", \"TI\"], \",\", StyleBox[\"labels\", \"TI\"]}], \"]\"}]\) gives one mode per label, named var with the label appended.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"labels\", \"TI\"]}], \"]\"}]\) uses \[FormalA].";

FieldVariables[] := FieldVariables[\[FormalA]]

FieldVariables[var_Symbol] := {{SuperDagger[var]}, {var}}

FieldVariables[vars : {__Symbol}] := {SuperDagger /@ vars, vars}

FieldVariables[labels_List] := FieldVariables[\[FormalA], labels]

FieldVariables[var_Symbol, labels_List] :=
    FieldVariables[Symbol[SymbolName[var] <> ToString[#]] & /@ labels]





LiftNCScalars[expr_, vars_List] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, Times[c_, op_], n_Integer ? Positive] /;
        FreeQ[c, Alternatives @@ vars] :>
            c^n GeneralizedPower[NonCommutativeMultiply, op, n],
    NonCommutativeMultiply[x___, Times[c_, y__], z___] /;
        FreeQ[c, Alternatives @@ vars] :>
            c NonCommutativeMultiply[x, y, z],
    NonCommutativeMultiply[x___, c_, z___] /;
        FreeQ[c, Alternatives @@ vars] && Length[{x, z}] > 0 :>
            c Replace[{x, z}, {{y_} :> y, ys_ :> NonCommutativeMultiply @@ ys}]
}


ToBosonicOperator::usage =
"\!\(\*RowBox[{\"ToBosonicOperator\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) realizes expr in the ladder operators vars as a truncated QuantumOperator.\n\!\(\*RowBox[{\"ToBosonicOperator\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) truncates each mode to dimension size.";

ToBosonicOperator[expr_] := (Message[ToBosonicOperator::argr, ToBosonicOperator, 2]; $Failed)


orderedPolynomial[(h : NormalOrdered | AntinormalOrdered)[f_, _], lv_, size_] :=
    Module[{x, t, k = Length[lv["Annihilators"]], xs, p, pow, block},
        xs = Array[x, 2 k];
        p = Normal @ Series[dequantize[f, Thread[lv["Generators"] -> t xs]], {t, 0, 2 k (size - 1)}] /. t -> 1;
        If[ ! PolynomialQ[p, xs], Return[$Failed, Module]];
        pow = NestList[SparseArray[Band[{1, 2}] -> Sqrt[Range[size - 1]], {size, size}] . # &, IdentityMatrix[size, SparseArray], size - 1];
        block = If[h === NormalOrdered, Transpose[pow[[#1 + 1]]] . pow[[#2 + 1]], pow[[#2 + 1]] . Transpose[pow[[#1 + 1]]]] &;
        QuantumOperator[
            Total[#2 Fold[KroneckerProduct, MapThread[block, TakeDrop[#1, k]]] & @@@ Select[CoefficientRules[p, xs], Max[First[#]] < size &]],
            Range[k], size]
    ]

ToBosonicOperator[expr_, vars_, size : _Integer?Positive : $FockSize] :=
    With[{lv = ResolveLadderVariables[ToBosonicOperator, vars]},
        If[ lv === $Failed,
            $Failed,
            With[{ops = Table[AnnihilationOperator[size, {k}], {k, Length[lv["Annihilators"]]}], g = Alternatives @@ lv["Generators"],
                    ex = expr /. o : (_NormalOrdered | _AntinormalOrdered) :> With[{p = orderedPolynomial[o, lv, size]}, p /; p =!= $Failed]},
                ex /. GeneralizedPower[NonCommutativeMultiply, b_, n_Integer] :> Power[b, n] /.
                    Power[base_, x_] /; FreeQ[base, g] && ! FreeQ[x, g] :> MatrixExp[Log[base] x] /. Join[
                    {NonCommutativeMultiply -> boseCompose},
                    Thread[lv["Creators"] -> (#["Dagger"] & /@ ops)],
                    Thread[lv["Annihilators"] -> ops]
                ]
            ]
        ]
    ]

boseCompose[] := 1

boseCompose[l___, c_, r___] /; FreeQ[c, _QuantumOperator | _boseCompose] := c boseCompose[l, r]

boseCompose[l___, HoldPattern[Times[x__]], r___] /; ! FreeQ[{x}, _QuantumOperator] := boseCompose[l, x, r]

boseCompose[ops__QuantumOperator] := Fold[#1 @ #2 &, {ops}]


G2Coherence::usage =
"\!\(\*RowBox[{\"G2Coherence\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"]}], \"]\"}]\) computes the second-order coherence function \!\(\*SuperscriptBox[\"g\", \"2\"]\)(0) for the given single-mode \!\(\*StyleBox[\"state\", \"TI\"]\).";

G2Coherence[\[Psi]_QuantumState] := Block[{aOp, nOp, a2Op, numerator, denominator},

    aOp  = AnnihilationOperator[\[Psi]["Dimension"]];
    
    nOp  = SuperDagger[aOp] @ aOp;
    
    a2Op = SuperDagger[aOp] @ SuperDagger[aOp] @ aOp @ aOp;
    
    If[ \[Psi]["VectorQ"],
        numerator   = (\[Psi]["Dagger"] @ a2Op @ \[Psi])["Scalar"];
        denominator = (\[Psi]["Dagger"] @ nOp  @ \[Psi])["Scalar"],

        numerator   = Tr[a2Op["Matrix"] . \[Psi]["DensityMatrix"]];
        denominator = Tr[nOp["Matrix"] . \[Psi]["DensityMatrix"]]
    ];
    numerator / denominator^2
]


G1Correlation::usage =
"\!\(\*RowBox[{\"G1Correlation\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", \"{{\", StyleBox[\"r1\", \"TI\"], \",\", StyleBox[\"t1\", \"TI\"], \"},{\", StyleBox[\"r2\", \"TI\"], \",\", StyleBox[\"t2\", \"TI\"], \"}}\"}], \"]\"}]\) computes the first-order correlation function \!\(\*SuperscriptBox[\"G\", \"1\"]\)(x1,x2) for the given single-mode state at space-time coordinates \!\(\*StyleBox[\"r1\", \"TI\"]\), \!\(\*StyleBox[\"t1\", \"TI\"]\) and \!\(\*StyleBox[\"r2\", \"TI\"]\), \!\(\*StyleBox[\"t2\", \"TI\"]\).";


G1Correlation[state_QuantumState, {{r1_, t1_}, {r2_, t2_}}] :=Block[{aOp, eMinus1, ePlus2},

    aOp    = AnnihilationOperator[state["Dimension"]];
    
    eMinus1 = -I SuperDagger[aOp] Exp[-I (\[FormalK] . r1 - \[FormalOmega] t1)];
    
    ePlus2  =  I aOp Exp[ I (\[FormalK] . r2 - \[FormalOmega] t2)];
    
    If[ state["PureStateQ"],
    
        (state["Dagger"] @ eMinus1 @ ePlus2 @ state)["Scalar"],
        
        Tr[(eMinus1 @ ePlus2)@state["Operator"]]
    ]
]


Options[CovarianceMatrix] = {"QuadratureScaling" -> 1/Sqrt[2], "ReturnMeans" -> False};

CovarianceMatrix::order = "Mode order `1` is not within the `2` mode(s) of the given state.";

CovarianceMatrix::norm = "The state has norm `1` rather than 1; moments were rescaled accordingly. In a truncated Fock space this usually means the state has leaked past the cutoff and the result is unreliable.";

CovarianceMatrix::usage =
"\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"]}], \"]\"}]\) computes the 2n\[Times]2n covariance matrix for all n modes of \!\(\*StyleBox[\"state\", \"TI\"]\) in the Serafini convention \!\(\*SubscriptBox[\(\[Sigma]\), \(vac\)]\) = 1/2 I.\n\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) restricts the computation to the modes specified by \!\(\*StyleBox[\"order\", \"TI\"]\).\n\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"\\\"QuadratureScaling\\\"->\", StyleBox[\"s\", \"TI\"]}], \"]\"}]\) uses \!\(\*OverscriptBox[\"X\", \"^\"]\) = s(a+\!\(\*SuperscriptBox[\"a\", \"\[Dagger]\"]\)), P = \[ImaginaryI] s(\!\(\*SuperscriptBox[\"a\", \"\[Dagger]\"]\)-a). Default s = 1/\!\(\*SqrtBox[\"2\"]\) (Serafini). Use s = 1/2 for \[HBar]=1/2 or s = 1 for Simon convention.\n\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"\\\"ReturnMeans\\\"->True\"}], \"]\"}]\) returns an Association with keys \"CovarianceMatrix\" and \"Means\", the latter being the vector of first moments \[LeftAngleBracket]\!\(\*SubscriptBox[OverscriptBox[\"R\", \"^\"], \"i\"]\)\[RightAngleBracket].";

CovarianceMatrix[state_QuantumState, opts : OptionsPattern[]] :=
    CovarianceMatrix[state, Range[state["Qudits"]], opts]

CovarianceMatrix[state_QuantumState, order_?orderQ, OptionsPattern[]] :=
    Block[{s = OptionValue["QuadratureScaling"], dims = state["Dimensions"], R, A, M, rho, norm, means, sigma},
        If[ ! AllTrue[order, 1 <= # <= Length[dims] &],
            Message[CovarianceMatrix::order, order, Length[dims]];
            Return[$Failed]
        ];
        R = Catenate[(2s * QuadratureOperators[dims[[#]], {#}]) & /@ order];
        If[ state["StateType"] === "Vector",
            A     = (# @ state)["StateVector"] & /@ R;
            norm  = Re[Conjugate[#] . #] & @ state["StateVector"];
            means = Re[Conjugate[A] . state["StateVector"]];
            sigma = Re[Conjugate[A] . Transpose[A]],
            rho   = state["Computational"];
            M     = #["Ordered", Range[Length[dims]], rho["Output"]]["Matrix"] & /@ R;
            rho   = rho["DensityMatrix"];
            A     = # . rho & /@ M;
            norm  = Re @ Tr @ rho;
            means = Re[Tr /@ A];
            sigma = Re @ Outer[Total[Transpose[#1] #2, 2] &, M, A, 1]
        ];
        If[ TrueQ[Abs[norm - 1] > 10^-8], Message[CovarianceMatrix::norm, norm]];
        means /= norm;
        sigma = sigma / norm - Outer[Times, means, means];
        If[ TrueQ[OptionValue["ReturnMeans"]],
            <|"CovarianceMatrix" -> sigma, "Means" -> means|>,
            sigma
        ]
    ]


SetFockSpaceSize[];
