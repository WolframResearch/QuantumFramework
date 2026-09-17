(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["$FockSize"]

PackageExport["SetFockSpaceSize"]

PackageExport["OperatorVariance"]

PackageExport["FieldVariables"]

PackageScope["OrderVariables"]

PackageExport["ToBosonicOperator"]

PackageScope["ExtractNCVars"]

PackageScope["LiftNCScalars"]

PackageScope["FormalSymbolQ"]

PackageScope["$nonuls"]

PackageScope["cNumberQ"]

PackageExport["G2Coherence"]

PackageExport["G1Correlation"]

PackageExport["CovarianceMatrix"]


$FockSize::usage = "Global variable holding the current Fock space truncation size (default 16).";


SetFockSpaceSize::usage = "\!\(\*RowBox[{\"SetFockSpaceSize\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) sets $FockSize to \!\(\*StyleBox[\"size\", \"TI\"]\). Default is 16.";

SetFockSpaceSize[size:_Integer?Positive:16]:= $FockSize = size;


OperatorVariance::usage =
"\!\(\*RowBox[{\"OperatorVariance\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", StyleBox[\"op\", \"TI\"]}], \"]\"}]\) computes \[LeftAngleBracket]\!\(\*SuperscriptBox[\"O\", \"2\"]\)\[RightAngleBracket] - \!\(\*SuperscriptBox[\"\[LeftAngleBracket]O\[RightAngleBracket]\", \"2\"]\) of \!\(\*StyleBox[\"op\", \"TI\"]\) in \!\(\*StyleBox[\"state\", \"TI\"]\).";

OperatorVariance[state_QuantumState, op_QuantumOperator]:= 
	Tr[state["Operator"]@(op^2)]-Tr[state["Operator"]@op]^2


FormalSymbolQ[var_Symbol]:=StringStartsQ[ToString[FullForm[var]],"\\[Formal"];

FormalSymbolQ[_]:=False;


$nonuls = {0. -> 0, 0. I -> 0, Complex[0., 0.] -> 0, Complex[x_, 0.] :> x,
    Complex[0., y_] :> I y};

cNumberQ[expr_, vars_] := FreeQ[expr /. $nonuls, Alternatives @@ vars]


FieldVariables::usage =
"\!\(\*RowBox[{\"FieldVariables\", \"[\", \"]\"}]\) returns the default field variables {\[FormalA], \!\(\*SuperscriptBox[\"\[FormalA]\", \"\[Dagger]\"])\)}.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"var\", \"TI\"]}], \"]\"}]\) returns {var, \!\(\*SuperscriptBox[StyleBox[\"var\", \"TI\"], \"\[Dagger]\"]\)} for a formal symbol \!\(\*StyleBox[\"var\", \"TI\"]\).\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"labels\", \"TI\"]}], \"]\"}]\) returns field variables indexed by \!\(\*StyleBox[\"labels\", \"TI\"]\), using \[FormalA] as base.\n\!\(\*RowBox[{\"FieldVariables\", \"[\", RowBox[{StyleBox[\"var\", \"TI\"], \",\", StyleBox[\"labels\", \"TI\"]}], \"]\"}]\) returns {\!\(\*SubscriptBox[StyleBox[\"var\", \"TI\"], \"L\"]\), \!\(\*SubsuperscriptBox[StyleBox[\"var\", \"TI\"], \"L\", \"\[Dagger]\"]\)} for each label L in \!\(\*StyleBox[\"labels\", \"TI\"]\).";

FieldVariables::notformal =
"The variable `1` is not a strict Mathematica Formal symbol. Please use \\[Formal`1`] instead."

FieldVariables[] := {\[FormalA], SuperDagger[\[FormalA]]}

FieldVariables[var_Symbol ?FormalSymbolQ] := {var, SuperDagger[var]}

FieldVariables[vars : {__Symbol ? FormalSymbolQ}] :=
    Flatten[{#, SuperDagger[#]} & /@ vars]

FieldVariables[labels_List] := FieldVariables[\[FormalA], labels]

FieldVariables[var_Symbol ?FormalSymbolQ, labels_List] :=
    Flatten[  
        Table[
            {Symbol[SymbolName[var] <> ToString[L]], SuperDagger[Symbol[SymbolName[var] <> ToString[L]]]},
            {L, labels}
        ]
    ]

FieldVariables[var_Symbol, args___] /; !FormalSymbolQ[var] :=
    (Message[FieldVariables::notformal, var]; $Failed)


ExtractNCVars[ops_List] :=
    DeleteDuplicates @ Cases[ops, (v : (_? FormalSymbolQ | SuperDagger[_? FormalSymbolQ])) :> v, {0, Infinity}]



LiftNCScalars[expr_, vars_List] := expr //. {
    GeneralizedPower[NonCommutativeMultiply, Times[c_, op_], n_Integer ? Positive] /;
        FreeQ[c, Alternatives @@ vars] :>
            c^n GeneralizedPower[NonCommutativeMultiply, op, n],
    NonCommutativeMultiply[x___, Times[c_, y__], z___] /;
        FreeQ[c, Alternatives @@ vars] :>
            c NonCommutativeMultiply[x, y, z]
}


(* Working order for default behavior of NCA in version 15. *)
OrderVariables[vars_List, direction_String : "Normal"] := Block[{annihilators, creators},

    annihilators = Select[vars, FreeQ[#, SuperDagger] &];

    creators     = Select[vars, !FreeQ[#, SuperDagger] &];

    If[ direction === "Antinormal",
        Join[Sort[annihilators],Sort[creators]],
        Join[Sort[creators],Sort[annihilators]]
    ]
]


ToBosonicOperator::usage =
"\!\(\*RowBox[{\"ToBosonicOperator\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"]}], \"]\"}]\) realizes a symbolic expression expr built from field variables (see FieldVariables) as a truncated \!\(\*RowBox[{\"QuantumOperator\", \"[\", \"]\"}]\), replacing each annihilation variable with \!\(\*RowBox[{\"AnnihilationOperator\", \"[\", \"]\"}]\) and each \!\(\*SuperscriptBox[\"\[Dagger]\", \"\"]\)-decorated variable with its adjoint, and realizing \[FormalCapitalX]**\[FormalCapitalY] as operator composition.\n\!\(\*RowBox[{\"ToBosonicOperator\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) uses vars as the ordered list of annihilation-type field variables, assigning mode k to vars[[k]].\n\!\(\*RowBox[{\"ToBosonicOperator\", \"[\", RowBox[{StyleBox[\"expr\", \"TI\"], \",\", StyleBox[\"vars\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) truncates each mode to Fock space dimension size (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).";

SetAttributes[ToBosonicOperator, HoldFirst];

ToBosonicOperator[expr_, vars_List : Automatic, size : _Integer?Positive : $FockSize] :=
    Block[{annihilators, rules, resolved},
        resolved = ReleaseHold[Hold[expr] /. HoldPattern[Exp[x_]] :> boseExpHold[x]];
        annihilators = If[vars === Automatic,
            DeleteDuplicates[
                Cases[resolved, (v : (_ ? FormalSymbolQ | SuperDagger[_ ? FormalSymbolQ])) :> v, {0, Infinity}] /.
                    SuperDagger[v_] :> v
            ],
            vars
        ];
        rules = Flatten @ MapIndexed[
            With[{op = AnnihilationOperator[size, {First[#2]}]},
                {#1 -> op, SuperDagger[#1] -> op["Dagger"]}
            ] &,
            annihilators
        ];
        (resolved /. {
            GeneralizedPower[NonCommutativeMultiply, b_, n_Integer] :> Power[b, n],
            Power[base_, x_] /; FreeQ[base, Alternatives @@ annihilators] && !FreeQ[x, Alternatives @@ annihilators] :>
                MatrixExp[Log[base] x]
        } /. rules //.
            NonCommutativeMultiply[a_, b__] :> Fold[#1 @ #2 &, a, {b}]) /.
            boseExpHold[qo_QuantumOperator] :> Exp[qo]
    ]


G2Coherence::usage =
"\!\(\*RowBox[{\"G2Coherence\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"]}], \"]\"}]\) computes the second-order coherence function \!\(\*SuperscriptBox[\"g\", \"2\"]\)(0) for the given single-mode \!\(\*StyleBox[\"state\", \"TI\"]\).";

G2Coherence[\[Psi]_QuantumState] := Block[{aOp, nOp, a2Op, numerator, denominator},

    aOp  = AnnihilationOperator[\[Psi]["Dimension"]];
    
    nOp  = SuperDagger[aOp] @ aOp;
    
    a2Op = SuperDagger[aOp] @ SuperDagger[aOp] @ aOp @ aOp;
    
    Switch[ TimeConstrained[\[Psi]["PureStateQ"], 3, $Aborted],
        True,
        (* pure state: inner product formula *)
        numerator   = (\[Psi]["Dagger"] @ a2Op @ \[Psi])["Scalar"];
        denominator = (\[Psi]["Dagger"] @ nOp  @ \[Psi])["Scalar"],

        False,
        (* mixed state: trace formula Tr(rho O) *)
        numerator   = Tr[a2Op @ \[Psi]["Operator"]];
        denominator = Tr[nOp @ \[Psi]["Operator"]],

        _,
        (* unclassified state: trace formula on the density matrix *)
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
    Block[{s = OptionValue["QuadratureScaling"], dims = state["Dimensions"], R, A, norm, means, sigma},
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
            A     = # @ state["Operator"] & /@ R;
            norm  = Re @ Tr @ state["DensityMatrix"];
            means = Re[Tr /@ A];
            sigma = Re @ Outer[Tr[#1 @ #2] &, R, A, 1]
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
