(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["$FockSize"]

PackageExport["SetFockSpaceSize"]

PackageExport["OperatorVariance"]

PackageExport["FieldVariables"]

PackageExport["OrderVariables"]

PackageScope["ExtractNCVars"]

PackageScope["FormalSymbolQ"]

PackageExport["G2Coherence"]

PackageExport["G1Correlation"]

PackageExport["CovarianceMatrix"]


$FockSize::usage = "Global variable holding the current Fock space truncation size (default 16).";


SetFockSpaceSize::usage = "\!\(\*RowBox[{\"SetFockSpaceSize\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) sets $FockSize to \!\(\*StyleBox[\"size\", \"TI\"]\). Default is 16.";

SetFockSpaceSize[size:_Integer?Positive:16]:= $FockSize = size;


OperatorVariance::usage =
"\!\(\*RowBox[{\"OperatorVariance\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", StyleBox[\"op\", \"TI\"]}], \"]\"}]\) computes \[LeftAngleBracket]\!\(\*SuperscriptBox[\"O\", \"2\"]\)\[RightAngleBracket] - \!\(\*SuperscriptBox[\"\[LeftAngleBracket]O\[RightAngleBracket]\", \"2\"]\) of \!\(\*StyleBox[\"op\", \"TI\"]\) in \!\(\*StyleBox[\"state\", \"TI\"]\).";

OperatorVariance[state_QuantumState, op_QuantumOperator]:= 
	(state["Dagger"]@ (op @ op)@ state)["Scalar"] - (state["Dagger"]@ op @ state)["Scalar"]^2


FormalSymbolQ[var_Symbol]:=StringStartsQ[ToString[FullForm[var]],"\\[Formal"];

FormalSymbolQ[_]:=False;


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


(* Working order for default behavior of NCA in version 15 *)
OrderVariables[vars_List] := Block[{annihilators, creators},

    annihilators = Select[vars, FreeQ[#, SuperDagger] &];
    
    creators     = Select[vars, !FreeQ[#, SuperDagger] &];
    
    Join[Sort[creators],Sort[annihilators]]
]


G2Coherence::usage =
"\!\(\*RowBox[{\"G2Coherence\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"]}], \"]\"}]\) computes the second-order coherence function \!\(\*SuperscriptBox[\"g\", \"2\"]\)(0) for the given single-mode \!\(\*StyleBox[\"state\", \"TI\"]\).";

G2Coherence[\[Psi]_QuantumState] := Block[{aOp, nOp, a2Op, numerator, denominator},

    aOp  = AnnihilationOperator[\[Psi]["Dimension"]];
    
    nOp  = SuperDagger[aOp] @ aOp;
    
    a2Op = SuperDagger[aOp] @ SuperDagger[aOp] @ aOp @ aOp;
    
    If[ \[Psi]["PureStateQ"],
        (* pure state: inner product formula *)
        numerator   = (\[Psi]["Dagger"] @ a2Op @ \[Psi])["Scalar"];
        denominator = (\[Psi]["Dagger"] @ nOp  @ \[Psi])["Scalar"],
        
        (* mixed state: trace formula Tr(rho O) *)
        numerator   = Tr[a2Op @ \[Psi]["Operator"]];
        denominator = Tr[nOp @ \[Psi]["Operator"]]
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


Options[CovarianceMatrix] = {"QuadratureScaling" -> 1/Sqrt[2]};

CovarianceMatrix::usage =
"\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"]}], \"]\"}]\) computes the 2\[Times]2 covariance matrix for the single-mode \!\(\*StyleBox[\"state\", \"TI\"]\) in the Serafini convention \!\(\*SubscriptBox[\(\[Sigma]\), \(vac\)]\) = 1/2 I.\n\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) computes the 2n\[Times]2n multi-mode covariance matrix for the n modes specified by \!\(\*StyleBox[\"order\", \"TI\"]\).\n\!\(\*RowBox[{\"CovarianceMatrix\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"\\\"QuadratureScaling\\\"->\", StyleBox[\"s\", \"TI\"]}], \"]\"}]\) uses \!\(\*OverscriptBox[\"X\", \"^\"]\) = s(a+\!\(\*SuperscriptBox[\"a\", \"\[Dagger]\"]\)), P = \[ImaginaryI] s(\!\(\*SuperscriptBox[\"a\", \"\[Dagger]\"]\)-a). Default s = 1/\!\(\*SqrtBox[\"2\"]\) (Serafini). Use s = 1/2 for \[HBar]=1/2 or s = 1 for Simon convention.";

CovarianceMatrix[state_QuantumState, opts : OptionsPattern[]] :=
    CovarianceMatrix[state, {1}, opts]

CovarianceMatrix[state_QuantumState, order_?orderQ, OptionsPattern[]] :=
    Block[{s = OptionValue["QuadratureScaling"], size, R, n, psi, Rpsi, rhoOp, Qops, means, upper},
        size = First[state["Dimensions"]];
        R    = Catenate[(2s * QuadratureOperators[size, {#}]) & /@ order];
        n    = Length[R];
        If[state["PureStateQ"],
            psi   = state["StateVector"];
            Rpsi  = (# @ state)["StateVector"] & /@ R;
            means = Re[Conjugate[Rpsi] . psi];
            Re[Conjugate[Rpsi] . Transpose[Rpsi]] - Outer[Times, means, means],
            rhoOp = state["Operator"];
            Qops  = # @ rhoOp & /@ R;
            means = Re[Tr /@ Qops];
            upper = Table[
                If[i <= j,
                    Re[Tr[R[[j]] @ Qops[[i]]]] - means[[i]] means[[j]],
                    0],
                {i, n}, {j, n}];
            upper + Transpose[upper] - DiagonalMatrix[Diagonal[upper]]
        ]
    ]


SetFockSpaceSize[];
