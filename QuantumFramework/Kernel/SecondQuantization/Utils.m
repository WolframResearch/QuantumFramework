(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["$FockSize"]

PackageExport["SetFockSpaceSize"]

PackageExport["OperatorVariance"]

PackageExport["FieldVariables"]

PackageExport["OrderVariables"]

PackageScope["ExtractNCVars"]

PackageScope["FormalSymbolQ"]

PackageExport["G2Coherence"]

PackageExport["G1Correlation"]


$FockSize::usage = "Global variable holding the current Fock space truncation size (default 16).";


SetFockSpaceSize::usage = "\!\(SetFockSpaceSize[size]\) Sets the global Fock space truncation size $FockSize. Default is 16.";

SetFockSpaceSize[size:_Integer?Positive:16]:=$FockSize = size;


OperatorVariance::usage = 
"\!\(OperatorVariance[state, op]\) Computes the variance \[LeftAngleBracket]O^2\[RightAngleBracket] - \[LeftAngleBracket]O\[RightAngleBracket]^2 of operator op for the given state.";

OperatorVariance[state_QuantumState, op_QuantumOperator]:= 
	(state["Dagger"]@ (op @ op)@ state)["Scalar"] - (state["Dagger"]@ op @ state)["Scalar"]^2


FormalSymbolQ[var_Symbol]:=StringStartsQ[ToString[FullForm[var]],"\\[Formal"];

FormalSymbolQ[_]:=False;


FieldVariables::usage =
"\!\(FieldVariables[]\) Returns the default field variables {\[FormalA], \!\(\*SuperscriptBox[\(\[FormalA]\), \(\[Dagger]\)]\)}.
\!\(FieldVariables[var]\) Returns {var, \!\(\*SuperscriptBox[\(var\), \(\[Dagger]\)]\)} for a formal symbol var.
\!\(FieldVariables[labels]\) Returns field variables indexed by labels, using \[FormalA] as the base symbol.
\!\(FieldVariables[var, labels]\) Returns {\!\(\*SubscriptBox[\(var\), \(L\)]\), \!\(\*SubsuperscriptBox[\(var\), \(L\), \(\[Dagger]\)]\)} for each label L in labels.";

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
"\!\(G2Coherence[state]\) Computes the second-order coherence function g^2(0) for the given single mode state.";

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
"\!\(G1Correlation[state,{{r1,t1},{r2,t2}}]\) Computes the first-order correlation function G^1(x1,x2) for the given single mode state for space-time coordinates r1,t1,r2,t2";


G1Correlation[state_QuantumState, {{r1_, t1_}, {r2_, t2_}}] :=Block[{aOp, eMinus1, ePlus2},

    aOp    = AnnihilationOperator[state["Dimension"]];
    
    eMinus1 = -I SuperDagger[aOp] Exp[-I (\[FormalK] . r1 - \[FormalOmega] t1)];
    
    ePlus2  =  I aOp Exp[ I (\[FormalK] . r2 - \[FormalOmega] t2)];
    
    If[ state["PureStateQ"],
    
        (state["Dagger"] @ eMinus1 @ ePlus2 @ state)["Scalar"],
        
        Tr[(eMinus1 @ ePlus2)@state["Operator"]]
    ]
]


SetFockSpaceSize[];
