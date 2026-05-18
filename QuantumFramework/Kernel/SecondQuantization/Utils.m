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
"\!\(FieldVariables[]\) Returns the default field variables {\[FormalA], SuperDagger[\[FormalA]]}.
\!\(FieldVariables[var]\) Returns {var, SuperDagger[var]} for a formal symbol var.
\!\(FieldVariables[labels]\) Returns field variables indexed by labels, using \[FormalA] as the base symbol.
\!\(FieldVariables[var, labels]\) Returns {var\!\(\*SubscriptBox[\()\), \(L\)]\), SuperDagger[var\!\(\*SubscriptBox[\()\), \(L\)]\)]} for each label L in labels.";

FieldVariables::notformal="The variable `1` is not a strict Mathematica Formal symbol. Please use \\[Formal`1`] instead.";

FieldVariables[]:={\[FormalA],SuperDagger[\[FormalA]]};

FieldVariables[var_Symbol?FormalSymbolQ]:={var,SuperDagger[var]};

FieldVariables[labels_List]:=FieldVariables[\[FormalA],labels];

FieldVariables[var_Symbol?FormalSymbolQ,labels_List]:=Flatten[Table[
{Symbol[SymbolName[var]<>ToString[L]],SuperDagger[Symbol[SymbolName[var]<>ToString[L]]]},
{L,labels}]];

FieldVariables[var_Symbol,args___]/;!FormalSymbolQ[var]:=Module[{},Message[FieldVariables::notformal,var];$Failed];


ExtractNCVars[ops_List]:=DeleteDuplicates@Cases[ops,(v:(_?FormalSymbolQ|SuperDagger[_?FormalSymbolQ])):>v,{0,Infinity}]


OrderVariables[vars_List]:=Module[{annihilators,creators},
annihilators=Select[vars,FreeQ[#,SuperDagger]&];
creators=Select[vars,!FreeQ[#,SuperDagger]&];
Join[ReverseSort[annihilators],ReverseSort[creators]]]


G2Coherence::usage = 
"\!\(G2Coherence[state]\) Computes the second-order coherence function g^2(0) for the given single mode state.";

G2Coherence[\[Psi]_QuantumState] := Module[{numerator, denominator, a2Op, aOp, nOp},
    aOp = AnnihilationOperator[\[Psi]["Dimension"]];
    nOp = SuperDagger[aOp] @ aOp;
    a2Op = SuperDagger[aOp] @ SuperDagger[aOp] @ aOp @ aOp;
    If[\[Psi]["PureStateQ"],
        (* Pure state: use inner product formula *)
        numerator = (SuperDagger[\[Psi]] @ a2Op @ \[Psi])["Scalar"];
        denominator = (SuperDagger[\[Psi]] @ nOp @ \[Psi])["Scalar"]
    ,
        (* Mixed state: use trace formula Tr(\[Rho]O) *)
        numerator = Tr[\[Psi]["DensityMatrix"] . a2Op["Matrix"]];
        denominator = Tr[\[Psi]["DensityMatrix"] . nOp["Matrix"]]
    ];
    numerator / denominator^2
]


G1Correlation::usage = 
"\!\(G1Correlation[state,{{r1,t1},{r2,t2}}]\) Computes the first-order correlation function G^1(x1,x2) for the given single mode state for space-time coordinates r1,t1,r2,t2";


G1Correlation[state_QuantumState,{{r1_,t1_},{r2_,t2_}}]:=Module[{aOp,Eminus1,Eplus2},
	aOp=AnnihilationOperator[state["Dimension"]];
	Eminus1=-I SuperDagger[aOp] Exp[-I(\[FormalK] . r1-\[FormalOmega] t1)];
	Eplus2=I aOp Exp[I(\[FormalK] . r2-\[FormalOmega] t2)];
	If[state["PureStateQ"],
		(SuperDagger[state]@Eminus1@Eplus2@state)["Scalar"],
		Tr[state["DensityMatrix"] . (Eminus1@Eplus2)["Matrix"]]
	]
]
