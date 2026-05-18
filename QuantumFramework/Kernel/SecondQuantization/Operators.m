(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["AnnihilationOperator"]

PackageExport["DisplacementOperator"]

PackageExport["SqueezeOperator"]

PackageExport["QuadratureOperators"]

PackageExport["BeamSplitterOperator"]
	
PackageExport["PhaseShiftOperator"]


AnnihilationOperator::usage = 
"\!\(AnnihilationOperator[]\) Creates the bosonic annihilation operator \[AHat].
\!\(AnnihilationOperator[size]\) Specifies the Fock space size (default: $FockSize).
\!\(AnnihilationOperator[size, order]\) Specifies which subsystem (order) for multi-mode systems.";

AnnihilationOperator[order_?orderQ]:= AnnihilationOperator[$FockSize, order]

AnnihilationOperator[size_:$FockSize, Optional[order_?orderQ, {1}]]:= AnnihilationOperator[size, order] = 
	QuantumOperator[
	SparseArray[Band[{1,2}]->Sqrt[Range[size-1]],{size,size}],
	order, size]


PhaseShiftOperator::usage = 
"\!\(PhaseShiftOperator[\[Theta]]\) Creates the phase shift operator e^{i\[Theta]n}.
\!\(PhaseShiftOperator[\[Theta], size]\) Specifies the Fock space size (default: $FockSize).
\!\(PhaseShiftOperator[\[Theta], size, order]\) Specifies the subsystem order.";

PhaseShiftOperator[\[Theta]_,order_?orderQ]:= PhaseShiftOperator[\[Theta], $FockSize, order]

PhaseShiftOperator[\[Theta]_,size_:$FockSize,Optional[order_?orderQ,{1}]]:= 
	QuantumOperator[SparseArray[Band[{1,1}]->Exp[I \[Theta] Range[0,size-1]],{size,size}],order,size]


DisplacementOperator::usage = 
"\!\(DisplacementOperator[\[Alpha]]\) Creates the displacement operator D(\[Alpha]) with complex amplitude \[Alpha].
\!\(DisplacementOperator[\[Alpha], size]\) Specifies the Fock space size (default: $FockSize).
\!\(DisplacementOperator[\[Alpha], size, order]\) Specifies the subsystem order.
\!\(DisplacementOperator[\[Ellipsis], \"Ordering\"\[Rule] ]\) \"Ordering\" accepts \"Normal\" | \"Weak\" | \"Antinormal\" for operator ordering.";

Options[DisplacementOperator] = {"Ordering" -> "Normal"};
DisplacementOperator::invalidorder = "The value for the 'Ordering' option, `1`, is invalid. Choose from 'Normal', 'Weak', or 'Antinormal'.";

DisplacementOperator[\[Alpha]_, opts : OptionsPattern[]] := DisplacementOperator[\[Alpha], $FockSize,{1}, opts];

DisplacementOperator[\[Alpha]_, size_, opts : OptionsPattern[]] := DisplacementOperator[\[Alpha], size,{1}, opts];

DisplacementOperator[\[Alpha]_, order_?orderQ, opts : OptionsPattern[]] := DisplacementOperator[\[Alpha],$FockSize,order, opts];

DisplacementOperator[\[Alpha]_, size_, order_?orderQ, OptionsPattern[]] :=
    Block[{a = AnnihilationOperator[size,order], ordering},
        
        ordering = OptionValue["Ordering"];
        
        Switch[ordering,
            "Normal",
                Exp[-\[Alpha] Conjugate[\[Alpha]]/2]  MatrixExp[\[Alpha] (a["Dagger"])] @ MatrixExp[-Conjugate[\[Alpha]] a],

            "Weak",
               MatrixExp[\[Alpha] a["Dagger"]-Conjugate[\[Alpha]] a],

            "Antinormal",
                Exp[\[Alpha] Conjugate[\[Alpha]]/2]  MatrixExp[-Conjugate[\[Alpha]] a] @ MatrixExp[\[Alpha] (a["Dagger"])] ,

            _, 
                Message[DisplacementOperator::invalidorder, ordering];
                Abort[]
        ]
    ]



SqueezeOperator::usage = 
"\!\(SqueezeOperator[\[Xi]]\) Creates the squeeze operator S(\[Xi]) with complex squeeze parameter \[Xi].
\!\(SqueezeOperator[\[Xi], size]\) Specifies the Fock space size (default: $FockSize).
\!\(SqueezeOperator[\[Xi], size, order]\) Specifies the subsystem order.
\!\(SqueezeOperator[\[Ellipsis], \"Ordering\" \[Rule]]\)\"Ordering\" accepts \"Normal\" | \"Weak\" | \"Antinormal\" for operator ordering.";

Options[SqueezeOperator] = {"Ordering" -> "Normal"};

SqueezeOperator::invalidorder = "The value for the 'Ordering' option, `1`, is invalid. Choose from 'Normal', 'Weak', or 'Antinormal'.";

SqueezeOperator[xi_, opts : OptionsPattern[]] := SqueezeOperator[xi, $FockSize,{1}, opts];

SqueezeOperator[xi_, size_, opts : OptionsPattern[]] := SqueezeOperator[xi, size,{1}, opts];

SqueezeOperator[xi_, order_?orderQ, opts : OptionsPattern[]] := SqueezeOperator[xi, $FockSize,order, opts];

SqueezeOperator[xi_, size_, order_?orderQ, OptionsPattern[]] :=
    Module[{tau, nu, a = AnnihilationOperator[size,order], ordering},
    
        ordering = OptionValue["Ordering"];
        
        tau = xi / Abs[xi] Tanh[Abs[xi]];
        
        nu = Log[Cosh[Abs[xi]]];
        
        
        Switch[ordering,
            "Normal",
                MatrixExp[-tau / 2 ((a["Dagger"]) @ (a["Dagger"]))] @ MatrixExp[-
                    nu ((a["Dagger"]) @ a + 1/2 )] @ MatrixExp[Conjugate[tau] / 2 (a 
                    @ a)]
            ,
            "Weak",
                MatrixExp[1/2 (Conjugate[xi] (a @ a) - xi (a["Dagger"] @ a["Dagger"]))]
            ,
            "Antinormal",
                MatrixExp[1/2 Conjugate[tau] (a @ a)] @ MatrixExp[-nu ((a["Dagger"
                    ]) @ a + 1/2 )] @ MatrixExp[-1/2 tau ((a["Dagger"]) @ (a["Dagger"
                    ]))]
            ,
            _,
                Message[SqueezeOperator::invalidorder, ordering];
                Abort[]
        ]
    ]


BeamsplitterMatrix[\[Theta]_,\[Phi]_,cutoff_]:=Module[{sqrt,ct,st,R,Z},
	sqrt=Sqrt[Range[0,cutoff-1]];
	ct=Cos[\[Theta]];
	st=Sin[\[Theta]] Exp[I \[Phi]];
	R={{0,0,ct,-Conjugate[st]},{0,0,st,ct},{ct,st,0,0},{-Conjugate[st],ct,0,0}};
	Z=ConstantArray[0,{cutoff,cutoff,cutoff,cutoff}];
	Z[[1,1,1,1]]=1;
	Do[
		If[0<p<cutoff,
			Z[[m+1,n+1,p+1,1]]=
			  If[m>0,R[[1,3]] sqrt[[m+1]]/sqrt[[p+1]] Z[[m,n+1,p,1]],0]+
			  If[n>0,R[[2,3]] sqrt[[n+1]]/sqrt[[p+1]] Z[[m+1,n,p,1]],0]
			],
			{m,0,cutoff-1},{n,0,cutoff-1},{p,m+n,m+n}];
	Do[
		With[{q=m+n-p},
			If[0<q<cutoff,
				Z[[m+1,n+1,p+1,q+1]]=
				  If[m>0,R[[1,4]] sqrt[[m+1]]/sqrt[[q+1]] Z[[m,n+1,p+1,q]],0]+
				  If[n>0,R[[2,4]] sqrt[[n+1]]/sqrt[[q+1]] Z[[m+1,n,p+1,q]],0]]
				],
			{m,0,cutoff-1},{n,0,cutoff-1},{p,0,cutoff-1}];
	Z]


BeamSplitterOperator::usage = 
"\!\(BeamSplitterOperator[{\[Theta], \[Phi]}]\) Creates a two-mode beam splitter operator with mixing angle \[Theta] and phase \[Phi].
\!\(BeamSplitterOperator[{\[Theta], \[Phi]}, size]\) Specifies the Fock space size (default: $FockSize).
\!\(BeamSplitterOperator[{\[Theta], \[Phi]}, size, order]\) Specifies the mode ordering.
\!\(BeamSplitterOperator[\[Ellipsis], Method \[Rule] ]\) Options: Method -> \"MatrixExp\" | \"Recurrence\".";

Options[BeamSplitterOperator] = {Method -> "MatrixExp"};

BeamSplitterOperator::badmethod = 
  "Method `1` not recognized. Supported methods: \"MatrixExp\", \"Recurrence\".";

BeamSplitterOperator[{\[Theta]_, \[Phi]_}, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{\[Theta], \[Phi]}, $FockSize, {1, 2}, opts]

BeamSplitterOperator[{\[Theta]_, \[Phi]_}, order_?orderQ, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{\[Theta], \[Phi]}, $FockSize, order, opts]
  
BeamSplitterOperator[{\[Theta]_, \[Phi]_}, size_Integer, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{\[Theta], \[Phi]}, size, {1,2}, opts]

BeamSplitterOperator[{\[Theta]_, \[CurlyPhi]_},size_Integer, order_?orderQ, OptionsPattern[]] :=
 Module[{method = OptionValue[Method], op1, op2},
  Switch[method,
   "MatrixExp",
    {op1, op2} = AnnihilationOperator[size, {#}] & /@ order;
    MatrixExp[
     \[Theta] (Exp[I \[CurlyPhi]] op1 @ SuperDagger[op2] - 
        Exp[-I \[CurlyPhi]] SuperDagger[op1] @ op2)
    ],
   "Recurrence",
    QuantumOperator[BeamsplitterMatrix[\[Theta], \[CurlyPhi], size], order, {size, size}],
   _,
    Message[BeamSplitterOperator::badmethod, method]; $Failed
  ]
]


QuadratureOperators::usage = 
"\!\(QuadratureOperators[]\) Returns {X, P} position and momentum quadrature operators.
\!\(QuadratureOperators[size]\) Specifies the Fock space size (default: $FockSize).
\!\(QuadratureOperators[size, order]\) Specifies the subsystem order.";

QuadratureOperators[order_?orderQ] := QuadratureOperators[$FockSize,order]

QuadratureOperators[size_:$FockSize, Optional[order_?orderQ, {1}]]:= Block[{a=AnnihilationOperator[size,order]},
							{1/2(a+a["Dagger"]),
							1/(2I)(a-a["Dagger"])}
						]
