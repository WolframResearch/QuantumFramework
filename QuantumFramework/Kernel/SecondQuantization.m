(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageImport["Wolfram`QuantumFramework`PackageScope`"]

PackageExport["$FockSize"]

PackageExport["SetFockSpaceSize"]

PackageExport["G2Coherence"]

PackageExport["OperatorVariance"]

PackageExport["CoherentState"]

PackageExport["FockState"]

PackageExport["AnnihilationOperator"]

PackageExport["DisplacementOperator"]

PackageExport["SqueezeOperator"]

PackageExport["ThermalState"]

PackageExport["CatState"]

PackageExport["QuadratureOperators"]

PackageExport["WignerRepresentation"]

PackageExport["HusimiQRepresentation"]

PackageExport["BeamSplitterOperator"]
	
PackageExport["PhaseShiftOperator"]


(* ::Section:: *)
(*General Definitions*)


SetFockSpaceSize[size:_Integer?Positive:16]:=$FockSize = size;


OperatorVariance[state_QuantumState, op_QuantumOperator]:= 
	(state["Dagger"]@ (op @ op)@ state)["Scalar"] - (state["Dagger"]@ op @ state)["Scalar"]^2


G2Coherence[\[Psi]_QuantumState,aOp_QuantumOperator]:= (SuperDagger[\[Psi]]@(SuperDagger[aOp]@SuperDagger[aOp]@aOp@aOp)@\[Psi])["Scalar"]/(SuperDagger[\[Psi]]@(SuperDagger[aOp]@aOp)@\[Psi])["Scalar"]^2


FockState::clip = "Index `1` was outside the valid range {0, `2`} and has been clipped.";

FockState[n_Integer, size_ : $FockSize] := Block[{nEff = n},
  If[n < 0 || n >= size,
    Message[FockState::clip, n, size - 1];
    nEff = Clip[n, {0, size - 1}];
  ];
  
  QuantumState[SparseArray[{nEff + 1 -> 1}, size], size]
];



FockState[vals_List, size_:$FockSize]:= 

If[!AllTrue[vals,IntegerQ[#]&&(0<=#<size)&], Message[FockVals::len,vals,size],

QuantumState[SparseArray[{FromDigits[vals,size]+1->1},size^Length[vals]],
			ConstantArray[size,Length[vals]]]
  ]

FockVals::len="Values of `1` must be non negative integers less that the desired size of the space: `2`";


AnnihilationOperator[order_?orderQ]:= AnnihilationOperator[$FockSize, order]

AnnihilationOperator[size_:$FockSize, Optional[order_?orderQ, {1}]]:= AnnihilationOperator[size, order] = 
	QuantumOperator[
	SparseArray[Band[{1,2}]->Sqrt[Range[size-1]],{size,size}],
	order, size]


CoherentState[size_:$FockSize] :=
    Block[{n=0},
       QuantumState[NestList[(n++; # \[FormalAlpha] / Sqrt[n])&, 1, size - 1], size, "Parameters" -> \[FormalAlpha]]
    ]


ThermalState[nbar_, size_:$FockSize] :=
    QuantumOperator[
        DiagonalMatrix[1 / (1 + nbar) 
            Table[(nbar / (1 + nbar)) ^ n, 
           {n, 0, size-1}
          ]
        ], 
    size]["MatrixQuantumState"]


PhaseShiftOperator[\[Theta]_,order_?orderQ]:= PhaseShiftOperator[\[Theta], $FockSize, order]

PhaseShiftOperator[\[Theta]_,size_:$FockSize,Optional[order_?orderQ,{1}]]:= 
	QuantumOperator[SparseArray[Band[{1,1}]->Exp[I \[Theta] Range[0,size-1]],{size,size}],order,size]


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


QuadratureOperators[order_?orderQ] := QuadratureOperators[$FockSize,order]

QuadratureOperators[size_:$FockSize, Optional[order_?orderQ, {1}]]:= Block[{a=AnnihilationOperator[size,order]},
							{1/2(a+a["Dagger"]),
							1/(2I)(a-a["Dagger"])}
						]


CatState[size_:$FockSize] :=
    Block[{amplitudes,n=0},
        amplitudes =
            Transpose[
                    NestList[
                        (
                            n++;
                            # {\[FormalAlpha] / Sqrt[n], -\[FormalAlpha] / Sqrt[n]}
                        )&
                        ,
                        {1, 1}
                        ,
                        size - 1
                    ]
            ];
        QuantumState[amplitudes[[1]] + E ^ (I \[FormalPhi]) amplitudes[[2]], size,
             "Parameters" -> {\[FormalAlpha], \[FormalPhi]}]["Normalize"]
    ]


SetFockSpaceSize[];


(* ::Section:: *)
(*Quasi-probability Distributions*)


(* ::Subsection:: *)
(*Wigner*)


(* ::Input::Initialization::Plain:: *)
Options[WignerRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};


(* ::Input::Initialization::Plain:: *)
WignerRepresentation[psi_QuantumState, {xmin_, xmax_}, {pmin_,pmax_}, OptionsPattern[]] := Module[{rho, M, X, Y, A2, B, w0, diag,
 
    g, xvec, pvec},

    rho = psi["DensityMatrix"];

    g = OptionValue["GaussianScaling"];

    M = Length[rho];

    xvec = N@Subdivide[xmin, xmax, OptionValue["GridSize"] - 1];

    pvec = N@Subdivide[pmin, pmax, OptionValue["GridSize"] - 1];

    {X, Y} = Transpose[Outer[List, xvec, pvec], {3, 2, 1}];

    A2 = g * (X + I Y);

    B = Abs[A2] ^ 2;

    w0 = ConstantArray[2 rho[[1, -1]], {Length[xvec], Length[pvec]}];
        
    While[
        M > 1
        ,
        M--;
        diag = Diagonal[rho, M - 1] If[M != 1,
            2
            ,
            1
        ];
        w0 = WigLaguerreVal[M - 1, B, diag] + w0 A2 * M ^ -0.5;
    ];

    Interpolation[MapThread[List, {Flatten[Outer[List, xvec, pvec], 1
        ], Flatten[Transpose[Re[w0] Exp[-B 0.5] (g^2 0.5 / \[Pi])]]}]]
]


(* ::Input::Initialization::Plain:: *)
WigLaguerreVal[L_, x_, c_] :=
    Module[{y0, y1, k, n},
        n = Length[c];
        Switch[n,
            1,
                y0 = c[[1]];
                y1 = 0.
            ,
            2,
                y0 = c[[1]];
                y1 = c[[2]]
            ,
            _,
                k = n;
                y0 = c[[-2]];
                y1 = c[[-1]];
                Do[
                    k--;
                    {y0, y1} = {c[[-i]] - y1 Sqrt[((k - 1.) (L + k - 
                        1.)) / ((L + k) k)], y0 - y1 (L + 2. k - 1 - x) Sqrt[1 / ((L + k) k)]
                        }
                    ,
                    {i, 3, n}
                ]
        ];
        y0 - y1 Sqrt[1 / (L + 1.)] (L + 1. - x)
    ]


(* ::Subsection:: *)
(*Husimi Q *)


Options[HusimiQRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};


(* ::Input::Initialization::Plain:: *)
HusimiQRepresentation[\[Psi]_QuantumState, {xmin_, xmax_}, {pmin_,pmax_}, OptionsPattern[]] :=
    Module[{X, Y, amat, qmat, d, v, qmatList, k, nonZeroEigenpairs, g,
         xvec, pvec},
        g = OptionValue["GaussianScaling"];

        xvec = Subdivide[xmin, xmax, OptionValue["GridSize"] - 1];

        pvec = Subdivide[pmin, pmax, OptionValue["GridSize"] - 1];

        {X, Y} = Transpose[Outer[List, xvec, pvec], {3, 2, 1}];

        amat = 0.5 g (X + I Y);

        qmat = ConstantArray[0, Dimensions[amat]];

        If[\[Psi]["PureStateQ"],
            qmat = HusimiPure[\[Psi], amat]
            ,
            nonZeroEigenpairs = Select[Transpose[Eigensystem[\[Psi]["DensityMatrix"
                ]]], First[#] != 0&];
            qmatList = Map[(#[[1]] HusimiPure[QuantumState[#[[2]], Length[
                #[[2]]]], amat])&, nonZeroEigenpairs];
            qmat = 0.25 Total[Re /@ qmatList] g^2
        ];

        Interpolation[MapThread[List, {Flatten[Outer[List, xvec, pvec
            ], 1], Flatten[Transpose@qmat]}]]
    ];


HusimiPure[psi_QuantumState, alphaMat_] :=
    Module[{n, psiVec, qmat,z},
        n = Times @@ psi["Dimensions"];
        psiVec = psi["StateVector"];
        qmat = Function[{q, x},
                    q . z ^ Range[Length[q] - 1, 0, -1] /. z -> x
                ][Reverse[psiVec / Sqrt[Factorial /@ Range[0, n - 1]]
                    ], Conjugate[alphaMat]] // Abs[#] ^ 2&;
        Re[qmat] Exp[-Abs[alphaMat] ^ 2] / Pi
    ]
