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
"\!\(\*RowBox[{\"AnnihilationOperator\", \"[\", \"]\"}]\) gives the bosonic annihilation operator \[AHat].\n\!\(\*RowBox[{\"AnnihilationOperator\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).\n\!\(\*RowBox[{\"AnnihilationOperator\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies the applied subsystem \!\(\*StyleBox[\"order\", \"TI\"]\) for multi-mode systems.";

AnnihilationOperator[order_?orderQ]:= AnnihilationOperator[$FockSize, order]

AnnihilationOperator[size_:$FockSize, Optional[order_?orderQ, {1}]]:= AnnihilationOperator[size, order] =
	QuantumOperator[
	SparseArray[Band[{1,2}]->Sqrt[Range[size-1]],{size,size}],
	order, size, "Label" -> OverHat["a"]]


PhaseShiftOperator::usage =
"\!\(\*RowBox[{\"PhaseShiftOperator\", \"[\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"]}], \"]\"}]\) gives the phase shift operator e^{i\[Theta]n}.\n\!\(\*RowBox[{\"PhaseShiftOperator\", \"[\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).\n\!\(\*RowBox[{\"PhaseShiftOperator\", \"[\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies the subsystem \!\(\*StyleBox[\"order\", \"TI\"]\).";

PhaseShiftOperator[\[Theta]_,order_?orderQ]:= PhaseShiftOperator[\[Theta], $FockSize, order]

PhaseShiftOperator[\[Theta]_,size_:$FockSize,Optional[order_?orderQ,{1}]]:=
	QuantumOperator[SparseArray[Band[{1,1}]->Exp[I \[Theta] Range[0,size-1]],{size,size}],order,size,
		"Label" -> "\[CapitalPhi]"[\[Theta]]]


DisplacementOperator::usage =
"\!\(\*RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]\) gives the displacement operator D(\[Alpha]) with complex amplitude \[Alpha].\n\!\(\*RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).\n\!\(\*RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies the subsystem \!\(\*StyleBox[\"order\", \"TI\"]\).\n\!\(\*RowBox[{\"DisplacementOperator\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"\\\"Ordering\\\"->\", StyleBox[\"ord\", \"TI\"]}], \"]\"}]\) \"Ordering\" accepts \"Normal\" | \"Weak\" | \"Antinormal\" for operator ordering definition; Automatic, the default, is \"Normal\" for a matrix and \"Weak\" for the symbolic form.\n\!\(\*RowBox[{\"DisplacementOperator\", \"[\", RowBox[{StyleBox[\"\[Alpha]\", \"TI\"], \",\", \"\[Infinity]\", \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives D(\[Alpha]) as a symbolic expression in the ladder operators vars.";


Options[DisplacementOperator] = {"Ordering" -> Automatic};
DisplacementOperator::invalidorder = "The value for the 'Ordering' option, `1`, is invalid. Choose from 'Normal', 'Weak', or 'Antinormal'.";

DisplacementOperator[\[Alpha]_, opts : OptionsPattern[]] := DisplacementOperator[\[Alpha], $FockSize,{1}, opts];

DisplacementOperator[\[Alpha]_, size : Except[\[Infinity]], opts : OptionsPattern[]] := DisplacementOperator[\[Alpha], size,{1}, opts];

DisplacementOperator[\[Alpha]_, order_?orderQ, opts : OptionsPattern[]] := DisplacementOperator[\[Alpha],$FockSize,order, opts];

symbolicMode[_, Automatic] := FieldVariables[]

symbolicMode[_, order_ ? orderQ] := FieldVariables[order]

symbolicMode[f_, vars_] :=
    Replace[ResolveLadderVariables[f, vars], {
        lv_Association /; Length[lv["Annihilators"]] == 1 :> {lv["Creators"], lv["Annihilators"]},
        _Association :> (Message[f::ladone, vars]; $Failed)
    }]

DisplacementOperator[\[Alpha]_, \[Infinity], vars : Except[_Rule | _RuleDelayed] : Automatic, OptionsPattern[]] :=
    Block[{
        ordering = Replace[OptionValue["Ordering"], Automatic -> "Weak"],
        lv = symbolicMode[DisplacementOperator, vars],
        a, adag
    },
        If[lv === $Failed, Return[$Failed, Block]];
        {{adag}, {a}} = lv;
        Switch[ordering,
            "Normal",
            Exp[-\[Alpha] Conjugate[\[Alpha]]/2] * Exp[\[Alpha] * adag] ** Exp[-Conjugate[\[Alpha]] * a],
        
            "Weak",
            Exp[\[Alpha] * adag - Conjugate[\[Alpha]] * a],
        
            "Antinormal",
            Exp[\[Alpha] Conjugate[\[Alpha]]/2] * Exp[-Conjugate[\[Alpha]] * a] ** Exp[\[Alpha] * adag],
        _,
            Message[DisplacementOperator::invalidorder, ordering];
            Abort[]
    ]
]


DisplacementOperator[\[Alpha]_, size_Integer ? Positive, order_?orderQ, OptionsPattern[]] :=
    Block[
        {a = AnnihilationOperator[size,order],
        ordering = Replace[OptionValue["Ordering"], Automatic -> "Normal"]
    },
        QuantumOperator[
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
            ],
            "Label" -> OverHat["D"][\[Alpha]]
        ]
    ]



SqueezeOperator::usage =
"\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"]}], \"]\"}]\) gives the squeeze operator S(\[Xi]) with complex squeeze parameter \[Xi].\n\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: $FockSize).\n\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"], \",\", StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies the subsystem \!\(\*StyleBox[\"order\", \"TI\"]\).\n\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"\\\"Ordering\\\"->\", StyleBox[\"ord\", \"TI\"]}], \"]\"}]\)\"Ordering\" accepts \"Normal\" | \"Weak\" | \"Antinormal\" for operator ordering; Automatic, the default, is \"Normal\" for a matrix and \"Weak\" for the symbolic form.\n\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{StyleBox[\"\[Xi]\", \"TI\"], \",\", \"\[Infinity]\", \",\", StyleBox[\"vars\", \"TI\"]}], \"]\"}]\) gives S(\[Xi]) as a symbolic expression in the ladder operators vars.\n\!\(\*RowBox[{\"SqueezeOperator\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"r\", \"TI\"], \",\", StyleBox[\"\[Theta]\", \"TI\"]}], \"}\"}], \",\", \"\[Ellipsis]\"}], \"]\"}]\) gives S(\!\(\*StyleBox[\"r\", \"TI\"]\) \!\(\*SuperscriptBox[\"e\", RowBox[{\"i\", StyleBox[\"\[Theta]\", \"TI\"]}]]\)) in polar form.";

Options[SqueezeOperator] = {"Ordering" -> Automatic};

SqueezeOperator::invalidorder = "The value for the 'Ordering' option, `1`, is invalid. Choose from 'Normal', 'Weak', or 'Antinormal'.";

SqueezeOperator[xi_, opts : OptionsPattern[]] := SqueezeOperator[xi, $FockSize,{1}, opts];

SqueezeOperator[xi_, size : Except[\[Infinity]], opts : OptionsPattern[]] := SqueezeOperator[xi, size,{1}, opts];

SqueezeOperator[xi_, order_?orderQ, opts : OptionsPattern[]] := SqueezeOperator[xi, $FockSize,order, opts];

SqueezeOperator[{r_, th_}, \[Infinity], vars : Except[_Rule | _RuleDelayed] : Automatic, OptionsPattern[]] :=
    With[{ordering = Replace[OptionValue["Ordering"], Automatic -> "Weak"], lv = symbolicMode[SqueezeOperator, vars]},
        If[ lv === $Failed, $Failed,
            With[{adag = lv[[1, 1]], a = lv[[2, 1]]},
                Replace[
                    Switch[ordering,
                        "Normal", Exp[-E^(I th) Tanh[r]/2 adag ** adag] ** Sech[r]^(adag ** a + 1/2) ** Exp[E^(-I th) Tanh[r]/2 a ** a],
                        "Weak", Exp[r (E^(-I th) a ** a - E^(I th) adag ** adag)/2],
                        "Antinormal", Exp[E^(-I th) Tanh[r]/2 a ** a] ** Cosh[r]^(adag ** a + 1/2) ** Exp[-E^(I th) Tanh[r]/2 adag ** adag],
                        _, Message[SqueezeOperator::invalidorder, ordering]; Abort[]
                    ],
                    NonCommutativeMultiply[1 ..] -> 1
                ]
            ]
        ]
    ]

SqueezeOperator[xi_, \[Infinity], vars : Except[_Rule | _RuleDelayed] : Automatic, opts : OptionsPattern[]] :=
    If[ Replace[OptionValue["Ordering"], Automatic -> "Weak"] === "Weak",
        Replace[symbolicMode[SqueezeOperator, vars], {{{adag_}, {a_}} :> Exp[(Conjugate[xi] a ** a - xi adag ** adag)/2], _ -> $Failed}],
        SqueezeOperator[{Abs[xi], Arg[xi]}, \[Infinity], vars, opts]
    ]

SqueezeOperator[xi_, size_Integer ? Positive, order_?orderQ, OptionsPattern[]] :=
    Module[{r, th, x, xb, tau, taub, nu, a = AnnihilationOperator[size,order], ordering},

        ordering = Replace[OptionValue["Ordering"], Automatic -> "Normal"];

        {r, th} = Replace[xi, {p : {_, _} :> p, _ :> {Abs[xi], Arg[xi]}}];

        {x, xb} = Replace[xi, {{_, _} :> r {E^(I th), E^(-I th)}, _ :> {xi, Conjugate[xi]}}];

        {tau, taub} = Tanh[r] {E^(I th), E^(-I th)};

        nu = Log[Cosh[r]];


        QuantumOperator[
            Switch[ordering,
                "Normal",
                    MatrixExp[-tau / 2 ((a["Dagger"]) @ (a["Dagger"]))] @ MatrixExp[-
                        nu ((a["Dagger"]) @ a + 1/2 )] @ MatrixExp[taub / 2 (a
                        @ a)]
                ,
                "Weak",
                    MatrixExp[1/2 (xb (a @ a) - x (a["Dagger"] @ a["Dagger"]))]
                ,
                "Antinormal",
                    MatrixExp[1/2 taub (a @ a)] @ MatrixExp[nu ((a["Dagger"
                        ]) @ a + 1/2 )] @ MatrixExp[-1/2 tau ((a["Dagger"]) @ (a["Dagger"
                        ]))]
                ,
                _,
                    Message[SqueezeOperator::invalidorder, ordering];
                    Abort[]
            ],
            "Label" -> OverHat["S"][xi]
        ]
    ]





BeamsplitterMatrix[th_,ph_,cutoff_]:= Block[{sqrt,ct,st,stconj,Z},
	sqrt=Sqrt[Range[0,cutoff-1]];
	ct=Cos[th];
	st=Sin[th] Exp[I ph];
    stconj = -Exp[-I ph] Sin[th];
	Z=ConstantArray[0,{cutoff,cutoff,cutoff,cutoff}];
	Z[[1,1,1,1]]=1;
	Do[
		If[0<p<cutoff,
			Z[[m+1,n+1,p+1,1]]=
			  If[m>0, ct sqrt[[m+1]]/sqrt[[p+1]] Z[[m,n+1,p,1]],0]+
			  If[n>0, st sqrt[[n+1]]/sqrt[[p+1]] Z[[m+1,n,p,1]],0]
			],
			{m,0,cutoff-1},{n,0,cutoff-1},{p,m+n,m+n}];
	Do[
		With[{q=m+n-p},
			If[0<q<cutoff,
				Z[[m+1,n+1,p+1,q+1]]=
				  If[m>0, stconj sqrt[[m+1]]/sqrt[[q+1]] Z[[m,n+1,p+1,q]],0]+
				  If[n>0, ct sqrt[[n+1]]/sqrt[[q+1]] Z[[m+1,n,p+1,q]],0]]
				],
			{m,0,cutoff-1},{n,0,cutoff-1},{p,0,cutoff-1}];
	Z]


BeamSplitterOperator::usage =
"\!\(\*RowBox[{\"BeamSplitterOperator\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"], \",\", StyleBox[\"\[Phi]\", \"TI\"]}], \"}\"}]}], \"]\"}]\) gives the two-mode beam splitter operator with mixing angle \[Theta] and phase \[Phi].\n\!\(\*RowBox[{\"BeamSplitterOperator\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"], \",\", StyleBox[\"\[Phi]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).\n\!\(\*RowBox[{\"BeamSplitterOperator\", \"[\", RowBox[{RowBox[{\"{\", RowBox[{StyleBox[\"\[Theta]\", \"TI\"], \",\", StyleBox[\"\[Phi]\", \"TI\"]}], \"}\"}], \",\", StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies subsystem \!\(\*StyleBox[\"order\", \"TI\"]\).\n\!\(\*RowBox[{\"BeamSplitterOperator\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"Method->\", StyleBox[\"\", \"TI\"]}], \"]\"}]\) option: Method -> \"MatrixExp\" | \"Recurrence\".";

Options[BeamSplitterOperator] = {Method -> "MatrixExp"};

BeamSplitterOperator::badmethod = 
  "Method `1` not recognized. Supported methods: \"MatrixExp\", \"Recurrence\".";

BeamSplitterOperator[{th_, ph_}, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{th, ph}, $FockSize, {1, 2}, opts]

BeamSplitterOperator[{th_, ph_}, order_?orderQ, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{th, ph}, $FockSize, order, opts]

BeamSplitterOperator[{th_, ph_}, size_Integer, opts:OptionsPattern[]] :=
  BeamSplitterOperator[{th, ph}, size, {1,2}, opts]

BeamSplitterOperator[{th_, ph_},size_Integer, order_?orderQ, OptionsPattern[]] :=
 Module[{method = OptionValue[Method], op1, op2},
  Switch[method,
   "MatrixExp",
    {op1, op2} = AnnihilationOperator[size, {#}] & /@ order;
    QuantumOperator[
     MatrixExp[
      th (Exp[I ph] op1 @ SuperDagger[op2] -
         Exp[-I ph] SuperDagger[op1] @ op2)
     ],
     "Label" -> OverHat["BS"][th, ph]
    ],
   "Recurrence",
    QuantumOperator[BeamsplitterMatrix[th, ph, size], order, {size, size},
     "Label" -> OverHat["BS"][th, ph]],
   _,
    Message[BeamSplitterOperator::badmethod, method]; $Failed
  ]
]


QuadratureOperators::usage =
"\!\(\*RowBox[{\"QuadratureOperators\", \"[\", \"]\"}]\) returns {X, P} position and momentum quadrature operators.\n\!\(\*RowBox[{\"QuadratureOperators\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"]}], \"]\"}]\) specifies the Fock space \!\(\*StyleBox[\"size\", \"TI\"]\) (default: \!\(\*StyleBox[\"$FockSize\", \"TI\"]\)).\n\!\(\*RowBox[{\"QuadratureOperators\", \"[\", RowBox[{StyleBox[\"size\", \"TI\"], \",\", StyleBox[\"order\", \"TI\"]}], \"]\"}]\) specifies the subsystem \!\(\*StyleBox[\"order\", \"TI\"]\).";

QuadratureOperators[order_?orderQ] := QuadratureOperators[$FockSize,order]

QuadratureOperators[size_:$FockSize, Optional[order_?orderQ, {1}]]:= Block[{a=AnnihilationOperator[size,order]},
							{QuantumOperator[1/2(a+a["Dagger"]), "Label" -> Subscript[OverHat["X"],1]],
							QuantumOperator[1/(2I)(a-a["Dagger"]), "Label" -> Subscript[OverHat["X"],2]]}
						]
