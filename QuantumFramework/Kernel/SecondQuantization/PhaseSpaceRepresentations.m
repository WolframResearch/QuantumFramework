(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["WignerRepresentation"]

PackageExport["HusimiQRepresentation"]

PackageExport["WignerFunction"]

PackageExport["HusimiQFunction"]


(* ::Input::Initialization::Plain:: *)
WignerRepresentation::usage =
"\!\(\*RowBox[{\"WignerRepresentation\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"xmin\",\"TI\"], \",\", StyleBox[\"xmax\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"pmin\", \"TI\"], \",\", StyleBox[\"pmax\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the numerical Wigner quasi-probability distribution W(x,p)for the single mode \!\(\*StyleBox[\"state\", \"TI\"]\). Returns an InterpolatingFunction over the specified x and p limits.\n\!\(\*RowBox[{\"WignerRepresentation\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", StyleBox[\"opts\", \"TI\"]}], \"]\"}]\) options: \"GaussianScaling\" \[Rule] \!\(\*SqrtBox[\(2\)]\)(default), \"GridSize\" \[Rule] 100 (default).";

WignerRepresentation::multimode = "WignerRepresentation is only defined for single-mode states. The provided state has `1` mode(s).";

Options[WignerRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};

WignerRepresentation[psi_QuantumState, {xmin_, xmax_}, {pmin_,pmax_}, OptionsPattern[]] := Module[
    {rho, M, X, Y, A2, B, w0, diag, g, xvec, pvec},

    If[psi["Qudits"] =!= 1, Message[WignerRepresentation::multimode, psi["Qudits"]]; Return[$Failed]];

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


(* Phase space Kernels *)
\[ScriptCapitalK]mn[\[Gamma]_, m_, n_] :=
 With[{\[Zeta] = 2 \[Gamma]},
  (-1)^n Exp[-(Abs[\[Zeta]]^2)/2] If[m >= n,
   Sqrt[n!/m!] \[Zeta]^(m - n) LaguerreL[n, m - n, Abs[\[Zeta]]^2],
   Sqrt[m!/n!] Conjugate[-\[Zeta]]^(n - m) LaguerreL[m, n - m, Abs[\[Zeta]]^2]]]
   
\[ScriptCapitalK]mnWirt[\[Gamma]_, \[Gamma]Star_, m_, n_] :=
 With[{\[Zeta] = 2 \[Gamma], \[Zeta]Star = 2 \[Gamma]Star},
  (-1)^n Exp[-(\[Zeta] \[Zeta]Star)/2] If[m >= n,
   Sqrt[n!/m!] \[Zeta]^(m - n) LaguerreL[n, m - n, \[Zeta] \[Zeta]Star],
   Sqrt[m!/n!] (-\[Zeta]Star)^(n - m) LaguerreL[m, n - m, \[Zeta] \[Zeta]Star]]]
   
\[ScriptCapitalK]mnWirtInactive[\[Gamma]_, \[Gamma]Star_, m_, n_] :=
 With[{\[Zeta] = 2 \[Gamma], \[Zeta]Star = 2 \[Gamma]Star},
  (-1)^n Exp[-(\[Zeta] \[Zeta]Star)/2] If[m >= n,
   Sqrt[n!/m!] \[Zeta]^(m - n) Inactive[LaguerreL][n, m - n, \[Zeta] \[Zeta]Star],
   Sqrt[m!/n!] (-\[Zeta]Star)^(n - m) Inactive[LaguerreL][m, n - m, \[Zeta] \[Zeta]Star]]]
   
\[ScriptCapitalK]HusimiMN[\[Alpha]_, m_, n_] := 
 Exp[-Abs[\[Alpha]]^2] \[Alpha]^n Conjugate[\[Alpha]]^m / Sqrt[m! n!]
   
SetAttributes[{\[ScriptCapitalK]mn, \[ScriptCapitalK]mnWirt, \[ScriptCapitalK]mnWirtInactive,\[ScriptCapitalK]HusimiMN}, Listable]



WignerFunction::usage =
"\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]\) gives the Wigner-function W(\[Alpha]) for the state \[Rho].\n\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"x\", \"TI\"], \",\", StyleBox[\"p\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Wigner quasi-probability distribution for quantum state \[Rho] using real quadrature variables x and p.\n\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"SymbolicForm->\", StyleBox[\"form\", \"TI\"]}], \"]\"}]\) \!\(\*StyleBox[\"'Wirtinger'\", \"TI\"]\) treats \[Alpha] and \[Alpha]\[Conjugate] as independent variables, \!\(\*StyleBox[\"'LaguerreForm'\", \"TI\"]\) holds Laguerre polynomials Inactive";


WignerFunction::multimode = "WignerFunction is only defined for single-mode states. The provided state has `1` mode(s).";

WignerFunction::badopt = "Unknown SymbolicForm `1`. Use Automatic, \"Wirtinger\", or \"LaguerreForm\".";

Options[WignerFunction] = {SymbolicForm -> Automatic};

WignerFunction[\[Rho]_QuantumState, {x_, p_}, opts:OptionsPattern[]] :=
  If[\[Rho]["Qudits"] =!= 1,
    (Message[WignerFunction::multimode, \[Rho]["Qudits"]]; $Failed),
    1/2 ComplexExpand[WignerFunction[\[Rho], (x + I p)/Sqrt[2], opts]]];

WignerFunction[\[Rho]_QuantumState, \[Alpha]_, opts:OptionsPattern[]] /; !ListQ[\[Alpha]] :=
  Block[{
    mat = \[Rho]["DensityMatrix"],
    form = OptionValue[SymbolicForm],
    pos,
    vals
  },
    If[\[Rho]["Qudits"] =!= 1, Message[WignerFunction::multimode, \[Rho]["Qudits"]]; Return[$Failed]];
    pos = mat["ExplicitPositions"];
    vals = mat["ExplicitValues"];
    Which[
      form === "Wirtinger",
       2/\[Pi] Dot[
           mat["ExplicitValues"],
           \[ScriptCapitalK]mnWirt[\[Alpha], Conjugate[\[Alpha]], pos[[All,2]] - 1, pos[[All,1]] - 1]
           ],
        
      form === "LaguerreForm",
       2/\[Pi] Dot[
           mat["ExplicitValues"],
           \[ScriptCapitalK]mnWirtInactive[\[Alpha], Conjugate[\[Alpha]], pos[[All,2]] - 1, pos[[All,1]] - 1] 
           ],
        
      form === Automatic,
       2/\[Pi] Dot[
           mat["ExplicitValues"],
           \[ScriptCapitalK]mn[\[Alpha], pos[[All,2]] - 1, pos[[All,1]]-1] 
           ],
        
      True,
       Message[WignerFunction::badopt, form]; $Failed
    ]
  ]


HusimiQFunction::usage =
"\!\(\*RowBox[{\"HusimiQFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]\) computes Q(\[Alpha]) for the state \[Rho] using the complex amplitude \[Alpha].\n\!\(\*RowBox[{\"HusimiQFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"x\", \"TI\"], \",\", StyleBox[\"p\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Husimi Q function for quantum state \[Rho] using real quadrature variables x and p."

HusimiQFunction::multimode = "HusimiQFunction is only defined for single-mode states. The provided state has `1` mode(s).";

HusimiQFunction[\[Rho]_QuantumState, {x_, p_}] :=
  If[\[Rho]["Qudits"] =!= 1,
    (Message[HusimiQFunction::multimode, \[Rho]["Qudits"]]; $Failed),
    1/2 ComplexExpand[HusimiQFunction[\[Rho], (x + I p)/Sqrt[2]]]];

HusimiQFunction[\[Rho]_QuantumState, \[Alpha]_] /; !ListQ[\[Alpha]] :=
 Block[{
   mat = \[Rho]["DensityMatrix"],
   pos,
   val},
   If[\[Rho]["Qudits"] =!= 1, Message[HusimiQFunction::multimode, \[Rho]["Qudits"]]; Return[$Failed]];
   pos = mat["ExplicitPositions"];
   val = mat["ExplicitValues"];
   1/\[Pi] Dot[val, \[ScriptCapitalK]HusimiMN[\[Alpha], pos[[All,2]] - 1, pos[[All,1]] - 1]]
   ]


(* ::Input::Initialization::Plain:: *)
HusimiQRepresentation::usage =
"\!\(\*RowBox[{\"HusimiQRepresentation\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"xmin\", \"TI\"], \",\", StyleBox[\"xmax\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"pmin\", \"TI\"], \",\", StyleBox[\"pmax\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Husimi Q quasi-probability distribution Q(x,p). Returns an InterpolatingFunction over the specified x and p limits.\n\!\(\*RowBox[{\"HusimiQRepresentation\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", StyleBox[\"opts\", \"TI\"]}], \"]\"}]\) Options: \"GaussianScaling\" \[Rule] \!\(\*SqrtBox[\(2\)]\) (default), \"GridSize\" -> 100 (default).";

HusimiQRepresentation::multimode = "HusimiQRepresentation is only defined for single-mode states. The provided state has `1` mode(s).";

Options[HusimiQRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};

HusimiQRepresentation[\[Psi]_QuantumState, {xmin_, xmax_}, {pmin_, pmax_},
     OptionsPattern[]] :=
    Module[{X, Y, amat, qmat, vals, vecs, g, xvec, pvec, outerList},
        
        If[\[Psi]["Qudits"] =!= 1, Message[HusimiQRepresentation::multimode, \[Psi]["Qudits"]]; Return[$Failed]];
        
        g = OptionValue["GaussianScaling"];

        xvec = Subdivide[xmin, xmax, OptionValue["GridSize"] - 1];

        pvec = Subdivide[pmin, pmax, OptionValue["GridSize"] - 1];

        outerList = Outer[List, xvec, pvec];

        {X, Y} = Transpose[outerList, {3, 2, 1}];

        amat = 0.5 g (X + I Y);

        qmat = ConstantArray[0, Dimensions[amat]];

        If[\[Psi]["PureStateQ"],
            qmat = HusimiPure[\[Psi], amat],
            {vals, vecs} = \[Psi]["Eigensystem"];
            {vals, vecs} =
                With[{mask = Unitize[vals]},
                    {Pick[vals, mask, 1], Pick[vecs, mask, 1]}
                ];
            qmat = Re[0.25 g^2 Dot[vals, Map[HusimiPure[QuantumState[
                #, Length @ #], amat]&, vecs]]]
        ];
        Interpolation[MapThread[List, {Flatten[outerList, 1], Flatten[
            Transpose @ qmat]}]]
    ]


HusimiPure[psi_QuantumState, alphaMat_] :=
    Module[{n,psiVecScal,qmat},
        n = Times @@ psi["Dimensions"];
        
        psiVecScal=psi["StateVector"]/Sqrt[Factorial @Range[0,n-1]];
        
		qmat=Abs[FromDigits[Reverse[psiVecScal], Conjugate[alphaMat]]]^2;
		
		Re[qmat] Exp[-Abs[alphaMat]^2]/Pi
    ]
