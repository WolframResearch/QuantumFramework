(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["WignerRepresentation"]

PackageExport["HusimiQRepresentation"]


(* ::Input::Initialization::Plain:: *)
WignerRepresentation::usage = 
"\!\(WignerRepresentation[state, {xmin, xmax}, {pmin, pmax}]\) Computes the Wigner quasi-probability distribution W(x,p). Returns an InterpolatingFunction over the specified phase space region.
\!\(WignerRepresentation[\[Ellipsis], opts]\) Options: \"GaussianScaling\" \[Rule] \!\(\*SqrtBox[\(2\)]\)(default), \"GridSize\" \[Rule] 100 (default).";

Options[WignerRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};

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


(* ::Input::Initialization::Plain:: *)
HusimiQRepresentation::usage = 
"\!\(HusimiQRepresentation[state, {xmin, xmax}, {pmin, pmax}]\) Computes the Husimi Q quasi-probability distribution Q(x,p). Returns an InterpolatingFunction over the specified phase space region.
\!\(HusimiQRepresentation[\[Ellipsis], opts]\) Options: \"GaussianScaling\" \[Rule] \!\(\*SqrtBox[\(2\)]\) (default), \"GridSize\" -> 100 (default).";

Options[HusimiQRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};

HusimiQRepresentation[\[Psi]_QuantumState, {xmin_, xmax_}, {pmin_, pmax_},
     OptionsPattern[]] :=
    Module[{X, Y, amat, qmat, vals, vecs, g, xvec, pvec, outerList},
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
