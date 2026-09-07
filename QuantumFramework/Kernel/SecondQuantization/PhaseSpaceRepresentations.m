(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageImport["Wolfram`QuantumFramework`"]

PackageExport["WignerRepresentation"]

PackageExport["HusimiQRepresentation"]

PackageExport["WignerFunction"]

PackageExport["HusimiQFunction"]

PackageExport["SOrderedFunction"]

PackageExport["SOrderedRepresentation"]


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
        ], Flatten[Transpose[Re[w0] Exp[-B 0.5] (g^2 0.5 / Pi)]]}]]
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

(* Power continued to the phase-space origin: base^0 = 1 (the m == n limit),
   so a diagonal element contributes its parity / vacuum value at the origin
   instead of raising 0^0. A positive exponent still gives 0 at base = 0. *)
ZeroLimitPower[_, 0] := 1
ZeroLimitPower[base_, exp_] := base ^ exp

WignerKernel[gamma_, m_, n_] :=
 With[{zeta = 2 gamma},
  (-1)^n Exp[-(Abs[zeta]^2)/2] If[m >= n,
   Sqrt[n!/m!] ZeroLimitPower[zeta, m - n] LaguerreL[n, m - n, Abs[zeta]^2],
   Sqrt[m!/n!] Conjugate[-zeta]^(n - m) LaguerreL[m, n - m, Abs[zeta]^2]]]
   
WignerKernelWirtinger[gamma_, gammaStar_, m_, n_] :=
 With[{zeta = 2 gamma, zetaStar = 2 gammaStar},
  (-1)^n Exp[-(zeta zetaStar)/2] If[m >= n,
   Sqrt[n!/m!] ZeroLimitPower[zeta, m - n] LaguerreL[n, m - n, zeta zetaStar],
   Sqrt[m!/n!] (-zetaStar)^(n - m) LaguerreL[m, n - m, zeta zetaStar]]]
   
WignerKernelWirtingerInactive[gamma_, gammaStar_, m_, n_] :=
 With[{zeta = 2 gamma, zetaStar = 2 gammaStar},
  (-1)^n Exp[-(zeta zetaStar)/2] If[m >= n,
   Sqrt[n!/m!] ZeroLimitPower[zeta, m - n] Inactive[LaguerreL][n, m - n, zeta zetaStar],
   Sqrt[m!/n!] (-zetaStar)^(n - m) Inactive[LaguerreL][m, n - m, zeta zetaStar]]]
   
(* Husimi kernel for the density-matrix element <n|rho|m> at position {row=n+1, col=m+1}:
   Q(alpha) = (1/pi) <alpha|rho|alpha> puts the bra amplitude Conjugate[alpha]^n on the
   row index and the ket amplitude alpha^m on the column index. *)
HusimiKernel[alpha_, m_, n_] :=
 Exp[-Abs[alpha]^2] ZeroLimitPower[Conjugate[alpha], n] ZeroLimitPower[alpha, m] / Sqrt[m! n!]
   
(* Cahill-Glauber s-ordered kernel <m|T(alpha,s)|n>, scaled so the caller applies the same
   outer 2/Pi as WignerFunction, so that s = 0 reproduces the Wirtinger kernel above term for
   term. Written with one branch in {lo, hi} = MinMax[{m, n}]: the per-branch signs spelled
   out above are exactly ((s+1)/(s-1))^lo here. *)
SOrderedKernel[alpha_, alphaStar_, m_, n_, s_] :=
 With[{zeta = 2 alpha, zetaStar = 2 alphaStar, u = 1 - s},
  With[{lo = Min[m, n], hi = Max[m, n]},
   Exp[-(zeta zetaStar)/(2 u)] Sqrt[lo!/hi!] ((s + 1)/(s - 1))^lo *
    ZeroLimitPower[If[m >= n, zeta, zetaStar]/u, hi - lo] *
    LaguerreL[lo, hi - lo, (zeta zetaStar)/(1 - s^2)] / u]]

(* At s = -1 the general form is 0*Infinity, not merely ill-conditioned: the weight
   ((s+1)/(s-1))^lo vanishes while the Laguerre argument (zeta zetaStar)/(1-s^2) divides by
   zero, and the product evaluates to Indeterminate. The limit is the Husimi kernel, written
   here in Wirtinger form rather than delegating so that it stays free of Abs and Conjugate
   and can be compiled over independent quadratures. *)
SOrderedKernel[alpha_, alphaStar_, m_, n_, s_ /; TrueQ[s == -1]] :=
 Exp[-alpha alphaStar] ZeroLimitPower[alphaStar, n] ZeroLimitPower[alpha, m] / (2 Sqrt[m! n!])

SetAttributes[{WignerKernel, WignerKernelWirtinger, WignerKernelWirtingerInactive,HusimiKernel, SOrderedKernel}, Listable]



WignerFunction::usage =
"\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]\) gives the Wigner-function W(\[Alpha]) for the state \[Rho].\n\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"x\", \"TI\"], \",\", StyleBox[\"p\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Wigner quasi-probability distribution for quantum state \[Rho] using real quadrature variables x and p.\n\!\(\*RowBox[{\"WignerFunction\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", \"SymbolicForm->\", StyleBox[\"form\", \"TI\"]}], \"]\"}]\) \!\(\*StyleBox[\"'Wirtinger'\", \"TI\"]\) treats \[Alpha] and \[Alpha]\[Conjugate] as independent variables, \!\(\*StyleBox[\"'LaguerreForm'\", \"TI\"]\) holds Laguerre polynomials Inactive";


WignerFunction::multimode = "WignerFunction is only defined for single-mode states. The provided state has `1` mode(s).";

WignerFunction::badopt = "Unknown SymbolicForm `1`. Use Automatic, \"Wirtinger\", or \"LaguerreForm\".";

Options[WignerFunction] = {SymbolicForm -> Automatic};

WignerFunction[state_QuantumState, {x_, p_}, opts:OptionsPattern[]] :=
  If[state["Qudits"] =!= 1,
    (Message[WignerFunction::multimode, state["Qudits"]]; $Failed),
    1/2 ComplexExpand[WignerFunction[state, (x + I p)/Sqrt[2], opts]]];

WignerFunction[state_QuantumState, alpha_, opts:OptionsPattern[]] /; !ListQ[alpha] :=
  Block[{
    mat = state["DensityMatrix"],
    form = OptionValue[SymbolicForm],
    pos,
    vals
  },
    If[state["Qudits"] =!= 1, Message[WignerFunction::multimode, state["Qudits"]]; Return[$Failed]];
    pos = mat["ExplicitPositions"];
    vals = mat["ExplicitValues"];
    Which[
      form === "Wirtinger",
       2/Pi Dot[
           mat["ExplicitValues"],
           WignerKernelWirtinger[alpha, Conjugate[alpha], pos[[All,2]] - 1, pos[[All,1]] - 1]
           ],
        
      form === "LaguerreForm",
       2/Pi Dot[
           mat["ExplicitValues"],
           WignerKernelWirtingerInactive[alpha, Conjugate[alpha], pos[[All,2]] - 1, pos[[All,1]] - 1] 
           ],
        
      form === Automatic,
       2/Pi Dot[
           mat["ExplicitValues"],
           WignerKernel[alpha, pos[[All,2]] - 1, pos[[All,1]]-1] 
           ],
        
      True,
       Message[WignerFunction::badopt, form]; $Failed
    ]
  ]


HusimiQFunction::usage =
"\!\(\*RowBox[{\"HusimiQFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", StyleBox[\"\[Alpha]\", \"TI\"]}], \"]\"}]\) computes Q(\[Alpha]) for the state \[Rho] using the complex amplitude \[Alpha].\n\!\(\*RowBox[{\"HusimiQFunction\", \"[\", RowBox[{StyleBox[\"\[Rho]\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"x\", \"TI\"], \",\", StyleBox[\"p\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Husimi Q function for quantum state \[Rho] using real quadrature variables x and p."

HusimiQFunction::multimode = "HusimiQFunction is only defined for single-mode states. The provided state has `1` mode(s).";

HusimiQFunction[state_QuantumState, {x_, p_}] :=
  If[state["Qudits"] =!= 1,
    (Message[HusimiQFunction::multimode, state["Qudits"]]; $Failed),
    1/2 ComplexExpand[HusimiQFunction[state, (x + I p)/Sqrt[2]]]];

HusimiQFunction[state_QuantumState, alpha_] /; !ListQ[alpha] :=
 Block[{
   mat = state["DensityMatrix"],
   pos,
   val},
   If[state["Qudits"] =!= 1, Message[HusimiQFunction::multimode, state["Qudits"]]; Return[$Failed]];
   pos = mat["ExplicitPositions"];
   val = mat["ExplicitValues"];
   1/Pi Dot[val, HusimiKernel[alpha, pos[[All,2]] - 1, pos[[All,1]] - 1]]
   ]


(* ::Input::Initialization::Plain:: *)
HusimiQRepresentation::usage =
"\!\(\*RowBox[{\"HusimiQRepresentation\", \"[\", RowBox[{StyleBox[\"state\", \"TI\"], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"xmin\", \"TI\"], \",\", StyleBox[\"xmax\", \"TI\"]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{StyleBox[\"pmin\", \"TI\"], \",\", StyleBox[\"pmax\", \"TI\"]}], \"}\"}]}], \"]\"}]\) computes the Husimi Q quasi-probability distribution Q(x,p). Returns an InterpolatingFunction over the specified x and p limits.\n\!\(\*RowBox[{\"HusimiQRepresentation\", \"[\", RowBox[{\"\[Ellipsis]\", \",\", StyleBox[\"opts\", \"TI\"]}], \"]\"}]\) Options: \"GaussianScaling\" \[Rule] \!\(\*SqrtBox[\(2\)]\) (default), \"GridSize\" -> 100 (default).";

HusimiQRepresentation::multimode = "HusimiQRepresentation is only defined for single-mode states. The provided state has `1` mode(s).";

Options[HusimiQRepresentation]={
"GaussianScaling"->Sqrt[2],
"GridSize"->100
};

HusimiQRepresentation[state_QuantumState, {xmin_, xmax_}, {pmin_, pmax_},
     OptionsPattern[]] :=
    Module[{X, Y, amat, qmat, vals, vecs, g, xvec, pvec, outerList},
        
        If[state["Qudits"] =!= 1, Message[HusimiQRepresentation::multimode, state["Qudits"]]; Return[$Failed]];
        
        g = OptionValue["GaussianScaling"];

        xvec = Subdivide[xmin, xmax, OptionValue["GridSize"] - 1];

        pvec = Subdivide[pmin, pmax, OptionValue["GridSize"] - 1];

        outerList = Outer[List, xvec, pvec];

        {X, Y} = Transpose[outerList, {3, 2, 1}];

        amat = 0.5 g (X + I Y);

        qmat = ConstantArray[0, Dimensions[amat]];

        (* HusimiPure returns Q(alpha) = (1/pi) |<alpha|psi>|^2, normalized to
           unit integral over d^2 alpha. With alpha = (g/2)(x + I p), the
           0.25 g^2 = |d^2 alpha / d(x,p)| Jacobian re-expresses it as a density
           over dx dp, so it integrates to 1 and agrees with the {x,p} closed
           form. The mixed branch below carries the same factor. *)
        If[state["PureStateQ"],
            qmat = 0.25 g^2 HusimiPure[state, amat],
            {vals, vecs} = state["Eigensystem"];
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


(* ::Input::Initialization::Plain:: *)
SOrderedFunction::usage =
"SOrderedFunction[state, alpha, s] gives the Cahill-Glauber s-ordered quasi-probability \
distribution W(alpha, s) = (1/Pi) Tr[rho T(alpha, s)]. s = 0 is the Wigner function, s = -1 \
the Husimi Q function, and s -> 1 the Glauber-Sudarshan P distribution.\n\
SOrderedFunction[state, {x, p}, s] uses real quadratures, alpha = (x + I p)/Sqrt[2].\n\
s may be symbolic; s = 1 is not a function and is rejected.";

SOrderedFunction::multimode = "SOrderedFunction is only defined for single-mode states. The provided state has `1` mode(s).";

SOrderedFunction::sunit = "s = 1 is the Glauber-Sudarshan P distribution, which is not a function: the kernel T(alpha, 1) is undefined and the s -> 1 limit is a distribution. Pair an s < 1 result with a test function instead.";

SOrderedFunction::snum = "s must be a real number less than 1; got `1`.";

SOrderedWirtinger[state_QuantumState, alpha_, alphaStar_, s_] :=
 With[{mat = state["DensityMatrix"]},

  With[{pos = mat["ExplicitPositions"]},
  
   (2/Pi) Dot[
    mat["ExplicitValues"],
    SOrderedKernel[alpha, alphaStar, pos[[All, 2]] - 1, pos[[All, 1]] - 1, s]]]]

SOrderedFunction[state_QuantumState, {x_, p_}, s_] :=
 Which[
  state["Qudits"] =!= 1, Message[SOrderedFunction::multimode, state["Qudits"]]; $Failed,
  TrueQ[s == 1], Message[SOrderedFunction::sunit]; $Failed,
  NumericQ[s] && ! TrueQ[Re[s] < 1], Message[SOrderedFunction::snum, s]; $Failed,
  True, 1/2 ComplexExpand @ SOrderedWirtinger[state, (x + I p)/Sqrt[2], (x - I p)/Sqrt[2], s]]

SOrderedFunction[state_QuantumState, alpha_, s_] /; ! ListQ[alpha] :=
 Which[
  state["Qudits"] =!= 1, Message[SOrderedFunction::multimode, state["Qudits"]]; $Failed,
  TrueQ[s == 1], Message[SOrderedFunction::sunit]; $Failed,
  NumericQ[s] && ! TrueQ[Re[s] < 1], Message[SOrderedFunction::snum, s]; $Failed,
  True, SOrderedWirtinger[state, alpha, Conjugate[alpha], s]]


(* ::Input::Initialization::Plain:: *)
SOrderedRepresentation::usage =
"SOrderedRepresentation[state, {xmin, xmax}, {pmin, pmax}, s] computes the s-ordered \
quasi-probability distribution over the given quadrature window. Returns an \
InterpolatingFunction. s must be numeric and less than 1.\n\
SOrderedRepresentation[..., opts] options: \"GaussianScaling\" -> Sqrt[2] (default), \
\"GridSize\" -> 100 (default).";

SOrderedRepresentation::multimode = "SOrderedRepresentation is only defined for single-mode states. The provided state has `1` mode(s).";

SOrderedRepresentation::sunit = "s = 1 is the Glauber-Sudarshan P distribution, which is not a function and cannot be sampled on a grid.";

SOrderedRepresentation::snum = "s must be a real number less than 1; got `1`.";

SOrderedRepresentation::coarse = "Grid spacing `1` cannot resolve structure of scale Sqrt[1-s] = `2`; increase \"GridSize\" or lower s.";

Options[SOrderedRepresentation] = {
"GaussianScaling" -> Sqrt[2],
"GridSize" -> 100
};

SOrderedRepresentation[state_QuantumState, {xmin_, xmax_}, {pmin_, pmax_}, s_?NumericQ,
    OptionsPattern[]] :=
 Module[{n, g, xvec, pvec, dx, xq, pq, cf, vals},

  Which[
   state["Qudits"] =!= 1,
    Message[SOrderedRepresentation::multimode, state["Qudits"]]; Return[$Failed],
   TrueQ[s == 1], Message[SOrderedRepresentation::sunit]; Return[$Failed],
   ! TrueQ[Re[s] < 1], Message[SOrderedRepresentation::snum, s]; Return[$Failed]];

  n = OptionValue["GridSize"];
  g = OptionValue["GaussianScaling"];
  xvec = N @ Subdivide[xmin, xmax, n - 1];
  pvec = N @ Subdivide[pmin, pmax, n - 1];

  (* Structure shrinks like Sqrt[1-s] while its amplitude grows, so a window that resolves
     the Wigner function can miss an s near 1 entirely: at s = 0.9 a hundred-point grid
     integrated to -122 instead of 1 while the pointwise values were still good to thirteen
     digits. *)
  dx = Min[xvec[[2]] - xvec[[1]], pvec[[2]] - pvec[[1]]];
  If[dx > Sqrt[1 - s]/2, Message[SOrderedRepresentation::coarse, dx, N @ Sqrt[1 - s]]];

  cf = Compile[{{xq, _Real}, {pq, _Real}},
   Evaluate @ SOrderedWirtinger[state, g (xq + I pq)/2, g (xq - I pq)/2, s],
   RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"];

  (* g^2/4 is the Jacobian from d^2 alpha to dx dp, as in WignerRepresentation. *)
  vals = (g^2/4) Re @ cf[ConstantArray[xvec, n], Transpose @ ConstantArray[pvec, n]];

  Interpolation @ MapThread[List, {Flatten[Outer[List, xvec, pvec], 1], Flatten[vals]}]]
