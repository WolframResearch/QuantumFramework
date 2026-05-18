(* ::Package:: *)

Package["Wolfram`QuantumFramework`SecondQuantization`"]

PackageScope["BlasiakNormalOrder"]

PackageScope["MultiModeBlasiakOrder"]


BlasiakBuildTerm[coeff_,p_,q_,b_,bdag_]:=Which[p==0&&q==0,coeff,
p==0,coeff GeneralizedPower[NonCommutativeMultiply,b,q],
q==0,coeff GeneralizedPower[NonCommutativeMultiply,bdag,p],
True,coeff GeneralizedPower[NonCommutativeMultiply,bdag,p]**GeneralizedPower[NonCommutativeMultiply,b,q]];

BlasiakBuildTerm[coeff_,1,q_,b_,bdag_]:=Which[q==0,coeff bdag,q==1,coeff bdag**b,True,coeff bdag**GeneralizedPower[NonCommutativeMultiply,b,q]];
BlasiakBuildTerm[coeff_,p_,1,b_,bdag_]:=Which[p==0,coeff b,p==1,coeff bdag**b,True,coeff GeneralizedPower[NonCommutativeMultiply,bdag,p]**b];
BlasiakBuildTerm[coeff_,1,1,b_,bdag_]:=coeff bdag**b;
BlasiakBuildTerm[coeff_,0,1,b_,bdag_]:=coeff b;
BlasiakBuildTerm[coeff_,1,0,b_,bdag_]:=coeff bdag;


BlasiakNormalOrder[rList_List, sList_List, var_] :=
  Module[{M = Length[rList], b = var, bdag = SuperDagger[var],
          cumD, SrsCache, Srs, dM},

   
    cumD = FoldList[Plus, 0, rList - sList];  
  
    SrsCache = <||>;
    Srs[k_] := Lookup[SrsCache, k,
      SrsCache[k] =
        (1/k!) * Sum[
          Binomial[k, j] * (-1)^(k - j) *
            Product[FactorialPower[cumD[[m]] + j, sList[[m]]], {m, 1, M}],
          {j, 0, k}]
    ];

    dM = cumD[[-1]];  
    If[dM >= 0,
      Sum[BlasiakBuildTerm[Srs[k], dM + k, k, b, bdag],
          {k, sList[[1]], Total[sList]}],

      Module[{rRev = Reverse[rList], sRev = Reverse[sList],
              cumDRev, SrsCacheRev, SrsRev},
        cumDRev    = FoldList[Plus, 0, sRev - rRev];
        SrsCacheRev = <||>;
        SrsRev[k_] := Lookup[SrsCacheRev, k,
          SrsCacheRev[k] =
            (1/k!) * Sum[
              Binomial[k, j] * (-1)^(k - j) *
                Product[FactorialPower[cumDRev[[m]] + j, rRev[[m]]], {m, 1, M}],
              {j, 0, k}]
        ];
        Sum[BlasiakBuildTerm[SrsRev[k], k, k - dM, b, bdag],
            {k, rList[[-1]], Total[rList]}]
      ]
    ]
  ];


ParseBlasiakMonomial[expr_, var_] :=
  Module[{bdag = SuperDagger[var], factors, blocks, cur, typ, pow,
          rList, sList},

    factors = If[Head[expr] === NonCommutativeMultiply, List @@ expr, {expr}];

    blocks = {};
    Do[
      {typ, pow} = Which[
        f === var,
          {"b", 1},
        f === bdag,
          {"bdag", 1},
        MatchQ[f, GeneralizedPower[NonCommutativeMultiply, var,  _]],
          {"b",    f[[3]]},
        MatchQ[f, GeneralizedPower[NonCommutativeMultiply, bdag, _]],
          {"bdag", f[[3]]},
        MatchQ[f, Power[var,  _]],
          {"b",    f[[2]]},
        MatchQ[f, Power[bdag, _]],
          {"bdag", f[[2]]},
        True,
          {"unknown", 0}
      ];
      If[typ === "unknown",
        Print["Error: unsupported non-ladder factor: ", f];
        Return[$Failed, Module]
      ];
      If[Length[blocks] > 0 && blocks[[-1, 1]] === typ,
        blocks[[-1, 2]] += pow,
        AppendTo[blocks, {typ, pow}]
      ],
      {f, factors}
    ];

    If[Length[blocks] > 0 && blocks[[-1, 1]] === "bdag",
      AppendTo[blocks, {"b", 0}]];
    If[Length[blocks] > 0 && blocks[[1, 1]] === "b",
      PrependTo[blocks, {"bdag", 0}]];
      
    sList = Reverse[blocks[[2 ;; ;; 2, 2]]];
    rList = Reverse[blocks[[1 ;; ;; 2, 2]]];
    {rList, sList}
  ];


MultiModeBlasiakOrder[exprIn_, vars_List, scalars_List] :=
  Module[{expr, factors, modeVars, modePolys, scalarFactors,
          finalExpr, getOrderWeight},

    expr = exprIn //. {
      GeneralizedPower[NonCommutativeMultiply,
        NonCommutativeMultiply[args__], n_Integer?Positive] :>
          NonCommutativeMultiply @@ Flatten[ConstantArray[{args}, n]],
      Power[NonCommutativeMultiply[args__], n_Integer?Positive] :>
          NonCommutativeMultiply @@ Flatten[ConstantArray[{args}, n]],
      NonCommutativeMultiply[a___, NonCommutativeMultiply[b__], c___] :>
          NonCommutativeMultiply[a, b, c]
    };

    factors = If[Head[expr] === NonCommutativeMultiply, List @@ expr, {expr}];

    scalarFactors = Select[factors, FreeQ[#, Alternatives @@ vars] &];
    factors       = Select[factors, !FreeQ[#, Alternatives @@ vars] &];

    modeVars = Partition[vars, 2];

    modePolys = Map[
      Function[modePair,
        Module[{b = modePair[[1]], bdag = modePair[[2]],
                modeFactors, parsed, rList, sList},
          modeFactors = Select[factors, !FreeQ[#, b] || !FreeQ[#, bdag] &];
          If[Length[modeFactors] == 0, Return[1, Module]];
          parsed = ParseBlasiakMonomial[
                     NonCommutativeMultiply @@ modeFactors, b];
          If[parsed === $Failed, Return[$Failed, Module]];
          {rList, sList} = parsed;
          BlasiakNormalOrder[rList, sList, b]
        ]
      ],
      modeVars
    ];

    If[MemberQ[modePolys, $Failed], Return[$Failed]];

    finalExpr = Fold[
      If[#1 === 1, #2, If[#2 === 1, #1, #1 ** #2]] &,
      1, modePolys
    ];

    finalExpr = finalExpr //. {
      NonCommutativeMultiply[a___, b_ + c_, d___] :>
        (NonCommutativeMultiply[a, b, d] + NonCommutativeMultiply[a, c, d]),
      NonCommutativeMultiply[a___, Times[c_, b__], d___] /;
        FreeQ[c, Alternatives @@ vars] :>
          c * NonCommutativeMultiply[a, b, d],
      NonCommutativeMultiply[a___, c_, d___] /;
        FreeQ[c, Alternatives @@ vars] :>
          c * NonCommutativeMultiply[a, d]
    };

    finalExpr = finalExpr /. {
      NonCommutativeMultiply[x_] :> x,
      NonCommutativeMultiply[]   -> 1
    };

    getOrderWeight[x_] :=
      Module[{isDagger, varIndex},
        isDagger = !FreeQ[x, SuperDagger];
        varIndex = SelectFirst[Range[Length[vars]],
                     !FreeQ[x, vars[[#]]] &, 999];
        If[isDagger, -1000 + varIndex, 1000 + varIndex]
      ];

    finalExpr = finalExpr /. NonCommutativeMultiply[args__] :>
      NonCommutativeMultiply @@ SortBy[{args}, getOrderWeight];

    If[Length[scalarFactors] > 0,
      (Times @@ scalarFactors) * finalExpr,
      finalExpr
    ]
  ];
