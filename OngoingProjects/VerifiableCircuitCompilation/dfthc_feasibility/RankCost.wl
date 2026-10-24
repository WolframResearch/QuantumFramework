BeginPackage["DFTHCFeasibility`"];
PSDTail::usage = "PSDTail[eigenvalues,k] gives the minimum Frobenius error of a PSD matrix of rank at most k. Values must be real numeric scalars; invalid input returns Failure.";
RankCap::usage = "RankCap[{R,B,C},d] gives Min[d,R Min[B,C]], a necessary DFTHC pair-rank cap.";
PreparationCost::usage = "PreparationCost[{R,B,C},beta] reproduces the official N=76, chi=10 state-preparation and Givens resource model.";
LookupCost::usage = "LookupCost[shape,beta,J,lambda] prices an ordinary unary angle loader plus the full official heuristic resources. Lambda must be a nonnegative real numeric scalar. It is a conditional estimate, not a compiled circuit score.";
SymbolicChecks::usage = "SymbolicChecks[] checks packing, the projector Gram identity, PSD truncation, limits, degeneracy, and invalid input.";

Begin["`Private`"];
positiveIntegerQ[x_] := IntegerQ[x] && x > 0;
shapeQ[s_] := MatchQ[s, {_?positiveIntegerQ, _?positiveIntegerQ, _?positiveIntegerQ}];
realNumberQ[x_] := NumericQ[x] && TrueQ[Im[x] == 0] && FreeQ[x, Indeterminate | _DirectedInfinity];

PSDTail[values_List, k_Integer?NonNegative] /; AllTrue[values, realNumberQ] :=
  With[{realValues = Re[values]},
    Sqrt[Total[Select[realValues, Negative]^2] +
      Total[Drop[Sort[Select[realValues, Positive], Greater], Min[k, Count[realValues, _?Positive]]]^2]]];
PSDTail[___] := Failure["InvalidSpectrum", <||>];
RankCap[s_?shapeQ, d_Integer?positiveIntegerQ] := Min[d, s[[1]] Min[s[[2]], s[[3]]]];
RankCap[___] := Failure["InvalidShape", <||>];

qrom[length_, width_] := With[{exponent = Max[0, Log[2, length/width]/2]},
  With[{k = First@MinimalBy[DeleteDuplicates[{Floor[exponent], Ceiling[exponent]}],
      length/2^# + width (2^# - 1) &]},
    {k, Ceiling[length/2^k + width (2^k - 1)]}]];

PreparationCost[s_?shapeQ, beta_Integer?positiveIntegerQ] := Module[
  {r, b, c, outerSize, selectBits, b1, b2, k1, k2, kToff, outerT, innerT,
   outerPersistent, innerPersistent, outerTemp, innerTemp, givensT},
  {r, b, c} = s;
  outerSize = 76 + r c; selectBits = Ceiling[Log[2, outerSize]];
  b1 = selectBits + 10; b2 = Ceiling[Log[2, b]] + 11;
  k1 = First@qrom[outerSize, b1];
  k2 = SelectFirst[Reverse@Range[0, First@qrom[r c b, b2]],
    Ceiling[Log[2, r c Ceiling[b/2^#]]] + (2^# - 1) b2 < beta 75 &];
  If[MissingQ[k2], Return[Failure["ScratchBudget", <|"Shape" -> s, "Precision" -> beta|>]]];
  kToff = Min[k2, Floor[Log[2, b]]];
  outerT = 8 selectBits + Last@qrom[outerSize, b1] + Last@qrom[outerSize, 1] + 2 b1 + 1;
  innerT = r c Ceiling[(b + 1)/2^kToff] + 2^kToff b2 +
    Last@qrom[r c, b + 1] + b + 1 + 76 + b2 + 8 Ceiling[Log[2, b + 1]];
  givensT = 150 (2 beta - 1);
  outerPersistent = b1 + 10 + selectBits + 2;
  innerPersistent = b2 + 10 + Ceiling[Log[2, r]] + 2 + Ceiling[Log[2, b + 1]] + 2;
  outerTemp = Ceiling[Log[2, outerSize/2^k1]] + 2^k1 b1 + Max[0, selectBits - 2];
  innerTemp = Max[0, Ceiling[Log[2, b + 1]] - 2] +
    Ceiling[Log[2, r c Ceiling[(b + 1)/2^k2]]] + (2^k2 - 1) b2;
  <|"outer_state_prep_toffoli" -> outerT, "inner_state_prep_toffoli" -> innerT,
    "givens_rotations_toffoli" -> givensT,
    "persistent_qubits" -> outerPersistent + innerPersistent,
    "temporary_qubits" -> Max[outerTemp, innerTemp],
    "toffoli" -> outerT + innerT + givensT,
    "qubits" -> outerPersistent + innerPersistent + Max[outerTemp, innerTemp]|>
];
PreparationCost[___] := Failure["InvalidShapeOrPrecision", <||>];

LookupCost[s_?shapeQ, beta_Integer?positiveIntegerQ, j_Integer?positiveIntegerQ, lambda_?realNumberQ] /; TrueQ[Re[lambda] >= 0] :=
  With[{a = 76 + s[[1]] s[[2]], h = PreparationCost[s, beta]},
    If[FailureQ[h], h, With[{selection = Ceiling[Log[2, a]], t = 2 Ceiling[75/j] (a - 2) + h["toffoli"]},
      With[{q = h["persistent_qubits"] + j beta + Max[h["temporary_qubits"], selection - 1]},
        <|"toffoli" -> t, "qubits" -> q, "angle_qubits" -> 2 selection - 1 + j beta,
          "within_cap" -> (2 selection - 1 + j beta <= 128),
          "conditional_cost" -> Re[lambda] t q|>]]]];
LookupCost[___] := Failure["InvalidCostInput", <||>];

SymbolicChecks[] := Module[{x, y, packedX, packedY, u, v, a, b, c, epsilon, diagonal, truncated,
  f0, f1, f2, f3, low0, low1, fourAddress, fourCases, signedWords},
  x = {{a, b}, {b, c}}; y = {{x11, x12}, {x12, x22}};
  packedX = {a, Sqrt[2] b, c}; packedY = {x11, Sqrt[2] x12, x22};
  u = {u1, u2, u3}; v = {v1, v2, v3};
  diagonal = DiagonalMatrix[{a, b, -c}]; truncated = DiagonalMatrix[{a, 0, 0}];
  fourAddress = f0 + low0 (f0+f1) + low1 (f0+f2) + low0 low1 (f0+f1+f2+f3);
  fourCases = MapThread[PolynomialMod[(fourAddress /. Thread[{low0,low1}->#1])-#2,2]&,
    {{{0,0},{1,0},{0,1},{1,1}},{f0,f1,f2,f3}}];
  signedWords = Mod[Range[-256,255],1024];
  <|"pair_inner_product" -> TrueQ[Expand[packedX.packedY - Tr[x.y]] == 0],
    "tensor_norm_preserved" -> TrueQ[Expand[(packedX.packedX)^2 - Total[Flatten[TensorProduct[x,x]]^2]] == 0],
    "projector_gram" -> TrueQ[Expand[Tr[Outer[Times,u,u].Outer[Times,v,v]] - (u.v)^2] == 0],
    "row_sign_preserves_projector" -> (Outer[Times,-u,-u] === Outer[Times,u,u]),
    "four_address_boolean_identity" -> AllTrue[fourCases,#===0&],
    "signed_word_redundant_bit" -> (BitXor[BitShiftRight[signedWords,9],BitAnd[BitShiftRight[signedWords,8],1]] === ConstantArray[0,512]),
    "symbolic_psd_tail" -> TrueQ[FullSimplify[
      Norm[diagonal-truncated,"Frobenius"]^2 == b^2+c^2,
      Assumptions -> a >= b >= 0 && c >= 0 && Element[{a,b,c},Reals]]],
    "decoupling_limit" -> TrueQ[Limit[Sqrt[epsilon^2 b^2+c^2],epsilon->0,Assumptions->c>0] == c],
    "mixed_sign_spectrum" -> (PSDTail[{5,3,-2},1] === Sqrt[13]),
    "negative_floor" -> (PSDTail[{5,3,-2},3] === 2),
    "zero_rank" -> (PSDTail[{5,3,-2},0] === Sqrt[38]),
    "degenerate_spectrum" -> (PSDTail[{2,2,0},1] === 2),
    "empty_spectrum" -> (PSDTail[{},0] === 0),
    "exact_radical_spectrum" -> (PSDTail[{Sqrt[2]},0] === Sqrt[2]),
    "exact_constant_spectrum" -> (PSDTail[{Pi},0] === Pi),
    "exact_mixed_ordering" -> (PSDTail[{Sqrt[2],3/2},1] === Sqrt[2]),
    "zero_imaginary_positive" -> TrueQ[PSDTail[{Complex[1.,0.]},0] == 1],
    "zero_imaginary_negative" -> TrueQ[PSDTail[{Complex[-1.,0.]},0] == 1],
    "high_precision_zero_imaginary" -> TrueQ[PSDTail[{Complex[SetPrecision[1.,40],SetPrecision[0.,40]]},0] == 1],
    "zero_imaginary_negative_floor" -> TrueQ[PSDTail[{Complex[SetPrecision[-1.,40],SetPrecision[0.,40]]},10] == 1],
    "zero_imaginary_cost" -> (Head[LookupCost[{1,1,1},10,1,Complex[1.,0.]]["conditional_cost"]] === Real),
    "malformed_spectrum" -> FailureQ[PSDTail[{1,I},1]],
    "negative_rank" -> FailureQ[PSDTail[{1},-1]],
    "invalid_shape" -> FailureQ[RankCap[{0,2,3},10]],
    "angle_cap_exact" -> With[{resource = LookupCost[{65,60,30},7,15,1]},
      resource["angle_qubits"] === 128 && resource["within_cap"] === True],
    "angle_cap_rejected" -> With[{resource = LookupCost[{65,60,30},10,11,1]},
      resource["angle_qubits"] === 133 && resource["within_cap"] === False],
    "negative_normalization_rejected" -> FailureQ[LookupCost[{1,1,1},10,1,-1]],
    "zero_normalization_limit" -> (LookupCost[{1,1,1},10,1,0]["conditional_cost"] === 0),
    "impossible_scratch_budget" -> FailureQ[PreparationCost[{2^75,1,1},1]],
    "propagate_scratch_failure" -> FailureQ[LookupCost[{2^75,1,1},1,1,1]],
    "span_rank_cap" -> (RankCap[{3,2,5},100] === 6),
    "power_two_B_cost" -> (PreparationCost[{3,32,4},16]["toffoli"] > 0)|>
];
End[];
EndPackage[];
