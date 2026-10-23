root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"Feasibility.wl"}]];
Needs["FunctionApproximations`"];
designs=Import[FileNameJoin[{root,"designs.json"}],"RawJSON"]["designs"];
phaseVariants=Import[FileNameJoin[{root,"phase-terms.json"}],"RawJSON"]["variants"];

(* Independent arithmetic formulations used by the checks below. *)
xorPolynomial[vars_List] := (1-Times@@(1-2 vars))/2;
parityProduct[vars_List] := Total[(-1)^(Length[#]-1) xorPolynomial[#]& /@
  Subsets[vars,{1,Length[vars]}]]/2^(Length[vars]-1);
squareTruncated[y_,w_,b_] := Total[Map[Function[k,
  If[2 k<=b,BitGet[y,w-k] 2^(b-2 k),0]],Range[w]]] +
  Total[Map[Function[kl,If[Total[kl]-1<=b,
    Times@@(BitGet[y,w-#]& /@ kl) 2^(b-Total[kl]+1),0]],Subsets[Range[w],{2}]]];
hornerInteger[coefs_List,delta_,r_,d_] := Fold[
  #2-Total[Map[Function[m,BitGet[delta,20-m] Floor[#1/2^m]],Range[r+1,Min[d,20]]]]&,
  Last[Abs[coefs]],Reverse[Most[Abs[coefs]]]];
hornerChecks[design_] := With[{r=design["interval_bits"],w=design["coef_bits"],d=design["delta_bits"]},
  And@@Map[Function[interval,With[{coefs=interval["coefs"],a=ToExpression[interval["a"]]},
    And@@Map[Function[delta,With[
      {observed=hornerInteger[coefs,delta,r,d]/2^w,
       exact=coefs . NestList[(delta/2^20) #&,1,design["degree"]]/2^w,
       allowance=ToExpression[design["horner_absolute_bound"]]+ToExpression[design["delta_absolute_bound"]]},
      Abs[observed-exact]<=allowance && 0<observed<1 &&
      N[Abs[Sqrt[a+delta/2^20] observed-1],40]<=ToExpression[design["seed_relative_bound"]]]],
      {0,2^(19-r),2^(20-r)-1}]]],design["intervals"]]];
phaseValue[block_,bits_] := Mod[Total[Map[
  ToExpression[#["turns"]] Mod[Total[bits[[#["bits"]+1]]],2]&,block["terms"]]],1];
pairValue[s_,bits_] := With[{y=FromDigits[Reverse[bits[[1;;17]]],2],
    z=FromDigits[Reverse[bits[[18;;46]]],2]},
  Mod[Total[Table[BitGet[y,17-k] Floor[z/2^(k-1)],{k,1,17}]] 2^(s-39),1]];
bitExamples=Map[IntegerDigits[#,2,47]&,{0,2^47-1,2^46+123456789,2^46+987654321,2^46+2^30+1}];

checks=<|
  "newton_symbolic" -> (InvSqrtFeasibility`NewtonRelativeError[e] === -3 e^2/2-e^3/2),
  "newton_exact_seed" -> (InvSqrtFeasibility`NewtonRelativeError[0]===0),
  "newton_quadratic_limit" -> (Limit[InvSqrtFeasibility`NewtonRelativeError[e]/e^2,e->0]===-3/2),
  "newton_uniform_bound" -> TrueQ[Resolve[ForAll[{e,h},0<=h<=1 && -h<=e<=h,
    Abs[-3 e^2/2-e^3/2]<=3 h^2/2+h^3/2],Reals]],
  "normalization_relative_error" -> TrueQ[FullSimplify[
    ((2^-k p-2^-k/Sqrt[m])/(2^-k/Sqrt[m])) == p Sqrt[m]-1,
    Assumptions->m>0&&Element[{k,p},Reals]]],
  "square_dropped_bits_exact_bound" -> And@@Flatten[Table[
    (y/2^w)^2-squareTruncated[y,w,b]/2^b<=InvSqrtFeasibility`SquareDropBound[w,b],
    {w,2,5},{b,2,10},{y,0,2^w-1}]],
  "square_bound_attained" -> And@@Flatten[Table[
    ((2^w-1)/2^w)^2-squareTruncated[2^w-1,w,b]/2^b === InvSqrtFeasibility`SquareDropBound[w,b],
    {w,2,5},{b,2,10}]],
  "row_truncation_is_not_whole_product" -> (Floor[3/2]+Floor[3/4]===1 && Floor[3 (1/2+1/4)]===2),
  "all_interval_certificates_true" -> And@@Flatten[Map[
    Flatten[Values /@ Lookup[#["intervals"],"certificates"]]&,designs]],
  "all_quantized_horner_examples" -> And@@(hornerChecks /@ designs),
  "all_error_budgets_positive" -> And@@Map[0<ToExpression[#["certified_relative_bound"]]<2^-24&,designs],
  "reject_inaccurate_polynomial" -> (InvSqrtFeasibility`CertifyPolynomial[1,t,1,1/4,2^-24]["relative"]===False),
  "reject_inexact_certificate" -> FailureQ[InvSqrtFeasibility`CertifyPolynomial[0.99-t/2,t,1,1/4,1/100]],
  "reject_empty_certificate_domain" -> FailureQ[InvSqrtFeasibility`CertifyPolynomial[1,t,1,0,1/100]],
  "two_bit_phase_identity" -> (Expand[parityProduct[{a,b}]-a b]===0),
  "three_bit_phase_identity" -> (Expand[parityProduct[{a,b,c}]-a b c]===0),
  "general_phase_identity" -> And@@Table[With[{vars=Array[v,n]},
    Expand[parityProduct[vars]-Times@@vars]===0],{n,1,6}],
  "boolean_anf_full_adder" -> TautologyQ[Equivalent[
    (a&&b)||(c&&Xor[a,b]),BooleanMinimize[(a&&b)||(c&&Xor[a,b]),"ANF"]]],
  "xor_real_polynomial" -> (Last[PolynomialReduce[
    (a+b-2 a b)^2-(a+b-2 a b),{a^2-a,b^2-b},{a,b}]]===0),
  "gf2_is_not_phase_arithmetic" -> (Exp[2 Pi I (1/8) (1+1)]===I && Exp[2 Pi I (1/8) Mod[1+1,2]]===1),
  "symbolic_diagonal_unitary" -> TrueQ[FullSimplify[
    ConjugateTranspose[DiagonalMatrix[Exp[2 Pi I theta {0,0,0,0,0,0,0,1}]]] .
      DiagonalMatrix[Exp[2 Pi I theta {0,0,0,0,0,0,0,1}]]==IdentityMatrix[8],
    Assumptions->Element[theta,Reals]]],
  "absolute_zero_state_phase" -> (Exp[2 Pi I theta parityProduct[{0,0,0}]]===1),
  "compiled_phase_terms_match_retained_rows" -> And@@Map[Function[variant,
    And@@Map[Function[block,And@@Map[Function[bits,
      phaseValue[block,bits] === Mod[If[variant["name"]=="all-exponents",Last[bits],1]
        pairValue[block["s"],bits],1]],bitExamples]],variant["exponent_blocks"]]],phaseVariants],
  "phase_modulo_one_turn" -> (Mod[3/2,1]===1/2 && Exp[2 Pi I 3/2]===Exp[2 Pi I 1/2]),
  "milp_documented_integer_signature" -> (q/.LinearOptimization[q,{q>=7/2,Element[q,Integers]},{q}])==4
|>;
result=<|"runtime"->$Version,"checks"->checks,"passed"->Count[Values[checks],True],
  "failed"->Keys[Select[checks,!TrueQ[#]&]],
  "source_sha256"->AssociationMap[FileHash[FileNameJoin[{root,#}],"SHA256","HexString"]&,
    {"Feasibility.wl","verify-feasibility.wls","run-feasibility.wls","designs.json","phase-terms.json"}]|>;
Export[FileNameJoin[{root,"wl-checks.json"}],result,"RawJSON"];
Echo[KeyTake[result,{"passed","failed"}]];
Exit[If[Length[result["failed"]]===0,0,1]];
