BeginPackage["CompensatedQuadratic`"];
CompensatedRows::usage="CompensatedRows[rows,widths] derives dyadic compensated coefficients from the fixed signed 36-bit table. Returns Failure for invalid input.";
CompensationDerivation::usage="CompensationDerivation[] solves the affine minimax problem exactly and checks the general phase/arithmetic identities.";
CoefficientCertificate::usage="CoefficientCertificate[integerRows,widths] checks all register bounds and a conservative continuous phase envelope. A rejected envelope is not a proof of integer-model failure.";
Begin["`Private`"];
validWidths[fs_]:=MatchQ[fs,{_Integer,_Integer,_Integer}] && And@@Thread[1<=fs<=36];
validRows[rows_]:=MatrixQ[rows,IntegerQ] && Dimensions[rows]=={192,3} &&
  And@@(MatchQ[Sign[#],{1,-1,1}]& /@ rows);

CompensatedRows[rows_,fs_] /; validRows[rows] && validWidths[fs]:=
  Map[Function[row,With[{original=Abs[row]/2^36,h=1/64},
    With[{coarse=Round[original[[3]] 2^fs[[3]]]/2^fs[[3]]},
      With[{delta=coarse-original[[3]]},
        Round[({original[[1]]+delta h^2/8,original[[2]]+delta h,coarse}) 2^fs]]]]],rows];
CompensatedRows[___]:=Failure["CoefficientDomain",<||>];

CompensationDerivation[]:=Module[{u,v,r,t,h,delta,a,b,c,e0,e1,x,y,s,phi},
  With[{solution=Solve[{v==r,1/4+u/2+v==-r,1+u+v==r},{u,v,r},Reals],
    feasible=Reduce[Exists[{u,v},ForAll[t,0<=t<=1,-r<=t^2+u t+v<=r]],r,Reals]},
    <|"affine_linear_and_constant"->(ToString[#,InputForm]& /@ ({u,v}/.First[solution])),
      "minimum_error_radius"->ToString[r/.First[solution],InputForm],
      "checks"-><|
        "solved_affine_compensation"->TrueQ[({u,v,r}/.First[solution])=={-1,1/8,1/8}],
        "global_minimax_radius"->TrueQ[feasible== (r>=1/8)],
        "scaled_polynomial_identity"->TrueQ[Expand[(a+delta h^2/8)-(b+delta h)t+(c+delta)t^2-(a-b t+c t^2)-delta(t^2-h t+h^2/8)]==0],
        "uniform_bound"->TrueQ[Resolve[ForAll[{t,h},h>0 && 0<=t<=h,-h^2/8<=t^2-h t+h^2/8<=h^2/8],Reals]],
        "zero_change_limit"->TrueQ[Limit[delta(t^2-h t+h^2/8),delta->0]==0],
        "signed_Horner_error"->TrueQ[Expand[(a-t(b-t c+e1)+e0)-(a-t(b-t c))-(e0-t e1)]==0],
        "normalization_preserves_relative_error"->TrueQ[FullSimplify[
          y 2^(s-10) Sqrt[x]==y Sqrt[x 4^s/2^20],x>0 && Element[{s,y},Reals]]],
        "unit_phase_modulus"->TrueQ[FullSimplify[Abs[Exp[2 Pi I phi]]^2,Element[phi,Reals]]==1],
        "ties_to_even"->(Round[{1/2,3/2,5/2,7/2}]=={0,2,2,4}),
        "malformed_rows_rejected"->FailureQ[CompensatedRows[{}, {30,23,10}]],
        "malformed_widths_rejected"->FailureQ[CompensatedRows[ConstantArray[{1,-1,1},192],{30,23,0}]]|>|>]];

CoefficientCertificate[rows_,fs_] /; MatrixQ[rows,IntegerQ] && Dimensions[rows]=={192,3} &&
    Min[rows]>0 && validWidths[fs]:=With[
  {lengths=IntegerLength[Max /@ Transpose[rows],2],h=1/64,
    outerShifts=Select[Range[7,20]+fs[[2]]-fs[[1]],Positive]},
  With[{outerError=Total[1-2^-outerShifts]/2^fs[[1]]},
    With[{certificates=MapIndexed[Function[{row,index},Module[{t},
      With[{origin=1+(First[index]-1)/64,p=Total[row (-t)^Range[0,2]/2^fs],
        innerError=Total[FractionalPart[row[[3]]/2^Select[Range[7,20]+fs[[3]]-fs[[2]],Positive]]]/2^fs[[2]]},
        <|"index"->First[index]-1,"inner_error_bound"->ToString[innerError,InputForm],
          "capacity"->TrueQ[And@@Thread[0<row<2^lengths] &&
            row[[2]]/2^fs[[2]]-h row[[3]]/2^fs[[3]]>0 &&
            row[[1]]/2^fs[[1]]-h row[[2]]/2^fs[[2]]-h innerError>0],
          "continuous_phase_envelope"->TimeConstrained[Resolve[ForAll[t,0<=t<=h,
            p-t innerError>0 && (origin+t)(p-t innerError)^2>=(1-2^-24)^2 &&
            (origin+t)(p+outerError)^2<=(1+2^-24)^2],Reals],10,"TIMEOUT"]|>]]],rows]},
      <|"physical_widths"->lengths,"outer_error_bound"->ToString[outerError,InputForm],
        "all_capacity_proved"->AllTrue[Lookup[certificates,"capacity"],TrueQ],
        "all_continuous_phase_envelopes_proved"->AllTrue[Lookup[certificates,"continuous_phase_envelope"],TrueQ],
        "intervals"->certificates|>]]];
CoefficientCertificate[___]:=Failure["CertificateDomain",<||>];
End[];
EndPackage[];
