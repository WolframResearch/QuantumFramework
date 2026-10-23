(* Unequal precision Horner arithmetic for inverse-square-root phase synthesis.
   All bounds are rational; approximate discovery is outside this module. *)
BeginPackage["InvSqrtWidths`"];
QuantizedRows::usage = "QuantizedRows[rows,widths] rounds the magnitudes of the frozen 36-bit coefficients to three fractional precisions.";
HornerError::usage = "HornerError[t,e0,e1] is the signed error from truncating each Horner product row.";
CertifyWidths::usage = "CertifyWidths[rows,widths,eta] certifies every interval and register capacity at relative tolerance eta.";
Begin["`Private`"];

validWidths[fs_] := MatchQ[fs,{_Integer,_Integer,_Integer}] && And@@Thread[1<=fs<=36];
validRows[rows_] := MatrixQ[rows,IntegerQ] && Dimensions[rows]=={192,3} &&
  And@@(MatchQ[Sign[#],{1,-1,1}]& /@ rows);
QuantizedRows[rows_,fs_] /; validRows[rows] && validWidths[fs] :=
  Round[# 2^(fs-36)]& /@ Abs[rows];
QuantizedRows[___] := Failure["WidthsDomain",<|"MessageTemplate"->"Expected the signed 192 by 3 integer table and three fractional widths in 1 through 36."|>];

HornerError[t_,e0_,e1_] := Expand[(c0-t (c1-t c2+e1)+e0)-(c0-t (c1-t c2))];

CertifyWidths[rows_,fs_,eta_] /; validRows[rows] && validWidths[fs] &&
    MatchQ[eta,_Integer|_Rational] && 0<eta<1 := With[
  {quantized=QuantizedRows[rows,fs],rowBounds=14/2^fs[[{1,2}]]},
  With[{lengths=IntegerLength[Max /@ Transpose[quantized],2]},
    With[{certificates=MapIndexed[Function[{row,index},With[
      {a=1+(First[index]-1)/64,p=Total[row (-t)^Range[0,2]/2^fs],
       firstStage=(row[[2]]-t row[[3]] 2^(fs[[2]]-fs[[3]]))/2^fs[[2]]},
      <|"index"->First[index]-1,
        "phase"->TimeConstrained[Resolve[ForAll[t,0<=t<=1/64,
          p-t rowBounds[[2]]>0 &&
          (a+t)(p-t rowBounds[[2]])^2 >= (1-eta)^2 &&
          (a+t)(p+rowBounds[[1]])^2 <= (1+eta)^2],Reals],15,"TIMEOUT"],
        "capacity"->TimeConstrained[Resolve[ForAll[t,0<=t<=1/64,
          firstStage>0 && And@@Thread[0<row<2^lengths] &&
          p-t rowBounds[[2]]>0],Reals],15,"TIMEOUT"]|>
      ]],quantized]},
    <|"fractional_bits"->fs,"register_bits"->lengths,"table"->quantized,
      "relative_bound"->ToString[eta,InputForm],
      "row_error_bounds"->(ToString[#,InputForm]& /@ rowBounds),
      "intervals"->certificates,
      "certified"->And@@(TrueQ /@ Catenate[Values[KeyDrop[#,"index"]]& /@ certificates])|>
    ]]];
CertifyWidths[___] := Failure["CertificateDomain",<|"MessageTemplate"->"Expected valid rows, fractional widths and a rational tolerance strictly between zero and one."|>];

End[];
EndPackage[];
