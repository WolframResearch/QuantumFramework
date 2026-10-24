(* Orthogonal rotations of SOS ensemble components; shifts and vectors stay fixed. *)
BeginPackage["OtherTasksMixing`"];
PairRotation::usage="PairRotation[t] is an exact rational orthogonal two-component rotation.";
MixingProofs::usage="MixingProofs[] proves the two-body, offset and energy invariants and pair objective.";
Begin["`Private`"];
PairRotation[t_] /; (Head[t]===Symbol && !MemberQ[{True,False,Indeterminate},t]) ||
   (NumericQ[t] && TrueQ[Element[t,Reals]]) := {{1-t^2,2 t},{-2 t,1-t^2}}/(1+t^2);
PairRotation[___] := Failure["RotationDomain",<|"MessageTemplate"->"Expected a real numeric parameter or a symbolic real parameter."|>];
MixingProofs[] := Module[{rotation=PairRotation[t],x,y,w,v,z,u,newX,newW,
  pairObjective,expectedPair},
 {newX,newW}={rotation.{x,y},rotation.{w,v}};
 pairObjective=TrigReduce[(a Cos[theta]+b Sin[theta])^2+
   (-c Sin[theta]+d Cos[theta])^2];
 expectedPair=(a^2+b^2+c^2+d^2)/2+
   (a^2+d^2-b^2-c^2) Cos[2 theta]/2+(a b-c d) Sin[2 theta];
 <|"orthogonal"->TrueQ[FullSimplify[rotation.Transpose[rotation]==IdentityMatrix[2],Element[t,Reals]]],
   "two_body_bilinear"->TrueQ[Together[
      (rotation.{x,y}).(rotation.{z,u})-(x z+y u)]===0],
   "offset_contraction"->TrueQ[Together[newX.newW-(x w+y v)]===0],
   "offset_squared_norm"->TrueQ[Together[newW.newW-(w^2+v^2)]===0],
   "identity_rotation"->TrueQ[PairRotation[0]===IdentityMatrix[2]],
   "large_parameter_limit"->TrueQ[Limit[rotation,t->Infinity]===-IdentityMatrix[2]],
   "pair_objective"->TrueQ[TrigExpand[pairObjective-expectedPair]===0],
   "pair_objective_expression"->ToString[expectedPair,InputForm],
   "complex_parameter_rejected"->FailureQ[PairRotation[I]],
   "exact_nontrivial_gain"->FullSimplify[
     Total[Total[Abs[{{1,1},{1,-1}}]]^2]/
     Total[Total[Abs[{{1,1},{1,-1}}.{{1,-1},{1,1}}/Sqrt[2]]]^2]],
   "physical_assumption"->"Real orthogonal transformations; rotate native offsets with weights."
 |>
];
End[];
EndPackage[];
