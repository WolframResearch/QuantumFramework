BeginPackage["ProjectiveDecoderRepresentations`"];
ProjectiveChecks::usage="ProjectiveChecks[] checks inverse-free projective locator synthesis and scale-invariant Forney values. It does not build a reversible decoder.";
Begin["`Private`"];

(* Inversionless BM. State is {locator,previous,order,gap,previousDiscrepancy}.
   Keeping the locator unnormalized avoids all divisions in this recurrence. *)
projectiveStep[syndromes_,z_,zero_][state_,n_]:=With[
  {locator=state[[1]],previous=state[[2]],order=state[[3]],gap=state[[4]],gamma=state[[5]]},
  With[{coefficients=CoefficientList[locator,z]},
    With[{delta=coefficients.PadRight[Reverse[Take[syndromes,n+1]],Length[coefficients],zero]},
      With[{next=Expand[gamma locator-delta z^gap previous]},
        If[delta!=zero && 2 order<=n,
          {next,locator,n+1-order,1,delta},
          {next,previous,order,gap+1,gamma}]]]]];

locatorCase[field_,locations_List,weights_List]:=Module[
  {z,syndromes,states,locator,expected,evaluator,one=field[1],zero=field[0]},
  syndromes=Table[Total[weights locations^n],{n,0,69}];
  states=FoldList[projectiveStep[syndromes,z,zero],{one,one,0,1,one},Range[0,69]];
  locator=First[Last[states]];
  expected=Expand[Times@@(one-# z& /@ locations)];
  evaluator=PolynomialRemainder[locator Total[syndromes z^Range[0,69]],z^70,z];
  <|"locator"->TrueQ[Expand[locator/(locator/.z->zero)-expected]==zero],
    "Forney"->And@@MapThread[Function[{location,weight},TrueQ[
      location (evaluator/.z->1/location)/(D[locator,z]/.z->1/location)==weight]],{locations,weights}],
    "degree_bound"->AllTrue[states,Exponent[First[#],z]<=35&],
    "nonzero_scale"->AllTrue[states,TrueQ[(First[#]/.z->zero)!=zero]&]|>
];

ProjectiveChecks[]:=Module[{a,b,c,p,q,d,z,g,gap,om,sd,field,full,single,zeroSyndromes,zeroState},
  field=FiniteField[2,(#^12+#^6+#^4+#+1)&];
  full=locatorCase[field,field /@ Range[2,36],field /@ Range[3,37]];
  single=locatorCase[field,{field[7]},{field[19]}];
  zeroSyndromes=ConstantArray[field[0],70];
  zeroState=Fold[projectiveStep[zeroSyndromes,z,field[0]],{field[1],field[1],0,1,field[1]},Range[0,69]];
  <|"projective_general_update_identity"->TrueQ[Expand[
      (c b)(a p)-(a d) z (c q)-(c b a)(p-d/b z q)]==0],
    "projective_Forney_scale_cancels"->TrueQ[Cancel[(a om)/(a sd)-om/sd]==0],
    "major_shifted_matrix_determinant"->TrueQ[Det[{{g,-d z^gap},{1,0}}]==d z^gap],
    "minor_shifted_matrix_determinant"->TrueQ[Det[{{g,-d z^gap},{0,1}}]==g],
    "major_formal_rational_inverse"->TrueQ[Cancel[(g p-(g p-d z^gap q))/(d z^gap)-q]==0],
    "truncated_major_coefficient_map_can_be_singular"->TrueQ[
      Det[{{g,0,0,0},{0,g,-d,0},{1,0,0,0},{0,1,0,0}}]==0],
    "native_field_35_errors_locator"->full["locator"],
    "native_field_35_errors_Forney"->full["Forney"],
    "native_field_35_errors_degree_bound"->full["degree_bound"],
    "native_field_35_errors_nonzero_scale"->full["nonzero_scale"],
    "projective_one_error_case"->And@@Values[single],
    "projective_zero_error_case"->TrueQ[First[zeroState]==field[1]],
    "zero_projective_scale_is_invalid"->TrueQ[(a sd/.a->0)==0],
    "scope"->"Classical algebraic recurrence; quantum branch history, inversion costs of reversible updates, and native syndrome adapter are not implemented."
  |>
];
End[];
EndPackage[];
