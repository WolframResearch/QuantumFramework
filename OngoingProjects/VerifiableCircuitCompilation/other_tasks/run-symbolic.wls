root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"FieldArithmetic.wl"}]];
Get[FileNameJoin[{root,"ComponentMixing.wl"}]];
plans=OtherTasksField`TowerPlan /@ {6,4,3};
six=OtherTasksField`SixProductProof[];
mix=OtherTasksMixing`MixingProofs[];
result=<|"kernel"->$Version,"tower_plans"->plans,"six_product"->six,"mixing"->mix,
 "source_sha256"->AssociationMap[FileHash[FileNameJoin[{root,#}],"SHA256","HexString"]&,
   {"FieldArithmetic.wl","ComponentMixing.wl","run-symbolic.wls"}]|>;
Export[FileNameJoin[{root,"symbolic.json"}],result,"RawJSON"];
Echo[<|"tower_certificates"->(KeyTake[#,{"subfield_degree","basis_certificate","interpolation_certificate","tensor_certificate"}]& /@ plans),
       "six_product"->six["certificate"],"mixing"->mix|>];
If[!And@@Join[Catenate[Lookup[plans,{"basis_certificate","interpolation_certificate","tensor_certificate"}]],
   {six["certificate"]},Select[Values[mix],BooleanQ]],Exit[1]];
