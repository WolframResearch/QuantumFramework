(* Exact forward-chain and Frobenius certificates consumed by the compiler. *)
root=DirectoryName[$InputFileName];
feasibility=FileNameJoin[{DirectoryName[root],"other_tasks"}];
Get[FileNameJoin[{feasibility,"FieldArithmetic.wl"}]];
saved=Import[FileNameJoin[{feasibility,"symbolic.json"}],"RawJSON"];
plan=OtherTasksField`TowerPlan[3];
field=FiniteField[2,(#^12+#^6+#^4+#+1)&];
bits[e_] := PadRight[e["Coefficients"],12];
frobenius[p_] := Transpose[bits /@ ((field /@ (2^Range[0,11]))^(2^p))];
powers={1,2,4,8,10,11};
maps=AssociationThread[ToString /@ powers,frobenius /@ powers];
beta[k_] := 2^k-1;
chainIdentity=Simplify[beta[u] 2^v+beta[v]==beta[u+v],
  Element[{u,v},Integers] && u>0 && v>0];
storedExponents={beta[2],beta[4],beta[8] 2^2,beta[10] 2};
inverseExponent=2 (beta[10] 2+beta[1]);
checks=<|
 "same_tower_plan"->MemberQ[saved["tower_plans"],plan],
 "same_six_product"->(saved["six_product"]===OtherTasksField`SixProductProof[]),
 "chain_composition"->TrueQ[chainIdentity],
 "native_scaffold_exponents"->(storedExponents==={3,15,1020,2046}),
 "inverse_exponent"->(inverseExponent===4094),
 "all_nonzero_inverses"->And@@(#^inverseExponent #==field[1]& /@ (field /@ Range[1,4095])),
 "zero_extension"->(field[0]^inverseExponent===field[0]),
 "frobenius_invertible"->And@@(Mod[Inverse[#,Modulus->2].#,2]===IdentityMatrix[12]& /@ Values[maps]),
 "frobenius_period"->(frobenius[12]===IdentityMatrix[12]),
 "frobenius_inverse"->And@@(Mod[maps[ToString[#]].maps[ToString[12-#]],2]===IdentityMatrix[12]& /@ powers),
 "frobenius_composition"->(Mod[maps["1"].maps["1"],2]===maps["2"] &&
   Mod[maps["2"].maps["2"],2]===maps["4"])
|>;
report=<|"status"->If[And@@Values[checks],"PASS","FAIL"],"checks"->checks,
 "frobenius_maps"->maps,"stored_scaffold_exponents"->storedExponents,
 "inverse_exponent"->inverseExponent,
 "source_sha256"->FileHash[$InputFileName,"SHA256","HexString"],
 "field_source_sha256"->FileHash[FileNameJoin[{feasibility,"FieldArithmetic.wl"}],"SHA256","HexString"],
 "symbolic_sha256"->FileHash[FileNameJoin[{feasibility,"symbolic.json"}],"SHA256","HexString"]|>;
Export[FileNameJoin[{root,"certificates.json"}],report,"RawJSON"];
Echo[<|"status"->report["status"],"checks"->checks|>];
If[report["status"]=!="PASS",Exit[1]];
