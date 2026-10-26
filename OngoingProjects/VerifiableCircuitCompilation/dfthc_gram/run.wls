root=DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"GramCores.wl"}]];
metadata=Import[FileNameJoin[{root,"inputs.json"}],"RawJSON"];
arrayFile=FileNameJoin[{root,"local","arrays.json"}];
If[FileHash[arrayFile,"SHA256","HexString"]=!=metadata["private_arrays_sha256"],Exit[1]];
arrays=Import[arrayFile,"RawJSON"];
symbolic=DFTHCGram`GramSymbolicChecks[];
offsetChecks=DFTHCGram`OffsetCompressionChecks[];
rescaleChecks=DFTHCGram`SchurRescaleChecks[];
If[!AllTrue[Values[symbolic["checks"]],TrueQ],Echo[symbolic];Exit[2]];
If[!AllTrue[Values[offsetChecks],TrueQ],Echo[offsetChecks];Exit[2]];
If[!AllTrue[Values[rescaleChecks],TrueQ],Echo[rescaleChecks];Exit[2]];
pairs=Select[Tuples[Range[76],2],First[#]>=Last[#]&];
pairRows=First/@pairs;pairCols=Last/@pairs;
scale=Sqrt[2-Boole[Equal@@#]]& /@ pairs;
coordinates[u_]:=Map[# scale&,u[[All,pairRows]] u[[All,pairCols]]];
records=MapThread[Function[{u,w,v},With[
 {a=Join[w,{v}],g=DFTHCGram`ProjectorGram[u,u],p=coordinates[u]},
 With[{x=Transpose[w].p,m=w.Transpose[w],q=w.v},
  <|"u_rank_mod_prime"->MatrixRank[SetPrecision[u,Infinity],Modulus->metadata["prime"]],
    "augmented_rank_mod_prime"->MatrixRank[SetPrecision[a,Infinity],Modulus->metadata["prime"]],
    "has_full_support_copy"->AnyTrue[Transpose[w],Function[column,AllTrue[column,#!=0&]]],
    "u_min_singular"->Min[SingularValueList[u,Tolerance->0]],
    "augmented_min_singular"->Min[SingularValueList[a,Tolerance->0]],
    "projector_gram_error"->Norm[p.Transpose[p]-g,"Frobenius"],
    "physical_core_norm_error"->Abs[Total[Flatten[(x.Transpose[x])^2]]-Tr[g.m.g.m]],
    "offset_vector_error"->Norm[Transpose[x].v-Transpose[p].q],
    "offset_norm_error"->Abs[Last[Last[DFTHCGram`AugmentedCore[w,v]]]-v.v]|>]]],
 {arrays["u"],arrays["w"],arrays["v"]}];
checks=<|"all_65_projector_frames_full_row_rank"->AllTrue[Lookup[records,"u_rank_mod_prime"],#==60&],
 "all_65_augmented_cores_rank_30"->AllTrue[Lookup[records,"augmented_rank_mod_prime"],#==30&],
 "all_65_have_a_factor_of_rank_60"->AllTrue[Lookup[records,"has_full_support_copy"],TrueQ],
 "all_65_reconstruction_checks"->AllTrue[records,Max[Lookup[#,
   {"projector_gram_error","physical_core_norm_error","offset_vector_error","offset_norm_error"}]]<10^-8&]|>;
(* A convex Schur epigraph fixture, with a separately derived exact optimum. *)
schurMatrix={{2,1,1},{1,3,-1},{1,-1,tau}};
exactMinimum=Minimize[{tau,And@@Thread[Flatten[Table[
  Det[schurMatrix[[#,#]]]& /@ Subsets[Range[3],{j}],{j,3}]]>=0]},tau,Reals];
sdpValue=SemidefiniteOptimization[tau,
 VectorGreaterEqual[{schurMatrix,0},{"SemidefiniteCone",3}],{tau},"PrimalMinimumValue"];
sdpCheck=NumericQ[sdpValue] && Abs[sdpValue-7/5]<10^-6;
gradientCase=<|"u"->{{3/5,4/5},{0,1}},"w"->{{1,2},{3,-1}},
 "u0"->{{1,0},{5/13,12/13}},"w0"->{{2,-1},{1,3}}|>;
gradient=DFTHCGram`ProjectorGradient@@Lookup[gradientCase,{"u","w","u0","w0"}];
report=<|"status"->If[AllTrue[Values[checks],TrueQ] && sdpCheck && First[exactMinimum]===7/5,"PASS","FAIL"],
 "kernel"->$Version,"symbolic"->symbolic,"offset_compression_checks"->offsetChecks,
 "schur_rescale_checks"->rescaleChecks,"data_checks"->checks,"core_records"->records,
 "rank_scope"->"Exact modular certificates for frozen Float64 arrays interpreted as rationals; not exact unrounded-trigonometric rank certificates.",
 "convex_schur_fixture"-><|"exact_minimum"->"7/5","numeric_sdp_minimum"->sdpValue,"passed"->sdpCheck|>,
 "gradient_fixture"->N[Append[gradientCase,"gradient"->gradient]],
 "source_sha256"->AssociationMap[FileHash[FileNameJoin[{root,#}],"SHA256","HexString"]&,
   {"GramCores.wl","run.wls","study.py","inputs.json"}]|>;
Export[FileNameJoin[{root,"symbolic.json"}],report,"RawJSON"];
Echo[KeyTake[report,{"status","data_checks","convex_schur_fixture"}]];
If[report["status"]=!="PASS",Exit[3]];
