root=DirectoryName[$InputFileName];
Get[FileNameJoin[{DirectoryName[root],"dfthc_feasibility","RankCost.wl"}]];
result=Import[FileNameJoin[{root,"results.json"}],"RawJSON"];
symbolic=Import[FileNameJoin[{root,"symbolic.json"}],"RawJSON"];
hashes=AllTrue[KeyValueMap[FileHash[FileNameJoin[{root,#1}],"SHA256","HexString"]===#2&,
 result["source_sha256"]],TrueQ];
wlHashes=AllTrue[KeyValueMap[FileHash[FileNameJoin[{root,#1}],"SHA256","HexString"]===#2&,
 symbolic["source_sha256"]],TrueQ];
rows=Join[Lookup[result["shape_estimates"],"reference_loader"],Lookup[result["candidates"],"reference_loader"]];
matches=Map[Function[row,With[
 {valid=Select[DFTHCFeasibility`LookupCost[row["shape"],row["beta"],#,row["assumed_lambda_eff"]]& /@ Range[75],#["within_cap"]===True&]},
 With[{best=First@MinimalBy[valid,#["conditional_cost"]&],
       actual=DFTHCFeasibility`LookupCost[row["shape"],row["beta"],row["batch"],row["assumed_lambda_eff"]]},
  AllTrue[{"toffoli","qubits","angle_qubits"},actual[#]===row[#]&] &&
  Abs[actual["conditional_cost"]/row["conditional_cost"]-1]<10^-12 &&
  Abs[best["conditional_cost"]/row["conditional_cost"]-1]<10^-12]]],rows];
incumbentFile=FileNameJoin[{DirectoryName[root],"dfthc_c30","verified.json"}];
buildFile=FileNameJoin[{DirectoryName[root],"dfthc_c30","build.json"}];
incumbent=Import[incumbentFile,"RawJSON"];build=Import[buildFile,"RawJSON"];
unchangedRows=Select[result["candidates"],KeyExistsQ[#,"unchanged_loader_estimate"]&];
unchangedMatches=Map[Function[candidate,With[
 {h=DFTHCFeasibility`PreparationCost[candidate["shape"],candidate["precision"]],
  row=candidate["unchanged_loader_estimate"],old=incumbent["resources"]},
 With[{t=2 Total[Lookup[build["profiles"],"static_toffoli"]]+h["toffoli"],
       q=h["persistent_qubits"]+old["angles_persistent_qubits"]+
         Max[h["temporary_qubits"],old["angles_temporary_qubits"]]},
  TrueQ[t==row["toffoli"]] && q===row["qubits"] &&
  Max[Lookup[build["profiles"],"num_qubits"]]===row["angle_qubits"] &&
  TrueQ[2 Total[Lookup[build["profiles"],"static_toffoli"]]==old["angles_toffoli"]] &&
  Abs[candidate["metrics"]["lambda_eff"] t q/row["cost"]-1]<10^-12]]],unchangedRows];
context=result["recorded_score_context"];
contextHashes=AllTrue[KeyValueMap[
 FileHash[FileNameJoin[{DirectoryName[root],#1}],"SHA256","HexString"]===#2&,
 context["report_sha256"]],TrueQ];
contextTracks=Map[Function[track,With[
 {reports=Import[FileNameJoin[{DirectoryName[root],#}],"RawJSON"]& /@ track["report_paths"]},
 AllTrue[reports,#["status"]==="PASS"&] &&
 reports[[1]]["score"]===track["baseline_cost"] && reports[[2]]["score"]===track["current_cost"] &&
 Abs[track["improvement"]/(reports[[1]]["score"]/reports[[2]]["score"])-1]<10^-12]],
 Values[context["tracks"]]];
combined=(Times@@Lookup[Values[context["tracks"]],"improvement"])^(1/3);
checks=<|"input_hashes_match"->hashes,"wl_hashes_match"->wlHashes,
 "all_whole_task_reference_costs_match"->AllTrue[matches,TrueQ],
 "unchanged_loader_whole_task_costs_match"->AllTrue[unchangedMatches,TrueQ],
 "incumbent_report_hash_matches"->(FileHash[incumbentFile,"SHA256","HexString"]===result["provenance"]["incumbent_report"]),
 "recorded_score_context_hashes_match"->contextHashes,
 "recorded_track_costs_match"->AllTrue[contextTracks,TrueQ],
 "combined_context_is_geometric_mean"->(Abs[combined/result["current_combined_local_estimate"]-1]<10^-12 &&
   result["current_combined_local_estimate"]===context["combined_local_estimate"]),
 "all_cost_rows_within_cap"->AllTrue[rows,#["angle_qubits"]<=128&],
 "no_circuit_score_claimed"->(result["circuit_implementation"]===False && result["new_competition_score"]===Null)|>;
report=<|"status"->If[AllTrue[Values[checks],TrueQ],"PASS","FAIL"],"checks"->checks,
 "cost_rows_checked"->Length[rows],"source_sha256"->FileHash[$InputFileName,"SHA256","HexString"],
 "unchanged_loader_rows_checked"->Length[unchangedRows],
 "static_loader_dependency_sha256"->FileHash[buildFile,"SHA256","HexString"],
 "resource_dependency_sha256"->FileHash[FileNameJoin[{DirectoryName[root],"dfthc_feasibility","RankCost.wl"}],"SHA256","HexString"],
 "results_sha256"->FileHash[FileNameJoin[{root,"results.json"}],"SHA256","HexString"]|>;
Export[FileNameJoin[{root,"wl-data-checks.json"}],report,"RawJSON"];
Echo[report];If[report["status"]=!="PASS",Exit[1]];
