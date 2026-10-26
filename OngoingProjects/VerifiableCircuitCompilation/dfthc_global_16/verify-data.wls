root=DirectoryName[$InputFileName];
Get[FileNameJoin[{DirectoryName[root],"dfthc_feasibility","RankCost.wl"}]];
result=Import[FileNameJoin[{root,"results.json"}],"RawJSON"];
symbolic=Import[FileNameJoin[{root,"symbolic.json"}],"RawJSON"];
inputs=Import[FileNameJoin[{root,"inputs.json"}],"RawJSON"];
hashes=AllTrue[KeyValueMap[FileHash[FileNameJoin[{root,#1}],"SHA256","HexString"]===#2&,
 Join[result["source_sha256"],symbolic["source_sha256"]]],TrueQ];
rows=Select[Lookup[result["candidates"],"whole_task_estimate"],AssociationQ];
matches=Map[Function[row,With[
 {actual=DFTHCFeasibility`LookupCost[row["shape"],row["beta"],row["batch"],row["assumed_lambda_eff"]],
  valid=Select[DFTHCFeasibility`LookupCost[row["shape"],row["beta"],#,row["assumed_lambda_eff"]]& /@ Range[75],
    #["within_cap"]===True&]},
  AllTrue[{"toffoli","qubits","angle_qubits"},actual[#]===row[#]&] &&
  Abs[actual["conditional_cost"]/row["conditional_cost"]-1]<10^-12 &&
  Abs[Min[Lookup[valid,"conditional_cost"]]/row["conditional_cost"]-1]<10^-12]],rows];
context=result["recorded_score_context"];
contextHashes=AllTrue[KeyValueMap[FileHash[FileNameJoin[{DirectoryName[root],#1}],"SHA256","HexString"]===#2&,
 context["report_sha256"]],TrueQ];
contextCosts=Map[Function[track,With[
 {reports=Import[FileNameJoin[{DirectoryName[root],#}],"RawJSON"]& /@ track["report_paths"]},
  AllTrue[reports,#["status"]==="PASS"&] &&
  reports[[1]]["score"]===track["baseline_cost"] && reports[[2]]["score"]===track["current_cost"] &&
  Abs[track["improvement"]/(reports[[1]]["score"]/reports[[2]]["score"])-1]<10^-12]],Values[context["tracks"]]];
combined=(Times@@Lookup[Values[context["tracks"]],"improvement"])^(1/3);
checks=<|"source_hashes_match"->hashes,
 "prepared_arrays_hash_matches"->(FileHash[FileNameJoin[{root,"local","prepared.npz"}],"SHA256","HexString"]===inputs["private_arrays_sha256"]),
 "symbolic_checks_pass"->(symbolic["status"]==="PASS" && AllTrue[Values[symbolic["checks"]],TrueQ]),
 "whole_task_costs_and_batches_match"->AllTrue[matches,TrueQ],
 "angle_caps_pass"->AllTrue[rows,#["angle_qubits"]<=128&],
 "actual_normalizations_are_priced"->AllTrue[result["candidates"],
   If[AssociationQ[#["whole_task_estimate"]],
    #["metrics"]["lambda_eff"]===#["whole_task_estimate"]["assumed_lambda_eff"],True]&],
 "recorded_context_hashes_match"->contextHashes,
 "recorded_track_costs_match"->AllTrue[contextCosts,TrueQ],
 "recorded_geometric_mean_matches"->(Abs[combined/context["combined_local_estimate"]-1]<10^-12),
 "floating_offsets_preserved"->AllTrue[result["floating_fits"],Max[Lookup[#,{"offset_matrix_drift","offset_norm_drift"}]]<10^-8&],
 "no_circuit_score_claimed"->(result["circuit_implementation"]===False && result["new_competition_score"]===Null)|>;
report=<|"status"->If[AllTrue[Values[checks],TrueQ],"PASS","FAIL"],"checks"->checks,
 "cost_rows_checked"->Length[rows],"source_sha256"->FileHash[$InputFileName,"SHA256","HexString"],
 "resource_dependency_sha256"->FileHash[FileNameJoin[{DirectoryName[root],"dfthc_feasibility","RankCost.wl"}],"SHA256","HexString"],
 "results_sha256"->FileHash[FileNameJoin[{root,"results.json"}],"SHA256","HexString"]|>;
Export[FileNameJoin[{root,"wl-data-checks.json"}],report,"RawJSON"];
Echo[report];If[report["status"]=!="PASS",Exit[1]];
