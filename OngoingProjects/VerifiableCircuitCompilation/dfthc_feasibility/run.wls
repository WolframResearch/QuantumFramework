root = DirectoryName[$InputFileName];
Get[FileNameJoin[{root,"RankCost.wl"}]];
inputs = Import[FileNameJoin[{root,"inputs.json"}],"RawJSON"];
checks = DFTHCFeasibility`SymbolicChecks[];
If[!AllTrue[Values[checks],TrueQ], Echo[checks,"Failed symbolic checks"]; Exit[1]];
matrixFile = FileNameJoin[{root,"local","shifted-pairs.f64"}];
If[FileHash[matrixFile,"SHA256","HexString"] =!= inputs["matrix_sha256"], Exit[2]];
dimension = inputs["pair_dimension"];
pairMatrix = ArrayReshape[BinaryReadList[matrixFile,"Real64",ByteOrdering->-1],{dimension,dimension}];
If[!SymmetricMatrixQ[pairMatrix], Exit[3]];
{spectrum, eigenvectors} = Eigensystem[pairMatrix];
eigenResidual = Norm[pairMatrix.Transpose[eigenvectors] - Transpose[eigenvectors] . DiagonalMatrix[spectrum],"Frobenius"];
orthogonality = Norm[eigenvectors.Transpose[eigenvectors]-IdentityMatrix[dimension],"Frobenius"];
scipySpectrum = Import[FileNameJoin[{root,"local","scipy-eigenvalues.json"}],"RawJSON"];
spectrumDifference = Norm[Sort[spectrum]-Sort[scipySpectrum]];
(* An explicit numerical guard, not a rigorous interval roundoff certificate. *)
guard = 10^-8;
If[Max[eigenResidual,orthogonality,spectrumDifference] >= guard/100, Exit[4]];
numericRank = SelectFirst[Range[0,dimension],DFTHCFeasibility`PSDTail[spectrum,#]<=.020 &];
costRows = inputs["cost_rows"];
costMatches = AllTrue[costRows,Function[row,
  With[{h = DFTHCFeasibility`PreparationCost[row["shape"],row["beta"]],
        model = DFTHCFeasibility`LookupCost[row["shape"],row["beta"],row["batch"],inputs["reference_lambda_eff"]]},
    model["within_cap"] === True && AllTrue[Keys[h],h[#] == row["heuristic"][#]&] &&
    AllTrue[{"toffoli","qubits","angle_qubits"},model[#] == row[#]&] &&
    Abs[model["conditional_cost"]/row["conditional_cost"]-1] < 10^-12]]];
If[!costMatches,Exit[5]];
comparison = Map[Function[shape,With[
  {rank = DFTHCFeasibility`RankCap[shape,dimension],
   rows = Select[costRows,#["shape"]===shape&]},
  With[{bound = DFTHCFeasibility`PSDTail[spectrum,rank]},
    <|"shape"->shape,"rank_cap"->rank,"psd_error_lower_bound"->bound,
      "rank_screen"->Which[bound-guard>.020,"REJECT",bound+guard<.020,"NOT_RULED_OUT",True,"REFINE"],
      "best_reference_loader"->First@MinimalBy[rows,#["conditional_cost"]&],
      "best_beta12_reference_loader"->First@MinimalBy[Select[rows,#["beta"]==12&],#["conditional_cost"]&]|>]]], inputs["shapes"]];
report = <|"status"->"FEASIBILITY_ONLY","kernel"->$Version,"symbolic_checks"->checks,
  "source_sha256"->AssociationMap[FileHash[FileNameJoin[{root,#}],"SHA256","HexString"]&,{"RankCost.wl","run.wls","prepare.py","inputs.json"}],
  "pair_dimension"->dimension,"numerical_guard"->guard,"eigen_residual_frobenius"->eigenResidual,
  "orthogonality_frobenius"->orthogonality,"scipy_spectrum_l2_difference"->spectrumDifference,
  "negative_spectrum_floor"->DFTHCFeasibility`PSDTail[spectrum,dimension],
  "minimum_unconstrained_psd_rank_at_20mHa"->numericRank,
  "bound_at_previous_rank"->DFTHCFeasibility`PSDTail[spectrum,numericRank-1],
  "bound_at_minimum_rank"->DFTHCFeasibility`PSDTail[spectrum,numericRank],
  "cost_rows_crosschecked"->Length[costRows],"cost_model_matches_official"->costMatches,
  "comparison"->comparison|>;
Export[FileNameJoin[{root,"comparison.json"}],report,"RawJSON"];
Echo[KeyDrop[report,{"comparison","source_sha256"}],"Feasibility results"];
report
