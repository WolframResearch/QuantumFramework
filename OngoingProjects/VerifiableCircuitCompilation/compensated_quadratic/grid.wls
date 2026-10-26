(* Close any conservative continuous-envelope gap on the complete binary grid.
   Integer arithmetic clears all denominators before checking relative error. *)
root=DirectoryName[$InputFileName];
saved=Import[FileNameJoin[{root,"coefficients.json"}],"RawJSON"];
If[saved["status"]!="DERIVATION_AND_CAPACITY_PASS" || saved["fractional_bits"]!={30,23,10},Exit[1]];
fs=saved["fractional_bits"];
unproved=Select[saved["certificate"]["intervals"],#["continuous_phase_envelope"]=!=True&];

gridCertificate[row_,index_]:=With[
  {offsets=Range[0,2^14-1],digits=IntegerDigits[Range[0,2^14-1],2,14]},
  With[{inner=row[[2]]-digits.Floor[row[[3]] 2^(fs[[2]]-fs[[3]]-Range[7,20])]},
    With[{outer=row[[1]]-Total[digits Transpose[
      (Floor[inner 2^#]& /@ (fs[[1]]-fs[[2]]-Range[7,20]))],{2}]},
      With[{squared=((64+index)2^14+offsets)outer^2},
        <|"index"->index,"normalized_offsets"->Length[offsets],
          "positive_accumulators"->TrueQ[Min[inner]>0 && Min[outer]>0],
          "relative_phase_bound"->TrueQ[And@@Thread[
            (2^24-1)^2 2^(20+2 fs[[1]]-48)<=squared<=(2^24+1)^2 2^(20+2 fs[[1]]-48)]]|>]]]];

gridCases=gridCertificate[saved["table"][[#["index"]+1]],#["index"]]& /@ unproved;
checks=<|"grid_indices_match_unproved"->(Lookup[gridCases,"index"]==Lookup[unproved,"index"]),
  "all_grid_cases_pass"->AllTrue[gridCases,TrueQ[#["positive_accumulators"]] && TrueQ[#["relative_phase_bound"]]&],
  "covers_all_192_intervals"->(Count[Lookup[saved["certificate"]["intervals"],"continuous_phase_envelope"],True]+Length[gridCases]==192)|>;
report=<|"kernel"->$Version,"coefficient_report_sha256"->FileHash[FileNameJoin[{root,"coefficients.json"}],"SHA256","HexString"],
  "source_sha256"->FileHash[$InputFileName,"SHA256","HexString"],
  "relative_tolerance"->"1/2^24","grid_denominator"->2^20,
  "continuous_intervals"->Count[Lookup[saved["certificate"]["intervals"],"continuous_phase_envelope"],True],
  "grid_cases"->gridCases,"checks"->checks,
  "scope"->"Exact relative error on every normalized 20-bit input in [1,4). Official uint64 phase rounding and synthesis margin are checked by the full circuit verifier.",
  "status"->If[AllTrue[Values[checks],TrueQ],"PASS","FAIL"]|>;
Export[FileNameJoin[{root,"grid-certificate.json"}],report,"RawJSON"];
Echo[KeyTake[report,{"continuous_intervals","grid_cases","checks","status"}]];
If[report["status"]!="PASS",Exit[1]];
