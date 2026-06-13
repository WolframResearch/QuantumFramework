Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

dir = "/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer";
files = If[Length[$ScriptCommandLine] > 1,
   {FileNameJoin[{dir, $ScriptCommandLine[[2]]}]},
   FileNames["*.wlt", dir]];

grand = {0, 0, 0};
Do[
  Module[{report, results, passed, failed, errors, total},
    report = TestReport[file];
    results = Values@report["TestResults"];
    passed = Count[results, _?(#["Outcome"] === "Success" &)];
    failed = Count[results, _?(#["Outcome"] === "Failure" &)];
    errors = Count[results, _?(MatchQ[#["Outcome"], "MessagesFailure" | "Error"] &)];
    total = Length[results];
    grand += {passed, failed + errors, total};
    Print[StringPadRight[FileNameTake[file], 36], " ", passed, "/", total,
      If[failed + errors > 0, "  *** FAIL " <> ToString[failed] <> " ERR " <> ToString[errors], "  ok"]];
    Do[
      If[results[[i]]["Outcome"] =!= "Success",
        Print["   [", results[[i]]["Outcome"], "] id=", results[[i]]["TestID"]];
        Print["      Expected: ", InputForm[results[[i]]["ExpectedOutput"]]];
        Print["      Actual:   ", InputForm[results[[i]]["ActualOutput"]]];
        If[results[[i]]["ActualMessages"] =!= {}, Print["      Msgs:     ", results[[i]]["ActualMessages"]]];
      ], {i, Length[results]}];
  ],
  {file, files}];
Print["======================================================="];
Print["TOTAL  passed=", grand[[1]], "  failed/err=", grand[[2]], "  total=", grand[[3]]];
