repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
report = TestReport[FileNameJoin[{repo, "Tests", "Stabilizer", "InversePhase.wlt"}]];
results = report["Results"];
failures = Select[results, #["Outcome"] =!= "Success" &];
Echo[<|"Count" -> Length[results], "Outcomes" -> Counts[#["Outcome"] & /@ results],
    "Failures" -> ( <|"ID" -> #["TestID"], "Expected" -> ToString[#["ExpectedOutput"], InputForm],
        "Actual" -> ToString[#["ActualOutput"], InputForm], "Messages" -> ToString[#["ActualMessages"], InputForm]|> & /@ failures),
    "OutsideMessages" -> ToString[$MessageList, InputForm]|>];
Exit[If[Length[results] > 0 && failures === {}, 0, 1]];
