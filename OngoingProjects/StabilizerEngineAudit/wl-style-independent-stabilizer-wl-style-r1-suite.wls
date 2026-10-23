repo = "/Users/mohammadb/Documents/GitHub/QuantumFramework";
auditDir = "/private/tmp";
If[! DirectoryQ[auditDir], CreateDirectory[auditDir]];
phase = "independent-style-r1";
$ProgressReporting = False;
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
SeedRandom[20261022];
rows = Map[Function[file, Module[{elapsed, report, results, row},
    report = TestReport[file]; elapsed = Null;
    results = report["Results"];
    row = <|"File" -> FileNameTake[file], "Count" -> Length[results], "Seconds" -> elapsed,
        "Outcomes" -> Counts[#["Outcome"] & /@ results],
        "Failures" -> Map[Function[r, <|"ID" -> ToString[r["TestID"], InputForm], "Outcome" -> r["Outcome"],
            "Expected" -> ToString[r["ExpectedOutput"], InputForm], "Actual" -> ToString[r["ActualOutput"], InputForm],
            "Messages" -> ToString[r["ActualMessages"], InputForm]|>], Select[results, #["Outcome"] =!= "Success" &]]|>;
    Echo[row]; row]], Sort[FileNames["*.wlt", FileNameJoin[{repo, "Tests", "Stabilizer"}]]]];
Export[FileNameJoin[{auditDir, "suite-" <> phase <> ".json"}],
    <|"Version" -> $Version, "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"], "Rows" -> rows,
        "Messages" -> ToString[$MessageList, InputForm]|>, "RawJSON"];
Exit[If[AllTrue[rows, #Count > 0 && #Failures === {} &], 0, 1]];
