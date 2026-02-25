#!/usr/bin/env wolframscript

(* Test Coverage Report for QuantumFramework *)
(* Usage: wolframscript -f Tests/RunTests.wls [pattern] *)
(* Examples: *)
(*   wolframscript -f Tests/RunTests.wls           — run all .wlt files *)
(*   wolframscript -f Tests/RunTests.wls Distance   — run only files matching "Distance" *)

$rootDir = DirectoryName[$InputFileName];
$repoDir = FileNameJoin[{$rootDir, ".."}];
$pacletDir = FileNameJoin[{$repoDir, "QuantumFramework"}];

(* Uninstall any repo-installed version to ensure local source takes priority *)
Quiet[PacletUninstall /@ PacletFind["Wolfram/QuantumFramework"]];
PacletDirectoryLoad[$pacletDir];
Needs["Wolfram`QuantumFramework`"];

(* Find test files *)
pattern = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], ""];
testFiles = FileNames["*.wlt", FileNameJoin[{$rootDir}]];
If[pattern =!= "",
    testFiles = Select[testFiles, StringContainsQ[FileNameTake[#], pattern, IgnoreCase -> True] &]
];

If[Length[testFiles] == 0,
    Print["No test files found" <> If[pattern =!= "", " matching \"" <> pattern <> "\"", ""] <> "."];
    Exit[1]
];

Print["======================================================="];
Print["  QuantumFramework Test Coverage Report"];
Print["======================================================="];
Print["  Date: ", DateString[]];
Print["  Paclet: ", $pacletDir];
Print["  Test files: ", Length[testFiles]];
Print[""];

(* Run each test file *)
allResults = {};
Do[
    Print["-------------------------------------------------------"];
    Print["  ", FileNameTake[file]];
    Print["-------------------------------------------------------"];
    report = TestReport[file];
    results = report["TestResults"];

    passed = Length[Select[results, #["Outcome"] === "Success" &]];
    failed = Length[Select[results, #["Outcome"] === "Failure" &]];
    errors = Length[Select[results, #["Outcome"] === "MessagesFailure" || #["Outcome"] === "Error" &]];
    total = Length[results];

    (* Print individual results *)
    Do[
        r = results[id];
        status = Switch[r["Outcome"],
            "Success", "  PASS ",
            "Failure", "  FAIL ",
            _, "  WARN "
        ];
        Print[status, r["TestID"]];
        If[r["Outcome"] =!= "Success",
            Print["      Expected: ", InputForm[r["ExpectedOutput"]]];
            Print["      Actual:   ", InputForm[r["ActualOutput"]]];
            If[r["ExpectedMessages"] =!= {} || r["ActualMessages"] =!= {},
                Print["      Messages: ", r["ActualMessages"]];
            ];
        ],
        {id, Keys[results]}
    ];

    Print[""];
    Print["  Results: ", passed, "/", total, " passed"];
    If[failed > 0, Print["  Failed:  ", failed]];
    If[errors > 0, Print["  Errors:  ", errors]];

    AppendTo[allResults, <|"File" -> FileNameTake[file], "Passed" -> passed, "Failed" -> failed, "Errors" -> errors, "Total" -> total|>],
    {file, testFiles}
];

(* Summary *)
Print[""];
Print["======================================================="];
Print["  Summary"];
Print["======================================================="];

totalPassed = Total[#["Passed"] & /@ allResults];
totalFailed = Total[#["Failed"] & /@ allResults];
totalErrors = Total[#["Errors"] & /@ allResults];
totalTests = Total[#["Total"] & /@ allResults];

Do[
    r = allResults[[i]];
    mark = If[r["Failed"] + r["Errors"] == 0, "PASS", "FAIL"];
    Print["  ", mark, " ", r["File"], ": ", r["Passed"], "/", r["Total"]],
    {i, Length[allResults]}
];

Print[""];
Print["  Total:  ", totalPassed, "/", totalTests, " passed (",
    ToString[Round[100. totalPassed / totalTests, 0.1]], "%)"];
If[totalFailed > 0, Print["  Failed: ", totalFailed]];
If[totalErrors > 0, Print["  Errors: ", totalErrors]];
Print["═══════════════════════════════════════════════════════"];

(* Exit with code *)
Exit[If[totalFailed + totalErrors > 0, 1, 0]]
