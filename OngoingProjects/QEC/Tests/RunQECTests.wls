#!/usr/bin/env wolframscript

(* Test runner for the QEC project.

   Usage: wolframscript -file OngoingProjects/QEC/Tests/RunQECTests.wls [pattern]
     ... RunQECTests.wls            run every .wlt beside this script
     ... RunQECTests.wls Noise      run only files whose name matches "Noise"

   Why this exists.  The repo's own Tests/RunTests.wls discovers test files with
   FileNames["*.wlt", <its own directory>, Infinity], so it only sees files under
   Tests/.  These live in the project folder instead, which keeps the QEC work in
   one place at the cost of falling outside that sweep.  This script is the local
   equivalent: same paclet loading, same per-file reporting, same green/red verdict.

   When the core migrates into QuantumFramework/Kernel/QEC/, these files should move
   back to Tests/QEC/ and this script goes away -- at that point the paclet's own
   test battery is the right home for them.

   Output carries the same machine-readable tags as the repo runner, so the two can
   be diffed the same way:

     FILE    <tab> file <tab> succeeded <tab> total
     TEST    <tab> file <tab> TestID <tab> outcome     (failures only)
     TOTAL   <tab> succeeded <tab> total
     STATUS  <tab> GREEN or RED
*)

Needs["MUnit`"];

$here = AbsoluteFileName[DirectoryName[$InputFileName]];
$repo = ParentDirectory[$here, 3];
$paclet = FileNameJoin[{$repo, "QuantumFramework"}];

If[ ! DirectoryQ[$paclet],
    Print["  ERROR: no QuantumFramework paclet directory at " <> $paclet];
    Exit[1]
];

(* The .wlt files resolve the core's location from the loaded paclet, so the repo
   checkout has to outrank the installed one before any of them is read. *)
PacletDirectoryLoad[$paclet];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

$pattern = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], ""];

$files = Select[
    Sort @ FileNames["*.wlt", $here],
    $pattern === "" || StringContainsQ[FileNameTake[#], $pattern, IgnoreCase -> True] &
];

If[ $files === {},
    Print["  No test files matched " <> $pattern <> " in " <> $here];
    Exit[1]
];

Print["  Paclet:     " <> $paclet];
Print["  Test files: " <> ToString[Length[$files]]];
Print[""];

row[tag_, cells___] := Print[StringRiffle[{tag, cells}, "\t"]];

runFile[file_] := Module[{name = FileNameTake[file], report, failures},
    report = TestReport[file];
    failures = Select[Values[report["TestResults"]], #["Outcome"] =!= "Success" &];
    row["FILE", name, ToString[report["TestsSucceededCount"]],
        ToString[report["TestsSucceededCount"] + report["TestsFailedCount"]]];
    Function[r, row["TEST", name, ToString[r["TestID"]], r["Outcome"]]] /@ failures;
    {report["TestsSucceededCount"], report["TestsFailedCount"]}
];

$results = runFile /@ $files;
$succeeded = Total[$results[[All, 1]]];
$failed = Total[$results[[All, 2]]];

Print[""];
row["TOTAL", ToString[$succeeded], ToString[$succeeded + $failed]];
row["STATUS", If[$failed === 0, "GREEN", "RED"]];

Exit[If[$failed === 0, 0, 1]]
