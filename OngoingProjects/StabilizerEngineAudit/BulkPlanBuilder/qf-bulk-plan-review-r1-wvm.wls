PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
$stabilizerCompileTarget = "WVM";
rvReport = TestReport["/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer/BulkPlanBuilder.wlt"];
Echo[{$Version, $stabilizerCompileTarget, rvReport["TestsSucceededCount"], rvReport["TestsFailedCount"]}, "WVM_RESULT"];
Exit[If[rvReport["TestsFailedCount"] === 0 && rvReport["TestsSucceededCount"] === 53, 0, 1]];
