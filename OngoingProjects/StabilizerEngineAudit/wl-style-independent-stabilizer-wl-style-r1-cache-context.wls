$ProgressReporting=False;
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
same[a_,b_]:=a["Signs"]===b["Signs"]&&a["Tableau"]===b["Tableau"];
Echo[Names["Wolfram`QuantumFramework`*`$stabilizerCircuitMemo"],"Cache symbols"];
Echo[Context[stabilizerExecutionPlan],"Plan context"];
Exit[];
