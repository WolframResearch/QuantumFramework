PacletDirectoryLoad["/private/tmp/qf-stabilizer-before/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
p=PauliStabilizer[<|"Signs"->{1,1,1,-1,1,1},"Tableau"->{{{1,1,0,1,1,1},{1,1,0,1,1,0},{1,0,0,1,0,1}},{{1,1,1,0,0,0},{0,0,1,1,0,1},{0,0,0,1,1,0}}}|>];
q=p["Dagger"]["Dagger"]; Echo[{FindFile["Wolfram`QuantumFramework`"],p["Signs"],q["Signs"],p["Tableau"]===q["Tableau"]}]; Echo[$MessageList]; Exit[];
