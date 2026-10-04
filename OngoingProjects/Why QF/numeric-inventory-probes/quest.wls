PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"];
Print["@@ QuEST backend available?                      ", TrueQ[Quiet[Wolfram`QuantumFramework`PackageScope`QuESTLinkQ[]]]];
