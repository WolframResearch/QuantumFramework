repo = Replace[Environment["STABILIZER_AUDIT_REPO"], $Failed -> "/Users/mohammadb/Documents/GitHub/QuantumFramework"];
auditDir = "/private/tmp/qf-inverse-review-timings-r1";
If[! DirectoryQ[auditDir], CreateDirectory[auditDir]];
phase = Environment["STABILIZER_AUDIT_PHASE"];
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
SeedRandom[20261022];
bag = Internal`Bag[];
If[FileExistsQ[FileNameJoin[{auditDir,"long-"<>phase<>".json"}]],
 Scan[Internal`StuffBag[bag,#]&, Import[FileNameJoin[{auditDir,"long-"<>phase<>".json"}],"RawJSON"]["Rows"]]];
measure[label_, f_, warm_:True] := Module[{trials, row},
    If[warm, f[]];
    trials = Table[ClearSystemCache[]; Block[{$MessageList = {}},
        With[{timing = AbsoluteTiming[f[]]},
            <|"Seconds" -> First[timing], "Head" -> ToString[Head[Last[timing]], InputForm],
                "Messages" -> ToString[$MessageList, InputForm]|>]], {5}];
    row = <|"Label" -> label, "Trials" -> trials, "MedianSeconds" -> Median[Lookup[trials, "Seconds"]]|>;
    Echo[row]; Internal`StuffBag[bag, row];
    Export[FileNameJoin[{auditDir, "long-" <> phase <> ".json"}],
        <|"Phase" -> phase, "Version" -> $Version, "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"],
            "Rows" -> Internal`BagPart[bag, All]|>, "RawJSON"];
];
makeSpecs[n_, m_] := Table[With[{j = Mod[37 i, n] + 1},
    Switch[Mod[i, 7], 0, "H" -> j, 1, "S" -> j,
        2, "CNOT" -> {j, Mod[j, n] + 1}, 3, "CZ" -> {j, Mod[j, n] + 1},
        4, "X" -> j, 5, "SWAP" -> {j, Mod[j, n] + 1}, 6, SuperDagger["S"] -> j]], {i, m}];

PauliStabilizer[2]["ApplyCircuit",{"H"->1}];
Scan[Function[n,With[{ps=PauliStabilizer[n],specs=Flatten[ConstantArray[makeSpecs[n,1000],10],1]},
 measure["DirectColdEncoding/n="<>ToString[n]<>"/g=10000",Function[stabilizerClearCaches[];ps["ApplyCircuit",specs]],False]]],{128,1024}];
SeedRandom[87923];
With[{p=PauliStabilizer[128],specs=Table[With[{w=RandomSample[Range[128],2],g=RandomChoice[{"H","S","X","Y","Z","CNOT","CZ","SWAP"}]},g->If[MemberQ[{"CNOT","CZ","SWAP"},g],w,First[w]]],{10000}]},
 measure["DirectRandomCold/n=128/g=10000",Function[stabilizerClearCaches[];p["ApplyCircuit",specs]],False]];
Echo[<|"Finished"->True,"Messages"->ToString[$MessageList,InputForm]|>]; Exit[];
