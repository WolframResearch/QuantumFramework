repo = Replace[Environment["STABILIZER_AUDIT_REPO"], $Failed -> "/Users/mohammadb/Documents/GitHub/QuantumFramework"];
auditDir = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/StabilizerEngineAudit";
If[! DirectoryQ[auditDir], CreateDirectory[auditDir]];
phase = Environment["STABILIZER_AUDIT_PHASE"];
PacletDirectoryLoad[FileNameJoin[{repo, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
SeedRandom[20261022];
bag = Internal`Bag[];
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
(* Gate objects are prepared outside timing. The aggregate circuit label is
   display metadata; omitting its large Composition avoids benchmarking drawing
   label construction instead of the stabilizer engine. Individual labels/matrices remain. *)
makeCircuit[specs_] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ specs), "Label" -> None|>];
PauliStabilizer[2]["ApplyCircuit", {"H" -> 1}];
Scan[Function[n, Module[{ps = PauliStabilizer[n], specs, qc, blocks, splitCircuit, evolved, base, baseCircuit},
    base = makeSpecs[n, 1000]; baseCircuit = makeCircuit[base];
    Scan[Function[m,
        Echo[{"Preparing", n, m}];
        specs = Flatten[ConstantArray[base, m/1000], 1];
        qc = If[m == 1000, baseCircuit, None];
        measure["Direct/n=" <> ToString[n] <> "/g=" <> ToString[m], Function[ps["ApplyCircuit", specs]]];
        If[m == 1000,
            measure["PublicWarm/n=" <> ToString[n] <> "/g=" <> ToString[m], Function[qc[ps]]];
            measure["PublicColdPlan/n=" <> ToString[n] <> "/g=" <> ToString[m], Function[stabilizerClearCaches[]; qc[ps]], False]];
        If[m == 10000,
            blocks = Partition[specs, 1000]; splitCircuit = ConstantArray[baseCircuit, 10];
            measure["TenDirectCalls/n=" <> ToString[n], Function[Fold[#1["ApplyCircuit", #2] &, ps, blocks]]];
            measure["TenPublicCalls/n=" <> ToString[n], Function[Fold[#2[#1] &, ps, splitCircuit]]];
            evolved = ps["ApplyCircuit", specs];
            measure["MaterializeTableau/n=" <> ToString[n], Function[evolved["Tableau"]]];
            Echo[<|"Storage" -> n, "Bytes" -> ByteCount[evolved], "Keys" -> Keys[First[evolved]]|>]
        ]
    ], {1000, 10000}]
]], {128, 1024}];
With[{n = 128, ps = PauliStabilizer[128], specs = makeSpecs[128, 10000]},
    With[{extended = Insert[specs, "V" -> 1, 5000]},
        measure["DirectOneV/n=128/g=10001", Function[ps["ApplyCircuit", extended]]]]
];
Echo[<|"Finished" -> True, "Messages" -> ToString[$MessageList, InputForm]|>];
Exit[];
