PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
SeedRandom[20261022];
out = Replace[Environment["STABILIZER_BULK_OUTPUT"], $Failed -> FileNameJoin[{DirectoryName[$InputFileName], "paired.json"}]];
currentDefinition = DownValues[stabilizerCircuitSpecs];
primitiveInitializer = OwnValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerPrimitiveDispatch];
With[{element = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerElementSpecs},
    beforeDefinition = {HoldPattern[stabilizerCircuitSpecs[qco_QuantumCircuitOperator]] :>
        Catenate[element /@ qco["Operators"]]}
];
makeSpecs[n_, m_] := Table[With[{j = Mod[37 i, n] + 1},
    Switch[Mod[i, 7], 0, "H" -> j, 1, "S" -> j,
        2, "CNOT" -> {j, Mod[j, n] + 1}, 3, "CZ" -> {j, Mod[j, n] + 1},
        4, "X" -> j, 5, "SWAP" -> {j, Mod[j, n] + 1}, 6, SuperDagger["S"] -> j]], {i, m}];
randomSpecs[n_, m_] := Table[With[{gate = RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"}], wires = RandomSample[Range[n], 2]},
    gate -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, gate], wires, First[wires]]], {m}];
makeCircuit[specs_] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ specs), "Label" -> None|>];
publicOnce[] := First[trialCircuits][trialState];
publicTen[] := Fold[#2[#1] &, trialState, trialCircuits];
directOnce[] := trialState["ApplyCircuit", trialSpecs];
rows = {};
measure[label_, fn_, cold_, initializePrimitives_:False] := Module[{trials, row},
    trials = Catenate @ Table[
        Map[Function[mode,
            DownValues[stabilizerCircuitSpecs] = If[mode === "Before", beforeDefinition, currentDefinition];
            stabilizerClearCaches[];
            If[!cold, fn[]];
            If[initializePrimitives && mode === "After",
                OwnValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerPrimitiveDispatch] = primitiveInitializer];
            ClearSystemCache[];
            Block[{$MessageList = {}}, With[{timing = AbsoluteTiming[fn[]]},
                <|"Mode" -> mode, "Trial" -> trial, "Seconds" -> First[timing],
                    "Signature" -> ToString[Hash[{Last[timing]["Tableau"], Last[timing]["Signs"]}, "SHA256"], InputForm],
                    "Head" -> ToString[Head[Last[timing]], InputForm], "Messages" -> ToString[$MessageList, InputForm]|>
            ]]
        ], If[OddQ[trial], {"Before", "After"}, {"After", "Before"}]], {trial, 7}];
    row = <|"Label" -> label, "Trials" -> trials,
        "BeforeMedianSeconds" -> Median[Lookup[Select[trials, #Mode === "Before" &], "Seconds"]],
        "AfterMedianSeconds" -> Median[Lookup[Select[trials, #Mode === "After" &], "Seconds"]]|>;
    rows = Append[rows, row]; Echo[KeyDrop[row, "Trials"]];
    Export[out, <|"Version" -> $Version, "CompileTarget" -> $stabilizerCompileTarget,
        "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"], "Rows" -> rows,
        "OutsideMessages" -> ToString[$MessageList, InputForm]|>, "RawJSON"]
];
PauliStabilizer[2]["ApplyCircuit", {"H" -> 1}];
stabilizerCircuitSpecs[QuantumCircuitOperator[{"H" -> 1}]];
Scan[Function[n,
    trialState = PauliStabilizer[n];
    trialSpecs = makeSpecs[n, 1000];
    trialCircuits = {makeCircuit[trialSpecs]};
    First[trialCircuits]["FullInputOrder"];
    measure["PublicCold/n=" <> ToString[n] <> "/g=1000", publicOnce, True];
    measure["PublicWarm/n=" <> ToString[n] <> "/g=1000", publicOnce, False];
    trialCircuits = ConstantArray[First[trialCircuits], 10];
    measure["TenPublicCold/n=" <> ToString[n] <> "/g=10000", publicTen, True];
    measure["TenPublicFirstPlan/n=" <> ToString[n] <> "/g=10000", publicTen, True, True];
    measure["TenPublicWarm/n=" <> ToString[n] <> "/g=10000", publicTen, False];
    trialSpecs = Catenate[ConstantArray[trialSpecs, 10]];
    measure["DirectCold/n=" <> ToString[n] <> "/g=10000", directOnce, True];
    measure["DirectWarm/n=" <> ToString[n] <> "/g=10000", directOnce, False];
    trialSpecs = randomSpecs[n, 1000];
    trialCircuits = {makeCircuit[trialSpecs]};
    First[trialCircuits]["FullInputOrder"];
    measure["RandomPublicCold/n=" <> ToString[n] <> "/g=1000", publicOnce, True];
    measure["RandomPublicWarm/n=" <> ToString[n] <> "/g=1000", publicOnce, False];
    trialCircuits = ConstantArray[First[trialCircuits], 10];
    measure["TenRandomPublicCold/n=" <> ToString[n] <> "/g=10000", publicTen, True];
    measure["TenRandomPublicWarm/n=" <> ToString[n] <> "/g=10000", publicTen, False];
], {128, 1024}];
DownValues[stabilizerCircuitSpecs] = currentDefinition;
Echo[<|"Finished" -> True, "OutsideMessages" -> ToString[$MessageList, InputForm]|>];
Exit[If[AllTrue[Catenate[Lookup[rows, "Trials"]], #Head === "PauliStabilizer" && #Messages === "{}" &] &&
    AllTrue[rows, Length[DeleteDuplicates[Lookup[#Trials, "Signature"]]] === 1 &], 0, 1]];
