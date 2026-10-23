PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
SeedRandom[20261022];
auditDir = DirectoryName[$InputFileName];
out = FileNameJoin[{auditDir, "paired.json"}];
currentSpecsDefinition = DownValues[stabilizerCircuitSpecs];
currentGroupDefinition = DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupSpecs];
Get[FileNameJoin[{auditDir, "BeforeGroupDefinition.wl"}]];
beforeGroupDefinition = DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupSpecs];
Echo[<|"DistinctGroupControls" -> (beforeGroupDefinition =!= currentGroupDefinition)|>];
If[beforeGroupDefinition === currentGroupDefinition, Exit[2]];
DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupSpecs] = currentGroupDefinition;
primitiveInitializer = OwnValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerPrimitiveDispatch];
With[{element = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerElementSpecs},
    legacySpecsDefinition = {HoldPattern[stabilizerCircuitSpecs[qco_QuantumCircuitOperator]] :>
        Catenate[element /@ qco["Operators"]]}
];
setMode[mode_] := (
    DownValues[stabilizerCircuitSpecs] = If[mode === "PreBulk", legacySpecsDefinition, currentSpecsDefinition];
    DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupSpecs] =
        If[mode === "BeforeHelper", beforeGroupDefinition, currentGroupDefinition]
);
makeSpecs[n_, m_] := Table[With[{j = Mod[37 i, n] + 1},
    Switch[Mod[i, 7], 0, "H" -> j, 1, "S" -> j,
        2, "CNOT" -> {j, Mod[j, n] + 1}, 3, "CZ" -> {j, Mod[j, n] + 1},
        4, "X" -> j, 5, "SWAP" -> {j, Mod[j, n] + 1}, 6, SuperDagger["S"] -> j]], {i, m}];
randomSpecs[n_, m_] := Table[With[{
    gate = RandomChoice[{"H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"}],
    wires = RandomSample[Range[n], 2]
}, gate -> If[MemberQ[{"CNOT", "CZ", "SWAP"}, gate], wires, First[wires]]], {m}];
makeCircuit[specs_] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ specs), "Label" -> None|>];
publicTen[] := Nest[trialCircuit[#] &, trialState, 10];
rows = {};
measure[label_, cold_, initializePrimitives_:False] := Module[{trials, row, modes, timings},
    trials = Catenate @ Table[
        modes = RotateLeft[{"PreBulk", "BeforeHelper", "AfterHelper"}, Mod[trial - 1, 3]];
        If[EvenQ[Quotient[trial - 1, 3]], modes = Reverse[modes]];
        Map[Function[mode,
            setMode[mode];
            stabilizerClearCaches[];
            If[!cold, publicTen[]];
            If[initializePrimitives && mode =!= "PreBulk",
                OwnValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerPrimitiveDispatch] = primitiveInitializer];
            ClearSystemCache[];
            Block[{$MessageList = {}}, With[{timing = AbsoluteTiming[publicTen[]]},
                <|"Mode" -> mode, "Trial" -> trial, "Seconds" -> First[timing],
                    "Signature" -> ToString[Hash[{Last[timing]["Tableau"], Last[timing]["Signs"]}, "SHA256"], InputForm],
                    "Head" -> ToString[Head[Last[timing]], InputForm], "Messages" -> ToString[$MessageList, InputForm]|>
            ]]
        ], modes], {trial, 9}];
    timings = GroupBy[trials, (Lookup[#, "Mode"] &) -> (Lookup[#, "Seconds"] &)];
    row = <|"Label" -> label, "Trials" -> trials,
        "MediansSeconds" -> Map[Median, timings],
        "RangesSeconds" -> Map[MinMax, timings],
        "MedianPairedAfterBefore" -> Median[Map[Function[t,
            With[{pair = AssociationThread[Lookup[Select[trials, #Trial === t &], "Mode"],
                Lookup[Select[trials, #Trial === t &], "Seconds"]]},
                pair["AfterHelper"] / pair["BeforeHelper"]]], Range[9]]]|>;
    rows = Append[rows, row]; Echo[KeyDrop[row, "Trials"]];
    Export[out, <|"Version" -> $Version, "CompileTarget" -> $stabilizerCompileTarget,
        "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"], "Rows" -> rows,
        "OutsideMessages" -> ToString[$MessageList, InputForm]|>, "RawJSON"]
];
PauliStabilizer[2]["ApplyCircuit", {"H" -> 1}];
stabilizerCircuitSpecs[QuantumCircuitOperator[{"H" -> 1}]];
Scan[Function[n,
    trialState = PauliStabilizer[n];
    trialCircuit = makeCircuit[makeSpecs[n, 1000]];
    trialCircuit["FullInputOrder"];
    measure["MixedCold/n=" <> ToString[n] <> "/g=10000", True];
    measure["MixedFirstPlan/n=" <> ToString[n] <> "/g=10000", True, True];
    measure["MixedWarm/n=" <> ToString[n] <> "/g=10000", False];
    trialCircuit = makeCircuit[randomSpecs[n, 1000]];
    trialCircuit["FullInputOrder"];
    measure["RandomCold/n=" <> ToString[n] <> "/g=10000", True];
    measure["RandomWarm/n=" <> ToString[n] <> "/g=10000", False];
], {128, 1024}];
setMode["AfterHelper"];
Echo[<|"Finished" -> True, "OutsideMessages" -> ToString[$MessageList, InputForm]|>];
Exit[If[AllTrue[Catenate[Lookup[rows, "Trials"]], #Head === "PauliStabilizer" && #Messages === "{}" &] &&
    AllTrue[rows, Length[DeleteDuplicates[Lookup[#Trials, "Signature"]]] === 1 &], 0, 1]];
