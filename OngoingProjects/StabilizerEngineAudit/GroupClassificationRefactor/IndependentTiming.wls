PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
SeedRandom[183947];
auditDir = DirectoryName[$InputFileName];
currentCircuitDefinition = DownValues[stabilizerCircuitSpecs];
Begin["Wolfram`QuantumFramework`GateMatrix`PackagePrivate`"];
afterDefinition = DownValues[stabilizerOperatorGroupSpecs];
Get["/private/tmp/qf-group-classification-refactor/BeforeGroupDefinition.wl"];
beforeDefinition = DownValues[stabilizerOperatorGroupSpecs];
DownValues[stabilizerOperatorGroupSpecs] = afterDefinition;
End[];
With[{element = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerElementSpecs},
  prebulkDefinition = {HoldPattern[stabilizerCircuitSpecs[qco_QuantumCircuitOperator]] :>
    Catenate[element /@ qco["Operators"]]}
];
setMode[mode_] := (
  DownValues[stabilizerCircuitSpecs] = If[mode === "PreBulk", prebulkDefinition, currentCircuitDefinition];
  DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupSpecs] =
    If[mode === "Before", Wolfram`QuantumFramework`GateMatrix`PackagePrivate`beforeDefinition,
      Wolfram`QuantumFramework`GateMatrix`PackagePrivate`afterDefinition]
);
independentCircuit[n_] := QuantumCircuitOperator[<|"Elements" -> Map[
  Function[i, FromCircuitOperatorShorthand @ With[{a = Mod[19 i, n] + 1, b = Mod[19 i + 7, n] + 1},
    Switch[Mod[i, 11], 0, "H" -> a, 1, "S" -> a, 2, "CNOT" -> {a, b},
      3, "CZ" -> {b, a}, 4, "V" -> a, 5, "Y" -> a, 6, "Z" -> a,
      7, "SWAP" -> {b, a}, 8, SuperDagger["S"] -> a, 9, "X" -> a,
      10, SuperDagger["V"] -> a]]], Range[1000]], "Label" -> None|>];
PauliStabilizer[2]["ApplyCircuit", {"H" -> 1}];
rows = Catenate @ Map[Function[n, With[{qc = independentCircuit[n], initial = PauliStabilizer[n]},
  qc["FullInputOrder"];
  Map[Function[warm, Module[{trials, medians},
    trials = Catenate @ Map[Function[t, Map[Function[mode,
      setMode[mode]; stabilizerClearCaches[];
      If[warm, qc[initial]];
      ClearSystemCache[];
      Block[{$MessageList = {}}, With[{timed = AbsoluteTiming[Nest[qc, initial, 10]]},
        <|"Trial" -> t, "Mode" -> mode, "Seconds" -> First[timed],
          "Valid" -> (Head[Last[timed]] === PauliStabilizer),
          "Signature" -> ToString[Hash[{Last[timed]["Tableau"], Last[timed]["Signs"]}, "SHA256"], InputForm],
          "Messages" -> ToString[$MessageList, InputForm]|>
      ]]
    ], RandomSample[{"PreBulk", "Before", "After"}]]], Range[7]];
    medians = Map[Median, GroupBy[trials, (#Mode &) -> (#Seconds &)]];
    With[{row = <|"Qubits" -> n, "Warm" -> warm, "Trials" -> trials,
      "MediansSeconds" -> medians,
      "RangesSeconds" -> Map[MinMax, GroupBy[trials, (#Mode &) -> (#Seconds &)]],
      "MedianPairedAfterBefore" -> Median @ Map[Function[t,
        With[{pair = AssociationThread[Lookup[Select[trials, #Trial === t &], "Mode"],
            Lookup[Select[trials, #Trial === t &], "Seconds"]]}, pair["After"]/pair["Before"]]
      ], Range[7]],
      "BaselineSpeedup" -> medians["PreBulk"]/medians["After"]|>},
      Echo[KeyDrop[row, "Trials"]]; row]
  ]], {False, True}]
]], {128, 1024}];
setMode["After"];
Export[FileNameJoin[{auditDir, "independent-timing.json"}], <|"Version" -> $Version,
  "CompileTarget" -> $stabilizerCompileTarget, "Rows" -> rows,
  "OutsideMessages" -> ToString[$MessageList, InputForm]|>, "RawJSON"];
pass = AllTrue[rows, AllTrue[#Trials, #Valid && #Messages === "{}" &] &&
  Length[DeleteDuplicates[Lookup[#Trials, "Signature"]]] === 1 &];
Echo[<|"Pass" -> pass|>]; Exit[If[pass, 0, 1]];
