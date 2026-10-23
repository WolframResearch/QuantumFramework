PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
SeedRandom[183947];
auditDir = DirectoryName[$InputFileName];
Echo[<|"Version" -> $Version, "LoadedFile" -> FindFile["Wolfram`QuantumFramework`"],
  "QFVersion" -> PacletFind["Wolfram/QuantumFramework"][[1]]["Version"]|>];
Begin["Wolfram`QuantumFramework`GateMatrix`PackagePrivate`"];
afterDefinition = DownValues[stabilizerOperatorGroupSpecs];
Get["/private/tmp/qf-group-classification-refactor/BeforeGroupDefinition.wl"];
beforeDefinition = DownValues[stabilizerOperatorGroupSpecs];
DownValues[stabilizerOperatorGroupSpecs] = afterDefinition;
End[];
beginningMessages = $MessageList;
bulkReport = TestReport["/Users/mohammadb/Documents/GitHub/QuantumFramework/Tests/Stabilizer/BulkPlanBuilder.wlt"];
Echo[<|"Suite" -> "BulkPlanBuilder", "Count" -> Length[bulkReport["Results"]],
  "Outcomes" -> Counts[#["Outcome"] & /@ bulkReport["Results"]],
  "Failures" -> Map[{#["TestID"], #["Outcome"], #["ActualOutput"], #["ActualMessages"]} &,
    Select[bulkReport["Results"], #["Outcome"] =!= "Success" &]]|>];
Begin["Wolfram`QuantumFramework`GateMatrix`PackagePrivate`"];
groupRun[definition_, operators_List] := (
  DownValues[stabilizerOperatorGroupSpecs] = definition;
  stabilizerClearCaches[];
  stabilizerOperatorGroupSpecs[operators, stabilizerOperatorData[First[operators]], Range[Length[operators]]]
);
groupComparisons = With[{groups = {
  {QuantumOperator["I", {9}], QuantumOperator["I", {3}]},
  {QuantumOperator["H", {9}, "Label" -> "T"], QuantumOperator["H", {3}, "Label" -> None]},
  {QuantumOperator["CNOT", {9, 2}], QuantumOperator["CNOT", {7, 1}]},
  {QuantumOperator["H", {9, 2, 7}], QuantumOperator["H", {3, 1, 8}]},
  {QuantumOperator[I {{1, 1}, {1, -1}}/Sqrt[2], {9}], QuantumOperator[I {{1, 1}, {1, -1}}/Sqrt[2], {3}]},
  {QuantumOperator["P"[ArcCos[reviewTheta]], {9}], QuantumOperator["P"[ArcCos[reviewTheta]], {3}, "Label" -> "H"]},
  {QuantumOperator[QuantumOperator["T"] @ QuantumOperator["H"], {9}], QuantumOperator[QuantumOperator["T"] @ QuantumOperator["H"], {3}, "Label" -> "H"]},
  {QuantumOperator[PauliMatrix[1], {{1}, {1}}], QuantumOperator[PauliMatrix[1], {{2}, {1}}]}
}}, MapIndexed[Function[{operators, index}, With[{
    before = groupRun[beforeDefinition, operators], after = groupRun[afterDefinition, operators]
  }, <|"ID" -> First[index], "Equal" -> SameQ[before, after], "Result" -> ToString[after, InputForm]|>]], groups]];
DownValues[stabilizerOperatorGroupSpecs] = afterDefinition;
Echo[groupComparisons, "Group equivalence: "];
lazyPrimitive = Module[{keyCalls = 0, classificationCalls = 0, op = QuantumOperator["H", {9}], result},
  stabilizerClearCaches[];
  result = Block[{
    stabilizerGateKey = Function[Null, keyCalls++; "poison"],
    stabilizerOperatorClassification = Function[Null, classificationCalls++; Missing["poison"]]
  }, stabilizerOperatorGroupClassification[op, stabilizerOperatorData[op], {{{9}, {9}}}]];
  <|"Pass" -> (keyCalls === 0 && classificationCalls === 0 && result === {"H" -> {1}}),
    "KeyCalls" -> keyCalls, "ClassificationCalls" -> classificationCalls, "Result" -> ToString[result, InputForm]|>
];
lazyRejected = Module[{classificationCalls = 0, dispatchCalls = 0, op = QuantumOperator["H", {9}], result},
  stabilizerClearCaches[];
  result = Block[{$stabilizerPrimitiveDispatch,
    stabilizerOperatorClassification = Function[Null, classificationCalls++; Missing["poison"]]},
    $stabilizerPrimitiveDispatch := (dispatchCalls++; Dispatch[{_ :> $Failed}]);
    stabilizerOperatorGroupClassification[op, stabilizerOperatorData[op], {{{9}, {9}}, {{3}, {2}}}]
  ];
  <|"Pass" -> (classificationCalls === 0 && dispatchCalls === 0 && result === $Failed),
    "DispatchCalls" -> dispatchCalls, "ClassificationCalls" -> classificationCalls, "Result" -> ToString[result, InputForm]|>
];
generalCached = Module[{calls = 0, op = QuantumOperator[I {{1, 1}, {1, -1}}/Sqrt[2], {9}], result},
  stabilizerClearCaches[];
  result = Block[{stabilizerOperatorClassification = Function[Null, calls++; {"H" -> {1}}]},
    {stabilizerOperatorGroupClassification[op, stabilizerOperatorData[op], {{{9}, {9}}}],
      stabilizerOperatorGroupClassification[QuantumOperator[op, {3}], stabilizerOperatorData[op], {{{3}, {3}}}]}];
  <|"Pass" -> (calls === 1 && result === {{"H" -> {1}}, {"H" -> {1}}}), "ClassificationCalls" -> calls|>
];
remapping = With[{
    wires = {{2, 7, 9}, {1, 3, 8}, {5, 11, 19}},
    specs = {"H" -> {3}, "CNOT" -> {3, 1}, "S" -> {2}, "SWAP" -> {1, 2}}
  }, {stabilizerRemapSpecs[specs, wires] ===
        Map[Function[w, Replace[specs, (g_ -> local_List) :> (g -> w[[local]]), {1}]], wires],
      stabilizerRemapSpecs[{}, wires] === {{}, {}, {}}, stabilizerRemapSpecs[{}, {}] === {}}];
Echo[<|"PrimitiveLazy" -> lazyPrimitive, "InvalidOccurrenceLazy" -> lazyRejected,
  "GeneralCached" -> generalCached, "Remapping" -> remapping|>];
End[];
extraPass = AllTrue[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`groupComparisons, TrueQ[#Equal] &] &&
  And @@ Lookup[{
      Wolfram`QuantumFramework`GateMatrix`PackagePrivate`lazyPrimitive,
      Wolfram`QuantumFramework`GateMatrix`PackagePrivate`lazyRejected,
      Wolfram`QuantumFramework`GateMatrix`PackagePrivate`generalCached}, "Pass"] &&
  And @@ Wolfram`QuantumFramework`GateMatrix`PackagePrivate`remapping;
Export[FileNameJoin[{auditDir, "correctness.json"}], <|"Version" -> $Version,
  "SuiteCount" -> Length[bulkReport["Results"]], "SuiteOutcomes" -> Counts[#["Outcome"] & /@ bulkReport["Results"]],
  "Groups" -> Wolfram`QuantumFramework`GateMatrix`PackagePrivate`groupComparisons,
  "ExtraPass" -> extraPass, "Messages" -> ToString[$MessageList, InputForm]|>, "RawJSON"];
Echo[<|"ExtraPass" -> extraPass, "Messages" -> $MessageList|>];
Exit[If[extraPass && AllTrue[bulkReport["Results"], #["Outcome"] === "Success" &], 0, 1]];
