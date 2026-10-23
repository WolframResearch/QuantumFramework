PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
rvCircuit[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];
rvGates = {"I", "H", "S", SuperDagger["S"], "X", "Y", "Z", "V", SuperDagger["V"], "CNOT", "CZ", "SWAP"};
rvOps = MapIndexed[QuantumOperator[#, If[MemberQ[{"CNOT", "CZ", "SWAP"}, #], {1, 3}, {2}], "Label" -> "T"] &, rvGates];
rvCalls = 0;
rvPlan = Internal`InheritedBlock[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification},
    Clear[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification];
    Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification[op_QuantumOperator] := (rvCalls++; Missing["Intercepted"]);
    stabilizerClearCaches[];
    stabilizerCircuitSpecs[rvCircuit[rvOps]]
];
Echo[{rvCalls, FreeQ[rvPlan, _Missing], rvPlan}, "EXACT_PRIMITIVES_WITH_GENERAL_CLASSIFIER_DISABLED"];
rvPhaseOps = {QuantumOperator[Exp[I Pi/9] {{1, 1}, {1, -1}}/Sqrt[2], {1}], QuantumOperator[N[{{1, 1}, {1, -1}}/Sqrt[2]], {3}]};
rvCalls = 0;
rvPhasePlan = Internal`InheritedBlock[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification},
    Clear[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification];
    Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification[op_QuantumOperator] := (rvCalls++; {"H" -> {1}});
    stabilizerClearCaches[];
    stabilizerCircuitSpecs[rvCircuit[rvPhaseOps]]
];
Echo[{rvCalls, rvPhasePlan}, "PHASE_AND_APPROXIMATE_USE_GENERAL_CLASSIFIER"];
rvCustom = {{1, I}, {I, 1}}/Sqrt[2];
rvCustomOps = QuantumOperator[rvCustom, {#}, "Label" -> ToString[#]] & /@ Range[1, 60];
rvCalls = 0;
rvCustomPlan = Internal`InheritedBlock[{Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification},
    Clear[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification];
    Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorClassification[op_QuantumOperator] := (rvCalls++; {"V" -> {1}});
    stabilizerClearCaches[];
    stabilizerCircuitSpecs[rvCircuit[rvCustomOps]]
];
Echo[{rvCalls, Length[rvCustomPlan], Last /@ rvCustomPlan === List /@ Range[60]}, "ONE_CUSTOM_CLASSIFICATION_PER_PHYSICAL_GROUP"];
rvMessages = TestReport[{
    VerificationTest[QuantumOperator["H", {1, 1}], $Failed, {QuantumOperator::invalidOutputOrder}, TestID -> "MalformedDuplicateWire"],
    VerificationTest[QuantumOperator["NoSuchBulkGate"], $Failed, {QuantumOperator::invalidName}, TestID -> "MalformedGateName"]
}];
Echo[{rvMessages["TestsSucceededCount"], rvMessages["TestsFailedCount"], ({#["TestID"], #["Outcome"], #["ActualOutput"], #["ActualMessages"]} & /@ Values[rvMessages["TestResults"]])}, "MALFORMED_CONSTRUCTORS"];
Exit[];
