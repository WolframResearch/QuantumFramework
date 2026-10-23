PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
rvCircuit[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];

rvA = QuantumOperator["H", {2}, "Label" -> "T"];
rvB = QuantumOperator[{{1, 1}, {1, -1}}/Sqrt[2]];
rvDA = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[rvA];
rvDB = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData[rvB];
rvDispatch = Wolfram`QuantumFramework`GateMatrix`PackagePrivate`$stabilizerPrimitiveDispatch;
Echo[{rvDA === rvDB, rvDA, rvDB}, "DATA_COMPARE"];
Echo[{Replace[rvDA, rvDispatch], Replace[rvDB, rvDispatch], Normal[rvDispatch]}, "DISPATCH_DIRECT"];
Echo[{stabilizerCircuitSpecs[rvCircuit[{rvA}]], stabilizerCircuitSpecs[rvCircuit[{rvB}]]}, "SPECS_DIRECT"];
Exit[];
