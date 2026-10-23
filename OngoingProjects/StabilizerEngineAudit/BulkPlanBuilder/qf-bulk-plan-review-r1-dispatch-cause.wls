PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
rvCircuit[ops_List] := QuantumCircuitOperator[<|"Elements" -> (FromCircuitOperatorShorthand /@ ops), "Label" -> None|>];


Echo[Block[{$ContextPath = {}}, ToString[DownValues[Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorData], InputForm]], "DEFINITION"];
rvA = QuantumOperator["H", {2}, "Label" -> "T"]; rvB = QuantumOperator[{{1, 1}, {1, -1}}/Sqrt[2]];
Echo[{MatchQ[rvA, QuantumOperator[QuantumState[a_, b_QuantumBasis], o_]], MatchQ[rvA[[1]], QuantumState[a_, b_QuantumBasis]], MatchQ[rvA[[1, 2]], b_QuantumBasis]}, "MATCHING"];
rvFixed[op_QuantumOperator] := {op[[1, 1]], Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerBasisData[op[[1, 2]]], Ordering /@ op[[2]]};
Echo[{rvFixed[rvA] === rvFixed[rvB], Head[rvFixed[rvA]], rvFixed[rvA][[2]]}, "PART_FIX"];
Exit[];
