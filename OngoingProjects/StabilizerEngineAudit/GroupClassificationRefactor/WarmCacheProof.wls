PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];
$ProgressReporting = False;
rows = Map[Function[n, Module[{qc, initial, calls = 0, helperCalls = 0, result},
  qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, n}, "S" -> n}];
  initial = PauliStabilizer[n]; stabilizerClearCaches[];
  qc[initial];
  result = Block[{
    stabilizerCircuitSpecs = Function[Null, calls++; $Failed],
    Wolfram`QuantumFramework`GateMatrix`PackagePrivate`stabilizerOperatorGroupClassification =
      Function[Null, helperCalls++; $Failed]
  }, Nest[qc, initial, 10]];
  <|"Qubits" -> n, "CircuitSpecsCalls" -> calls, "HelperCalls" -> helperCalls,
    "Valid" -> (Head[result] === PauliStabilizer)|>
]], {128, 1024}];
Echo[<|"Rows" -> rows, "Messages" -> $MessageList|>];
Exit[If[AllTrue[rows, #CircuitSpecsCalls === 0 && #HelperCalls === 0 && #Valid &] && $MessageList === {}, 0, 1]];
