PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

VerificationTest[
	With[{psi = QuantumCircuitOperator[{"RY"[\[Theta]] -> 1, "CNOT" -> {1, 2}}][]},
		FullSimplify[QuantumEntanglementMonotone[psi, "Concurrence"], 0 < \[Theta] < Pi]
	],
	Sin[\[Theta]],
	SameTest -> SameQ,
	TestID -> "RY-CNOT concurrence is Sin[theta]"
]

VerificationTest[
	QuantumEntangledQ[QuantumState["00"]],
	False,
	SameTest -> SameQ,
	TestID -> "Product state |00> is not entangled"
]
