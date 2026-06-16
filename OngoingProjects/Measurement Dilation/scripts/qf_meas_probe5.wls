PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "RY"[0.5] -> 2}];
psi = qc[];
(* (1) single-outcome prob of |00> via circuit /* bra, then ["Number"] *)
amp00 = (qc /* {QuantumState["00"]["Dagger"]})[]["Number"];
(* (2) expectation <ZZ> via qs^dagger @ op @ qs  ["Scalar"] *)
zz = QuantumOperator["ZZ"];
exp1 = (psi["Dagger"][zz[psi]])["Scalar"];
(* (3) literal sandwich circuit form bra ... ket *)
exp2 = QuantumCircuitOperator[{psi, zz, psi["Dagger"]}][]["Scalar"];
(* references *)
pl = psi["ProbabilitiesList"];
trace = Tr[zz["MatrixRepresentation"] . psi["DensityMatrix"]];
Print["amp00^2 vs pl[[1]] : ", {Abs[amp00]^2, First[pl], Abs[Abs[amp00]^2 - First[pl]]} // InputForm];
Print["exp1 vs Tr : ", {Re[exp1], Re[trace], Abs[exp1 - trace]} // InputForm];
Print["exp2 (sandwich) vs Tr : ", {Re[exp2], Re[trace], Abs[exp2 - trace]} // InputForm];
