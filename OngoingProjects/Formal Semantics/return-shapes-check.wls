PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

emit[label_, val_] := WriteString[$Output, label, " -> ", ToString[val, InputForm], "\n"];

(* internal wrapped forms *)
qoX = QuantumOperator["X"];
emit["QuantumOperator wraps head", Head[qoX[[1]]]];                 (* QuantumState *)
emit["QuantumOperator order {out,in}", {qoX["OutputOrder"], qoX["InputOrder"]}];
chan = QuantumChannel["BitFlip"[1/10]];
emit["Channel outQ>=inQ", {chan["QuantumOperator"]["OutputQudits"], chan["QuantumOperator"]["InputQudits"]}];
qmo = QuantumMeasurementOperator[{1}];
emit["QMO wraps head", Head[qmo[[1]]]];

(* interaction return heads *)
emit["qo[qs]", Head[qoX[QuantumState["0"]]]];
emit["qo1[qo2]", Head[QuantumOperator["Z"][QuantumOperator["X"]]]];
emit["qo[qmo]", Head[qoX[qmo]]];
emit["qo[qc(hannel)]", Head[qoX[chan]]];
emit["qmo[qs]", Head[qmo[QuantumState["0"]]]];
emit["qmo[qo]", Head[qmo[qoX]]];
emit["chan[qs]", Head[chan[QuantumState["0"]]]];
emit["qco[qs]", Head[QuantumCircuitOperator["Bell"][QuantumState["00"]]]];
emit["qco[]", Head[QuantumCircuitOperator["Bell"][]]];
emit["qco[op] (extend)", Head[QuantumCircuitOperator["Bell"][qoX]]];
emit["qc1 /* qc2", Head[QuantumCircuitOperator["Bell"] /* QuantumCircuitOperator[{"X" -> 1}]]];
emit["PartialTrace[qs]", Head[QuantumPartialTrace[QuantumState["Bell"], {1}]]];
emit["PartialTrace[qo]", Head[QuantumPartialTrace[QuantumOperator["CNOT"], {{1, 1}}]]];
emit["PartialTrace[qc]", Head[QuantumPartialTrace[QuantumCircuitOperator[{"CNOT"}], {{1, 1}}]]];
emit["SuperDagger[qo]", Head[SuperDagger[qoX]]];
emit["qo1 + qo2", Head[QuantumOperator["X"] + QuantumOperator["Z"]]];
emit["TensorProduct[qo,qo]", Head[QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]]];
emit["QuantumEvolve[H,psi0,t]", Head[QuantumEvolve[QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 1}]]];

(* measurement bundle: classical distribution + post-states *)
qmZ = QuantumMeasurementOperator[QuantumOperator["Z"], {1}][QuantumState["+"]];
emit["qm Head", Head[qmZ]];
emit["qm Distribution head", Head[qmZ["Distribution"]]];
emit["qm PostMeasurementState head", Head[qmZ["PostMeasurementState"]]];

WriteString[$Output, "DONE\n"];
