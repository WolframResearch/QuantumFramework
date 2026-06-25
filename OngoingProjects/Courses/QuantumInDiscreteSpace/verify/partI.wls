PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
SetAttributes[rec, HoldRest];
rec[id_, e_] := Module[{v = e}, WriteString[$Output, id <> " = " <> ToString[InputForm[v]] <> "\n"]; v];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* I1 stabilizer state + tableau *)
ps0 = PauliStabilizer[2];
rec["I1-tableau", ps0["Tableau"] // Normal // safe];
rec["I1-gens", ps0["Stabilizers"] // safe];

(* I2 apply Clifford gates -> Bell; read generators *)
psBell = PauliStabilizer[2][{"H", "CNOT"}] // safe;
rec["I2-bell-gens", psBell["Stabilizers"] // safe];
rec["I2-bell-state", psBell["QuantumState"]["StateVector"] // Normal // safe];

(* I3 large Clifford circuit, tableau-efficient *)
n = 60;
circ = Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}]];
psBig = PauliStabilizer[n][circ] // safe;
rec["I3-qubits", psBig["Qubits"] // safe];
rec["I3-stabQ", StabilizerStateQ[psBig] // safe];

(* I4 graph / cluster state on path 1-2-3 (|+++> then CZ on edges) *)
gcirc = {"H" -> 1, "H" -> 2, "H" -> 3, "CZ" -> {1, 2}, "CZ" -> {2, 3}};
rec["I4-graph-gens", PauliStabilizer[3][gcirc]["Stabilizers"] // safe];

(* I5 three-qubit bit-flip code: encode, inject X on qubit 2, measure syndrome *)
enc = QuantumCircuitOperator[{"CNOT" -> {1, 2}, "CNOT" -> {1, 3}}];
rec["I5-enc|000>", enc[QuantumState["000"]]["StateVector"] // Normal // safe];
logical1 = enc[QuantumState["100"]];
rec["I5-enc|100>", logical1["StateVector"] // Normal // safe];
errored = QuantumCircuitOperator[{"X" -> 2}][logical1];
rec["I5-synZ12", QuantumMeasurementOperator[QuantumOperator["ZZI"]][errored]["Mean"] // safe];
rec["I5-synZ23", QuantumMeasurementOperator[QuantumOperator["IZZ"]][errored]["Mean"] // safe];

(* I6 five-qubit code stabilizers commute *)
s1 = QuantumOperator["XZZXI"]; s2 = QuantumOperator["IXZZX"];
s3 = QuantumOperator["XIXZZ"]; s4 = QuantumOperator["ZXIXZ"];
commute[a_, b_] := (a @ b)["MatrixRepresentation"] == (b @ a)["MatrixRepresentation"];
rec["I6-5qubit-commute", And @@ (commute @@@ Subsets[{s1, s2, s3, s4}, {2}]) // safe];

(* I7 magic / T gate breaks the tableau: |+> is stabilizer, T|+> is not *)
rec["I7-+isStab", MatchQ[PauliStabilizer[QuantumState["+"]], _PauliStabilizer]];
rec["I7-T+isStab", MatchQ[PauliStabilizer[QuantumOperator["T"][QuantumState["+"]]], _PauliStabilizer]];
rec["I7-frame-components", Length[(PauliStabilizer[1]["T", 1])["Components"]]];
