(* Validates prototype.wl against the live paclet + qiskit 2.4.1.
   Uses the chosen CamelCase WL option surface (algorithmically mapped to qiskit
   snake_case). Run: wolframscript -file validate.wl *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "prototype.wl"}]];

qc = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, {1}, {2}}]; qk = qc["Qiskit"];
cmap = {{1, 0}, {1, 3}, {2, 0}, {2, 3}};
results = {};
add[id_, cond_] := AppendTo[results, {If[TrueQ[cond], "PASS", "FAIL"], id}];

(* CamelCase keys: BasisGates / CouplingMap / OptimizationLevel -> snake_case *)
t = qTranspile[qk, <|"BasisGates" -> {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "OptimizationLevel" -> 3|>];
add["transpile-returns-QiskitCircuit", Head[t] === QiskitCircuit];
add["native-count-ops", Sort@Keys@t["Ops"] === Sort@{"cz", "measure", "rz", "sx"}];
nq = qDumpQASM[t, 3];
add["faithful-native-qasm3", StringContainsQ[nq, "cz $1, $0"] && ! StringContainsQ[nq, "u2"] && ! StringContainsQ[nq, "u3"]];
gq = qDumpQASM[qk, 3];
add["faithful-generic-qasm3", StringContainsQ[gq, "h q[0]"] && StringContainsQ[gq, "cx q[0], q[1]"]];

(* raw snake_case keys still pass through (resilience: unmapped/new params reachable) *)
t2 = qTranspile[qk, <|"basis_gates" -> {"x", "sx", "rz", "cz"}, "seed_transpiler" -> 7, "OptimizationLevel" -> 1|>];
add["snake-case-passthrough", Head[t2] === QiskitCircuit];

(* unknown option -> Failure naming the (normalized) bad key *)
f = qTranspile[qk, <|"BasisGates" -> {"x"}, "NoSuchParam" -> 5|>];
add["unknown-option-failure", Head[f] === Failure && MemberQ[f[[2]]["Unknown"], "no_such_param"]];

(* non-universal basis -> Failure, never fake QASM *)
f2 = qTranspile[qk, <|"BasisGates" -> {"x", "cz"}, "OptimizationLevel" -> 1|>];
add["non-universal-basis-failure", Head[f2] === Failure];

(* Target object via CamelCase spec *)
tgt = QiskitTarget[<|"BasisGates" -> {"x", "sx", "rz", "cz"}, "NumQubits" -> 4, "CouplingMap" -> cmap, "CustomNameMapping" -> Null|>];
add["target-builds", Head[tgt] === QiskitTargetObject];

(* end-to-end one-liner QF -> native QASM3 *)
ol = qDumpQASM[qTranspile[qc["Qiskit"], <|"BasisGates" -> {"x", "sx", "rz", "cz"}, "CouplingMap" -> cmap, "OptimizationLevel" -> 3|>], 3];
add["end-to-end-oneliner", StringStartsQ[ol, "OPENQASM 3.0"]];

Scan[Echo[StringRiffle[#, " "]] &, results];
Echo[StringForm["`` / `` passed", Count[results, {"PASS", _}], Length[results]]];
