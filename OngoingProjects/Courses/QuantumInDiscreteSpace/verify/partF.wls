PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* F1 two-qubit and three-qubit gates *)
f1cnot = QuantumOperator["CNOT"]["MatrixRepresentation"] // Normal;
f1cnot10 = QuantumOperator["CNOT"][QuantumState["10"]]["StateVector"] // Normal;     (* -> 11 *)
f1cz = QuantumOperator["CZ"]["MatrixRepresentation"] // Normal;
f1swap10 = QuantumOperator["SWAP"][QuantumState["10"]]["StateVector"] // Normal;       (* -> 01 *)
f1tof = QuantumOperator["Toffoli"][QuantumState["110"]]["StateVector"] // Normal // safe;  (* -> 111 *)

(* F2 compose a circuit -> Bell *)
bell = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}];
f2 = bell[QuantumState["00"]]["StateVector"] // Normal;
f2elements = bell["Operators"] // Length // safe;

(* F3 depth / layers *)
f3depth = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 1, "X" -> 2}]["Depth"] // safe;
f3layers = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}, "H" -> 1, "X" -> 2}]["Layers"] // Length // safe;

(* F4 overall unitary *)
f4 = bell["QuantumOperator"]["MatrixRepresentation"] // Normal // safe;

(* F5 KAK decomposition of CNOT: circuit of local U gates + controlled-X entanglers *)
f5kak = QuantumOperator["CNOT"]["KAK"];
f5head = Head[f5kak];
f5num = Length[f5kak["Operators"]];
f5tally = Tally[Replace[#["Label"],
      {"U"[__] -> "U", Subscript["C", "X"][__] -> "CX", "GlobalPhase"[_] -> "GlobalPhase"}] & /@
    f5kak["Operators"]];

(* F6 universality: decompose a random 2-qubit unitary into CNOT+locals, recompose, check *)
SeedRandom[7];
u = QuantumOperator["RandomUnitary"[{2, 2}]];
dec = u["Decompose"] // safe;
f6faithful = Chop[Norm[(u["MatrixRepresentation"] // Normal) -
     (dec["MatrixRepresentation"] // Normal), "Frobenius"]] // safe;
f6dist = QuantumDistance[u["QuantumOperator"], dec["QuantumOperator"]] // safe;

res = <|
 "F1-cnot" -> f1cnot, "F1-cnot|10>" -> f1cnot10, "F1-cz" -> f1cz, "F1-swap|10>" -> f1swap10, "F1-toffoli|110>" -> f1tof,
 "F2-bell-from-circuit" -> f2, "F2-numops" -> f2elements,
 "F3-depth" -> f3depth, "F3-numlayers" -> f3layers,
 "F4-unitary" -> f4,
 "F5-kak-head" -> f5head, "F5-kak-numgates" -> f5num, "F5-kak-tally" -> f5tally,
 "F6-recompose-frobnorm" -> f6faithful, "F6-distance" -> f6dist
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
