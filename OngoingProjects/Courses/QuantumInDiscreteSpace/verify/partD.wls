PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* D1 tensor products *)
d1state = QuantumTensorProduct[QuantumState["0"], QuantumState["1"]]["StateVector"] // Normal // safe;
d1op = QuantumTensorProduct[QuantumOperator["X"], QuantumOperator["Z"]]["MatrixRepresentation"] // Normal // safe;

(* D2 Bell / GHZ / W *)
d2bell = QuantumState["PhiPlus"]["StateVector"] // Normal;
d2ghz = QuantumState["GHZ"[3]]["StateVector"] // Normal;
d2w = QuantumState["W"[3]]["StateVector"] // Normal;

(* D3 partial trace -> reduced states *)
d3bell = QuantumPartialTrace[QuantumState["PhiPlus"], {2}]["DensityMatrix"] // Normal;
d3ghz = QuantumPartialTrace[QuantumState["GHZ"[3]], {2, 3}]["DensityMatrix"] // Normal;

(* D4 entanglement test *)
d4bell = QuantumEntangledQ[QuantumState["PhiPlus"]] // safe;
d4bellpart = QuantumEntangledQ[QuantumState["PhiPlus"], {1}] // safe;
d4prod = QuantumEntangledQ[QuantumTensorProduct[QuantumState["0"], QuantumState["1"]], {1}] // safe;

(* D5 quantify: entropy of entanglement, concurrence, negativity *)
d5ent = QuantumPartialTrace[QuantumState["PhiPlus"], {2}]["VonNeumannEntropy"] // safe;
d5conc = QuantumEntanglementMonotone[QuantumState["PhiPlus"], {1}, "Concurrence"] // safe;
d5neg = QuantumEntanglementMonotone[QuantumState["PhiPlus"], {1}, "Negativity"] // safe;
d5concProd = QuantumEntanglementMonotone[QuantumTensorProduct[QuantumState["0"], QuantumState["1"]], {1}, "Concurrence"] // safe;

(* D6 Schmidt decomposition *)
d6 = QuantumState["PhiPlus"]["SchmidtDecompose"] // safe;
d6coeffs = QuantumState["PhiPlus"]["SchmidtCoefficients"] // safe;
d6props = Select[QuantumState["PhiPlus"]["Properties"], StringContainsQ[ToLowerCase[#], "schmidt"] &] // safe;

(* D7 no-signalling: local op on qubit 1 leaves reduced state of qubit 2 unchanged *)
bellX1 = QuantumCircuitOperator[{"X" -> 1}][QuantumState["PhiPlus"]] // safe;
d7a = QuantumPartialTrace[QuantumState["PhiPlus"], {1}]["DensityMatrix"] // Normal // safe;
d7b = QuantumPartialTrace[bellX1, {1}]["DensityMatrix"] // Normal // safe;
d7equal = (d7a === d7b) // safe;

res = <|
 "D1-state-01" -> d1state, "D1-op-XZ" -> d1op,
 "D2-bell" -> d2bell, "D2-ghz3" -> d2ghz, "D2-w3" -> d2w,
 "D3-bell-reduced" -> d3bell, "D3-ghz-reduced" -> d3ghz,
 "D4-bell-ent" -> d4bell, "D4-bell-ent-part" -> d4bellpart, "D4-prod-ent" -> d4prod,
 "D5-entropy" -> d5ent, "D5-concurrence" -> d5conc, "D5-negativity" -> d5neg, "D5-conc-prod" -> d5concProd,
 "D6-schmidt" -> d6, "D6-coeffs" -> d6coeffs, "D6-props" -> d6props,
 "D7-red1" -> d7a, "D7-red1-afterX" -> d7b, "D7-equal" -> d7equal
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
