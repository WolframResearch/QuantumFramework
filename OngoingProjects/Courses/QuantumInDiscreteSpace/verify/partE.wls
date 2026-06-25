PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* E1 ensemble -> density matrix + purity *)
rho = 1/2 QuantumState["0"]["DensityMatrix"] + 1/2 QuantumState["+"]["DensityMatrix"];
mixed = QuantumState[rho];
e1rho = mixed["DensityMatrix"] // Normal;
e1purity = mixed["Purity"];
e1pure0 = QuantumState["0"]["Purity"];
e1maxmix = QuantumState[{{1/2, 0}, {0, 1/2}}]["Purity"];

(* E2 von Neumann entropy + bipartition workaround *)
e2pure = QuantumState["0"]["VonNeumannEntropy"];
e2maxmix = QuantumState[{{1/2, 0}, {0, 1/2}}]["VonNeumannEntropy"];
e2ent = QuantumPartialTrace[QuantumState["PhiPlus"], {2}]["VonNeumannEntropy"];

(* E3 Werner state threshold *)
e3werner = QuantumState["Werner"[p, 2]]["DensityMatrix"] // Normal // safe;
e3ent5 = QuantumEntangledQ[QuantumState["Werner"[0.5, 2]], {1}] // safe;
e3ent2 = QuantumEntangledQ[QuantumState["Werner"[0.2, 2]], {1}] // safe;
e3neg = QuantumEntanglementMonotone[QuantumState["Werner"[p, 2]], {1}, "Negativity"] // safe;

(* E4 fidelity + trace distance *)
e4fid = QuantumDistance[QuantumState["0"], QuantumState["+"], "Fidelity"] // safe;
e4td = QuantumDistance[QuantumState["0"], QuantumState["+"], "TraceDistance"] // safe;
e4bures = QuantumDistance[QuantumState["0"], QuantumState["+"], "BuresDistance"] // safe;

(* E5 purification *)
e5 = QuantumState[{{1/2, 0}, {0, 1/2}}]["Purify"] // safe;
e5vec = (QuantumState[{{1/2, 0}, {0, 1/2}}]["Purify"]["StateVector"] // Normal) // safe;
e5check = (QuantumPartialTrace[QuantumState[{{1/2, 0}, {0, 1/2}}]["Purify"], {2}]["DensityMatrix"] // Normal) // safe;
e5props = Select[QuantumState[{{1/2, 0}, {0, 1/2}}]["Properties"], StringContainsQ[ToLowerCase[#], "purif"] &] // safe;

res = <|
 "E1-rho" -> e1rho, "E1-purity-mixed" -> e1purity, "E1-purity-pure" -> e1pure0, "E1-purity-maxmix" -> e1maxmix,
 "E2-vne-pure" -> e2pure, "E2-vne-maxmix" -> e2maxmix, "E2-vne-entanglement" -> e2ent,
 "E3-werner-rho" -> e3werner, "E3-ent-p0.5" -> e3ent5, "E3-ent-p0.2" -> e3ent2, "E3-negativity(p)" -> e3neg,
 "E4-fidelity" -> e4fid, "E4-tracedist" -> e4td, "E4-bures" -> e4bures,
 "E5-purify-vec" -> e5vec, "E5-purify-reduced" -> e5check, "E5-props" -> e5props
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
