(* derived-functionals-check.wls: the facts of section 10 (distances and entanglement),
   run against the working tree. Each line is the call followed by what the kernel returns. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
out[label_, expr_] := WriteString[$Output, label, ": ", ToString[expr, InputForm], "\n"];

(* QuantumDistance default measure is "Fidelity", which is the DISTANCE 1 - Re Tr sqrt(rho1 rho2) *)
rho1 = QuantumState[{{3/5, 0}, {0, 2/5}}];
rho2 = QuantumState[{{1/2, 1/5}, {1/5, 1/2}}];
hand = 1 - Re[Tr[MatrixPower[Normal[rho1["DensityMatrix"]] . Normal[rho2["DensityMatrix"]], 1/2]]];
out["E1 QuantumDistance default == explicit Fidelity == 1 - Re Tr sqrt(rho1 rho2)",
  {Simplify[QuantumDistance[rho1, rho2] - QuantumDistance[rho1, rho2, "Fidelity"]] === 0,
   Simplify[QuantumDistance[rho1, rho2] - hand] === 0}];
out["E1b QuantumDistance vanishes on identical states, and default != Uhlmann fidelity",
  {QuantumDistance[QuantumState["0"], QuantumState["0"]] === 0,
   QuantumDistance[QuantumState["0"], QuantumState["1"]] === 1}];
out["E1c QuantumSimilarity default is 1 - d for the Fidelity distance",
  Simplify[QuantumSimilarity[rho1, rho2] - (1 - QuantumDistance[rho1, rho2])] === 0];

(* QuantumEntanglement defaults: monotone Concurrence, QuantumEntangledQ Realignment *)
bell = QuantumState["PhiPlus"];
prod = QuantumState["00"];
out["E2 QuantumEntanglementMonotone default (Concurrence): Bell 1, product 0",
  {Simplify[QuantumEntanglementMonotone[bell] - 1] === 0, Simplify[QuantumEntanglementMonotone[prod]] === 0}];
out["E2b QuantumEntangledQ (default Realignment): Bell True, product False",
  {QuantumEntangledQ[bell], QuantumEntangledQ[prod]}];
