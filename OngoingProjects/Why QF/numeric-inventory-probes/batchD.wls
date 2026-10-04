PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];
b[e_, t_:240] := TimeConstrained[MemoryConstrained[e, 6*^9, "OOM"], t, "TIMEOUT"];
P[l_, x_] := Print["@@ ", StringPadRight[l, 50], " ", ToString[Short[x, 2], InputForm]];
tm[e_, t_:240] := Module[{r = AbsoluteTiming[b[e, t]]}, If[StringQ[r[[2]]] , r[[2]], {Round[r[[1]], 0.001], r[[2]]}]];
SetAttributes[tm, HoldFirst];

(* stabilizer: the showcase's own claim, then 10x the qubits and gates *)
Do[SeedRandom[1];
  Module[{n = c[[1]], g = c[[2]], gates},
   gates = Flatten @ Table[{"H" -> RandomInteger[{1, n}], "CNOT" -> RandomSample[Range[n], 2], "S" -> RandomInteger[{1, n}]}, {Floor[g/3]}];
   P["stabilizer n=" <> ToString[n] <> ", " <> ToString[3 Floor[g/3]] <> " gates (s, entropy)",
     tm[Module[{st = PauliStabilizer[n]["ApplyCircuit", gates]}, st["Entropy", Range[n/2]]]]]],
  {c, {{1000, 20000}, {10000, 200000}}}];

Print[];
(* tensor network: the SAME depth-4 Trotter quench on the FULL chain, <Z_center>, greedy path *)
d = 4; aZZ = 0.4; aX = 0.6;
quench[m_] := QuantumCircuitOperator @ Join[{StringJoin @ ConstantArray["0", m]},
   Flatten @ Table[Join[("R"[aZZ, "ZZ"] -> # &) /@ Partition[Range[m], 2, 1], ("RX"[aX] -> # &) /@ Range[m]], {d}]];
meanZ[m_, c_] := quench[m] /* QuantumOperator["Z", {c}] /* quench[m]["Dagger"];
Do[P["TN full chain n=" <> ToString[n] <> " <Z_center> (s, value)",
     tm[Chop @ meanZ[n, Ceiling[n/2]][Method -> {"TensorNetwork", "Path" -> "Greedy"}]["Scalar"], 180]],
  {n, {9, 16, 24, 32}}];

Print[];
(* dense state vector *)
Do[P["Schrodinger dense quench n=" <> ToString[n] <> " (s)",
     First @ tm[quench[n][Method -> "Schrodinger"]["Dimension"], 240]], {n, {12, 16, 20, 22}}];

Print[];
P["QuEST backend available locally?", b[TrueQ[Wolfram`QuantumFramework`PackageScope`QuESTLinkQ[]] , 30]];
