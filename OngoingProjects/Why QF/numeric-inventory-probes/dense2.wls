PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 56], " ", ToString[x, InputForm]]);
(* The 4-step Trotter quench of dense.wls, applied with QF's default tensor-network method and with
   "Schrodinger", one bounded run each: seconds, or "timeout" / "memory". *)
d = 4; aZZ = 0.4; aX = 0.6;
quench[m_] := QuantumCircuitOperator @ Join[{StringJoin @ ConstantArray["0", m]},
   Flatten @ Table[Join[("R"[aZZ, "ZZ"] -> # &) /@ Partition[Range[m], 2, 1], ("RX"[aX] -> # &) /@ Range[m]], {d}]];
SetAttributes[sec, HoldFirst];
sec[e_, limit_] := MemoryConstrained[TimeConstrained[With[{t = AbsoluteTiming[e]}, {N[Round[First[t], 1 / 100]], Last[t]}], limit, {"timeout", None}], 4 10^9, {"memory", None}];
zc[psi_, n_] := If[MatchQ[psi, _QuantumState], Chop @ (psi["Dagger"] @ QuantumOperator["Z", {Ceiling[n/2]}] @ psi)["Scalar"], psi];
P["TensorNetworks", First[PacletFind["Wolfram/TensorNetworks"]]["Version"] <> " from " <> First[PacletFind["Wolfram/TensorNetworks"]]["Location"]];
Do[Module[{qc = quench[n], t, s},
   t = sec[qc[], 240];
   P["default method n=" <> ToString[n] <> " (s), <Z_center>", {t[[1]], zc[t[[2]], n]}];
   If[n <= 16, s = sec[qc[Method -> "Schrodinger"], 240]; P["Schrodinger n=" <> ToString[n] <> " (s), <Z_center>", {s[[1]], zc[s[[2]], n]}]]],
  {n, {12, 16, 20}}];
