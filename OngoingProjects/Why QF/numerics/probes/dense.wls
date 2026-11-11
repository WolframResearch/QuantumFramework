PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 52], " ", ToString[x, InputForm]]);
d = 4; aZZ = 0.4; aX = 0.6;
quench[m_] := QuantumCircuitOperator @ Join[{StringJoin @ ConstantArray["0", m]},
   Flatten @ Table[Join[("R"[aZZ, "ZZ"] -> # &) /@ Partition[Range[m], 2, 1], ("RX"[aX] -> # &) /@ Range[m]], {d}]];
Do[Module[{t, psi}, {t, psi} = AbsoluteTiming[quench[n][Method -> "Schrodinger"]];
   P["dense Schrodinger quench n=" <> ToString[n] <> " (s)", t];
   P["  <Z_center> from the state vector", Chop @ (psi["Dagger"] @ QuantumOperator["Z", {Ceiling[n/2]}] @ psi)["Scalar"]]],
  {n, {12, 16, 20, 22}}];
P["QuEST backend available?", TrueQ[Quiet[Wolfram`QuantumFramework`PackageScope`QuESTLinkQ[]]]];
