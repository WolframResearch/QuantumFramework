PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 52], " ", ToString[Short[x, 2], InputForm]]);
SeedRandom[1]; n = 10000;
gates = Flatten @ Table[{"H" -> RandomInteger[{1, n}], "CNOT" -> RandomSample[Range[n], 2], "S" -> RandomInteger[{1, n}]}, {6666}];
{tg, st} = AbsoluteTiming[PauliStabilizer[n]["ApplyCircuit", gates]];
P["stabilizer n=10000, 19998 gates: ApplyCircuit (s)", tg];
P["  per-gate cost (microseconds)", 10.^6 tg/19998];
{te, ent} = AbsoluteTiming[st["Entropy", Range[n/2]]];
P["  half-chain entropy (s), value", {te, ent}];
