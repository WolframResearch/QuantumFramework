PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[expr_] := Module[{t, best = Infinity}, Do[ClearSystemCache[]; t = First@AbsoluteTiming[expr]; best = Min[best, t], {3}]; best];

runOps[ps0_, ops_] := Fold[
   Function[{ps, op},
     Switch[op[[1]],
       "H", ps["H", op[[2]] + 1],
       "S", ps["S", op[[2]] + 1],
       "CNOT", ps["CNOT", op[[2]] + 1, op[[3]] + 1]]],
   ps0, ops];

Do[
  Module[{ops = Import["/tmp/stab_ops_" <> ToString[n] <> ".json"], ms},
   ms = 1000*tm[runOps[PauliStabilizer[n], ops]["Signs"]];
   pr["stab_n" <> ToString[n] <> "_m" <> ToString[Length[ops]] <> "_ms", ms];
   pr["stab_n" <> ToString[n] <> "_ms_per_gate", ms/Length[ops]];
  ], {n, {100, 500, 1000}}];
Print["DONE"];
