PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "..", "QuantumFramework"}]]; Needs["Wolfram`QuantumFramework`"]; Needs["Wolfram`QuantumFramework`SecondQuantization`"];
P[l_, x_] := (Print["@@ ", StringPadRight[l, 52], " ", ToString[Short[x, 2], InputForm]]);
(* qubit, H = Omega(t) X/2, a pi pulse of width w centred at t0 = 50 in [0, 100]. Captured: P1 = 1. Skipped: P1 = 0. *)
p1[om_, opts___] := Module[{psi = QuantumEvolve[(om/2) QuantumOperator["X"], QuantumState["0"], {\[FormalT], 0, 100}, opts]},
   Chop @ Abs[Normal[psi["StateVector"] /. \[FormalT] -> 100][[2]]]^2];
Do[Module[{sq, ga},
   sq = Piecewise[{{Pi/w, Abs[\[FormalT] - 50] < w/2}}, 0];
   ga = (Pi/(Sqrt[2 Pi] (w/4))) Exp[-(\[FormalT] - 50)^2/(2 (w/4)^2)];
   P["w=" <> ToString[w] <> " square (Piecewise) pulse: P1", p1[sq]];
   P["w=" <> ToString[w] <> " smooth Gaussian pulse:   P1", p1[ga]]],
  {w, {1., 0.1, 0.01}}];
