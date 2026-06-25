PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* J1 discrete Wigner of a qutrit stabilizer state (odd d -> d x d, nonnegative) *)
w0 = QuantumState["0", 3]["PhaseSpace"] // Normal // safe;
w0qp = QuantumState["0", 3]["QuasiProbability"] // Normal // safe;
j1total = Total[Flatten[QuantumState["0", 3]["PhaseSpace"] // Normal]] // safe;
j1min = Min[QuantumState["0", 3]["PhaseSpace"] // Normal] // safe;

(* J2 Wigner negativity of a qutrit magic ("strange") state *)
strange = QuantumState[{0, 1, -1}/Sqrt[2], 3];
wS = strange["PhaseSpace"] // Normal // safe;
j2min = Min[strange["PhaseSpace"] // Normal] // safe;
j2neg = Total[Select[Flatten[strange["PhaseSpace"] // Normal], # < 0 &]] // safe;

(* J3 mana = Log of the 1-norm of the Wigner function *)
mana[qs_] := Log[Total[Abs[Flatten[qs["PhaseSpace"] // Normal]]]];
j3strange = FullSimplify[mana[strange]] // safe;
j3strangeN = N[mana[strange]] // safe;
j3log53 = N[Log[5/3]];
j3stab = FullSimplify[mana[QuantumState["0", 3]]] // safe;     (* stabilizer -> 0 *)

res = <|
 "J1-wigner|0>3" -> w0, "J1-total" -> j1total, "J1-min" -> j1min,
 "J2-wigner-strange" -> wS, "J2-min" -> j2min, "J2-negsum" -> j2neg,
 "J3-mana-strange" -> j3strange, "J3-mana-N" -> j3strangeN, "J3-log(5/3)" -> j3log53, "J3-mana-stab" -> j3stab
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
