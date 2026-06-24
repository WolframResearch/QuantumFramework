(* Path-sum simulator + [HH] elimination, self-contained, verified vs QF.      *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
w = Exp[I Pi/4];

idem[p_] := Expand[p] /. Power[x[i_], _] :> x[i];
f2[p_] := PolynomialMod[idem[p], 2];
boolToInt[p_] := idem @ Fold[#1 + #2 - 2 #1 #2 &, 0, MonomialList[f2[p]]];

initPS[n_] := <|"val" -> Association[Table[q -> 0, {q, n}]], "phase" -> 0, "vars" -> {}, "h" -> 0, "scale" -> 1|>;
applyGate[ps_, g_] := Module[{val = ps["val"], phase = ps["phase"], vars = ps["vars"], h = ps["h"], v},
  Switch[g,
    {"H", _}, v = x[Length[vars] + 1]; AppendTo[vars, v]; h++; phase += 4 boolToInt[val[g[[2]]]] v; val[g[[2]]] = v,
    {"T", _}, phase += boolToInt[val[g[[2]]]],
    {"Z", _}, phase += 4 boolToInt[val[g[[2]]]],
    {"CZ", _, _}, phase += 4 boolToInt[val[g[[2]]]] boolToInt[val[g[[3]]]],
    {"CCZ", _, _, _}, phase += 4 boolToInt[val[g[[2]]]] boolToInt[val[g[[3]]]] boolToInt[val[g[[4]]]],
    {"CNOT", _, _}, val[g[[3]]] = f2[val[g[[2]]] + val[g[[3]]]]];
  <|ps, "val" -> val, "phase" -> idem[phase], "vars" -> vars, "h" -> h|>];
buildPS[n_, gates_] := Fold[applyGate, initPS[n], gates];

cofactor[phase_, v_] := PolynomialMod[idem[(phase /. v -> 1) - (phase /. v -> 0)], 8];
elimQ[phase_, v_] := Module[{A = cofactor[phase, v], Q},
  If[PolynomialMod[A, 4] =!= 0, Return[False]];
  Q = f2[A/4];
  Length[Variables[Q]] >= 1 && AllTrue[MonomialList[Q], Length[Variables[#]] <= 1 &]];
internalQ[ps_, v_] := AllTrue[Values[ps["val"]], FreeQ[#, v] &];
elimStep[ps_] := Module[{phase = ps["phase"], v, A, Q, u, sub, val},
  v = SelectFirst[ps["vars"], elimQ[phase, #] && internalQ[ps, #] &, Missing[]];
  If[MissingQ[v], Return[ps]];
  A = cofactor[phase, v]; Q = f2[A/4]; u = First[Variables[Q]]; sub = u -> f2[u + Q];
  phase = idem[(phase /. v -> 0) /. sub];
  val = f2[#] & /@ (ps["val"] /. sub);
  <|ps, "phase" -> phase, "val" -> val, "vars" -> DeleteCases[ps["vars"], v | u], "scale" -> ps["scale"]*2|>];
reduce[ps_] := FixedPoint[elimStep, ps];

amp[ps_, n_, y_] := Module[{vars = ps["vars"]},
  (ps["scale"]/Sqrt[2]^ps["h"]) Total[
    (w^Mod[ps["phase"] /. Thread[vars -> #], 8]) Boole[
       And @@ Table[Mod[ps["val"][q] /. Thread[vars -> #], 2] == y[[q]], {q, n}]] & /@ Tuples[{0, 1}, Length[vars]]]];
qfAmp[n_, gates_, y_] := Module[{qco},
  qco = QuantumCircuitOperator[gates /. {{"H", q_} :> "H" -> q, {"T", q_} :> "T" -> q, {"Z", q_} :> "Z" -> q,
     {"CZ", a_, b_} :> "CZ" -> {a, b}, {"CCZ", a_, b_, c_} :> QuantumOperator[DiagonalMatrix[{1,1,1,1,1,1,1,-1}], {a, b, c}],
     {"CNOT", a_, b_} :> "CNOT" -> {a, b}}];
  (Normal[qco[QuantumState[StringRepeat["0", n]], Method -> "Schrodinger"]["StateVector"]])[[1 + FromDigits[y, 2]]]];

matches[n_, gates_] := AllTrue[Tuples[{0, 1}, n], Chop[N[amp[reduce[buildPS[n, gates]], n, #] - qfAmp[n, gates, #]]] == 0 &];

Echo["==== builder correctness (no elimination) vs QF ====", ""];
Echo[AllTrue[Tuples[{0,1}, 1], Chop[N[amp[buildPS[1,{{"H",1},{"T",1},{"H",1},{"T",1},{"H",1}}], 1, #] - qfAmp[1,{{"H",1},{"T",1},{"H",1},{"T",1},{"H",1}}, #]]]==0 &], "B1 interleaved 1q:"];

Echo["==== elimination: #vars before/after, exact vs QF ====", ""];
demo[label_, n_, gates_] := With[{ps = buildPS[n, gates], r = reduce[buildPS[n, gates]]},
  Echo[{Length[ps["vars"]], Length[r["vars"]], matches[n, gates]}, label]];
(* Internal Hadamards (a wire re-Hadamarded later) are the eliminable ones.     *)
cliffChain[k_] := Join[{{"H", 1}}, Catenate[Table[{{"Z", 1}, {"H", 1}}, k]]];  (* H (Z H)^k *)
magicChain[k_] := Join[{{"H", 1}}, Catenate[Table[{{"T", 1}, {"H", 1}}, k]]];  (* H (T H)^k *)
deep2q[k_] := Join[{{"H", 1}, {"H", 2}}, Catenate[Table[{{"CZ", 1, 2}, {"H", 1}, {"H", 2}}, k]]];
demo["Clifford H(ZH)^8  1q:", 1, cliffChain[8]];                 (* 9 H's, internal ones eliminate *)
demo["Magic    H(TH)^8  1q:", 1, magicChain[8]];                 (* 9 H's, T's block elimination *)
demo["Clifford 2q deep   :", 2, deep2q[4]];
demo["Mixed: Clifford+1T :", 2, {{"H",1},{"Z",1},{"H",1},{"CZ",1,2},{"H",1},{"T",1},{"H",1},{"H",2}}];

Echo["pathsum-final done", "STATUS:"];
