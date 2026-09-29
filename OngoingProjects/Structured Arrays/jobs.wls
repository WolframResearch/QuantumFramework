PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

logFile = "/private/tmp/claude-10838/-Users-mohammadb-Documents-GitHub-QuantumFramework/36bcbe1e-22b4-45d6-a9b2-99d6625ab2c6/scratchpad/jobs.log";
log[s__] := With[{st = OpenAppend[logFile]}, WriteString[st, StringJoin[ToString[#, InputForm, PageWidth -> 300] & /@ {s}], "\n"]; Close[st]];
log["START ", $Version, " ", DateString[]];

SetAttributes[tc, HoldFirst];
tc[e_] := Block[{r, t}, ClearSystemCache[]; {t, r} = AbsoluteTiming[TimeConstrained[e, 300, $TimedOut]]; {If[r === $TimedOut, "timeout>300s", t], r}];
vec[st_] := Normal[st["StateVector"]];
ts = TargetStructure -> "Structured";

(* inverse QFT on the first n qubits of a state with 2^n x r amplitudes, as an FFT along that axis *)
invQFT[v_, n_] := With[{m = ArrayReshape[v, {2^n, Length[v]/2^n}]}, Flatten[Transpose[InverseFourier /@ Transpose[m]]]];

zzJ[J_, n_] := With[{z = 1 - 2 Tuples[{0, 1}, n]}, Total[Flatten[Table[J[[i, j]] z[[All, i]] z[[All, j]], {i, n}, {j, i + 1, n}], 1]]];
qaoaGates[J_, n_] := Join[
    Table[QuantumOperator["H", {i}], {i, n}],
    Flatten @ Table[QuantumOperator["R"[2 \[FormalGamma] J[[i, j]], "ZZ"], {i, j}], {i, n}, {j, i + 1, n}],
    Table[QuantumOperator["RX"[2 \[FormalBeta]], {i}], {i, n}]];
mixer[n_] := QuantumCircuitOperator[Table[QuantumOperator["RX"[2 \[FormalBeta]], {i}], {i, n}], "Parameters" -> {\[FormalBeta]}];

(* ---------- convention checks at small size ---------- *)
SeedRandom[3]; u = QuantumOperator["RandomUnitary"[{2, 2}]];
Module[{n = 3, qpe, els, body, front, ref, hyb},
  qpe = QuantumCircuitOperator["PhaseEstimation"[u, n]]; els = qpe["Elements"];
  body = QuantumCircuitOperator[Drop[els, -n]]; front = QuantumCircuitOperator[Drop[els, -(n + 1)]];
  log["CHECK QPE last element of body: ", Head[Last[body["Elements"]]], " ", Last[body["Elements"]]["Label"]];
  ref = vec[body[Method -> "Schrodinger"]];
  hyb = invQFT[vec[front[]], n];
  log["CHECK QPE FFT-replacement vs QF body (Schrodinger): ", Norm[ref - hyb]];
];
Module[{n = 3, J, h, cost, g0 = 0.37, U},
  SeedRandom[5]; J = RandomReal[{-1, 1}, {n, n}]; h = zzJ[J, n];
  cost = QuantumCircuitOperator[Flatten @ Table[QuantumOperator["R"[2 g0 J[[i, j]], "ZZ"], {i, j}], {i, n}, {j, i + 1, n}]];
  U = Normal[cost["QuantumOperator"]["Matrix"]];
  log["CHECK QAOA cost layer == Exp[-I g h] up to phase: ", With[{r = U . ConjugateTranspose[DiagonalMatrix[Exp[-I g0 h]]]}, Max[Abs[r - r[[1, 1]] IdentityMatrix[2^n]]]]];
];

(* ---------- QPE jobs ---------- *)
qpeRun[n_] := Module[{qpe, els, body, front, res = <||>, ref, t, st, tf, sf, tfft, hyb},
  qpe = QuantumCircuitOperator["PhaseEstimation"[u, n]]; els = qpe["Elements"];
  body = QuantumCircuitOperator[Drop[els, -n]]; front = QuantumCircuitOperator[Drop[els, -(n + 1)]];
  Scan[Function[m,
     {t, st} = tc[body[Method -> m]];
     res[m] = t;
     If[NumberQ[t], If[! ValueQ[ref], ref = vec[st]]; log["QPE n=", n, " QF full circuit, Method -> ", m, ": ", t, " s, diff vs first: ", Norm[vec[st] - ref]], log["QPE n=", n, " QF full circuit, Method -> ", m, ": ", t]]],
    If[n <= 8, {Automatic, {"TensorNetwork", "Path" -> "Greedy"}, {"TensorNetwork", "Path" -> "Optimal"}, "Schrodinger"}, {Automatic, {"TensorNetwork", "Path" -> "Greedy"}, "Schrodinger"}]];
  Scan[Function[m,
     {tf, sf} = tc[front[Method -> m]];
     If[NumberQ[tf],
      tfft = First[AbsoluteTiming[hyb = invQFT[vec[sf], n];]];
      log["QPE n=", n, " QF without inverse QFT (Method -> ", m, "): ", tf, " s + FFT: ", tfft, " s = ", tf + tfft, " s, diff vs QF: ", If[ValueQ[ref], Norm[hyb - ref], "no ref"]],
      log["QPE n=", n, " QF without inverse QFT (Method -> ", m, "): ", tf]]],
    {Automatic, {"TensorNetwork", "Path" -> "Greedy"}}];
];
Scan[qpeRun, {8, 10, 12}];

(* ---------- QAOA jobs: <C> on a 5 x 5 grid of (gamma, beta), p = 1 ---------- *)
qaoaRun[n_] := Module[{J, h, qc, grid, v0, mix, t, vals, ref, th, hv},
  SeedRandom[5]; J = RandomReal[{-1, 1}, {n, n}]; h = zzJ[J, n];
  qc = QuantumCircuitOperator[qaoaGates[J, n], "Parameters" -> {\[FormalGamma], \[FormalBeta]}];
  grid = Tuples[{Subdivide[0.1, 1.0, 4], Subdivide[0.1, 1.0, 4]}];
  Scan[Function[m,
     {t, vals} = tc[Map[With[{v = vec[qc[#[[1]], #[[2]]][Method -> m]]}, Re[Conjugate[v] . (h v)]] &, grid]];
     If[NumberQ[t] && ! ValueQ[ref], ref = vals];
     log["QAOA n=", n, " QF full circuit, Method -> ", m, ": ", t, If[NumberQ[t], {" s, max diff vs first: ", Max[Abs[vals - ref]]}, ""]]],
    {Automatic, {"TensorNetwork", "Path" -> "Greedy"}, "Schrodinger"}];
  mix = mixer[n]; v0 = ConstantArray[1/Sqrt[2.^n], 2^n];
  Scan[Function[m,
     {th, hv} = tc[Map[Function[gb, With[{v1 = DiagonalMatrix[Exp[-I gb[[1]] h], ts] . v0},
          With[{v = vec[mix[gb[[2]]][QuantumState[v1], Method -> m]]}, Re[Conjugate[v] . (h v)]]]], grid]];
     log["QAOA n=", n, " structured cost layer + QF mixer (Method -> ", m, "): ", th, If[NumberQ[th] && ValueQ[ref], {" s, max diff vs QF: ", Max[Abs[hv - ref]]}, ""]]],
    {Automatic, {"TensorNetwork", "Path" -> "Greedy"}}];
];
Scan[qaoaRun, {10, 12, 14}];
log["DONE ", DateString[]];
