PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Get["/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Structured Arrays/prototype/fusion.wl"];

stage = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], "check"];
logFile = "/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Structured Arrays/prototype/bench-" <> stage <> ".log";
log[s__] := With[{st = OpenAppend[logFile]}, WriteString[st, StringJoin[ToString[#, InputForm, PageWidth -> 300] & /@ {s}], "\n"]; Close[st]];
log["START ", stage, " ", $Version, " ", DateString[]];

SetAttributes[tc, HoldFirst];
tc[e_] := Block[{r, t}, ClearSystemCache[]; {t, r} = AbsoluteTiming[TimeConstrained[e, 300, $TimedOut]]; {If[r === $TimedOut, "timeout>300s", t], r}];

(* ---------- benchmark circuits ---------- *)
qftC[n_] := QuantumCircuitOperator["Fourier"[n]];
qpeC[n_] := Module[{u, els}, SeedRandom[3]; u = QuantumOperator["RandomUnitary"[{2, 2}]];
    els = QuantumCircuitOperator["PhaseEstimation"[u, n]]["Elements"]; QuantumCircuitOperator[Drop[els, -n]]];
bvC[n_] := (SeedRandom[4]; QuantumCircuitOperator["BernsteinVazirani"[StringJoin[ToString /@ RandomInteger[1, n]]]]);
qaoaC[n_, p_] := Module[{J}, SeedRandom[5]; J = RandomReal[{-1, 1}, {n, n}];
    QuantumCircuitOperator[Join[
        Table[QuantumOperator["H", {i}], {i, n}],
        Catenate @ Table[Join[
            Flatten @ Table[QuantumOperator["R"[2 0.3 l J[[i, j]], "ZZ"], {i, j}], {i, n}, {j, i + 1, n}],
            Table[QuantumOperator["RX"[2 0.2 l], {i}], {i, n}]], {l, p}]]]];
cliffordTC[n_, g_] := (SeedRandom[6]; QuantumCircuitOperator[Table[
    With[{r = RandomReal[], q = RandomSample[Range[n], 2]},
        Which[r < 0.2, QuantumOperator["H", {q[[1]]}], r < 0.4, QuantumOperator["T", {q[[1]]}], r < 0.5, QuantumOperator["S", {q[[1]]}],
            r < 0.6, QuantumOperator["X", {q[[1]]}], r < 0.85, QuantumOperator["CNOT", q], True, QuantumOperator["CZ", q]]], g]]);
denseC[n_, g_] := (SeedRandom[8]; QuantumCircuitOperator[Table[QuantumOperator[QuantumOperator["RandomUnitary"[{2, 2}]], RandomSample[Range[n], 2]], g]]);

width[qc_] := Max[Flatten[{#[[2]]} & /@ circuitItems[qc, True]]];

(* ---------- one circuit, all arms ---------- *)
runCircuit[name_, qc0_, qfMethods_] := Module[{qc = withoutMeasurements[qc0], n, v0, qs, ref = None, t, st, steps, tcomp, tapp, out, share},
    n = width[qc]; SeedRandom[1]; v0 = Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]]; qs = QuantumState[v0];
    share = fusedShare[qc];
    log[name, " n=", n, " gates=", Length[circuitItems[qc, False]], " share in monomial runs=", share];
    Scan[Function[m,
        {t, st} = tc[qc[qs, Method -> m]];
        If[NumberQ[t], If[ref === None, ref = Normal[st["StateVector"]]]; log["  QF ", m, ": ", t, " s, diff ", Norm[Normal[st["StateVector"]] - ref]], log["  QF ", m, ": ", t]]],
        qfMethods];
    {t, out} = tc[applyGateByGate[qc, v0, n]];
    log["  B gate-by-gate state vector: ", t, " s", If[ref =!= None && NumberQ[t], {" diff ", Norm[out - ref]}, ""]];
    Scan[Function[{arm},
        {tcomp, steps} = tc[compileSteps[qc, n, arm[[2]], arm[[3]]]];
        {tapp, out} = tc[applySteps[steps, v0, n, arm[[2]]]];
        log["  ", arm[[1]], ": compile ", tcomp, " s + apply ", tapp, " s = ", If[NumberQ[tcomp] && NumberQ[tapp], tcomp + tapp, "n/a"], " s, steps ", Length[steps],
            If[ref =!= None && NumberQ[tapp], {" diff ", Norm[out - ref]}, ""]]],
        {{"C fusion as structured arrays (Dot)", "Structured", False}, {"D fusion as vectors (Part, Times)", "Vector", False}, {"E D + QFT blocks as FFT", "Vector", True}}];
];

greedy = {"TensorNetwork", "Path" -> "Greedy"};
Switch[stage,
    "check",
        Scan[runCircuit[#[[1]], #[[2]], {"Schrodinger", Automatic}] &, {
            {"QFT", qftC[4]}, {"QPE", qpeC[3]}, {"BV", bvC[4]}, {"QAOA p=2", qaoaC[4, 2]}, {"Clifford+T", cliffordTC[4, 40]}, {"dense", denseC[4, 12]}}],
    "bench",
        Scan[runCircuit[#[[1]], #[[2]], {Automatic, greedy}] &, {
            {"QFT", qftC[14]}, {"QFT", qftC[16]},
            {"QPE", qpeC[12]}, {"QPE", qpeC[14]},
            {"BV", bvC[16]},
            {"QAOA p=2", qaoaC[12, 2]}, {"QAOA p=2", qaoaC[14, 2]},
            {"Clifford+T", cliffordTC[12, 500]}, {"Clifford+T", cliffordTC[14, 500]},
            {"dense", denseC[12, 150]}, {"dense", denseC[14, 150]}}]
];
log["DONE ", DateString[]];
