PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

pr["names_prop", Names["*`QuantumOperatorProp"]];
pr["names_tnc", Names["*`TensorNetworkCompile"]];
pr["names_fos", Names["*`FromOperatorShorthand"]];
pr["n_qf_names", Length[Names["Wolfram`QuantumFramework`*"]]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

n = 12;
psi = QuantumState[Table[0, {n}], 2];
qc = circ[n, 0.3];

mem[l_] := pr[l <> "_MB", Round[MemoryInUse[]/2.^20, 1.]];

mem["m0"];
t = First @ AbsoluteTiming[res = MemoryConstrained[qc[psi], 8 * 2^30]];
pr["apply_ms", Round[1000. t, 1.]]; pr["apply_aborted", res === $Aborted];
mem["m1_after_apply"];

t = First @ AbsoluteTiming[composed = qc /* QuantumCircuitOperator[{}];
    composed["InputDimensions"][[Ordering[composed["InputOrder"]]]] == psi["OutputDimensions"]];
pr["s0_ms", Round[1000. t, 1.]]; mem["m2_after_s0"];

t = First @ AbsoluteTiming[qcoF = composed["Flatten"]];
pr["s1_ms", Round[1000. t, 1.]]; mem["m3_after_s1"];

t = First @ AbsoluteTiming[sorted = qcoF["Sort"]];
pr["s2_ms", Round[1000. t, 1.]]; mem["m4_after_s2"];

res = MemoryConstrained[
    First @ AbsoluteTiming[circuit2 = QuantumCircuitOperator[{psi -> sorted["FullInputOrder"]} /* sorted, "Label" -> None]],
    8 * 2^30];
pr["s3_ms_or_abort", If[NumberQ[res], Round[1000. res, 1.], "ABORT"]]; mem["m5_after_s3"];

res = MemoryConstrained[First @ AbsoluteTiming[circN = circuit2["Normal"]], 8 * 2^30];
pr["s4_normal_ms_or_abort", If[NumberQ[res], Round[1000. res, 1.], "ABORT"]]; mem["m6_after_s4"];
