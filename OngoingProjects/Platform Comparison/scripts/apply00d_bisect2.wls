(* Bisect level 2: which exact call inside FromOperatorShorthand[psi -> Range[12]] allocates GBs,
   and what the real apply's memory high-water mark is. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
probe[l_, expr_] := Block[{r},
    r = MemoryConstrained[{First[AbsoluteTiming[expr]], MaxMemoryUsed[]}, 4 * 2^30];
    If[r === $Aborted, pr[l, "ABORT_4GB"],
        pr[l, {Round[1000. r[[1]], 1.], Round[r[[2]]/2.^20, 1.]}]];
];
SetAttributes[probe, HoldRest];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

n = 12;
psi = QuantumState[Table[0, {n}], 2];
qc = circ[n, 0.3];

(* real apply high-water mark, fresh kernel *)
pr["maxmem_before_apply_MB", Round[MaxMemoryUsed[]/2.^20, 1.]];
t = First @ AbsoluteTiming[qc[psi]];
pr["apply_ms", Round[1000. t, 1.]];
pr["maxmem_after_apply_MB", Round[MaxMemoryUsed[]/2.^20, 1.]];

op0 = QuantumOperator[psi];
pr["op0_order", op0["Order"]];

probe["p2_qo_psi_range", QuantumOperator[psi, Range[12]]];
probe["p3_qo_op0_range", QuantumOperator[op0, Range[12]]];
probe["p4_reorder", op0["Reorder", {Range[12], Range[12]}, True]];
probe["p5_qmo_rule", QuantumMeasurementOperator[psi -> Range[12]]];
probe["p6_fos_rule", Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand[psi -> Range[12]]];
