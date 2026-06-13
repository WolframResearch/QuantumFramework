(* Bisect the n=12 memory blowup in QuantumCircuitOperator[{psi -> order} /* circuit]. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

probe[l_, expr_] := Block[{r},
    r = MemoryConstrained[{First[AbsoluteTiming[expr]], MaxMemoryUsed[]}, 6 * 2^30];
    If[r === $Aborted,
        pr[l, "ABORT_6GB"],
        pr[l, {Round[1000. r[[1]], 1.], Round[r[[2]]/2.^20, 1.]}]];
];
SetAttributes[probe, HoldRest];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

n = 12;
psi = QuantumState[Table[0, {n}], 2];
qc = circ[n, 0.3];
sorted = (qc /* QuantumCircuitOperator[{}])["Flatten"]["Sort"];

probe["a_state_to_op", op0 = QuantumOperator[psi]];
probe["b_fos_rule", fos = Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand[psi -> Range[12]]];
probe["c_qco_of_list", x1 = QuantumCircuitOperator[{psi -> Range[12]}]];
probe["d_rightcomp", x2 = {psi -> Range[12]} /* sorted];
probe["e_label_none", x3 = QuantumCircuitOperator[x2, "Label" -> None]];
probe["f_normal", x4 = x3["Normal"]];
pr["fos_dv_count", Length[DownValues[Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand]]];
