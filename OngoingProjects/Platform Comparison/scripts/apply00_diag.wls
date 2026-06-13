(* Diagnostic: resolve internal symbol contexts; probe the state-prepend memory blowup at small n. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

pr["ctx_QuantumOperatorProp", Names["Wolfram`QuantumFramework`*QuantumOperatorProp"]];
pr["ctx_TensorNetworkCompile", Names["Wolfram`QuantumFramework`*TensorNetworkCompile*"]];
pr["ctx_FromOperatorShorthand", Names["Wolfram`QuantumFramework`*FromOperatorShorthand*"]];
pr["ctx_SparseArrayFlatten", Names["Wolfram`QuantumFramework`*SparseArrayFlatten*"]];

qop = Symbol[First @ Names["Wolfram`QuantumFramework`*QuantumOperatorProp"]];
pr["dv_count_qop", Length[DownValues[qop]]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

(* probe the prepend-state blowup with growing n under a 6 GB ceiling *)
Do[
    Module[{qc = circ[n, 0.3], psi = QuantumState[Table[0, {n}], 2], sorted, t, m0, m1, res},
        sorted = qc["Flatten"]["Sort"];
        m0 = MemoryInUse[];
        res = MemoryConstrained[
            t = AbsoluteTiming[QuantumCircuitOperator[{psi -> sorted["FullInputOrder"]} /* sorted, "Label" -> None]];
            {First[t], MaxMemoryUsed[]},
            6 * 2^30
        ];
        m1 = MemoryInUse[];
        If[res === $Aborted,
            pr["prepend_n" <> ToString[n], "MEMORY_ABORT"],
            pr["prepend_n" <> ToString[n] <> "_ms_maxmem", {Round[1000. res[[1]], 0.1], Round[res[[2]]/2.^20, 1.]}]
        ];
        pr["meminuse_after_n" <> ToString[n] <> "_MB", Round[(m1 - m0)/2.^20, 1.]];
    ],
    {n, {4, 6, 8, 10}}
]
