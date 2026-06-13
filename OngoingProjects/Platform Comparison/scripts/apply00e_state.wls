(* What is QuantumState[Table[0,{n}],2] really, and does the benchmark contract a zero state? *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

n = 12;
psiBad = QuantumState[Table[0, {n}], 2];
pr["bad_dims", psiBad["Dimensions"]];
pr["bad_qudits", psiBad["OutputQudits"]];
pr["bad_norm", psiBad["Norm"]];
pr["bad_vec_explicit", psiBad["State"]["ExplicitLength"]];

psiReg = QuantumState["Register"[n]];
pr["reg_dims", psiReg["Dimensions"]];
pr["reg_norm", psiReg["Norm"]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];
qc = circ[n, 0.3];

(* timings cold-in-this-kernel, one shot each *)
t1 = First @ AbsoluteTiming[out1 = qc[psiBad]];
pr["apply_bad_ms", Round[1000. t1, 1.]];
pr["out_bad_norm", out1["Norm"]];
pr["out_bad_qudits", out1["OutputQudits"]];

t2 = First @ AbsoluteTiming[out2 = qc[psiReg]];
pr["apply_reg_ms", Round[1000. t2, 1.]];
pr["out_reg_norm", N[out2["Norm"]]];
pr["out_reg_explicit", out2["State"]["ExplicitLength"]];

(* second application, warm *)
t3 = First @ AbsoluteTiming[qc[psiReg]];
pr["apply_reg_warm_ms", Round[1000. t3, 1.]];
t4 = First @ AbsoluteTiming[qc[psiBad]];
pr["apply_bad_warm_ms", Round[1000. t4, 1.]];

(* zero-state artifact check: dense random state *)
psiRnd = QuantumState["RandomPure"[n]];
t5 = First @ AbsoluteTiming[qc[psiRnd]];
pr["apply_rnd_ms", Round[1000. t5, 1.]];
