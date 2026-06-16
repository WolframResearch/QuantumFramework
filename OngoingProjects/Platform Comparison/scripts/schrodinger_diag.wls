(* Diagnosis of the dense state-vector path, special attention to Method -> "Schrodinger".
   Reproduces TN-default vs Schrodinger across (n, depth), finds the crossover, and
   decomposes the Schrodinger fold into (per-gate TN apply) + (FullSimplify tax).
   Working tree at HEAD f9dc1cdc. HoldFirst timing, ClearSystemCache per rep, min-of-reps. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[expr_, reps_ : 3] := Block[{t, best = Infinity},
   Do[ClearSystemCache[]; t = First @ AbsoluteTiming[expr]; best = Min[best, t], {reps}];
   Round[1000. best, 0.001]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_, L_ : 3] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {L}], 1]];
reg[n_] := QuantumState["Register"[n]];

Print["=== PART 1: CROSSOVER SWEEP (1 layer, numeric) TN vs Schrodinger ==="];
Do[Block[{qc = circ[n, 0.3, 1], psi = reg[n], g, tTN, tSch},
   g = qc["Gates"];
   tTN = tm[qc[psi]["StateVector"], 4];
   tSch = tm[qc[psi, Method -> "Schrodinger"]["StateVector"], 3];
   pr["sweep1_n" <> ToString[n], {"gates", g, "TN_ms", tTN, "Sch_ms", tSch, "ratio_Sch/TN", Round[tSch/tTN, 0.01]}]
   ], {n, {1, 2, 3, 4, 5, 6, 8}}];

Print["=== PART 1b: tiny fixed circuits (Bell, 1-3 gates) ==="];
Block[{bell = QuantumCircuitOperator[{"H" -> 1, "CNOT" -> {1, 2}}], psi = reg[2]},
   pr["bell_TN_ms", tm[bell[psi]["StateVector"], 6]];
   pr["bell_Sch_ms", tm[bell[psi, Method -> "Schrodinger"]["StateVector"], 6]]];
Block[{one = QuantumCircuitOperator[{"H" -> 1}], psi = reg[1]},
   pr["1gate_TN_ms", tm[one[psi]["StateVector"], 8]];
   pr["1gate_Sch_ms", tm[one[psi, Method -> "Schrodinger"]["StateVector"], 8]]];

Print["=== PART 2: DEPTH SWEEP at fixed n=4 ==="];
Do[Block[{qc = circ[4, 0.3, L], psi = reg[4], g, tTN, tSch},
   g = qc["Gates"];
   tTN = tm[qc[psi]["StateVector"], 4];
   tSch = tm[qc[psi, Method -> "Schrodinger"]["StateVector"], 3];
   pr["depth_n4_L" <> ToString[L], {"gates", g, "TN_ms", tTN, "Sch_ms", tSch}]
   ], {L, {1, 2, 4, 8, 16}}];

Print["=== PART 3: PROFILE the Schrodinger fold (decompose FullSimplify tax) ==="];
Do[Block[{qc = circ[n, 0.3, 3], psi = reg[n], ops, full, nofs, fsOnly, perGate, oneApply, st},
   ops = qc["NormalOperators"];
   (* (a) full Schrodinger semantics, exactly as the kernel does it *)
   full = tm[Fold[(FullSimplify @ #2[#1, Method -> {"TensorNetwork", "Computational" -> False}]) &, psi, ops], 2];
   (* (b) same fold WITHOUT the FullSimplify wrapper *)
   nofs = tm[Fold[#2[#1, Method -> {"TensorNetwork", "Computational" -> False}] &, psi, ops], 2];
   (* (c) cost of FullSimplify alone on the final numeric state object, x (#gates) *)
   st = Fold[#2[#1, Method -> {"TensorNetwork", "Computational" -> False}] &, psi, ops];
   fsOnly = tm[FullSimplify[st], 3];
   (* (d) one single-gate TN apply (an H on qubit 1) *)
   oneApply = tm[QuantumOperator["H", {1}][psi, Method -> {"TensorNetwork", "Computational" -> False}], 5];
   pr["profile_n" <> ToString[n], {"gates", Length[ops],
      "full_ms", full, "noFS_ms", nofs, "FStax_ms", full - nofs,
      "FS_on_final_state_ms", fsOnly, "single_gate_apply_ms", oneApply,
      "bulk_TN_ms", tm[qc[psi]["StateVector"], 3]}]
   ], {n, {6, 8, 10}}];

Print["=== PART 4: is the state packed or sparse on the dense path? ==="];
Block[{qc = circ[10, 0.3, 3], psi = reg[10], out, sv},
   out = qc[psi];
   sv = out["StateVector"];
   pr["n10_out_head", Head[sv]];
   pr["n10_out_packedQ", Developer`PackedArrayQ[sv]];
   pr["n10_out_explicitLen_vs_dim", {out["State"]["ExplicitLength"], 2^10}];
   pr["n10_normalSV_packedQ", Developer`PackedArrayQ[Normal[sv]]]];

Print["DONE"];
