(* Does Method -> "Schrodinger" have a SYMBOLIC niche? The per-gate FullSimplify could
   keep symbolic amplitudes compact where the bulk TN path returns un-simplified ones,
   and TN can OOM on QPE-style circuits. Test compactness + time on symbolic circuits. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[expr_, reps_ : 3] := Block[{t, best = Infinity},
   Do[ClearSystemCache[]; t = First @ AbsoluteTiming[expr]; best = Min[best, t], {reps}];
   Round[1000. best, 0.001]];

slayer[n_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[\[Theta]] -> q, {q, n}]];
scirc[n_, L_ : 1] := QuantumCircuitOperator[Flatten[Table[slayer[n], {L}], 1]];
reg[n_] := QuantumState["Register"[n]];

Print["=== SYMBOLIC circuits: TN-default vs Schrodinger (time + LeafCount of statevector) ==="];
Do[Block[{qc = scirc[n, 1], psi = reg[n], tnOut, schOut, tnT, schT},
   tnT = tm[qc[psi]["StateVector"], 2];
   schT = tm[qc[psi, Method -> "Schrodinger"]["StateVector"], 2];
   tnOut = qc[psi]["StateVector"];
   schOut = qc[psi, Method -> "Schrodinger"]["StateVector"];
   pr["sym_n" <> ToString[n], {"gates", qc["Gates"],
      "TN_ms", tnT, "Sch_ms", schT,
      "TN_leafcount", LeafCount[Normal[tnOut]], "Sch_leafcount", LeafCount[Normal[schOut]],
      "amplitudes_equal_afterFS", Simplify[Normal[tnOut] - Normal[schOut]] === Table[0, {2^n}]}]
   ], {n, {1, 2, 3}}];

Print["=== Single symbolic RZ amplitude shape (what FullSimplify buys per gate) ==="];
Block[{qc = QuantumCircuitOperator[{"H" -> 1, "RZ"[\[Theta]] -> 1}], psi = reg[1]},
   pr["sym1_TN_sv", Normal[qc[psi]["StateVector"]]];
   pr["sym1_Sch_sv", Normal[qc[psi, Method -> "Schrodinger"]["StateVector"]]]];

Print["DONE"];
