(* Corrected state-vector rows for QF-Master-Platform-Comparison.md section 3.1.
   Input state is the true n-qubit register QuantumState["Register"[n]], not the
   misparsed QuantumState[Table[0,{n}],2] (a 4-qubit norm-0 state; see
   QF-Apply-Path-Deep-Audit.md section 0.1 and apply00e_state.wls).
   Readout is the statevector, matching what Aer/Cirq read; a density-matrix
   readout is a 2^16 x 2^16 dense object at n=16 on a real state. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[expr_, reps_ : 3] := Block[{t, best = Infinity},
    Do[ClearSystemCache[]; t = First @ AbsoluteTiming[expr]; best = Min[best, t], {reps}]; best];

svCircuit[n_] := QuantumCircuitOperator @ Flatten @ Table[
    Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[0.3] -> q, {q, n}]], {3}];

pr["sv_circuit_n12_gatecount", svCircuit[12]["GateCount"]];

(* default (TensorNetwork) method, min of 3 *)
Do[Block[{qc = svCircuit[n], psi = QuantumState["Register"[n]], out},
    pr["reg_n" <> ToString[n] <> "_norm", N[psi["Norm"]]];
    out = qc[psi];
    pr["out_n" <> ToString[n] <> "_norm", N[out["Norm"]]];
    pr["out_n" <> ToString[n] <> "_explicit", out["State"]["ExplicitLength"]];
    pr["sv_n" <> ToString[n] <> "_default_TN_ms", 1000 * tm[qc[psi]["StateVector"]]]
  ], {n, {12, 16}}];

(* Schrodinger method, min of 2 (slow) *)
Do[Block[{qc = svCircuit[n], psi = QuantumState["Register"[n]]},
    pr["sv_n" <> ToString[n] <> "_schrodinger_ms", 1000 * tm[qc[psi, Method -> "Schrodinger"]["StateVector"], 2]]
  ], {n, {12, 16}}];
