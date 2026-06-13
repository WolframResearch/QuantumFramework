(* Packing diagnostics + ablations: cache off, pre-normalized circuit, raw-matrix gates,
   depth sweep (superlinearity check), width sweep. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
ms[t_] := Round[1000. t, 0.1];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_, L_ : 3] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {L}], 1]];
reg[n_] := QuantumState["Register"[n]];

(* ---------- representation audit ---------- *)
n = 12; psi = reg[n]; qc = circ[n, 0.3];
pr["psi_state_head", Head[psi["State"]]];
pr["psi_packedQ", Developer`PackedArrayQ[psi["State"]]];
h = QuantumOperator["H", {1}];
pr["h_tensor_head", Head[h["Tensor"]]];
out = qc[psi];
pr["out_state_head", Head[out["State"]]];
pr["out_explicit", out["State"]["ExplicitLength"]];

(* unpacking diagnostics on a warm apply *)
On["Packing"];
qc[psi];
Off["Packing"];
Print["PACKING_SECTION_DONE"];

(* ---------- ablation: property cache off ---------- *)
t = First @ AbsoluteTiming[qc[psi]];
pr["warm_cacheON_ms", ms[t]];
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
t = First @ AbsoluteTiming[qc[psi]];
pr["warm_cacheOFF_ms", ms[t]];
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = True;

(* ---------- ablation: qudit-basis construction cache off ---------- *)
Wolfram`QuantumFramework`PackageScope`$QuditBasisCaching = False;
t = First @ AbsoluteTiming[qc[psi]];
pr["warm_quditcacheOFF_ms", ms[t]];
Wolfram`QuantumFramework`PackageScope`$QuditBasisCaching = True;

(* ---------- ablation: pre-normalized circuit ---------- *)
qcN = qc["Normal"];
t = First @ AbsoluteTiming[qcN[psi]];
pr["prenormal_first_ms", ms[t]];
t = First @ AbsoluteTiming[qcN[psi]];
pr["prenormal_second_ms", ms[t]];

(* ---------- ablation: raw numeric matrices instead of named gates ---------- *)
rawGate[m_, qs_] := QuantumOperator[m, qs];
hM = N @ Normal @ QuantumOperator["H"]["MatrixRepresentation"];
cxM = N @ Normal @ QuantumOperator["CNOT"]["MatrixRepresentation"];
rzM = N @ Normal @ QuantumOperator["RZ"[0.3]]["MatrixRepresentation"];
qcRaw = QuantumCircuitOperator[Flatten[Table[Join[
    Table[rawGate[hM, {q}], {q, n}],
    Table[rawGate[cxM, {q, q + 1}], {q, n - 1}],
    Table[rawGate[rzM, {q}], {q, n}]], {3}], 1]];
t = First @ AbsoluteTiming[qcRaw[psi]];
pr["raw_first_ms", ms[t]];
t = First @ AbsoluteTiming[qcRaw[psi]];
pr["raw_second_ms", ms[t]];

(* ---------- depth sweep: superlinearity check (warm, 2nd apply reported) ---------- *)
Do[
    Module[{qcL = circ[12, 0.3 + L 10.^-6, L], t1, t2},
        t1 = First @ AbsoluteTiming[qcL[psi]];
        t2 = First @ AbsoluteTiming[qcL[psi]];
        pr["depth_L" <> ToString[L] <> "_gates" <> ToString[qcL["GateCount"]] <> "_cold_warm_ms", {ms[t1], ms[t2]}]
    ],
    {L, {1, 2, 4, 8, 16, 24}}
]

(* ---------- width sweep ---------- *)
Do[
    Module[{qcW = circ[w, 0.3 + w 10.^-6], psiW = reg[w], t1, t2},
        t1 = First @ AbsoluteTiming[qcW[psiW]];
        t2 = First @ AbsoluteTiming[qcW[psiW]];
        pr["width_n" <> ToString[w] <> "_gates" <> ToString[qcW["GateCount"]] <> "_cold_warm_ms", {ms[t1], ms[t2]}]
    ],
    {w, {8, 12, 16}}
]

(* ---------- parametric reuse idiom ---------- *)
qcP = QuantumCircuitOperator[Flatten[Table[layer[12, \[Theta]], {3}], 1], "Parameters" -> {\[Theta]}];
t = First @ AbsoluteTiming[qcP[<|\[Theta] -> 0.37|>][psi]];
pr["param_subst_apply_first_ms", ms[t]];
t = First @ AbsoluteTiming[qcP[<|\[Theta] -> 0.38|>][psi]];
pr["param_subst_apply_second_ms", ms[t]];
