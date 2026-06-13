(* Prototype validation of the three candidate fixes from GitHub issue #5:
   A. numeric fast path (packed statevector fold, bypassing the TN object pipeline)
   B. memoized gate construction (name + numeric parameter + order)
   C. compile-once / reuse for repeated application
   Each prototype is validated against the default path's amplitudes. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
ms[t_] := Round[1000. t, 0.01];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_, L_ : 3] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {L}], 1]];
reg[n_] := QuantumState["Register"[n]];

n = 12; psi = reg[n]; qc = circ[n, 0.3];

(* reference *)
refT = First @ AbsoluteTiming[ref = qc[psi]];   (* cold *)
refT2 = First @ AbsoluteTiming[qc[psi]];        (* warm *)
pr["ref_cold_ms", ms[refT]]; pr["ref_warm_ms", ms[refT2]];
refV = Normal[N[ref["StateVector"]]];

(* ---------- Prototype A: numeric fast path ---------- *)
(* compile: dedup gate matrices by label, keep (dense packed matrix, qubit list) per gate *)
ClearAll[matOf];
matOf[op_] := matOf[op["Label"]] = Developer`ToPackedArray @ N @ Normal @ op["MatrixRepresentation"];
compileT = First @ AbsoluteTiming[
    specs = {matOf[#], #["InputOrder"]} & /@ qc["Flatten"]["Operators"];
];
pr["A_compile_ms", ms[compileT]];
pr["A_specs_packedQ", AllTrue[specs[[All, 1]], Developer`PackedArrayQ]];

applyGate[v_, {m_, qs_}, n_] := Block[{k = Length[qs], perm, inv, t},
    perm = Join[qs, Complement[Range[n], qs]];
    inv = Ordering[perm];
    t = Transpose[ArrayReshape[v, ConstantArray[2, n]], inv];
    t = m . ArrayReshape[t, {2^k, 2^(n - k)}];
    Flatten @ Transpose[ArrayReshape[t, ConstantArray[2, n]], perm]
];
fastApply[specs_, v0_, n_] := Fold[applyGate[#1, #2, n] &, v0, specs];

v0 = Developer`ToPackedArray @ N @ Normal @ SparseArray[{1} -> 1. + 0. I, {2^n}];
(* warm up + verify *)
outV = fastApply[specs, v0, n];
pr["A_correct_maxdiff", Max[Abs[outV - refV]]];
pr["A_out_packedQ", Developer`PackedArrayQ[outV]];
Do[
    ClearSystemCache[];
    pr["A_fold_run" <> ToString[k] <> "_ms", ms[First @ AbsoluteTiming[fastApply[specs, v0, n]]]],
    {k, 3}
];
(* end-to-end second call: compile is memoized via matOf, specs rebuild + fold *)
t = First @ AbsoluteTiming[
    specs2 = {matOf[#], #["InputOrder"]} & /@ qc["Flatten"]["Operators"];
    fastApply[specs2, v0, n]];
pr["A_end2end_warm_ms", ms[t]];

(* scaling: n = 16 *)
n16 = 16; qc16 = circ[n16, 0.3]; psi16 = reg[n16];
t = First @ AbsoluteTiming[qc16[psi16]]; pr["ref16_cold_ms", ms[t]];
t = First @ AbsoluteTiming[qc16[psi16]]; pr["ref16_warm_ms", ms[t]];
ClearAll[matOf16]; matOf16[op_] := matOf16[op["Label"]] = Developer`ToPackedArray @ N @ Normal @ op["MatrixRepresentation"];
specs16 = {matOf16[#], #["InputOrder"]} & /@ qc16["Flatten"]["Operators"];
v016 = Developer`ToPackedArray @ N @ Normal @ SparseArray[{1} -> 1. + 0. I, {2^n16}];
ref16V = Normal[N[qc16[psi16]["StateVector"]]];
out16 = fastApply[specs16, v016, n16];
pr["A16_correct_maxdiff", Max[Abs[out16 - ref16V]]];
Do[
    ClearSystemCache[];
    pr["A16_fold_run" <> ToString[k] <> "_ms", ms[First @ AbsoluteTiming[fastApply[specs16, v016, n16]]]],
    {k, 3}
];

(* ---------- Prototype B: memoized parametric-gate construction ---------- *)
t = First @ AbsoluteTiming[circ[12, 0.301]];
pr["B_build_named_cold_ms", ms[t]];
t = First @ AbsoluteTiming[circ[12, 0.301]];
pr["B_build_named_again_ms", ms[t]];
ClearAll[memoGate];
memoGate[spec_, ord_] := memoGate[spec, ord] = Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand[spec -> ord];
memoLayer[n_, a_] := Join[Table[memoGate["H", q], {q, n}], Table[memoGate["CNOT", {q, q + 1}], {q, n - 1}], Table[memoGate["RZ"[a], q], {q, n}]];
t = First @ AbsoluteTiming[qcM = QuantumCircuitOperator[Flatten[Table[memoLayer[12, 0.302], {3}], 1]]];
pr["B_build_memo_cold_ms", ms[t]];
t = First @ AbsoluteTiming[QuantumCircuitOperator[Flatten[Table[memoLayer[12, 0.302], {3}], 1]]];
pr["B_build_memo_warm_ms", ms[t]];
(* sanity: memoized circuit applies identically *)
outM = qcM[psi];
pr["B_correct_maxdiff", Max[Abs[Normal[N[outM["StateVector"]]] - Normal[N[circ[12, 0.302][psi]["StateVector"]]]]]];

(* ---------- Prototype C: compile-once for repeated application ---------- *)
(* C1: materialized operator at n=8 (2^8 x 2^8) *)
n8 = 8; qc8 = circ[n8, 0.3]; psi8 = reg[n8];
t = First @ AbsoluteTiming[qop8 = qc8["QuantumOperator"]];
pr["C_compile_qop8_ms", ms[t]];
t = First @ AbsoluteTiming[qop8[psi8]]; pr["C_qop8_apply1_ms", ms[t]];
t = First @ AbsoluteTiming[qop8[psi8]]; pr["C_qop8_apply2_ms", ms[t]];
t = First @ AbsoluteTiming[qc8[psi8]]; pr["C_ref8_warm_ms", ms[t]];

(* C2: dense matrix-vector at n=8 after extracting the matrix once *)
m8 = Developer`ToPackedArray @ N @ Normal @ qop8["MatrixRepresentation"];
v08 = Developer`ToPackedArray @ N @ Normal @ SparseArray[{1} -> 1. + 0. I, {2^n8}];
t = First @ AbsoluteTiming[m8 . v08]; pr["C_matvec8_ms", ms[t]];
