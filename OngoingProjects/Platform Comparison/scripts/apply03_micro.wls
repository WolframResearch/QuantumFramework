(* Micro-costs of the WL mechanisms in the apply path: property dispatch, cache guard,
   delegation, Sort/Tensor/Matrix chains, constructor costs, RZ decomposition. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

pr[l_, v_] := Print["RESULT|", l, "|", v];
SetAttributes[tm, HoldFirst];
tm[expr_] := Module[{t, best = Infinity}, Do[ClearSystemCache[]; t = First @ AbsoluteTiming[expr]; best = Min[best, t], {5}]; best];
SetAttributes[tmr, HoldFirst];  (* RepeatedTiming-style: average over k reps, min of 3 batches *)
tmr[expr_, k_ : 100] := Module[{t, best = Infinity}, Do[ClearSystemCache[]; t = First[AbsoluteTiming[Do[expr, {k}]]]/k; best = Min[best, t], {3}]; best];
us[t_] := Round[10.^6 t, 0.1];
msr[t_] := Round[1000. t, 0.001];

resolve[name_] := With[{hits = Names["*`" <> name]}, If[hits === {}, $Failed, Symbol[First[hits]]]];
$cacheVar = resolve["$QuantumFrameworkPropCache"];
pr["resolved_cache_var", $cacheVar =!= $Failed];

cacheOff := (Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False);
cacheOn := (Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = True);

(* ---------- 1. dispatch + cache-guard overhead ---------- *)
h = QuantumOperator["H", {3}];
h["Order"];  (* warm *)
pr["op_Order_warm_us", us @ tmr[h["Order"], 1000]];
cacheOff;
pr["op_Order_nocache_us", us @ tmr[h["Order"], 1000]];
pr["op_BasisParameterArity_nocache_us", us @ tmr[h["Basis"]["ParameterArity"], 200]];
pr["op_State_us", us @ tmr[h["State"], 1000]];
cacheOn;
pr["op_BasisParameterArity_warm_us", us @ tmr[h["Basis"]["ParameterArity"], 1000]];

(* ---------- 2. the Sort/Tensor/Matrix chains (recompute cost, cache off) ---------- *)
cacheOff;
pr["op_SortedQ_us", us @ tmr[h["SortedQ"], 200]];
pr["op_Sort_us", us @ tmr[h["Sort"], 200]];
pr["op_StateTensor_us", us @ tmr[h["StateTensor"], 200]];
pr["op_Tensor_us", us @ tmr[h["Tensor"], 200]];
pr["op_Matrix_us", us @ tmr[h["Matrix"], 200]];
pr["op_Computational_us", us @ tmr[h["Computational"], 200]];
pr["op_MatrixRepresentation_us", us @ tmr[h["MatrixRepresentation"], 200]];
pr["op_Reorder_us", us @ tmr[h["Reorder", {{3}, {3}}, False], 200]];
cnot = QuantumOperator["CNOT", {3, 4}];
pr["cnot_Sort_us", us @ tmr[cnot["Sort"], 200]];
pr["cnot_Tensor_us", us @ tmr[cnot["Tensor"], 200]];
pr["cnot_MatrixRepresentation_us", us @ tmr[cnot["MatrixRepresentation"], 100]];

(* ---------- 3. delegation cost (state property reached through the operator) ---------- *)
pr["op_AllProperties_us", us @ tmr[h["AllProperties"], 100]];
pr["state_AllProperties_us", us @ tmr[h["State"]["AllProperties"], 100]];
pr["op_delegated_AmplitudesList_us", us @ tmr[h["AmplitudesList"], 100]];
cacheOn;
pr["op_delegated_AmplitudesList_warm_us", us @ tmr[h["AmplitudesList"], 200]];

(* ---------- 4. constructor costs ---------- *)
(* named, fresh order each time to dodge any expression-level caches *)
pr["ctor_H_us", us @ tmr[QuantumOperator["H", {RandomInteger[{1, 10^9}]}], 50]];
pr["ctor_CNOT_us", us @ tmr[QuantumOperator["CNOT", {1, 2}], 50]];
pr["ctor_X_us", us @ tmr[QuantumOperator["X", {1}], 50]];
pr["ctor_state12_us", us @ tmr[QuantumState[Table[0, {12}], 2], 50]];
pr["ctor_QuantumBasis22_us", us @ tmr[QuantumBasis[2, 2], 100]];
pr["ctor_QuditBasis2_us", us @ tmr[QuditBasis[2], 200]];
Wolfram`QuantumFramework`PackageScope`$QuditBasisCaching = False;
pr["ctor_QuditBasis2_nocache_us", us @ tmr[QuditBasis[2], 100]];
pr["ctor_H_nobasiscache_us", us @ tmr[QuantumOperator["H", {1}], 50]];
Wolfram`QuantumFramework`PackageScope`$QuditBasisCaching = True;

(* ---------- 5. RZ cold decomposition ---------- *)
(* whole constructor, distinct numeric angles (cold) *)
angles = 0.3 + Range[20] 10.^-7;
t = First[AbsoluteTiming[QuantumOperator["RZ"[#], {1}] & /@ angles]]/Length[angles];
pr["ctor_RZ_numeric_cold_ms", msr[t]];
t = First[AbsoluteTiming[QuantumOperator["RZ"[0.3], {1}] & /@ Range[20]]]/20;
pr["ctor_RZ_numeric_samearg_ms", msr[t]];

(* decompose: PauliZ basis op, Exp, FullSimplify, outer wrap *)
pr["rz_pauliZ_ctor_ms", msr @ tmr[QuantumOperator["PauliZ"[2]], 20]];
pz = QuantumOperator["PauliZ"[2]];
pr["rz_exp_ms", msr @ tmr[Exp[-I 0.3/2 pz], 20]];
ex = Exp[-I 0.3/2 pz];
pr["rz_fullsimplify_ms", msr @ tmr[FullSimplify[ex], 20]];
fs = FullSimplify[ex];
pr["rz_label_wrap_ms", msr @ tmr[QuantumOperator[fs, "Label" -> Subscript["R", "Z"][0.3]], 20]];
lw = QuantumOperator[fs, "Label" -> Subscript["R", "Z"][0.3]];
pr["rz_outer_order_wrap_ms", msr @ tmr[QuantumOperator[lw, {1}], 20]];

(* symbolic angle for comparison *)
t = First[AbsoluteTiming[QuantumOperator["RZ"[\[Theta]], {1}]]];
pr["ctor_RZ_symbolic_cold_ms", msr[t]];

(* memoization prototype: name+numeric key *)
ClearAll[memoGate];
memoGate[name_, a_?NumericQ, order_] := memoGate[name, a, order] = QuantumOperator[name[a], order];
memoGate["RZ", 0.7, {1}];
pr["rz_memo_hit_us", us @ tmr[memoGate["RZ", 0.7, {1}], 1000]];

(* ---------- 6. what does op["Normal"] even do ---------- *)
res = Quiet @ Check[h["Normal"], "MESSAGE"];
pr["op_Normal_head", Head[res]];
t = tm[QuantumOperator["H", {1}]["Normal"]];
pr["op_Normal_ms", msr[t]];
(* circuit-level single-gate Normal *)
qc1 = QuantumCircuitOperator[{"H" -> 1}];
pr["qco1_Normal_ms", msr @ tm[qc1["Normal"]]];
