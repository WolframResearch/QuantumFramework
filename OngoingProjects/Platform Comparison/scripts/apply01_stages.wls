(* Deep audit of the QF circuit-apply hot path: stage decomposition of TensorNetworkApply.
   Output: RESULT|label|value lines, times in ms.
   Input state: QuantumState["Register"[n]] (the true |0...0> register).
   NOTE: the first-pass benchmark state QuantumState[Table[0,{n}],2] is a 4-qubit ZERO vector
   (norm 0); see apply00e_state.wls. We time both.
   Cache regimes: run1 = cold kernel, run2 = perturbed circuit (RZ angle differs, H/CNOT warm),
   run3 = identical circuit (fully warm). *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`TensorNetworks`"];

pr[l_, v_] := Print["RESULT|", l, "|", v];
ms[t_] := Round[1000. t, 0.01];

resolve[name_] := With[{hits = Names["*`" <> name]}, If[hits === {}, $Failed, Symbol[First[hits]]]];
QTN = resolve["QuantumTensorNetwork"];
saFlatten = resolve["SparseArrayFlatten"];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

n = 12;
psi = QuantumState["Register"[n]];

props = resolve /@ {"QuantumOperatorProp", "QuantumStateProp", "QuantumBasisProp", "QuantumCircuitOperatorProp", "QuditBasisProp"};
dvCount[] := Length[DownValues[#]] & /@ props;
dvBytes[] := ByteCount[DownValues[#]] & /@ props;

pr["dv_count_after_load", dvCount[]];

stageRun[qc_, label_] := Block[{t, composed, qcoF, sorted, circuit2, circN, basis, info, computationalQ, net, raw, res, state, total, sum = 0.},
    total = First @ AbsoluteTiming[qc[psi]];
    pr[label <> "_total_ms", ms[total]];

    (* s0: SubValue wrapper at QuantumCircuitOperator.m:120-137 + dim guard at :97 *)
    t = First @ AbsoluteTiming[
        composed = qc /* QuantumCircuitOperator[{}];
        composed["InputDimensions"][[Ordering[composed["InputOrder"]]]] == psi["OutputDimensions"]];
    sum += t; pr[label <> "_s0_wrapper_dimguard_ms", ms[t]];

    (* s1: qco["Flatten"] (QuantumCircuitOperator.m:106) *)
    t = First @ AbsoluteTiming[qcoF = composed["Flatten"]];
    sum += t; pr[label <> "_s1_flatten_ms", ms[t]];

    (* s2: qco["Sort"] (TensorNetwork.m:252) *)
    t = First @ AbsoluteTiming[sorted = qcoF["Sort"]];
    sum += t; pr[label <> "_s2_sort_ms", ms[t]];

    (* s3: prepend state as a gate (TensorNetwork.m:255) *)
    t = First @ AbsoluteTiming[circuit2 = QuantumCircuitOperator[{psi -> sorted["FullInputOrder"]} /* sorted, "Label" -> None]];
    sum += t; pr[label <> "_s3_prepend_state_ms", ms[t]];

    (* s4: qco["Normal"] (TensorNetwork.m:185) *)
    t = First @ AbsoluteTiming[circN = circuit2["Normal"]];
    sum += t; pr[label <> "_s4_normal_ms", ms[t]];

    (* s5: TensorNetworkBasis incl TensorNetworkInfo (TensorNetwork.m:189) *)
    t = First @ AbsoluteTiming[circN["Width"]; basis = circN["TensorNetworkBasis"]];
    sum += t; pr[label <> "_s5_tnbasis_ms", ms[t]];

    (* s6: computational-basis test (TensorNetwork.m:191-199); info is now cached *)
    t = First @ AbsoluteTiming[
        info = circN["TensorNetworkInfo"];
        computationalQ = TrueQ[Or @@ Unequal @@@ Partition[Lookup[info["QuditBases"], Catenate[info["ContractionIndices"]]], 2]]];
    sum += t; pr[label <> "_s6_computationalq_ms", ms[t]];

    (* s7: trace/eigen/order/bend scans (TensorNetwork.m:200-203) *)
    t = First @ AbsoluteTiming[
        circN["TraceOrder"]; circN["Eigenorder"]; Sort /@ circN["Order"];
        AnyTrue[circN["Operators"], #["MatrixQ"] &]];
    sum += t; pr[label <> "_s7_trace_eigen_order_ms", ms[t]];

    (* s8: build the tensor network (TensorNetwork.m:213) *)
    t = First @ AbsoluteTiming[net = QTN[circN, "Computational" -> True, "PrependInitial" -> False]];
    sum += t; pr[label <> "_s8_build_network_ms", ms[t]];

    (* s9: contraction (TensorNetwork.m:215) *)
    t = First @ AbsoluteTiming[raw = TensorNetworkContract[net, Automatic]];
    sum += t; pr[label <> "_s9_contract_ms", ms[t]];

    (* s10: result wrap (TensorNetwork.m:216-246 + Apply unwrap :258-265) *)
    t = First @ AbsoluteTiming[
        res = QuantumState[saFlatten[raw],
            QuantumBasis[QuditBasis[basis["OutputDimensions"]], QuditBasis[basis["InputDimensions"]]]];
        res = QuantumState[res, basis];
        state = QuantumOperator[res, Sort /@ circN["Order"]]["State"]];
    sum += t; pr[label <> "_s10_wrap_ms", ms[t]];

    pr[label <> "_stage_sum_ms", ms[sum]];
    pr[label <> "_residual_ms", ms[total - sum]];
];

(* run 1: cold kernel *)
qc1 = circ[n, 0.3];
pr["gatecount", qc1["GateCount"]];
stageRun[qc1, "run1_cold"];
pr["dv_count_after_run1", dvCount[]];
pr["dv_bytes_after_run1_MB", Round[Total[dvBytes[]]/2.^20, 0.1]];

(* run 2: perturbed RZ angle -> fresh RZ ops, H/CNOT identical (semi-warm) *)
qc2 = circ[n, 0.3000001];
stageRun[qc2, "run2_perturbed"];

(* run 3: same circuit as run2 reapplied -> fully warm *)
stageRun[qc2, "run3_warm"];
pr["dv_count_after_run3", dvCount[]];
pr["dv_bytes_after_run3_MB", Round[Total[dvBytes[]]/2.^20, 0.1]];

(* warm totals: same circuit+state applied repeatedly *)
Do[pr["warm_apply_" <> ToString[k] <> "_ms", ms[First @ AbsoluteTiming[qc1[psi]]]], {k, 3}];

(* correctness anchor *)
out = qc1[psi];
pr["norm_check", N[out["Norm"]]];
pr["explicit_amplitudes", out["State"]["ExplicitLength"]];
