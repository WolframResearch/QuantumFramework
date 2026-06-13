(* Validate the cache-off fast apply: correct state, and decompose 226 ms into contraction + rest. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`TensorNetworks`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
ms[t_] := Round[1000. t, 0.1];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
qc = QuantumCircuitOperator[Flatten[Table[layer[12, 0.3], {3}], 1]];
psi = QuantumState["Register"[12]];

(* warm up with cache ON *)
ref = qc[psi];
refV = Normal[N[ref["StateVector"]]];
t = First @ AbsoluteTiming[qc[psi]]; pr["warm_ON_ms", ms[t]];

Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
t = First @ AbsoluteTiming[outOFF = qc[psi]];
pr["warm_OFF_ms", ms[t]];
pr["OFF_is_state", QuantumStateQ[outOFF]];
pr["OFF_maxdiff", Max[Abs[Normal[N[outOFF["StateVector"]]] - refV]]];
t = First @ AbsoluteTiming[qc[psi]]; pr["warm_OFF_repeat_ms", ms[t]];
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = True;

(* contraction share inside the OFF apply: time the contract alone on the same net *)
resolveTNC = Symbol[First @ Names["*`TensorNetworkCompile"]];
circuit = QuantumCircuitOperator[{psi -> qc["Flatten"]["Sort"]["FullInputOrder"]} /* qc["Flatten"]["Sort"], "Label" -> None];
net = resolveTNC[circuit, "ReturnTensorNetwork" -> True];
t = First @ AbsoluteTiming[TensorNetworkContract[net, Automatic]];
pr["contract_same_net_ms", ms[t]];
t = First @ AbsoluteTiming[TensorNetworkContract[net, Automatic]];
pr["contract_same_net_repeat_ms", ms[t]];
