(* Verify three sharp claims:
   1. Memoize (Utilities.m:382) appends memo entries AFTER its catch-all trampoline,
      so they are never matched and DownValues grow on every call (incl. duplicates).
   2. The property-dispatch cache re-executes Set on every HIT (write-on-hit waste).
   3. qco["InputDimensions"] (the apply dim guard) goes through ResourceFunction["LookupPart"]. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

fos = Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand;
n0 = Length[DownValues[fos]];
fos["H" -> 7];
n1 = Length[DownValues[fos]];
fos["H" -> 7];
n2 = Length[DownValues[fos]];
fos["H" -> 7];
n3 = Length[DownValues[fos]];
pr["memoize_dv_growth_same_arg", {n0, n1, n2, n3}];
pr["memoize_first_rule_is_catchall", MatchQ[DownValues[fos][[1, 1]], Verbatim[HoldPattern][fos[___]] | HoldPattern[HoldPattern[fos[Pattern[_, ___]]]]]];
pr["memoize_first_rule_lhs", DownValues[fos][[1, 1]] /. s_Symbol /; Context[s] =!= "System`" :> ToString[Unevaluated[s]]];
t = First[AbsoluteTiming[Do[fos["H" -> 7], {100}]]]/100;
pr["memoize_samecall_us", Round[10.^6 t, 1.]];
pr["memoize_dv_after_100", Length[DownValues[fos]]];

(* 2. write-on-hit: count Set executions via a hook on the Prop symbol's memo *)
qop = Symbol[First @ Names["*`QuantumOperatorProp"]];
h = QuantumOperator["H", {5}];
h["Order"]; (* populate memo *)
c0 = Length[DownValues[qop]];
Do[h["Order"], {10}];
c1 = Length[DownValues[qop]];
pr["prop_dv_stable_on_hits", {c0, c1}];
(* time the Quiet[Set] alone on a representative memo rule *)
t = First[AbsoluteTiming[Do[Quiet[qop[h, "Order"] = {{5}, {5}}, Rule::rhs], {1000}]]]/1000;
pr["quiet_set_us", Round[10.^6 t, 0.1]];
t = First[AbsoluteTiming[Do[qop[h, "Order"], {1000}]]]/1000;
pr["bare_prop_lookup_us", Round[10.^6 t, 0.1]];

(* 3. InputDimensions cost on the standard circuit *)
layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
qc = QuantumCircuitOperator[Flatten[Table[layer[12, 0.3], {3}], 1]];
t = First @ AbsoluteTiming[qc["InputDimensions"]];
pr["inputdims_first_ms", Round[1000. t, 0.1]];
t = First @ AbsoluteTiming[qc["InputDimensions"]];
pr["inputdims_warm_ms", Round[1000. t, 0.1]];

(* 4. what H["Normal"] returns *)
res = Quiet[QuantumOperator["H", {1}]["Normal"]];
pr["h_normal_result", ToString[Head[res]] <> "/" <> ToString[res === QuantumOperator["H", {1}]["Normal"] // Quiet]];
