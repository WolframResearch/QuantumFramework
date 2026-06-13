(* Call-count instrumentation of one circuit apply.
   A counting rule with an always-False condition is prepended to each internal Prop
   dispatch table: the matcher tries it first, the condition bumps a tally and returns
   False, and evaluation falls through to the real rules. With the property cache ON,
   memoized values are stored as more-specific (literal) rules that match before the
   counter, so the counter sees cache MISSES (actual computes) only. With the cache OFF,
   it sees every property call. The difference is the cache hit count. *)

PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];

pr[l_, v_] := Print["RESULT|", l, "|", v];
resolve[name_] := With[{hits = Names["*`" <> name]}, If[hits === {}, $Failed, Symbol[First[hits]]]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {3}], 1]];

n = 12;
psi = QuantumState["Register"[n]];

$counts = <||>;
bump[k_] := ($counts[k] = Lookup[$counts, Key[k], 0] + 1;);
pName[p_] := Replace[p, {s_String :> s, {s_String, ___} :> s, _ :> "OTHER"}];

addPropCounter[name_String] := With[{sym = resolve[name], tag = name},
    If[sym === $Failed, pr["MISSING_SYMBOL", name],
        PrependTo[DownValues[sym], HoldPattern[sym[_, p_, ___]] :> Null /; (bump[{tag, pName[p]}]; False)]]
];

heads = {"QuantumOperatorProp", "QuantumStateProp", "QuantumBasisProp", "QuantumCircuitOperatorProp", "QuditBasisProp"};
Scan[addPropCounter, heads];

(* constructor touch counters (counts every evaluator pass over the head, incl. valid objects) *)
addCtorCounter[sym_Symbol, tag_String] := PrependTo[DownValues[sym], HoldPattern[sym[___]] :> Null /; (bump[{"CTOR", tag}]; False)];
addCtorCounter[QuantumOperator, "QuantumOperator"];
addCtorCounter[QuantumState, "QuantumState"];
addCtorCounter[QuantumBasis, "QuantumBasis"];
addCtorCounter[QuditBasis, "QuditBasis"];
addCtorCounter[QuantumCircuitOperator, "QuantumCircuitOperator"];

report[label_] := Block[{byHead, top},
    byHead = GroupBy[Normal[$counts], #[[1, 1]] &, Total @* Map[Last]];
    pr[label <> "_totals_by_head", byHead];
    top = TakeLargest[$counts, UpTo[40]];
    KeyValueMap[pr[label <> "_count_" <> #1[[1]] <> "_" <> #1[[2]], #2] &, top];
];

(* run A: cache ON, fresh circuit -> counts are computes (cache misses) *)
qcA = circ[n, 0.31];
$counts = <||>;
tA = First @ AbsoluteTiming[qcA[psi]];
report["A_cacheON"];
pr["A_time_ms", Round[1000. tA, 0.1]];

(* run B: cache OFF, fresh circuit -> counts are ALL property calls *)
qcB = circ[n, 0.32];
$counts = <||>;
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
tB = First @ AbsoluteTiming[qcB[psi]];
Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = True;
report["B_cacheOFF"];
pr["B_time_ms", Round[1000. tB, 0.1]];
