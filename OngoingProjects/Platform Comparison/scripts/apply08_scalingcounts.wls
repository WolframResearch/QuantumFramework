(* Does the number of circuit-level O(G) property computes grow with G?
   Counts NormalOrders / Flatten / Operators computes per apply at two depths, cache OFF. *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];
resolve[name_] := Symbol[First[Names["*`" <> name]]];

layer[n_, a_] := Join[Table["H" -> q, {q, n}], Table["CNOT" -> {q, q + 1}, {q, n - 1}], Table["RZ"[a] -> q, {q, n}]];
circ[n_, a_, L_] := QuantumCircuitOperator[Flatten[Table[layer[n, a], {L}], 1]];
psi = QuantumState["Register"[12]];

$counts = <||>;
bump[k_] := ($counts[k] = Lookup[$counts, Key[k], 0] + 1;);
pName[p_] := Replace[p, {s_String :> s, {s_String, ___} :> s, _ :> "OTHER"}];
With[{qcop = resolve["QuantumCircuitOperatorProp"]},
    PrependTo[DownValues[qcop], HoldPattern[qcop[_, p_, ___]] :> Null /; (bump[pName[p]]; False)]];

Wolfram`QuantumFramework`PackageScope`$QuantumFrameworkPropCache = False;
Do[
    Module[{qcL = circ[12, 0.3 + L 10.^-6, L], t},
        $counts = <||>;
        t = First @ AbsoluteTiming[qcL[psi]];
        pr["L" <> ToString[L] <> "_ms", Round[1000. t, 1.]];
        pr["L" <> ToString[L] <> "_counts", Normal @ TakeLargest[$counts, UpTo[12]]];
    ],
    {L, {1, 4, 16}}
]
