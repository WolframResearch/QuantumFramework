(* Correct DownValues inspection of the Memoize'd FromOperatorShorthand (DownValues holds its arg). *)
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
pr[l_, v_] := Print["RESULT|", l, "|", v];

With[{f = Wolfram`QuantumFramework`PackageScope`FromOperatorShorthand},
    pr["dv0", Length[DownValues[f]]];
    f["H" -> 7]; pr["dv1", Length[DownValues[f]]];
    f["H" -> 7]; pr["dv2", Length[DownValues[f]]];
    f["H" -> 7]; pr["dv3", Length[DownValues[f]]];
    pr["first_lhs_shape", DownValues[f][[1, 1]] // Hold // ToString];
    pr["last_lhs_shape", DownValues[f][[-1, 1]] // Hold // ToString];
    t = First[AbsoluteTiming[Do[f["H" -> 7], {100}]]]/100;
    pr["repeat_call_us", Round[10.^6 t, 1.]];
    pr["dv_after_100_more", Length[DownValues[f]]];
    (* does the memo entry ever match? compare timing against the raw un-memoized path *)
]
