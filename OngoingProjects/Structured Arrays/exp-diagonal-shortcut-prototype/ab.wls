(* The transverse-field control of section 4.5, measured alternately on the two copies:
     wolframscript -file ab.wls tree
   the smallest of seven runs of Exp[-I 0.3 sum_k X_k] on 10 qubits after
   ClearSystemCache[], and the smallest of seven runs of the strict diagonal test on
   the matrix it exponentiates. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

SetAttributes[minTime, HoldFirst];
minTime[e_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[e;]], 7]];

With[{transverse = Total[QuantumOperator["X", {#}] & /@ Range[10]]},
    report["Exp, smallest of 7", minTime[Exp[-I 0.3 transverse]]];
    With[{m = (-I 0.3 transverse)["Sort"]["Matrix"]},
        report["strict test on its matrix, smallest of 7", minTime[DiagonalMatrixQ[m, Tolerance -> 0]]]
    ]
]
