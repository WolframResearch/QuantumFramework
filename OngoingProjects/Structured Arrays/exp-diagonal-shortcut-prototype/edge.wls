(* Edge cases quoted in the report's risks, on one paclet copy:
     wolframscript -file edge.wls tree *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

report["MatrixExp of diag(-720., 0.)", Normal[MatrixExp[DiagonalMatrix[{-720., 0.}]]]];
report["Exp[I t PauliZ[3]], its diagonal and the difference from Exp of the clock phases", With[{d = Normal[Diagonal[Exp[I t QuantumOperator["PauliZ"[3]]]["Matrix"]]]},
    {d, Simplify[d - Exp[I t Normal[Diagonal[QuantumOperator["PauliZ"[3]]["Matrix"]]]]]}
]];
report["Exp[-2000. diag(1., 1.)] entries", Normal[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 1.}]]]["Matrix"]]];
report["Exp with gamma declared at gamma = 0.", Normal[Exp[QuantumOperator["Diagonal"[-I \[Gamma] {0.5, -0.25}], {1}, "Parameters" -> {\[Gamma]}]][0.]["Matrix"]]];
report["e^800: relative error at 40 digits, precision", With[{m = Normal[Exp[-800. QuantumOperator[DiagonalMatrix[{-1., 1.}]]]["Matrix"]]}, {N[Abs[SetPrecision[m[[1, 1]], 40] / Exp[800] - 1]], Precision[m[[1, 1]]]}]];
report["diag(1/a, b) declared, at a = 0: head", Head[Exp[QuantumOperator[DiagonalMatrix[{1 / \[Alpha], \[Beta]}], "Parameters" -> {\[Alpha], \[Beta]}]][0, 1]]];
