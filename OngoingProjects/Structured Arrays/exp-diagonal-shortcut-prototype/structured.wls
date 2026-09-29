(* Two WL 15.0.1 results on structured diagonal matrices, independent of QF:
     wolframscript -file structured.wls *)
report[label_, value_] := WriteString[$Output, $VersionNumber, "\t", label, "\t", ToString[value, InputForm], "\n"];
report["MatrixExp of a structured diagonal with an underflowing complex entry", MatrixExp[DiagonalMatrix[{-1000. + 3. I, 1.}, TargetStructure -> "Structured"]]];
report["MatrixExp of the same diagonal as a SparseArray", Normal[MatrixExp[DiagonalMatrix[{-1000. + 3. I, 1.}, TargetStructure -> "Sparse"]]]];
report["Eigensystem of a structured diagonal", Eigensystem[DiagonalMatrix[{1., 3., 4.}, TargetStructure -> "Structured"]]];
