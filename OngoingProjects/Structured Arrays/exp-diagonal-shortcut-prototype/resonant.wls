(* Section 3's resonant example on one paclet copy:
     wolframscript -file resonant.wls tree
   Levels 2 and 3 are degenerate at energy 0 with a coupling g = 10^-9 between them,
   beside a level at 10^6; over t = 10^6 the pair rotates by g t = 10^-3. The exact
   2-3 block of e^(-i t M) is cos(g t) 1 - i sin(g t) X. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

resonant = {{1.*^6, 0., 0.}, {0., 0., 1.*^-9}, {0., 1.*^-9, 0.}};
t = 1.*^6;
exactBlock = Cos[10.^-3] IdentityMatrix[2] - I Sin[10.^-3] {{0, 1}, {1, 0}};
blockError[u_] := Max[Abs[Normal[u][[2 ;;, 2 ;;]] - exactBlock]];

report["default-tolerance and strict DiagonalMatrixQ", {DiagonalMatrixQ[resonant], DiagonalMatrixQ[resonant, Tolerance -> 0]}];
report["2-3 block error, MatrixExp of the full matrix", blockError[MatrixExp[-I t resonant]]];
report["2-3 block error, entry by entry", blockError[DiagonalMatrix[Exp[-I t Diagonal[resonant]]]]];
report["2-3 block error, QF Exp[-I t qo]", blockError[Exp[-I t QuantumOperator[resonant, 3]]["Matrix"]]];
report["2-3 block error, MatrixExp of the block alone", blockError[ArrayFlatten[{{1, 0}, {0, MatrixExp[-I t resonant[[2 ;;, 2 ;;]]]}}]]];
