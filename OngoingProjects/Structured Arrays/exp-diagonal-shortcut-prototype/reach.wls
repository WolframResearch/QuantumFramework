(* What section 5.4 reaches besides the exponential, on one paclet copy:
     wolframscript -file reach.wls tree
   the head of "Matrix" and "StateMatrix" for sparse mixed states, and the exponential
   of an amplitude-damping Liouvillian against MatrixExp of its dense matrix. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

With[{rho = QuantumState[SparseArray[{{1, 1} -> 1/2, {1, 4} -> 1/2, {4, 1} -> 1/2, {4, 4} -> 1/2}, {4, 4}], QuantumBasis[{2, 2}]]},
    report["Bell density matrix: heads of Matrix, StateMatrix", Head /@ {rho["Matrix"], rho["StateMatrix"]}];
    report["its partial trace over qudit 1: head of Matrix, the matrix", {Head[QuantumPartialTrace[rho, {1}]["Matrix"]], Normal[QuantumPartialTrace[rho, {1}]["Matrix"]]}]
];
With[{l = QuantumOperator["Liouvillian"[QuantumOperator["Z"], {QuantumOperator[{{0, 1}, {0, 0}}]}, {0.3}]]},
    report["amplitude damping: head of Matrix", Head[l["Matrix"]]];
    report["amplitude damping: Exp[0.7 L] against MatrixExp of the dense matrix", Max[Abs[Normal[Exp[0.7 l]["Matrix"]] - MatrixExp[Normal[(0.7 l)["Sort"]["Matrix"]]]]]]
];
