(* A nested map with a declared parameter, e^(sin(gamma H)) for a diagonal 10-qubit
   Ising H, on one paclet copy:
     wolframscript -file nested.wls tree
   the time of one substitution (best of three after ClearSystemCache[]) and the
   largest distance to the diagonal of e^(sin(gamma E_x)) at gamma = 0.3. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

SetAttributes[bestTime, HoldFirst];
bestTime[e_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[e;]], 3]];

isingEnergies[n_, seed_] := BlockRandom[
    With[{spins = 1 - 2 Tuples[{0, 1}, n]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {n, n}], 1]) spins, {2}]],
    RandomSeeding -> seed
];

With[{h = isingEnergies[10, 7]},
    With[{u = Exp[Sin[QuantumOperator["Diagonal"[\[Gamma] h], Range[10], "Parameters" -> {\[Gamma]}, "Label" -> None]]]},
        report["one substitution, 10 qubits", bestTime[u[0.3]["Matrix"]]];
        report["largest distance to e^(sin(0.3 E_x))", Max[Abs[Normal[u[0.3]["Matrix"]] - DiagonalMatrix[Exp[Sin[0.3 h]]]]]]
    ]
]
