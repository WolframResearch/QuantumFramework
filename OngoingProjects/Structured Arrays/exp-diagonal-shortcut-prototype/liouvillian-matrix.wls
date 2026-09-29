(* The matrix the exponential sees for the n = 6 pure-dephasing Liouvillian:
     wolframscript -file liouvillian-matrix.wls tree
   its head, whether a List is packed, and the time of the strict diagonal test on it
   (best of three after ClearSystemCache[]). *)
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
pauliOn[n_, k_, p_] := QuantumOperator[StringReplacePart[StringRepeat["I", n], p, {k, k}]];

With[{n = 6},
    With[{m = (0.7 QuantumOperator["Liouvillian"[
            QuantumOperator["Diagonal"[isingEnergies[n, 7]], Range[n], "Label" -> None],
            pauliOn[n, #, "Z"] & /@ Range[n], ConstantArray[0.2, n]]])["Sort"]["Matrix"]},
        report["head, dimensions, packed", {Head[m], Dimensions[m], Developer`PackedArrayQ[m]}];
        report["strict test and its time", {DiagonalMatrixQ[m, Tolerance -> 0], bestTime[DiagonalMatrixQ[m, Tolerance -> 0]]}]
    ]
]
