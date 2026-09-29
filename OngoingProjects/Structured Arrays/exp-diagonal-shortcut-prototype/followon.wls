(* Follow-ons on a diagonal operator, measured on one paclet copy:
     wolframscript -file followon.wls tree
   non-integer powers, the spectrum, and QuantumEvolve, each against its closed form.
   A time is the best of three runs after ClearSystemCache[] when one run takes
   under 5 s, else that one run; the first call of a property is timed on a fresh
   operator, since QF caches property values. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

SetAttributes[{timeOnce, bestTime}, HoldFirst];
timeOnce[e_, limit_] := (ClearSystemCache[]; TimeConstrained[Round[First[AbsoluteTiming[e;]], 0.0001], limit, "TIMEOUT"]);
bestTime[e_, limit_ : 300] := bestOf[timeOnce[e, limit], timeOnce[e, limit] &];
bestOf[t : "TIMEOUT", _] := t;
bestOf[t_ ? (# > 5 &), _] := t;
bestOf[t_, again_] := Min[t, again[], again[]];

isingEnergies[n_, seed_] := BlockRandom[
    With[{spins = 1 - 2 Tuples[{0, 1}, n]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {n, n}], 1]) spins, {2}]],
    RandomSeeding -> seed
];
diagonalOperator[d_, n_] := QuantumOperator["Diagonal"[d], Range[n], "Label" -> None];

(* A non-integer power goes through MatrixPower; on a diagonal it is the principal
   branch entry by entry, negative entries included. *)
powerCase[n_] := With[{h = isingEnergies[n, 7]},
    With[{qo = diagonalOperator[h, n]},
        report["qo^0.5 n=" <> ToString[n], bestTime[qo ^ 0.5]];
        report["  qo^0.5 against h^0.5 entry by entry n=" <> ToString[n], Max[Abs[Normal[Diagonal[(qo ^ 0.5)["Matrix"]]] - h ^ 0.5]]];
        report["  h^0.5 as a sparse diagonal n=" <> ToString[n], bestTime[DiagonalMatrix[h ^ 0.5, TargetStructure -> "Sparse"]]];
        report["qo^3 n=" <> ToString[n], bestTime[qo ^ 3]]
    ]
];
Scan[powerCase, {10, 12}];

(* The spectrum: the first call on a fresh operator, the order Eigensystem gives, and
   the structured-array route. *)
spectrumCase[n_] := With[{h = isingEnergies[n, 7]},
    report["first-call Eigenvalues n=" <> ToString[n], timeOnce[diagonalOperator[h, n]["Eigenvalues"], 600]];
    report["first-call Eigensystem n=" <> ToString[n], timeOnce[diagonalOperator[h, n]["Eigensystem"], 600]];
    With[{values = diagonalOperator[h, n]["Eigenvalues"]},
        report["  order is decreasing absolute value n=" <> ToString[n], values === h[[Reverse[OrderingBy[h, Abs]]]]];
        report["  Eigenvalues of the structured diagonal n=" <> ToString[n], bestTime[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]]];
        report["  same values n=" <> ToString[n], Sort[Eigenvalues[DiagonalMatrix[h, TargetStructure -> "Structured"]]] == Sort[values]]
    ]
];
Scan[spectrumCase, {8, 10, 12}];
report["Eigensystem of a structured DiagonalMatrix", Eigensystem[DiagonalMatrix[{1., 3., 4.}, TargetStructure -> "Structured"]]];
report["Eigensystem of an exact degenerate diagonal", Eigensystem[DiagonalMatrix[{2, 1, 2}]]];

(* QuantumEvolve with a time-independent diagonal H, against e^(-i h t) psi(0). *)
evolveCase[n_] := With[{h = isingEnergies[n, 7], v = BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]], RandomSeeding -> 2]},
    With[{psi = QuantumEvolve[diagonalOperator[h, n], QuantumState[v], {\[FormalT], 0, 1}]},
        report["QuantumEvolve over t in [0, 1] n=" <> ToString[n], bestTime[QuantumEvolve[diagonalOperator[h, n], QuantumState[v], {\[FormalT], 0, 1}], 600]];
        report["  error at t = 1 against Exp[-I h] psi(0) n=" <> ToString[n], Max[Abs[Normal[psi[1]["StateVector"]] - Exp[-I h] v]]]
    ]
];
Scan[evolveCase, {6, 8, 10}];
report["QuantumEvolve with symbolic t, 2 qubits", Normal[QuantumEvolve[diagonalOperator[{1, -1, 2, 0}, 2], QuantumState["++"], \[FormalT]]["StateVector"]]];

(* f[qo] of a diagonal operator: the diagonal branch builds its result with Band and
   an exact zero background, so Normal of a machine result is not packed. *)
With[{qo = diagonalOperator[isingEnergies[16, 7], 16]},
    report["Sin of a 16-qubit diagonal: time, background, Normal packed", {bestTime[Sin[qo]], Sin[qo]["Matrix"]["Background"], Developer`PackedArrayQ[Normal[Sin[qo]["Matrix"]]]}]
];
