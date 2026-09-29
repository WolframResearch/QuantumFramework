(* Accuracy of the exponential of an 8-qubit Ising Hamiltonian, and two costs around it:
     wolframscript -file accuracy.wls tree
   with tree a paclet copy beside this script (base, proto). Unitarity is measured as
   the largest deviation of |U_kk| from 1 in units of $MachineEpsilon; the reference
   error is the largest distance to the 40-digit exponential of the exact product of
   theta and each stored entry, in units of (|theta h_k|/2 + 1) $MachineEpsilon, one
   rounding of that product plus one unit of the result. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", label, "\t", ToString[value, InputForm], "\n"];

SetAttributes[bestTime, HoldFirst];
bestTime[e_] := Min[Table[ClearSystemCache[]; First[AbsoluteTiming[e;]], 3]];

(* sum_(i<j) J_ij Z_i Z_j from QF's two-qubit ZZ operators, J uniform in [-1, 1]. *)
zzOperator[n_, seed_] := With[{couplings = BlockRandom[RandomReal[{-1, 1}, {n, n}], RandomSeeding -> seed]},
    Total[Apply[couplings[[#1, #2]] QuantumOperator["ZZ", {#1, #2}] &, Subsets[Range[n], {2}], {1}]]
];

hamiltonian = zzOperator[8, 7];
stored = Normal[Diagonal[hamiltonian["Sort"]["Matrix"]]];
diagonalAt[theta_] := Normal[Diagonal[Exp[-I theta hamiltonian]["Matrix"]]];

unitarityUlps[theta_] := Max[Abs[Abs[diagonalAt[theta]] - 1]] / $MachineEpsilon;
referenceError[theta_] := With[{phases = SetPrecision[theta, 40] SetPrecision[stored, 40]},
    Max[Abs[SetPrecision[diagonalAt[theta], 40] - Exp[-I phases]] / ((Abs[phases] / 2 + 1) $MachineEpsilon)]
];

report["unitarity, ulps, theta = 0.3, 1000, 10^6", N[unitarityUlps /@ {0.3, 1000., 1.*^6}, 4]];
report["reference error, theta = 0.3, 7.1, 1000, 10^6", N[referenceError /@ {0.3, 7.1, 1000., 1.*^6}, 4]];

(* QF's "Diagonal" constructor labels the operator with its whole diagonal; at 16
   qubits the label arithmetic costs more than the exponential on the prototype. On
   the base copy the exponential itself takes minutes there, so it is skipped. *)
labelCost["proto"] := With[
    {energies = BlockRandom[With[{spins = 1 - 2 Tuples[{0, 1}, 16]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {16, 16}], 1]) spins, {2}]], RandomSeeding -> 7]},
    With[{labelled = QuantumOperator["Diagonal"[energies], Range[16]], unlabelled = QuantumOperator["Diagonal"[energies], Range[16], "Label" -> None]},
        report["16 qubits, Exp time with the list label, without a label", {bestTime[Exp[-I 0.3 labelled]], bestTime[Exp[-I 0.3 unlabelled]]}]
    ]
];
labelCost[_] := Null;
labelCost[tree];

(* The QAOA mixer e^(-i beta sum_k X_k) is a product of single-qubit rotations, so it
   is the Kronecker product of n 2 x 2 exponentials; QF exponentiates the 2^n-dimensional
   sum with MatrixExp. *)
With[{n = 9, beta = 0.4},
    With[{mixer = Total[QuantumOperator["X", {#}] & /@ Range[n]], single = MatrixExp[-I beta {{0, 1}, {1, 0}}]},
        report["mixer, 9 qubits: QF Exp time, Kronecker product time", {bestTime[Exp[-I beta mixer]], bestTime[KroneckerProduct @@ ConstantArray[single, n]]}];
        report["mixer, 9 qubits: largest difference", Max[Abs[Normal[Exp[-I beta mixer]["Matrix"]] - KroneckerProduct @@ ConstantArray[single, n]]]]
    ]
];
