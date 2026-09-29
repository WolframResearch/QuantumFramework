(* Timing harness for the diagonal exponential:
     wolframscript -file bench.wls tree group
   tree is a paclet copy beside this script (base, proto, protos) and group one of
   "ising", "liouvillian", "qaoa", "controls", "detect". A time is the best of three
   runs after ClearSystemCache[] when one run takes under 5 s, else that one run,
   each within TimeConstrained; "TIMEOUT" marks a case that did not finish. *)
{tree, group} = Take[$ScriptCommandLine, -2];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[label_, value_] := WriteString[$Output, tree, "\t", group, "\t", label, "\t", ToString[value, InputForm], "\n"];

SetAttributes[{timeOnce, bestTime, timedValue}, HoldFirst];
timeOnce[e_, limit_] := (ClearSystemCache[]; TimeConstrained[Round[First[AbsoluteTiming[e;]], 0.0001], limit, "TIMEOUT"]);
bestTime[e_, limit_ : 300] := bestOf[timeOnce[e, limit], timeOnce[e, limit] &];
bestOf[t : "TIMEOUT", _] := t;
bestOf[t_ ? (# > 5 &), _] := t;
bestOf[t_, again_] := Min[t, again[], again[]];

(* One timed run and its value; "TIMEOUT" and a Missing value when it does not finish. *)
timedValue[e_, limit_] := (ClearSystemCache[]; TimeConstrained[MapAt[Round[#, 0.0001] &, AbsoluteTiming[e], 1], limit, {"TIMEOUT", Missing["TimedOut"]}]);

(* The Ising energies E_x = sum_(i<j) J_ij z_i z_j of every basis state, J uniform in [-1, 1]. *)
isingEnergies[n_, seed_] := BlockRandom[
    With[{spins = 1 - 2 Tuples[{0, 1}, n]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {n, n}], 1]) spins, {2}]],
    RandomSeeding -> seed
];

(* A diagonal operator on qubits 1..n with the given diagonal. QF's "Diagonal"
   constructor labels the operator with the whole list, and arithmetic on that label
   would dominate the timings at 16 qubits, so the label is dropped. *)
diagonalOperator[d_, n_, opts___] := QuantumOperator["Diagonal"[d], Range[n], opts, "Label" -> None];

(* The Pauli string with p on qubit k and the identity on the other n - 1. *)
pauliOn[n_, k_, p_] := QuantumOperator[StringReplacePart[StringRepeat["I", n], p, {k, k}]];

diagonalError[u_QuantumOperator, reference_] := Max[Abs[Normal[Diagonal[u["Matrix"]]] - reference]];
diagonalError[_, _] := Missing["TimedOut"];

run["ising"] := Scan[isingCase, {12, 13, 14, 15, 16}];

isingCase[n_] := With[{
    h = isingEnergies[n, 7],
    state = QuantumState[BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 2^n]], RandomSeeding -> 1]]
},
    With[{qo = diagonalOperator[h, n]},
        With[{exponential = timedValue[Exp[-I 0.3 qo], 400]},
            report["Exp[-I 0.3 qo] n=" <> ToString[n], First[exponential]];
            report["  error against Exp of the energies n=" <> ToString[n], diagonalError[Last[exponential], Exp[-I 0.3 h]]]
        ];
        report["MatrixExp[-I 0.3 qo, qs] n=" <> ToString[n], bestTime[MatrixExp[-I 0.3 qo, state], 400]]
    ]
];

run["liouvillian"] := Scan[liouvillianCase, If[tree === "protos", {4, 5, 6, 7}, {4, 5, 6}]];

(* Pure dephasing of every qubit at rate 0.2 under a diagonal Ising Hamiltonian. *)
liouvillianCase[n_] := With[{
    build = timedValue[
        QuantumOperator["Liouvillian"[diagonalOperator[isingEnergies[n, 7], n], pauliOn[n, #, "Z"] & /@ Range[n], ConstantArray[0.2, n]]],
        600
    ]
},
    report["build Liouvillian n=" <> ToString[n], First[build]];
    liouvillianExponential[n, Last[build]]
];

liouvillianExponential[n_, l_QuantumOperator] := With[{m = (0.7 l)["Sort"]["Matrix"]},
    report["  head of the stored state, time and head of Matrix n=" <> ToString[n], {Head[l["State"]["State"]], bestTime[(0.7 l)["Sort"]["Matrix"]], Head[m]}];
    With[{exponential = timedValue[Exp[0.7 l], 400]},
        report["Exp[0.7 L] n=" <> ToString[n], First[exponential]];
        report["  error against Exp of the diagonal n=" <> ToString[n], diagonalError[Last[exponential], Exp[Normal[Diagonal[m]]]]]
    ]
];
liouvillianExponential[_, _] := Null;

run["qaoa"] := Scan[qaoaCase, If[tree === "base", {8, 10, 12}, {8, 10, 12, 14, 16}]];

(* The cost layer e^(-i gamma H) with gamma declared, read at 25 and 100 values. *)
qaoaCase[n_] := With[{h = isingEnergies[n, 7]},
    With[{build = timedValue[Exp[diagonalOperator[-I \[Gamma] h, n, "Parameters" -> {\[Gamma]}]], 600]},
        report["build Exp with gamma declared n=" <> ToString[n], First[build]];
        qaoaScan[n, h, Last[build]]
    ]
];

qaoaScan[n_, h_, u_QuantumOperator] := With[{one = bestTime[u[0.3]["Matrix"], 600]},
    report["  ByteCount n=" <> ToString[n], ByteCount[u]];
    report["  one value n=" <> ToString[n], one];
    report["  error at 0.3 against Exp of the energies n=" <> ToString[n], diagonalError[u[0.3], Exp[-I 0.3 h]]];
    If[ NumberQ[one] && one < 2,
        report["  25 values n=" <> ToString[n], bestTime[u[#]["Matrix"] & /@ Subdivide[0.05, 1.25, 24], 600]];
        report["  100 values n=" <> ToString[n], bestTime[u[#]["Matrix"] & /@ Subdivide[0.05, 1.25, 99], 600]]
    ]
];
qaoaScan[_, _, _] := Null;

run["controls"] := With[{
    transverse = Total[QuantumOperator["X", {#}] & /@ Range[10]],
    hermitian = BlockRandom[QuantumOperator[With[{m = RandomComplex[{-1 - I, 1 + I}, {2^#, 2^#}]}, (m + ConjugateTranspose[m]) / 2], Range[#]] & /@ {8, 10}, RandomSeeding -> 3],
    nearlyDiagonal = diagonalOperator[isingEnergies[10, 7], 10] + QuantumOperator[SparseArray[{{1, 1024} -> 1.*^-12, {1024, 1} -> 1.*^-12}, {1024, 1024}], Range[10]],
    symbolicBlocks = QuantumOperator[{{a, b, 0, 0}, {b, a, 0, 0}, {0, 0, c, d}, {0, 0, d, c}}, {1, 2}]
},
    report["transverse field, 10 qubits, Exp", bestTime[Exp[-I 0.3 transverse]]];
    report["random Hermitian, 8 qubits, Exp", bestTime[Exp[-I 0.3 First[hermitian]]]];
    report["random Hermitian, 10 qubits, Exp", bestTime[Exp[-I 0.3 Last[hermitian]]]];
    report["diagonal plus a 1e-12 corner coupling, 10 qubits, Exp", bestTime[Exp[-I 0.3 nearlyDiagonal]]];
    report["symbolic 4 x 4 blocks, Exp", bestTime[Exp[symbolicBlocks]]];
    report["build RX[theta]", bestTime[QuantumOperator["RX"[\[Theta]]]]];
    report["build R[theta, XX]", bestTime[QuantumOperator["R"[\[Theta], "XX"]]]];
    report["build RZ[theta]", bestTime[QuantumOperator["RZ"[\[Theta]]]]];
    report["build R[theta, ZZ]", bestTime[QuantumOperator["R"[\[Theta], "ZZ"]]]];
    report["build R[0.3, ZZZZ]", bestTime[QuantumOperator["R"[0.3, "ZZZZ"]]]];
    report["build a 4-qubit order-2 Trotter circuit, 5 steps", bestTime[QuantumCircuitOperator["Trotterization"[
        {QuantumOperator["ZZ", {1, 2}], QuantumOperator["ZZ", {2, 3}], QuantumOperator["ZZ", {3, 4}], Total[QuantumOperator["X", {#}] & /@ Range[4]]}, 2, 5]]]]
];

run["detect"] := With[{
    dense = BlockRandom[RandomComplex[{-1 - I, 1 + I}, {1024, 1024}], RandomSeeding -> 3],
    sparseDiagonal = SparseArray[Band[{1, 1}] -> BlockRandom[RandomReal[1, 2^16], RandomSeeding -> 3]],
    transverse = Total[QuantumOperator["X", {#}] & /@ Range[12]]["Matrix"],
    symbolic = Array[Subscript[x, #1, #2] + Subscript[y, #1, #2] &, {64, 64}],
    denseDiagonal = Normal[SparseArray[Band[{1, 1}] -> BlockRandom[RandomComplex[1, 4096], RandomSeeding -> 3]]],
    hiddenZeros = DiagonalMatrix[Array[Subscript[x, #] &, 64]] + (Sin[y] ^ 2 + Cos[y] ^ 2 - 1) (ConstantArray[1, {64, 64}] - IdentityMatrix[64]),
    values = BlockRandom[RandomReal[1, 2 ^ 16], RandomSeeding -> 3]
},
    report["strict test, dense random 1024 (not diagonal)", bestTime[DiagonalMatrixQ[dense, Tolerance -> 0]]];
    report["strict test, SparseArray diagonal 2^16", bestTime[DiagonalMatrixQ[sparseDiagonal, Tolerance -> 0]]];
    report["strict test, sparse transverse field, 12 qubits (not diagonal)", bestTime[DiagonalMatrixQ[transverse, Tolerance -> 0]]];
    report["strict test, symbolic dense 64 (not diagonal)", bestTime[DiagonalMatrixQ[symbolic, Tolerance -> 0]]];
    report["strict test, dense diagonal 4096", bestTime[DiagonalMatrixQ[denseDiagonal, Tolerance -> 0]]];
    report["strict test, symbolic 64 whose off-diagonal entries are all Sin[y]^2 + Cos[y]^2 - 1 (diagonal)", {DiagonalMatrixQ[hiddenZeros, Tolerance -> 0], bestTime[DiagonalMatrixQ[hiddenZeros, Tolerance -> 0]]}];
    report["sparse diagonal of 2^16 entries: DiagonalMatrix Sparse, index pairs, Band", {
        bestTime[DiagonalMatrix[values, TargetStructure -> "Sparse"]],
        bestTime[SparseArray[Transpose[{Range[2 ^ 16], Range[2 ^ 16]}] -> values, {2 ^ 16, 2 ^ 16}]],
        bestTime[SparseArray[Band[{1, 1}] -> values, {2 ^ 16, 2 ^ 16}]]
    }]
];

run[other_] := report["unknown group", other];

run[group]
