(* Correctness battery for the diagonal exponential, run on two paclet copies and diffed:
     wolframscript -file correct.wls tree
   with tree a copy beside this script (base, proto). Lines tagged SAME print a hash
   of the binary serialization of the full result, so that equal lines mean results
   identical to the last bit in the two copies (Hash alone gives machine numbers that
   differ in their last bits the same value); ERR lines print an error against an
   independent reference; VAL lines print a value; MSG lines print how many messages
   a case showed the user. *)
tree = Last[$ScriptCommandLine];
PacletDirectoryLoad[FileNameJoin[{DirectoryName[$InputFileName], tree, "QuantumFramework"}]];
Needs["Wolfram`QuantumFramework`"];

report[tag_, label_, value_] := WriteString[$Output, tag, "\t", label, "\t", ToString[value, InputForm], "\n"];
SetAttributes[{same, messageCount}, HoldRest];
same[label_, e_] := report["SAME", label, Hash[BinarySerialize[Developer`FromPackedArray[Normal[e]]]]];
messageCount[label_, e_] := report["MSG", label, Length[Block[{$MessageList = {}}, e; $MessageList]]];
matrixOf[qo_] := Normal[qo["Matrix"]];
maxError[a_, b_] := Max[Abs[N[Normal[a] - Normal[b]]]];

(* The Ising energies E_x = sum_(i<j) J_ij z_i z_j of every basis state, J uniform in [-1, 1]. *)
isingEnergies[n_, seed_] := BlockRandom[
    With[{spins = 1 - 2 Tuples[{0, 1}, n]}, Total[(spins . UpperTriangularize[RandomReal[{-1, 1}, {n, n}], 1]) spins, {2}]],
    RandomSeeding -> seed
];
pauliOn[n_, k_, p_] := QuantumOperator[StringReplacePart[StringRepeat["I", n], p, {k, k}]];

h8 = isingEnergies[8, 11];
ising8 = QuantumOperator["Diagonal"[h8], Range[8]];
u8 = Exp[-I 0.3 ising8];

(* 1. a numeric diagonal: the spellings agree, against MatrixExp and Exp of the entries *)
report["ERR", "Ising8 Exp vs dense MatrixExp", maxError[matrixOf[u8], MatrixExp[-I 0.3 DiagonalMatrix[h8]]]];
report["ERR", "Ising8 Exp vs Exp of the entries", maxError[Diagonal[matrixOf[u8]], Exp[-I 0.3 h8]]];
report["ERR", "Ising8 U U^dagger - 1", maxError[matrixOf[u8] . ConjugateTranspose[matrixOf[u8]], IdentityMatrix[256]]];
report["ERR", "Ising8 E^ and MatrixExp vs Exp", maxError[matrixOf[#], matrixOf[u8]] & /@ {E ^ (-I 0.3 ising8), MatrixExp[-I 0.3 ising8]}];
report["ERR", "Ising8 2^qo vs MatrixExp of Log[2] qo", maxError[matrixOf[2 ^ (-I 0.3 ising8)], MatrixExp[Log[2] (-I 0.3) DiagonalMatrix[h8]]]];
report["VAL", "Ising8 result matrix head", Head[u8["Matrix"]]];

(* 2. exact and symbolic diagonals *)
With[{exact = QuantumOperator[DiagonalMatrix[{I Pi, Log[2], 0, -I Pi / 2}], {1, 2}], symbolic = QuantumOperator[DiagonalMatrix[{a, b, c, d}], {1, 2}]},
    same["exact diag Exp", Exp[exact]["Matrix"]];
    report["VAL", "exact diag Exp", matrixOf[Exp[exact]]];
    same["symbolic diag Exp", Exp[symbolic]["Matrix"]];
    report["VAL", "symbolic diag Exp", matrixOf[Exp[symbolic]]];
    same["symbolic diag x^qo", (x ^ symbolic)["Matrix"]];
    same["symbolic MatrixExp[qo]", MatrixExp[symbolic]["Matrix"]];
    same["symbolic Exp label", Exp[symbolic]["Label"]]
];

(* 3. named operators built through the exponential *)
Scan[same[ToString[#, InputForm], QuantumOperator[#]["Matrix"]] &,
    {"RZ"[\[Theta]], "RZ"[0.4], "RZ"[\[Theta], 3], "R"[\[Theta], "ZZ"], "R"[\[Theta], "ZZZ"], "R"[0.7, "ZIZ"], "RX"[\[Theta]], "R"[\[Theta], "XX"]}];
same["RZ label", QuantumOperator["RZ"[\[Theta]]]["Label"]];

(* 4. declared parameters: substitutions and collisions *)
With[{u = Exp[QuantumOperator["Diagonal"[-I \[Gamma] h8], Range[8], "Parameters" -> {\[Gamma]}]]},
    report["VAL", "parametrized Exp amplitudes head", Head[u["State"]["State"]]];
    report["ERR", "parametrized Exp at 0.3 vs Exp[-I 0.3 qo]", maxError[matrixOf[u[0.3]], matrixOf[u8]]];
    report["ERR", "parametrized Exp at 0 vs identity", maxError[matrixOf[u[0]], IdentityMatrix[256]]];
    same["parametrized Exp closed form", u["Matrix"]]
];
With[{collision = QuantumOperator[DiagonalMatrix[{\[Alpha], \[Beta]}], "Parameters" -> {\[Alpha], \[Beta]}]},
    report["VAL", "collision Exp at a = b = 1", matrixOf[Exp[collision][1, 1]]];
    report["VAL", "collision Sin at a = b = 0", matrixOf[Sin[collision][0, 0]]];
    report["VAL", "collision Exp[Sin] at a = b = 0", matrixOf[Exp[Sin[collision]][0, 0]]];
    report["VAL", "collision MatrixExp at a = b = 2", matrixOf[MatrixExp[collision][2, 2]]]
];
With[{offDiagonal = QuantumOperator[{{\[Alpha], \[Beta]}, {\[Beta], -\[Alpha]}}, "Parameters" -> {\[Alpha], \[Beta]}]},
    report["VAL", "non-diagonal parametrized Exp at 0, 0", matrixOf[Exp[offDiagonal][0, 0]]];
    report["ERR", "non-diagonal parametrized Exp at 1/3, 1/5", maxError[matrixOf[Exp[offDiagonal][1/3, 1/5]], MatrixExp[{{1/3, 1/5}, {1/5, -1/3}}]]]
];
With[{u = Exp[QuantumOperator[{{\[Alpha], \[Beta] \[Alpha]}, {0, 1}}, "Parameters" -> {\[Alpha], \[Beta]}]]},
    report["ERR", "coupling vanishing at b = 0: value at b = 0 vs Exp of the entries", maxError[matrixOf[u[0.5, 0.]], DiagonalMatrix[Exp[{0.5, 1.}]]]];
    report["ERR", "coupling vanishing at b = 0: value at b = 0.2 vs MatrixExp", maxError[matrixOf[u[0.5, 0.2]], MatrixExp[{{0.5, 0.1}, {0, 1}}]]]
];

(* 5. action on states *)
With[{v = BlockRandom[Normalize[RandomComplex[{-1 - I, 1 + I}, 256]], RandomSeeding -> 3]},
    report["ERR", "MatrixExp[qo, ket] vs MatrixExp of the matrix on the vector", maxError[MatrixExp[-I 0.3 ising8, QuantumState[v]]["StateVector"], MatrixExp[-I 0.3 DiagonalMatrix[h8], v]]]
];
With[{rho = BlockRandom[QuantumState["RandomMixed"[2]], RandomSeeding -> 4]},
    report["ERR", "MatrixExp[-I 0.3 XI, rho] vs Exp[-I 0.3 XI][rho]", maxError[MatrixExp[-I 0.3 QuantumOperator["XI"], rho]["DensityMatrix"], Exp[-I 0.3 QuantumOperator["XI"]][rho]["DensityMatrix"]]];
    report["ERR", "MatrixExp[-I 0.3 ZZ, rho] vs Exp[-I 0.3 ZZ][rho]", maxError[MatrixExp[-I 0.3 QuantumOperator["ZZ"], rho]["DensityMatrix"], Exp[-I 0.3 QuantumOperator["ZZ"]][rho]["DensityMatrix"]]];
    report["VAL", "trace of MatrixExp[-I 0.3 XI, rho]", Tr[MatrixExp[-I 0.3 QuantumOperator["XI"], rho]["DensityMatrix"]]]
];

(* 6. a matrix-type operator: pure dephasing of three qubits *)
With[{h3 = isingEnergies[3, 5]},
    With[{l3 = QuantumOperator["Liouvillian"[QuantumOperator["Diagonal"[h3], Range[3]], pauliOn[3, #, "Z"] & /@ Range[3], ConstantArray[0.2, 3]]],
          rho0 = BlockRandom[QuantumState["RandomMixed"[3]], RandomSeeding -> 3],
          bits = Tuples[{0, 1}, 3]},
        With[{rhoT = Normal[Exp[0.7 l3][rho0]["DensityMatrix"]], r0 = Normal[rho0["DensityMatrix"]]},
            report["ERR", "Liouvillian Exp vs dense MatrixExp", maxError[matrixOf[Exp[0.7 l3]], MatrixExp[0.7 Normal[l3["Matrix"]]]]];
            report["ERR", "Liouvillian trace change", Abs[Tr[rhoT] - Tr[r0]]];
            report["ERR", "Liouvillian populations change", maxError[Diagonal[rhoT], Diagonal[r0]]];
            report["ERR", "Liouvillian vs closed form, all entries", maxError[rhoT, r0 Exp[-0.7 I Outer[Subtract, h3, h3] - 2 0.2 0.7 Outer[HammingDistance, bits, bits, 1]]]];
            report["ERR", "Liouvillian MatrixExp[qo, qs] vs Exp[qo][qs]", maxError[MatrixExp[0.7 l3, rho0]["DensityMatrix"], rhoT]]
        ]
    ]
];

(* 7. non-diagonal controls: the results must be identical in the two copies *)
same["transverse field Exp", Exp[-I 0.3 Total[QuantumOperator["X", {#}] & /@ Range[4]]]["Matrix"]];
same["random Hermitian Exp", Exp[-I 0.3 QuantumOperator[BlockRandom[With[{m = RandomComplex[{-1 - I, 1 + I}, {16, 16}]}, (m + ConjugateTranspose[m]) / 2], RandomSeeding -> 9], Range[4]]]["Matrix"]];
same["exact non-diagonal Exp", Exp[QuantumOperator[{{0, 1}, {1, 0}}]]["Matrix"]];
same["symbolic non-diagonal Exp", Exp[QuantumOperator[{{a, b}, {b, a}}]]["Matrix"]];
same["coupling 1e-12 Exp", Exp[-I QuantumOperator[{{1., 1.*^-12}, {1.*^-12, 2.}}]]["Matrix"]];
same["coupling 1e-15 Exp", Exp[-I QuantumOperator[{{1., 1.*^-15}, {1.*^-15, 2.}}]]["Matrix"]];

(* 8. a resonant coupling at roundoff beside a large entry *)
With[{resonant = {{1.*^6, 0., 0.}, {0., 0., 1.*^-9}, {0., 1.*^-9, 0.}}, t = 1.*^6},
    report["VAL", "default-tolerance DiagonalMatrixQ calls it diagonal", DiagonalMatrixQ[resonant]];
    report["ERR", "resonant Exp, 2-3 block vs MatrixExp of the block alone", maxError[matrixOf[Exp[-I t QuantumOperator[resonant, 3]]][[2 ;;, 2 ;;]], MatrixExp[-I t resonant[[2 ;;, 2 ;;]]]]];
    report["ERR", "resonant Cos[qo], 2-3 block vs Cos of the block alone", maxError[matrixOf[Cos[t QuantumOperator[resonant, 3]]][[2 ;;, 2 ;;]], {{Cos[10.^-3], 0}, {0, Cos[10.^-3]}}]]
];

(* 9. edge cases *)
messageCount["underflow Exp[-2000 qo]", report["VAL", "underflow Exp[-2000 qo]", matrixOf[Exp[-2000. QuantumOperator[DiagonalMatrix[{1., 0.}]]]]]];
messageCount["complex underflow", report["VAL", "complex underflow", matrixOf[Exp[QuantumOperator[DiagonalMatrix[{-2000. + 3. I, 1. I}]]]]]];
report["ERR", "e^0 beside e^-800 vs 1", Abs[matrixOf[Exp[QuantumOperator[DiagonalMatrix[{-800., 0}]]]][[2, 2]] - 1]];
report["VAL", "modulus of the phase e^(i 10^20)", Abs[matrixOf[Exp[QuantumOperator[DiagonalMatrix[{1.*^20 I, 1.}]]]][[1, 1]]]];
messageCount["Ising8 Exp", Exp[-I 0.3 ising8]];
report["VAL", "zero operator Exp", matrixOf[Exp[QuantumOperator[SparseArray[{}, {4, 4}], {1, 2}]]]];
same["non-square E^qo", E ^ QuantumOperator[{{1, 0}, {0, 1}, {0, 0}, {0, 0}}]];
report["VAL", "30-digit diagonal: precision of the result", Precision[matrixOf[Exp[QuantumOperator[DiagonalMatrix[N[{1, 2}, 30]]]]]]];
report["VAL", "hidden zero proved: strict test, result", {DiagonalMatrixQ[{{a, Sin[y]^2 + Cos[y]^2 - 1}, {0, b}}, Tolerance -> 0], matrixOf[Exp[QuantumOperator[{{a, Sin[y]^2 + Cos[y]^2 - 1}, {0, b}}]]]}];
report["VAL", "Normal of a machine exponential is packed", Developer`PackedArrayQ[matrixOf[u8]]];
