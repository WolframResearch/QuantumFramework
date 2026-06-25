p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
ket[bits__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[2, # + 1] & /@ {bits}];

(* F1 controlled gates via projector construction *)
P0 = {{1, 0}, {0, 0}}; P1 = {{0, 0}, {0, 1}};
cnot = KroneckerProduct[P0, IdentityMatrix[2]] + KroneckerProduct[P1, PauliMatrix[1]];
swap = ArrayFlatten[{{P0, {{0, 0}, {0, 0}}}, {{{0, 0}, {0, 0}}, P1}}]; (* placeholder, test below *)
toffoli = KroneckerProduct[IdentityMatrix[4] - KroneckerProduct[P1, P1], IdentityMatrix[2]] +
   KroneckerProduct[P1, P1, PauliMatrix[1]];
p["F1-cnot", cnot];
p["F1-cnot|10>", cnot . ket[1, 0]];
p["F1-toffoli|110>", toffoli . ket[1, 1, 0]];
swapM = 1/2 Sum[KroneckerProduct[s, s], {s, Table[PauliMatrix[j], {j, 0, 3}]}];
p["F1-swap|10>", swapM . ket[1, 0]];

(* F2 compose a circuit (matrix product) and apply *)
H = {{1, 1}, {1, -1}}/Sqrt[2];
bellU = cnot . KroneckerProduct[H, IdentityMatrix[2]];
p["F2-bell", bellU . ket[0, 0]];

(* F4 overall unitary *)
p["F4-unitary", bellU];
