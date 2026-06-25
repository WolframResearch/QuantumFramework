p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
ket[bits__] := Fold[Flatten@KroneckerProduct[#1, #2] &, UnitVector[2, # + 1] & /@ {bits}];
H = {{1, 1}, {1, -1}}/Sqrt[2];

(* K1 Deutsch-Jozsa *)
H3 = KroneckerProduct[H, H, H];
dj[oracle_] := H3 . oracle . H3 . ket[0, 0, 0];
p["K1-constant", Abs[dj[IdentityMatrix[8]][[1]]]^2];
p["K1-balanced", Abs[dj[KroneckerProduct[PauliMatrix[3], IdentityMatrix[2], IdentityMatrix[2]]][[1]]]^2];

(* K2 Bernstein-Vazirani, s = 101 (Z on qubits 1 and 3) *)
bvOracle = KroneckerProduct[PauliMatrix[3], IdentityMatrix[2], PauliMatrix[3]];
p["K2-output", H3 . bvOracle . H3 . ket[0, 0, 0] // Chop];

(* K5 QFT *)
p["K5-qft|000>", Abs[FourierMatrix[8] . UnitVector[8, 1]]^2];

(* K6 phase-estimation readout (inverse QFT on a phase-kicked register) *)
phi = 3/8;
kicked = Table[Exp[2 Pi I phi k], {k, 0, 7}]/Sqrt[8];
p["K6-readout", Abs[ConjugateTranspose[FourierMatrix[8]] . kicked]^2 // Chop];

(* K7 HHL: Pauli decomposition (block-encoding ingredient) and A.x *)
A = {{1, 1/2}, {1/2, 1}}; x = {3/5, 4/5};
c = Table[Tr[s . A]/2, {s, Table[PauliMatrix[j], {j, 0, 3}]}];
p["K7-pauli", c];
p["K7-Ax", Sum[c[[j + 1]] PauliMatrix[j] . x, {j, 0, 3}]];

(* K8 VQE: ground energy of Z via RY ansatz *)
ry[\[Theta]_] := MatrixExp[-I \[Theta]/2 PauliMatrix[2]];
cost[\[Theta]_] := Re[Conjugate[#] . PauliMatrix[3] . #] &[ry[\[Theta]] . {1, 0}];
p["K8-VQE", NMinimize[{cost[\[Theta]], -\[Pi] <= \[Theta] <= \[Pi]}, \[Theta]][[1]]];

(* K9 QAOA / Max-Cut: max eigenvalue of the cost Hamiltonian (triangle) *)
Zat[positions_, n_] := KroneckerProduct @@ Table[If[MemberQ[positions, k], PauliMatrix[3], IdentityMatrix[2]], {k, n}];
Hc = Total[1/2 (IdentityMatrix[8] - Zat[#, 3]) & /@ {{1, 2}, {2, 3}, {1, 3}}];
p["K9-maxcut", Max[Eigenvalues[Hc]]];
