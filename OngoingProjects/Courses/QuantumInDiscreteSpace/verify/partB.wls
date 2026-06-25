PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
show[id_, e_] := id <> " = " <> ToString[InputForm[e]];
SetAttributes[safe, HoldFirst];
safe[e_] := Check[e, "ERR"];
X = QuantumOperator["X"]; Y = QuantumOperator["Y"]; Z = QuantumOperator["Z"];
mat[o_] := o["MatrixRepresentation"] // Normal;

(* B1 Pauli algebra *)
b1a = mat[X];
b1b = (X @ X)["MatrixRepresentation"] // Normal;            (* = I *)
b1c = ((X @ Y)["MatrixRepresentation"] // Normal) == (I (Z["MatrixRepresentation"] // Normal));  (* XY = iZ *)
b1d = ((X @ Y) + (Y @ X))["MatrixRepresentation"] // Normal;  (* anticommutator = 0 *)

(* B2 expectation four ways on psi = (Sqrt3/2)|0> + (1/2)|1> -> <X> = Sqrt3/2 *)
psi = QuantumState[{Sqrt[3]/2, 1/2}];
e1 = Tr[X["Matrix"] . psi["DensityMatrix"]] // FullSimplify;
e2 = QuantumMeasurementOperator[X][psi]["Mean"] // FullSimplify;
e3 = (SuperDagger[psi] @ X @ psi)["Scalar"] // FullSimplify;
e4 = QuantumPartialTrace[X[psi["Operator"]]]["Scalar"] // FullSimplify // safe;
b2 = {e1, e2, e3, e4};

(* B3 eigen-decomposition of Y + outcomes *)
b3a = Y["Eigenvalues"];
b3b = QuantumMeasurementOperator[Y][QuantumState["0"]]["ProbabilitiesList"] // safe;

(* B4 rotations + U3 *)
b4rx = QuantumOperator["RX"[th]]["MatrixRepresentation"] // Normal;
b4rz = QuantumOperator["RZ"[th]]["MatrixRepresentation"] // Normal;
b4u3 = QuantumOperator["U3"[th, ph, lam]]["MatrixRepresentation"] // Normal // FullSimplify;

(* B5 unitary / hermitian checks (symbolic caveat) *)
ry = QuantumOperator["RY"[th]];
b5naive = ry["UnitaryQ"];                                  (* symbolic -> may be False *)
b5proper = Assuming[th \[Element] Reals,
   FullSimplify[(ry["Matrix"] . ry["Dagger"]["Matrix"] // Normal) == IdentityMatrix[2]]];
b5herm = X["HermitianQ"];

(* B6 ZYZ / Euler decomposition of H -> U3 angles in the Label *)
b6angles = QuantumOperator["H"]["ZYZ"]["Label"] // safe;
b6faithful = (QuantumOperator["H"]["ZYZ"]["MatrixRepresentation"] // Normal // FullSimplify) ==
   (QuantumOperator["H"]["MatrixRepresentation"] // Normal);

(* B7 commutator + uncertainty *)
comm = (X @ Y) - (Y @ X);
b7a = comm["MatrixRepresentation"] // Normal;             (* = 2 i Z *)
b7b = b7a == 2 I (Z["MatrixRepresentation"] // Normal);

res = <|
 "B1a-X" -> b1a, "B1b-XX=I" -> b1b, "B1c-XY=iZ" -> b1c, "B1d-anticomm=0" -> b1d,
 "B2-expectX-ways{Tr,Mean,SuperDag,PT}" -> b2,
 "B3a-Yeigen" -> b3a, "B3b-Yon0-probs" -> b3b,
 "B4rx" -> b4rx, "B4rz" -> b4rz, "B4u3" -> b4u3,
 "B5naive-unitaryQ" -> b5naive, "B5proper-simplify" -> b5proper, "B5herm-X" -> b5herm,
 "B6-ZYZ-H-angles" -> b6angles, "B6-faithful" -> b6faithful,
 "B7a-comm-XY" -> b7a, "B7b-=2iZ" -> b7b
|>;
WriteString[$Output, StringRiffle[show @@@ Normal[res], "\n"] <> "\n"];
