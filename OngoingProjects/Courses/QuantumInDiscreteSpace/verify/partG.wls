PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
SetAttributes[rec, HoldRest];
rec[id_, e_] := Module[{v = e}, WriteString[$Output, id <> " = " <> ToString[InputForm[v]] <> "\n"]; v];
SetAttributes[safe, HoldFirst]; safe[e_] := Check[e, "ERR"];

(* G1 clock / shift Weyl-Heisenberg for d=3 *)
Xd = QuantumOperator["X"[3]]; Zd = QuantumOperator["Z"[3]];
rec["G1-X3", Xd["MatrixRepresentation"] // Normal];
rec["G1-Z3", Zd["MatrixRepresentation"] // Normal];
om = Exp[2 Pi I/3];
rec["G1-ZX=wXZ", Simplify[(Zd @ Xd)["MatrixRepresentation"] == om (Xd @ Zd)["MatrixRepresentation"]] // safe];
rec["G1-ZX=w^-1XZ", Simplify[(Zd @ Xd)["MatrixRepresentation"] == Conjugate[om] (Xd @ Zd)["MatrixRepresentation"]] // safe];

(* G2 qutrit SUM gate (qudit CNOT) *)
rec["G2-sum3", QuantumOperator["SUM"[3]]["MatrixRepresentation"] // Normal // safe];
rec["G2-sum|10>", QuantumOperator["SUM"[3]][QuantumTensorProduct[QuantumState["1", 3], QuantumState["0", 3]]]["StateVector"] // Normal // safe];

(* G3 spin-1 operators (named JZ in computational basis; explicit Jx,Jy,Jz for the algebra) *)
rec["G3-jz1-named", QuantumOperator["JZ"[1]]["MatrixRepresentation"] // Normal // safe];
rec["G3-jz-eig", QuantumOperator["JZ"[1]]["Eigenvalues"] // safe];
Jz = QuantumOperator[{{1, 0, 0}, {0, 0, 0}, {0, 0, -1}}];
Jx = QuantumOperator[{{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}/Sqrt[2]];
Jy = QuantumOperator[{{0, -I, 0}, {I, 0, -I}, {0, I, 0}}/Sqrt[2]];
rec["G3-[Jx,Jy]=iJz", Simplify[((Jx @ Jy) - (Jy @ Jx))["MatrixRepresentation"] ==
   I Jz["MatrixRepresentation"]] // safe];
rec["G3-J^2", Simplify[((Jx @ Jx) + (Jy @ Jy) + (Jz @ Jz))["MatrixRepresentation"] // Normal]];

(* G4 qutrit Fourier / mutually unbiased basis *)
w = Exp[2 Pi I/3];
F = (1/Sqrt[3]) Table[w^(j k), {j, 0, 2}, {k, 0, 2}];
rec["G4-fourier", F];
rec["G4-MUB-overlaps", Table[Simplify[Abs[F[[1, k]]]^2], {k, 1, 3}]];      (* all 1/3 *)
rec["G4-unitary", Simplify[ConjugateTranspose[F] . F] == IdentityMatrix[3]];

(* G5 KCBS qutrit contextuality: S = sum_k |<psi|l_k>|^2 = Sqrt5 > 2 *)
ct2 = Cos[Pi/5]/(1 + Cos[Pi/5]); ct = Sqrt[ct2]; st = Sqrt[1 - ct2];
lvec[k_] := {ct, st Sin[4 Pi k/5], st Cos[4 Pi k/5]};
psiK = {1, 0, 0};
rec["G5-ortho", FullSimplify[lvec[0] . lvec[1]]];
gS = FullSimplify[Sum[Abs[Conjugate[psiK] . lvec[k]]^2, {k, 0, 4}]];
rec["G5-S", gS];
rec["G5-N", N[gS]];
rec["G5-=Sqrt5", FullSimplify[gS == Sqrt[5]] // safe];
