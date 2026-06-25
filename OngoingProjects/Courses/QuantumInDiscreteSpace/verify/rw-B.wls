p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
blochVector[\[Psi]_?VectorQ] := Re@Table[Conjugate[\[Psi]] . s . \[Psi], {s, \[Sigma]}];

(* B1 Pauli algebra: product identity + commutator (the example) *)
p["B1-product", Outer[Dot, \[Sigma], \[Sigma], 1] ==
   TensorProduct[IdentityMatrix[3], IdentityMatrix[2]] +
    I TensorContract[TensorProduct[LeviCivitaTensor[3], \[Sigma]], {{3, 4}}]];
p["B1-commutator", Outer[Commutator[#1, #2, Dot] &, \[Sigma], \[Sigma], 1] ==
   2 I TensorContract[TensorProduct[LeviCivitaTensor[3], \[Sigma]], {{3, 4}}]];

(* B2 expectation two ways + Bloch component *)
\[Psi] = {Sqrt[3]/2, 1/2};
p["B2-braket", Conjugate[\[Psi]] . PauliMatrix[1] . \[Psi]];
p["B2-trace", Tr[PauliMatrix[1] . TensorProduct[\[Psi], Conjugate@\[Psi]]]];
p["B2-bloch", blochVector[\[Psi]]];

(* B3 eigenvalues + outcome probabilities (measuring Y on |0>) *)
p["B3-eigs", Eigenvalues[PauliMatrix[2]]];
p["B3-probs", Abs[Conjugate[Normalize[#]] . {1, 0}]^2 & /@ Eigenvectors[PauliMatrix[2]]];

(* B4 rotations + closed-form axis-angle exponential *)
p["B4-Rx", MatrixExp[-I \[Theta]/2 PauliMatrix[1]]];
p["B4-Rz", MatrixExp[-I \[Theta]/2 PauliMatrix[3]]];
p["B4-closedform", With[{n = {Cos[\[Phi]] Sin[\[Theta]], Sin[\[Phi]] Sin[\[Theta]], Cos[\[Theta]]}},
   FullSimplify[MatrixExp[-I \[Alpha]/2 n . \[Sigma]] ==
     Cos[\[Alpha]/2] IdentityMatrix[2] - I Sin[\[Alpha]/2] n . \[Sigma],
    (\[Theta] | \[Phi] | \[Alpha]) \[Element] Reals]]];

(* B5 unitary / hermitian checks (symbolic trap) *)
ry = MatrixExp[-I \[Theta]/2 PauliMatrix[2]];
p["B5-naiveUnitaryQ", UnitaryMatrixQ[ry]];
p["B5-proper", FullSimplify[ry . ConjugateTranspose[ry] == IdentityMatrix[2], \[Theta] \[Element] Reals]];
p["B5-hermX", HermitianMatrixQ[PauliMatrix[1]]];

(* B6 Euler (ZYZ) decomposition *)
{\[Alpha]e, \[Beta]e, \[Gamma]e} = EulerAngles[RotationMatrix[\[Pi]/3, {Cos[\[Pi]/8], Sin[\[Pi]/8], 0}]];
p["B6-angles", {\[Alpha]e, \[Beta]e, \[Gamma]e}];
p["B6-su2", FullSimplify[
   MatrixExp[-I \[Alpha]e/2 PauliMatrix[3]] . MatrixExp[-I \[Beta]e/2 PauliMatrix[2]] .
     MatrixExp[-I \[Gamma]e/2 PauliMatrix[3]] ==
    MatrixExp[-I (\[Pi]/3)/2 {Cos[\[Pi]/8], Sin[\[Pi]/8], 0} . \[Sigma]]]];

(* B7 commutator + Robertson uncertainty (saturated on |0>) *)
p["B7-comm", Commutator[PauliMatrix[1], PauliMatrix[2], Dot] == 2 I PauliMatrix[3]];
\[Psi]0 = {1, 0};
var[A_] := Conjugate[\[Psi]0] . A . A . \[Psi]0 - (Conjugate[\[Psi]0] . A . \[Psi]0)^2;
p["B7-uncert", {Sqrt[var[PauliMatrix[1]] var[PauliMatrix[2]]], Abs[Conjugate[\[Psi]0] . PauliMatrix[3] . \[Psi]0]}];
