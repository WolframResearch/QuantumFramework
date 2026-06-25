Needs["Wolfram`QuantumFramework`"];
p[l_, e_] := WriteString[$Output, l <> " = " <> ToString[InputForm[e]] <> "\n"];
(* setup toolkit (reused across parts) *)
\[Sigma] = Table[PauliMatrix[j], {j, 3}];
blochVector[\[Psi]_?VectorQ] := Re@Table[Conjugate[\[Psi]] . s . \[Psi], {s, \[Sigma]}];
blochVector[\[Rho]_?MatrixQ] := Re@Table[Tr[\[Rho] . s], {s, \[Sigma]}];
densityMatrix[r_] := 1/2 (IdentityMatrix[2] + r . \[Sigma]);

(* A1 *)
\[Psi] = {Cos[\[Pi]/6], Sin[\[Pi]/6] Exp[I \[Pi]/4]};
p["A1-vec", \[Psi]];
p["A1-born", Abs[\[Psi]]^2];

(* A2 *)
p["A2-norm", Norm[\[Psi]]];
p["A2-gphase", FullSimplify[Abs[{\[Alpha], \[Beta]}]^2 == Abs[Exp[I \[CurlyPhi]] {\[Alpha], \[Beta]}]^2, \[CurlyPhi] \[Element] Reals]];

(* A3 *)
p["A3-plus", blochVector[{1, 1}/Sqrt[2]]];
p["A3-R", blochVector[{1, I}/Sqrt[2]]];

(* A4 *)
p["A4-angles", FullSimplify[blochVector[{Cos[\[Theta]/2], Exp[I \[Phi]] Sin[\[Theta]/2]}], {\[Theta], \[Phi]} \[Element] Reals]];

(* A5 *)
amp = Conjugate[{1, 0}] . ({1, 1}/Sqrt[2]);
p["A5-amp", amp];
p["A5-prob", Abs[amp]^2];

(* A6 *)
p["A6-qutrit0", UnitVector[3, 1]];
p["A6-flat3", Normalize@ConstantArray[1, 3]];
p["A6-flat4", Normalize@ConstantArray[1, 4]];
