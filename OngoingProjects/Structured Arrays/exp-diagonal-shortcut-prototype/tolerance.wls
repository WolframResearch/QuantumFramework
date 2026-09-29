(* What DiagonalMatrixQ counts as diagonal, with its default tolerance and with
   Tolerance -> 0 (WL, independent of QF):
     wolframscript -file tolerance.wls *)
report[label_, value_] := WriteString[$Output, label, "\t", ToString[value, InputForm], "\n"];
report["Options[DiagonalMatrixQ]", Options[DiagonalMatrixQ]];
report["default, strict: off-diagonal 5.6e-14 beside a largest entry 2", {DiagonalMatrixQ[#], DiagonalMatrixQ[#, Tolerance -> 0]} & @ {{2., 5.6*^-14}, {0., 1.}}];
report["default, strict: off-diagonal 1e-12 beside a largest entry 2", {DiagonalMatrixQ[#], DiagonalMatrixQ[#, Tolerance -> 0]} & @ {{2., 1.*^-12}, {0., 1.}}];
report["default, strict: 1e-9 between two entries 1, an entry 1e6 elsewhere", {DiagonalMatrixQ[#], DiagonalMatrixQ[#, Tolerance -> 0]} & @ {{1.*^6, 0., 0.}, {0., 1., 1.*^-9}, {0., 1.*^-9, 1.}}];
report["default, strict: the same coupling without the entry 1e6", {DiagonalMatrixQ[#], DiagonalMatrixQ[#, Tolerance -> 0]} & @ {{1., 1.*^-9}, {1.*^-9, 1.}}];
report["strict: machine zeros 0., -0., 0. + 0. I off the diagonal", DiagonalMatrixQ[{{1., 0., -0.}, {0. + 0. I, 1., 0.}, {0., 0., 1.}}, Tolerance -> 0]];
report["strict: Sin[y]^2 + Cos[y]^2 - 1 off the diagonal", DiagonalMatrixQ[{{a, Sin[y]^2 + Cos[y]^2 - 1}, {0, b}}, Tolerance -> 0]];
report["strict: Gamma[g + 1] - g Gamma[g] off the diagonal, and FullSimplify of it", {DiagonalMatrixQ[{{a, Gamma[g + 1] - g Gamma[g]}, {0, b}}, Tolerance -> 0], FullSimplify[Gamma[g + 1] - g Gamma[g]]}];
report["strict: a 2 x 3 matrix with zeros off the diagonal", DiagonalMatrixQ[{{1, 0, 0}, {0, 2, 0}}, Tolerance -> 0]];
