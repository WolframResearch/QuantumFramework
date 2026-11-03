{a, ad} = FieldVariables[{1}];
ov[x_, y_] := BosonicMatrixElement[{x, y}, 1, "Basis" -> "Coherent"];
{Simplify[ComplexExpand[ov[al, be] ov[be, al] - Exp[-Abs[al - be]^2] /. {al -> x1 + I y1, be -> x2 + I y2}]],
 Simplify[ComplexExpand[BosonicMatrixElement[{d, d}, Exp[I Pi ad ** a], "Basis" -> "Coherent"] - Exp[-2 Abs[d]^2] /. d -> x + I y]],
 Table[Simplify[Abs[BosonicMatrixElement[{m, 0}, DisplacementOperator[al]]]^2/Pi - Exp[-Abs[al]^2] Abs[al]^(2 m)/(Pi m!)], {m, 0, 5}],
 Table[Simplify[2/Pi (-1)^m BosonicMatrixElement[{m, m}, DisplacementOperator[2 al]] - 2/Pi (-1)^m LaguerreL[m, 4 Abs[al]^2] Exp[-2 Abs[al]^2]], {m, 0, 5}]}
